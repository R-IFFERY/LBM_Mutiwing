subroutine init_motion_runtime()
    use mpi
    use com_date
    implicit none

    integer :: env_status, tmp_int
    character(len=256) :: env_value

    motion_source_kind = motion_source_file
    motion_source_name = 'file'
    motion_file = ''
    flow_save_k = 0.0d0
    zmq_endpoint = 'tcp://*:5555'
    zmq_flow_output_interval = 0
    motion_stream_ready = .false.
    motion_stream_start_row = 0.0d0
    motion_stream_target_row = 0.0d0
    latest_force_row = 0.0d0
    latest_force_sum_row = 0.0d0
    latest_force_row_len = 0
    latest_force_sample_count = 0
    episode_index = 0

    if (myid.eq.0) then
        env_value = ''
        call get_environment_variable('LBM_MOTION_SOURCE', env_value, status=env_status)
        if (env_status.eq.0) then
            call uppercase_ascii_text(env_value)
            if (trim(env_value).eq.'ZMQ') then
                motion_source_kind = motion_source_zmq
                motion_source_name = 'zmq'
            endif
        endif

        env_value = ''
        call get_environment_variable('LBM_ZMQ_ENDPOINT', env_value, status=env_status)
        if ((env_status.eq.0).and.(len_trim(env_value).gt.0)) then
            zmq_endpoint = trim(env_value)
        endif

    endif

    call MPI_BCAST(motion_source_kind, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(zmq_endpoint, len(zmq_endpoint), MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)

    if (motion_source_kind.eq.motion_source_zmq) then
        motion_source_name = 'zmq'
    else
        motion_source_name = 'file'
    endif

    if (myid.eq.0) then
        print*, 'Motion source:', trim(motion_source_name)
        if (motion_source_kind.eq.motion_source_zmq) then
            print*, 'ZeroMQ endpoint:', trim(zmq_endpoint)
            print*, 'ZeroMQ flow output mode: explicit SAVE command'
        endif
    endif
end subroutine

subroutine run_zmq_driver(mesh, h)
    use mpi
    use cublas
    use cusolverDn
    use com_date
    use meshtype
    use zmq_comm_mod
    implicit none

    type(mesh_dat) mesh(max_meshid)
    type(cusolverDnHandle) :: h

    integer :: n, nforce, flag2, cmd
    logical :: ok, stop_requested, have_pending_reset, pending_flow_output
    real(kind=ps) :: row(max_motion_columns), pending_reset_row(max_motion_columns), next_row(max_motion_columns)
    character(len=256) :: error_text

    n = 1
    nforce = 101
    flag2 = 1
    stop_requested = .false.
    have_pending_reset = .false.
    pending_flow_output = .false.
    pending_reset_row = 0.0d0

    call zmq_comm_initialize()
    if (myid.eq.0) print*, 'ZeroMQ 仿真进程已启动，等待 RESET / STEP / STOP'

    do while (.not. stop_requested)
        if (.not. have_pending_reset) then
            call zmq_receive_command(cmd, row, ok)
            if (.not. ok) cycle
        else
            cmd = zmq_cmd_reset
            row = pending_reset_row
            have_pending_reset = .false.
        endif

        if (cmd.eq.zmq_cmd_stop) then
            call zmq_send_status_reply('STOPPING')
            stop_requested = .true.
            exit
        endif

        if (cmd.eq.zmq_cmd_output) then
            call zmq_send_error_reply('online mode requires RESET before SAVE')
            cycle
        endif

        if (cmd.ne.zmq_cmd_reset) then
            call zmq_send_error_reply('online mode requires RESET before STEP')
            cycle
        endif

        if (abs(row(1)).gt.1.0d-12) then
            call zmq_send_error_reply('RESET command must provide initial time 0.0')
            cycle
        endif

        episode_index = episode_index + 1
        episode_flow_output_dir = ''
        motion_stream_start_row = 0.0d0
        motion_stream_target_row = 0.0d0
        motion_stream_start_row(1) = row(1)
        motion_stream_target_row(1) = row(1)
        motion_stream_ready = .true.
        motion_end_time = row(1)
        latest_force_row = 0.0d0
        latest_force_sum_row = 0.0d0
        latest_force_row_len = 0
        latest_force_sample_count = 0

        call open_zmq_episode_outputs()
        call reset_zmq_episode(mesh, n, flag2, pending_flow_output)
        call reset_force_reply_average()
        call advance_zmq_segment(mesh, n, nforce, h, flag2, pending_flow_output, row(1), ok, error_text)
        if (.not. ok) then
            call zmq_send_error_reply(trim(error_text))
            call close_zmq_episode_outputs()
            cycle
        endif
        call finalize_force_reply_average()
        call zmq_send_force_reply()

        do
            call zmq_receive_command(cmd, row, ok)
            if (.not. ok) cycle

            if (cmd.eq.zmq_cmd_stop) then
                call zmq_send_status_reply('STOPPING')
                stop_requested = .true.
                exit
            endif

            if (cmd.eq.zmq_cmd_reset) then
                pending_reset_row = row
                have_pending_reset = .true.
                exit
            endif

            if (cmd.eq.zmq_cmd_output) then
                call save_current_flow_output(mesh, n, flag2, ok, error_text)
                if (.not. ok) then
                    call zmq_send_error_reply(trim(error_text))
                else
                    call zmq_send_status_reply(trim(error_text))
                endif
                cycle
            endif

            if (cmd.ne.zmq_cmd_step) then
                call zmq_send_error_reply('expected STEP, SAVE, RESET or STOP')
                cycle
            endif

            motion_stream_start_row = motion_stream_target_row
            next_row = 0.0d0
            call build_stream_target_row(motion_stream_start_row(1:motion_column_count), row(1), row(1:motion_stream_value_count), &
                next_row(1:motion_column_count), ok, error_text)
            if (.not. ok) then
                call zmq_send_error_reply(trim(error_text))
                cycle
            endif
            motion_stream_target_row = 0.0d0
            motion_stream_target_row(1:motion_column_count) = next_row(1:motion_column_count)
            motion_stream_ready = .true.
            motion_end_time = row(1)

            call reset_force_reply_average()
            call advance_zmq_segment(mesh, n, nforce, h, flag2, pending_flow_output, row(1), ok, error_text)
            if (.not. ok) then
                motion_stream_target_row = motion_stream_start_row
                call zmq_send_error_reply(trim(error_text))
                cycle
            endif
            call finalize_force_reply_average()
            call zmq_send_force_reply()
        enddo

        call close_zmq_episode_outputs()
    enddo

    call close_zmq_episode_outputs()
    call zmq_comm_finalize()
end subroutine

subroutine reset_zmq_episode(mesh, n, flag2, pending_flow_output)
    use com_date
    use meshtype
    implicit none

    type(mesh_dat) mesh(max_meshid)
    integer, intent(out) :: n, flag2
    logical, intent(out) :: pending_flow_output

    n = 1
    call init(mesh)
    call read_boundary(mesh(max_meshid))

    body_rel_a = 0.0d0
    body_rel_a_prev = 0.0d0
    body_rel_v = 0.0d0
    body_rel_v_prev = 0.0d0
    body_rel_s = 0.0d0
    body_rel_s_prev = 0.0d0
    body_rel_pos = 0.0d0
    body_rel_theta = 0.0d0
    body_rel_theta_prev = 0.0d0
    body_rel_temp_theta = 0.0d0
    body_rel_temp_theta_prev = 0.0d0
    mom = 0.0d0
    mom2 = 0.0d0
    pm1 = 0.0d0
    pm2 = 0.0d0
    body_extra_v = 0.0d0
    body_extra_s = 0.0d0
    body_extra_pos = 0.0d0
    wing_force = 0.0d0
    wing_moment = 0.0d0
    latest_force_row = 0.0d0
    latest_force_sum_row = 0.0d0
    latest_force_row_len = 0
    latest_force_sample_count = 0
    call update_body_extra_from_rel()

    movenum = 0
    flag2 = 1
    pending_flow_output = .false.
    comj = 0.0d0

    !$acc update device (mesh)
    !$acc update device (myid,numprocs,leftpro,rightpro,layer1,matm,inversematm,e,w,gforce,comj,wing_force,wing_moment)
    !$acc update device (b_point1,b_point,b_segment_count,wing_point_start,wing_segment_start,adjoint_wing_id)
    !$acc update device (sb_arr,sbv_arr,sbac_arr,sa_arr,sau_arr,saac_arr,angle_arr,omega_arr,omgac_arr)
    !$acc update device (body_rel_a,body_rel_v,body_rel_s,body_rel_pos,body_rel_theta,body_rel_temp_theta)
    !$acc update device (body_rel_a_prev,body_rel_v_prev,body_rel_s_prev,body_rel_theta_prev,body_rel_temp_theta_prev)
    !$acc update device (body_extra_v,body_extra_s,body_extra_pos)
    !$acc update device (xb0,yb0,xb1,yb1,xb2,yb2,arc_l,adjoint_n,adjoint_p)
    !$acc update device (mom,mom2,movenum,pm1,pm2)
end subroutine

subroutine reset_force_reply_average()
    use com_date
    implicit none

    latest_force_row = 0.0d0
    latest_force_sum_row = 0.0d0
    latest_force_row_len = 0
    latest_force_sample_count = 0
end subroutine

subroutine finalize_force_reply_average()
    use com_date
    implicit none

    integer :: i, avg_end_col

    if (latest_force_sample_count.le.0) then
        latest_force_row = 0.0d0
        latest_force_row_len = 0
        return
    endif

    avg_end_col = min(latest_force_row_len, 1 + 3 * wing_count + 3)
    do i = 2, avg_end_col
        latest_force_row(i) = latest_force_sum_row(i) / dble(latest_force_sample_count)
    enddo

    ! Keep the segment end time and the last-step motion columns.
    latest_force_row(1) = motion_end_time
end subroutine

subroutine save_current_flow_output(mesh, n, flag2, ok, status_text)
    use com_date
    use meshtype
    implicit none

    type(mesh_dat) mesh(max_meshid)
    integer, intent(in) :: n, flag2
    logical, intent(out) :: ok
    character(len=*), intent(out) :: status_text

    integer :: current_cycle
    character(len=32) :: cycle_text

    ok = .false.
    status_text = ''

	    if (flag2.ne.0) then
	        status_text = '当前时刻处于移网格后的过渡阶段，请在下一次 STEP 后再保存流场'
	        return
	    endif
	
	    current_cycle = max(1, n - 1)
	    call ensure_episode_flow_output_dir(episode_index)
	    !$acc update host(mesh(1)%u,mesh(1)%v,mesh(1)%pre,mesh(1)%body,mesh(1)%dvordt,xb1,yb1)
	    call output(current_cycle, mesh)

    write(cycle_text, '(i0)') current_cycle
    status_text = 'FLOW_SAVED F' // trim(cycle_text) // '.dat'
    ok = .true.
end subroutine

subroutine advance_zmq_segment(mesh, n, nforce, h, flag2, pending_flow_output, target_time, ok, error_text)
    use cublas
    use cusolverDn
    use com_date
    use meshtype
    implicit none

    type(mesh_dat) mesh(max_meshid)
    integer, intent(inout) :: n, flag2
    integer, intent(in) :: nforce
    logical, intent(inout) :: pending_flow_output
    type(cusolverDnHandle) :: h
    real(kind=ps), intent(in) :: target_time
    logical, intent(out) :: ok
    character(len=*), intent(out) :: error_text

    integer :: target_cycle
    real(kind=ps) :: solver_dt, aligned_time

    ok = .false.
    error_text = ''
    solver_dt = mesh(1)%scale * UU
    if (solver_dt.le.0.0d0) then
        error_text = '求解器时间步无效'
        return
    endif

    target_cycle = nint(target_time / solver_dt) + 1
    aligned_time = dble(target_cycle - 1) * solver_dt

    if (abs(aligned_time - target_time).gt.1.0d-10) then
        error_text = 'STEP 时间必须与 LBM 时间步严格对齐'
        return
    endif

    if (target_cycle.lt.n) then
        error_text = 'STEP 时间必须大于等于当前待推进时刻'
        return
    endif

    do while (n.le.target_cycle)
        call run_zmq_time_step(mesh, n, nforce, h, flag2, pending_flow_output)
    enddo

    ok = .true.
end subroutine

subroutine run_zmq_time_step(mesh, n, nforce, h, flag2, pending_flow_output)
    use mpi
    use cublas
    use cusolverDn
    use com_date
    use meshtype
    implicit none

    type(mesh_dat) mesh(max_meshid)
    integer, intent(inout) :: n, flag2
    integer, intent(in) :: nforce
    logical, intent(inout) :: pending_flow_output
    type(cusolverDnHandle) :: h

    time = dble(n - 1) * (mesh(1)%scale * UU)
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    if (n.ge.2) call erranalysis(mesh, n)
    call evolution(mesh, n, nforce, h)

    if (body_rel_pos(1).ge.1.0d0) then
        call move_mesh(1, mesh)
        body_rel_pos(1) = body_rel_pos(1) - 1.0d0
        call update_body_extra_from_rel()
        movenum = movenum + 1
        flag2 = 3
    elseif (body_rel_pos(1).le.-1.0d0) then
        call move_mesh(-1, mesh)
        body_rel_pos(1) = body_rel_pos(1) + 1.0d0
        call update_body_extra_from_rel()
        movenum = movenum - 1
        flag2 = 3
    endif

    if (flag2.gt.0) flag2 = flag2 - 1
    if ((myid.eq.0).and.(mod(n,100).eq.0)) then
        print*, '在线模式已迭代', n, '步，误差为：', err
    endif

    n = n + 1
end subroutine

subroutine open_zmq_episode_outputs()
    use com_date
    implicit none

    character(len=128) :: force_name, motion_name
    character(len=512) :: out_path

    call close_zmq_episode_outputs()

    write(force_name, '("force_episode_",i0,".dat")') episode_index
    write(motion_name, '("airfoil_motion_episode_",i0,".dat")') episode_index

    if (myid.eq.0) then
        call build_output_path(trim(force_name), out_path)
        open(101, file=trim(out_path), status='replace')
        call build_output_path(trim(motion_name), out_path)
        open(102, file=trim(out_path), status='replace')
    endif
end subroutine

subroutine close_zmq_episode_outputs()
    use com_date
    implicit none

    logical :: is_open

    if (myid.ne.0) return

    inquire(unit=101, opened=is_open)
    if (is_open) close(101)
    inquire(unit=102, opened=is_open)
    if (is_open) close(102)
end subroutine

subroutine uppercase_ascii_text(text)
    implicit none

    character(len=*), intent(inout) :: text
    integer :: i, code

    do i = 1, len_trim(text)
        code = iachar(text(i:i))
        if ((code.ge.iachar('a')).and.(code.le.iachar('z'))) then
            text(i:i) = achar(code - 32)
        endif
    enddo
end subroutine
