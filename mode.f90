subroutine init_output_dir
    use com_date
    implicit none

    integer :: values(8), ios, unit_id
    logical :: exists
    character(len=256) :: latest_file, case_name
    character(len=512) :: cmd

    output_root = 'output'
    output_dir = ''
    episode_flow_output_dir = ''
    restart_input_dir = ''
    latest_file = trim(output_root) // '/latest_run.txt'

    cmd = 'mkdir -p ' // trim(output_root)
    call system(trim(cmd))

    inquire(file=trim(latest_file), exist=exists)
    if (exists) then
        unit_id = 902
        open(unit_id, file=trim(latest_file), status='old', action='read', iostat=ios)
        if (ios .eq. 0) then
            read(unit_id, '(A)', iostat=ios) restart_input_dir
            if (ios .ne. 0) restart_input_dir = ''
            close(unit_id)
        endif
    endif

    if (motion_source_kind.eq.motion_source_file) then
        call resolve_offline_motion_file(case_name)
        output_dir = trim(output_root) // '/' // trim(case_name)
    else
        call date_and_time(values=values)
        write(output_dir,'(a,"/run_",i4.4,i2.2,i2.2,"_",i2.2,i2.2,i2.2,"_",i3.3)') &
            trim(output_root), values(1), values(2), values(3), values(5), values(6), values(7), values(8)
    endif

    cmd = 'mkdir -p ' // trim(output_dir)
    call system(trim(cmd))

    unit_id = 903
    open(unit_id, file=trim(latest_file), status='replace', action='write', iostat=ios)
    if (ios .eq. 0) then
        write(unit_id, '(A)') trim(output_dir)
        close(unit_id)
    endif
end subroutine

subroutine resolve_offline_motion_file(case_name)
    use com_date
    implicit none

    character(len=*), intent(out) :: case_name
    integer :: values(8), unit_id, ios, count, dot_pos, i, k_pos, st_pos, slash_pos, env_status
    logical :: exists
    character(len=64) :: k_text
    character(len=256) :: filename
    character(len=512) :: line, list_file, cmd, selected_file

    selected_file = ''
    if (command_argument_count() .ge. 1) then
        call get_command_argument(1, selected_file)
    endif
    if (len_trim(selected_file) .eq. 0) then
        call get_environment_variable('LBM_MOTION_FILE', selected_file, status=env_status)
        if (env_status .ne. 0) selected_file = ''
    endif

    if (len_trim(selected_file) .gt. 0) then
        if (index(selected_file, '/') .eq. 0) then
            motion_file = 'input/' // trim(selected_file)
        else
            motion_file = trim(selected_file)
        endif

        inquire(file=trim(motion_file), exist=exists)
        if (.not. exists) then
            write(*,*) '指定的离线运动输入文件不存在: ', trim(motion_file)
            stop 1
        endif

        slash_pos = 0
        do i = len_trim(motion_file), 1, -1
            if (motion_file(i:i) .eq. '/') then
                slash_pos = i
                exit
            endif
        enddo
        filename = motion_file(slash_pos+1:len_trim(motion_file))
    else
        call date_and_time(values=values)
        write(list_file,'("/tmp/lbm_input_files_",i4.4,i2.2,i2.2,"_",i2.2,i2.2,i2.2,"_",i3.3,".txt")') &
            values(1), values(2), values(3), values(5), values(6), values(7), values(8)

        cmd = 'find input -maxdepth 1 -type f -printf "%f\n" > ' // trim(list_file)
        call system(trim(cmd))

        unit_id = 904
        open(unit_id, file=trim(list_file), status='old', action='read', iostat=ios)
        if (ios .ne. 0) then
            write(*,*) '无法读取 input 文件列表，请确认当前目录下存在 input 文件夹'
            stop 1
        endif

        count = 0
        filename = ''
        do
            read(unit_id, '(A)', iostat=ios) line
            if (ios .ne. 0) exit
            if (len_trim(line) .eq. 0) cycle
            count = count + 1
            filename = trim(line)
        enddo
        close(unit_id)

        cmd = 'rm -f ' // trim(list_file)
        call system(trim(cmd))

        if (count .lt. 1) then
            write(*,*) 'input 文件夹中没有离线运动输入文件'
            stop 1
        elseif (count .gt. 1) then
            write(*,*) 'input 文件夹中只能包含一个离线运动输入文件，当前数量为: ', count
            stop 1
        endif

        motion_file = 'input/' // trim(filename)
    endif
    k_pos = 0
    st_pos = 0
    do i = 1, len_trim(filename)
        if ((filename(i:i) .eq. 'k') .or. (filename(i:i) .eq. 'K')) then
            k_pos = i
            exit
        endif
    enddo
    if (k_pos .gt. 0) then
        do i = k_pos + 1, len_trim(filename) - 1
            if (((filename(i:i) .eq. 's') .or. (filename(i:i) .eq. 'S')) .and. &
                ((filename(i+1:i+1) .eq. 't') .or. (filename(i+1:i+1) .eq. 'T'))) then
                st_pos = i
                exit
            endif
        enddo
    endif
    if ((k_pos .lt. 1) .or. (st_pos .le. k_pos + 1)) then
        write(*,*) '无法从输入文件名解析流场保存用的 k 值，文件名应类似 k2.0st0.2.dat: ', trim(filename)
        stop 1
    endif

    k_text = filename(k_pos+1:st_pos-1)
    read(k_text, *, iostat=ios) flow_save_k
    if ((ios .ne. 0) .or. (flow_save_k .le. 0.0d0)) then
        write(*,*) '输入文件名中的 k 值无效: ', trim(k_text), ' file=', trim(filename)
        stop 1
    endif

    case_name = trim(filename)
    dot_pos = 0
    do i = len_trim(case_name), 1, -1
        if (case_name(i:i) .eq. '.') then
            dot_pos = i
            exit
        endif
    enddo
    if (dot_pos .gt. 1) case_name = case_name(1:dot_pos-1)
end subroutine

subroutine build_output_path(filename, path)
    use com_date
    implicit none

    character(len=*), intent(in) :: filename
    character(len=*), intent(out) :: path

    path = trim(output_dir) // '/' // trim(filename)
end subroutine

subroutine build_flow_output_path(filename, path)
    use com_date
    implicit none

    character(len=*), intent(in) :: filename
    character(len=*), intent(out) :: path

    if (len_trim(episode_flow_output_dir) .gt. 0) then
        path = trim(episode_flow_output_dir) // '/' // trim(filename)
    else
        path = trim(output_dir) // '/' // trim(filename)
    endif
end subroutine

subroutine ensure_episode_flow_output_dir(episode_id)
    use mpi
    use com_date
    implicit none

    integer, intent(in) :: episode_id
    character(len=64) :: dirname
    character(len=512) :: cmd

    if (myid.eq.0) then
        if (len_trim(episode_flow_output_dir).eq.0) then
            write(dirname, '("episode_",i0)') episode_id
            episode_flow_output_dir = trim(output_dir) // '/' // trim(dirname)
            cmd = 'mkdir -p ' // trim(episode_flow_output_dir)
            call system(trim(cmd))
        endif
    endif

    call MPI_BCAST(episode_flow_output_dir, len(episode_flow_output_dir), MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
end subroutine

subroutine build_restart_input_path(filename, path)
    use com_date
    implicit none

    character(len=*), intent(in) :: filename
    character(len=*), intent(out) :: path

    if (len_trim(restart_input_dir) .gt. 0) then
        path = trim(restart_input_dir) // '/' // trim(filename)
    else
        path = trim(filename)
    endif
end subroutine

subroutine load_motion_table
    use com_date
    implicit none

    integer :: unit_id, ios, count, i
    real(kind=ps) :: row(max_motion_columns)
    character(len=512) :: line

    if (motion_loaded) return

    if (len_trim(motion_file) .eq. 0) then
        call resolve_offline_motion_file(line)
    endif
    unit_id = 901
    count = 0

    open(unit_id, file=trim(motion_file), status='old', action='read', iostat=ios)
    if (ios .ne. 0) then
        write(*,*) '无法打开运动参数文件: ', trim(motion_file)
        stop 1
    endif

    do
        read(unit_id, '(A)', iostat=ios) line
        if (ios .ne. 0) exit
        if (len_trim(line) .eq. 0) cycle
        row = 0.0d0
        read(line, *, iostat=ios) row(1:motion_column_count)
        if (ios .ne. 0) then
            write(*,*) '运动参数文件格式错误，列数应为: ', motion_column_count, trim(motion_file)
            stop 1
        endif
        count = count + 1
    enddo

    if (count .lt. 1) then
        write(*,*) '运动参数文件为空: ', trim(motion_file)
        stop 1
    endif

    rewind(unit_id)
    allocate(motion_table(count, motion_column_count))

    i = 0
    do
        read(unit_id, '(A)', iostat=ios) line
        if (ios .ne. 0) exit
        if (len_trim(line) .eq. 0) cycle
        i = i + 1
        row = 0.0d0
        read(line, *, iostat=ios) row(1:motion_column_count)
        motion_table(i, :) = row(1:motion_column_count)
        if (ios .ne. 0) then
            write(*,*) '运动参数文件格式错误，列数应为: ', motion_column_count, trim(motion_file)
            stop 1
        endif
    enddo
    close(unit_id)

    motion_points = count
    motion_index = 1
    motion_end_time = motion_table(motion_points, 1)
    motion_loaded = .true.

    do i = 2, motion_points
        if (motion_table(i, 1) .lt. motion_table(i-1, 1)) then
            write(*,*) '运动参数文件时间列必须非递减: ', trim(motion_file)
            stop 1
        endif
    enddo

    if (motion_points .gt. 1) then
        do i = 2, motion_points
            if (motion_table(i, 1) .gt. motion_table(i-1, 1)) exit
        enddo
    endif
end subroutine

subroutine apply_motion_row(row)
    use com_date
    implicit none

    real(kind=ps), intent(in) :: row(motion_column_count)
    integer :: wing, col
    real(kind=ps) :: x_base, y_base

    angle_arr = 0.0d0
    omega_arr = 0.0d0
    omgac_arr = 0.0d0
    sa_arr = 0.0d0
    sau_arr = 0.0d0
    saac_arr = 0.0d0
    sb_arr = 0.0d0
    sbv_arr = 0.0d0
    sbac_arr = 0.0d0

    do wing = 1, wing_count
        col = 2 + (wing - 1) * motion_columns_per_wing
        angle_arr(wing) = row(col)
        omega_arr(wing) = row(col + 1)
        omgac_arr(wing) = row(col + 2)
        sa_arr(wing) = row(col + 3)
        sau_arr(wing) = row(col + 4)
        saac_arr(wing) = row(col + 5)
        sb_arr(wing) = row(col + 6)
        sbv_arr(wing) = row(col + 7)
        sbac_arr(wing) = row(col + 8)

        if (motion_source_kind.eq.motion_source_file) then
            x_base = wing_initial_x_origin + dble(wing - 1) * wing_initial_x_spacing
            y_base = wing_initial_y_origin
            sa_arr(wing) = sa_arr(wing) - x_base
            sb_arr(wing) = sb_arr(wing) - y_base
        endif
    enddo

end subroutine

real(kind=ps) function smooth_r(x)
    use com_date
    implicit none

    real(kind=ps), intent(in) :: x

    if (x .lt. 0.0d0) then
        smooth_r = 0.0d0
    elseif (x .le. 1.0d0) then
        smooth_r = 0.5d0 * x * x + (cos(2.0d0 * pai * x) - 1.0d0) / (4.0d0 * pai * pai)
    else
        smooth_r = x - 0.5d0
    endif
end function

real(kind=ps) function smooth_h(x)
    use com_date
    implicit none

    real(kind=ps), intent(in) :: x

    if (x .lt. 0.0d0) then
        smooth_h = 0.0d0
    elseif (x .le. 1.0d0) then
        smooth_h = x - sin(2.0d0 * pai * x) / (2.0d0 * pai)
    else
        smooth_h = 1.0d0
    endif
end function

real(kind=ps) function smooth_d(x)
    use com_date
    implicit none

    real(kind=ps), intent(in) :: x

    if ((x .ge. 0.0d0) .and. (x .le. 1.0d0)) then
        smooth_d = 1.0d0 - cos(2.0d0 * pai * x)
    else
        smooth_d = 0.0d0
    endif
end function

subroutine interpolate_motion_triplet(s0, v0, v1, dt, local_t, s, v, a)
    use com_date
    implicit none

    real(kind=ps), intent(in) :: s0, v0, v1, dt, local_t
    real(kind=ps), intent(out) :: s, v, a
    real(kind=ps) :: delta_v, x
    real(kind=ps), external :: smooth_r, smooth_h, smooth_d

    delta_v = v1 - v0

    if (dt .le. 1.0d-14) then
        s = (v0 + delta_v) * local_t + s0
        v = v0 + delta_v
        a = 0.0d0
        return
    endif

    x = local_t / dt
    s = delta_v * dt * smooth_r(x) + v0 * local_t + s0
    v = delta_v * smooth_h(x) + v0
    a = delta_v / dt * smooth_d(x)
end subroutine

subroutine build_stream_target_row(start_row, target_time, delta_row, target_row, ok, error_text)
    use com_date
    implicit none

    real(kind=ps), intent(in) :: start_row(motion_column_count)
    real(kind=ps), intent(in) :: target_time
    real(kind=ps), intent(in) :: delta_row(*)
    real(kind=ps), intent(out) :: target_row(motion_column_count)
    logical, intent(out) :: ok
    character(len=*), intent(out) :: error_text

    real(kind=ps) :: dt
    integer :: wing, col, delta_col

    ok = .false.
    error_text = ''
    dt = target_time - start_row(1)
    if (dt .lt. -1.0d-12) then
        error_text = 'STEP 时间必须大于等于当前段起始时刻'
        return
    endif

    target_row = start_row
    target_row(1) = target_time

    do wing = 1, wing_count
        col = 2 + (wing - 1) * motion_columns_per_wing
        delta_col = 2 + (wing - 1) * 3

        call interpolate_motion_triplet(start_row(col), start_row(col + 1), start_row(col + 1) + delta_row(delta_col), &
            dt, max(dt, 0.0d0), target_row(col), target_row(col + 1), target_row(col + 2))
        call interpolate_motion_triplet(start_row(col + 3), start_row(col + 4), start_row(col + 4) + delta_row(delta_col + 1), &
            dt, max(dt, 0.0d0), target_row(col + 3), target_row(col + 4), target_row(col + 5))
        call interpolate_motion_triplet(start_row(col + 6), start_row(col + 7), start_row(col + 7) + delta_row(delta_col + 2), &
            dt, max(dt, 0.0d0), target_row(col + 6), target_row(col + 7), target_row(col + 8))
    enddo

    ok = .true.
end subroutine

subroutine interpolate_motion_row(i1, i2, target_time, row)
    use com_date
    implicit none

    integer, intent(in) :: i1, i2
    real(kind=ps), intent(in) :: target_time
    real(kind=ps), intent(out) :: row(motion_column_count)
    real(kind=ps) :: t1, t2, dt, ratio
    real(kind=ps), parameter :: time_tol = 1.0d-14

    t1 = motion_table(i1, 1)
    t2 = motion_table(i2, 1)
    dt = t2 - t1

    if (abs(dt) .le. time_tol) then
        row = motion_table(i2, :)
        row(1) = target_time
        return
    endif

    ratio = (target_time - t1) / dt
    ratio = min(max(ratio, 0.0d0), 1.0d0)
    row = motion_table(i1, :) + ratio * (motion_table(i2, :) - motion_table(i1, :))
    row(1) = target_time
end subroutine

subroutine lookup_motion_row(target_time, row)
    use com_date
    implicit none

    real(kind=ps), intent(in) :: target_time
    real(kind=ps), intent(out) :: row(motion_column_count)
    real(kind=ps), parameter :: time_tol = 1.0d-12

    if (target_time .le. motion_table(1, 1) + time_tol) then
        motion_index = 1
        row = motion_table(1, :)
        row(1) = target_time
        return
    endif

    if (target_time .ge. motion_table(motion_points, 1) - time_tol) then
        motion_index = motion_points
        row = motion_table(motion_points, :)
        row(1) = target_time
        return
    endif

    motion_index = max(1, min(motion_index, motion_points - 1))

    do while ((motion_index .gt. 1) .and. (target_time .le. motion_table(motion_index, 1) + time_tol))
        motion_index = motion_index - 1
    enddo

    do while ((motion_index .lt. motion_points - 1) .and. (target_time .gt. motion_table(motion_index + 1, 1) + time_tol))
        motion_index = motion_index + 1
    enddo

    call interpolate_motion_row(motion_index, motion_index + 1, target_time, row)
end subroutine

subroutine build_motion_step_table(dt)
    use com_date
    implicit none

    real(kind=ps), intent(in) :: dt
    integer :: i
    real(kind=ps) :: target_time
    real(kind=ps) :: row(max_motion_columns)

    if (.not. motion_loaded) call load_motion_table()
    if (dt .le. 0.0d0) stop '仿真时间步必须大于0'

    motion_dt = dt
    motion_step_count = int((motion_end_time + 1.0d-12) / motion_dt) + 1

    if (allocated(motion_step_table)) deallocate(motion_step_table)
    allocate(motion_step_table(motion_step_count, motion_column_count))

    motion_index = 1
    do i = 1, motion_step_count
        target_time = dble(i - 1) * motion_dt
        call lookup_motion_row(target_time, row)
        motion_step_table(i, :) = row(1:motion_column_count)
    enddo
end subroutine

subroutine lookup_stream_motion_row(target_time, row)
    use com_date
    implicit none

    real(kind=ps), intent(in) :: target_time
    real(kind=ps), intent(out) :: row(motion_column_count)
    real(kind=ps) :: t1, t2, dt, local_t
    integer :: wing, col
    real(kind=ps), parameter :: time_tol = 1.0d-12

    if (.not. motion_stream_ready) stop '在线运动参数尚未初始化'

    t1 = motion_stream_start_row(1)
    t2 = motion_stream_target_row(1)
    dt = t2 - t1

    if (target_time .le. t1 + time_tol) then
        row = motion_stream_start_row(1:motion_column_count)
        row(1) = target_time
        return
    endif

    local_t = min(max(target_time - t1, 0.0d0), max(dt, 0.0d0))
    row = 0.0d0
    row(1) = target_time

    do wing = 1, wing_count
        col = 2 + (wing - 1) * motion_columns_per_wing
        call interpolate_motion_triplet(motion_stream_start_row(col), motion_stream_start_row(col + 1), motion_stream_target_row(col + 1), &
            dt, local_t, row(col), row(col + 1), row(col + 2))
        call interpolate_motion_triplet(motion_stream_start_row(col + 3), motion_stream_start_row(col + 4), motion_stream_target_row(col + 4), &
            dt, local_t, row(col + 3), row(col + 4), row(col + 5))
        call interpolate_motion_triplet(motion_stream_start_row(col + 6), motion_stream_start_row(col + 7), motion_stream_target_row(col + 7), &
            dt, local_t, row(col + 6), row(col + 7), row(col + 8))
    enddo
end subroutine

subroutine mode
    use com_date
    implicit none

    real(kind=ps) :: row(max_motion_columns)

    if (motion_source_kind.eq.motion_source_zmq) then
        call lookup_stream_motion_row(time, row(1:motion_column_count))
        call apply_motion_row(row(1:motion_column_count))
        return
    endif

    if (.not. motion_loaded) call load_motion_table()
    if (.not. allocated(motion_step_table)) stop '运动参数步表尚未初始化'
    if (motion_dt .le. 0.0d0) stop '运动参数时间步无效'

    motion_step_count = size(motion_step_table, 1)
    motion_index = nint(time / motion_dt) + 1
    motion_index = max(1, min(motion_index, motion_step_count))
    call apply_motion_row(motion_step_table(motion_index, :))
end subroutine

subroutine update_body_extra_from_rel
    use com_date
    implicit none

    integer :: wing
    real(kind=ps) :: accum_v, accum_s, accum_pos

    body_extra_v = 0.0d0
    body_extra_s = 0.0d0
    body_extra_pos = 0.0d0
    accum_v = 0.0d0
    accum_s = 0.0d0
    accum_pos = 0.0d0

    do wing = 1, wing_count
        accum_v = accum_v + body_rel_v(wing)
        accum_s = accum_s + body_rel_s(wing)
        accum_pos = accum_pos + body_rel_pos(wing)
        body_extra_v(wing) = accum_v
        body_extra_s(wing) = accum_s
        body_extra_pos(wing) = accum_pos
    enddo
end subroutine
