module zmq_comm_mod
    use iso_c_binding, only: c_ptr, c_null_ptr, c_int, c_size_t, c_char, c_null_char, c_loc, c_associated, c_sizeof
    use mpi
    use com_date
    implicit none

    private
    public :: zmq_comm_initialize
    public :: zmq_comm_finalize
    public :: zmq_receive_command
    public :: zmq_send_force_reply
    public :: zmq_send_error_reply
    public :: zmq_send_status_reply

    integer, parameter :: zmq_buffer_bytes = 8192
    integer(c_int), parameter :: zmq_rep_socket = 4_c_int
    integer(c_int), parameter :: zmq_linger_option = 17_c_int

    type(c_ptr), save :: zmq_context = c_null_ptr
    type(c_ptr), save :: zmq_reply_socket = c_null_ptr
    logical, save :: zmq_ready = .false.

    interface
        function zmq_ctx_new_c() bind(C, name='zmq_ctx_new') result(res)
            import :: c_ptr
            type(c_ptr) :: res
        end function

        function zmq_ctx_term_c(context) bind(C, name='zmq_ctx_term') result(res)
            import :: c_ptr, c_int
            type(c_ptr), value :: context
            integer(c_int) :: res
        end function

        function zmq_socket_c(context, socket_type) bind(C, name='zmq_socket') result(res)
            import :: c_ptr, c_int
            type(c_ptr), value :: context
            integer(c_int), value :: socket_type
            type(c_ptr) :: res
        end function

        function zmq_bind_c(socket, endpoint) bind(C, name='zmq_bind') result(res)
            import :: c_ptr, c_char, c_int
            type(c_ptr), value :: socket
            character(kind=c_char), dimension(*) :: endpoint
            integer(c_int) :: res
        end function

        function zmq_setsockopt_c(socket, option_name, option_value, option_len) bind(C, name='zmq_setsockopt') result(res)
            import :: c_ptr, c_int, c_size_t
            type(c_ptr), value :: socket
            integer(c_int), value :: option_name
            type(c_ptr), value :: option_value
            integer(c_size_t), value :: option_len
            integer(c_int) :: res
        end function

        function zmq_send_c(socket, buffer, buffer_len, flags) bind(C, name='zmq_send') result(res)
            import :: c_ptr, c_char, c_size_t, c_int
            type(c_ptr), value :: socket
            character(kind=c_char), dimension(*) :: buffer
            integer(c_size_t), value :: buffer_len
            integer(c_int), value :: flags
            integer(c_int) :: res
        end function

        function zmq_recv_c(socket, buffer, buffer_len, flags) bind(C, name='zmq_recv') result(res)
            import :: c_ptr, c_char, c_size_t, c_int
            type(c_ptr), value :: socket
            character(kind=c_char), dimension(*) :: buffer
            integer(c_size_t), value :: buffer_len
            integer(c_int), value :: flags
            integer(c_int) :: res
        end function

        function zmq_close_c(socket) bind(C, name='zmq_close') result(res)
            import :: c_ptr, c_int
            type(c_ptr), value :: socket
            integer(c_int) :: res
        end function
    end interface

contains

subroutine zmq_comm_initialize()
    implicit none

    integer :: root_status
    integer(c_int) :: rc
    integer(c_int), target :: linger_value
    character(kind=c_char), allocatable :: endpoint_c(:)

    root_status = 0

    if (myid.eq.0) then
        zmq_context = zmq_ctx_new_c()
        if (.not. c_associated(zmq_context)) then
            root_status = 1
        else
            zmq_reply_socket = zmq_socket_c(zmq_context, zmq_rep_socket)
            if (.not. c_associated(zmq_reply_socket)) then
                root_status = 2
            endif
        endif

        if (root_status.eq.0) then
            linger_value = 0_c_int
            rc = zmq_setsockopt_c(zmq_reply_socket, zmq_linger_option, c_loc(linger_value), c_sizeof(linger_value))
            if (rc.ne.0_c_int) root_status = 3
        endif

        if (root_status.eq.0) then
            call to_c_string(trim(zmq_endpoint), endpoint_c)
            rc = zmq_bind_c(zmq_reply_socket, endpoint_c)
            deallocate(endpoint_c)
            if (rc.ne.0_c_int) root_status = 4
        endif
    endif

    call MPI_BCAST(root_status, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    if (root_status.ne.0) then
        if (myid.eq.0) write(*,*) 'ZeroMQ 初始化失败，状态码: ', root_status
        call MPI_FINALIZE(ierr)
        stop 1
    endif

    zmq_ready = .true.
end subroutine

subroutine zmq_comm_finalize()
    implicit none

    integer(c_int) :: rc

    if (.not. zmq_ready) return

    if (myid.eq.0) then
        if (c_associated(zmq_reply_socket)) then
            rc = zmq_close_c(zmq_reply_socket)
            zmq_reply_socket = c_null_ptr
        endif
        if (c_associated(zmq_context)) then
            rc = zmq_ctx_term_c(zmq_context)
            zmq_context = c_null_ptr
        endif
    endif

    zmq_ready = .false.
end subroutine

subroutine zmq_receive_command(cmd, row, ok)
    implicit none

    integer, intent(out) :: cmd
    real(kind=ps), intent(out) :: row(max_motion_columns)
    logical, intent(out) :: ok

    character(len=zmq_buffer_bytes) :: message
    character(len=256) :: error_text
    integer :: ok_int

    cmd = zmq_cmd_invalid
    row = 0.0d0
    ok = .false.
    ok_int = 0
    error_text = ''

    if (myid.eq.0) then
        call zmq_recv_text(message, ok)
        if (ok) then
            call parse_command_message(message, cmd, row, ok, error_text)
            if (.not. ok) call zmq_send_error_reply(trim(error_text))
        endif
        if (ok) ok_int = 1
    endif

    call MPI_BCAST(cmd, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(row, max_motion_columns, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(ok_int, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    ok = (ok_int.eq.1)
end subroutine

subroutine zmq_send_force_reply()
    implicit none

    character(len=zmq_buffer_bytes) :: message
    character(len=64) :: token
    integer :: i

    if (myid.ne.0) return
    if (.not. zmq_ready) return

    message = 'FORCE'
    write(token, '(i0)') latest_force_row_len
    message = trim(message) // ' ' // trim(token)

    do i = 1, latest_force_row_len
        write(token, '(es24.16)') latest_force_row(i)
        message = trim(message) // ' ' // trim(adjustl(token))
    enddo

    call zmq_send_text(trim(message))
end subroutine

subroutine zmq_send_error_reply(message_text)
    implicit none

    character(len=*), intent(in) :: message_text

    if (myid.ne.0) return
    if (.not. zmq_ready) return
    call zmq_send_text('ERROR ' // trim(message_text))
end subroutine

subroutine zmq_send_status_reply(message_text)
    implicit none

    character(len=*), intent(in) :: message_text

    if (myid.ne.0) return
    if (.not. zmq_ready) return
    call zmq_send_text('STATUS ' // trim(message_text))
end subroutine

subroutine parse_command_message(message, cmd, row, ok, error_text)
    implicit none

    character(len=*), intent(in) :: message
    integer, intent(out) :: cmd
    real(kind=ps), intent(out) :: row(max_motion_columns)
    logical, intent(out) :: ok
    character(len=*), intent(out) :: error_text

    character(len=32) :: keyword
    integer :: ios

    cmd = zmq_cmd_invalid
    row = 0.0d0
    ok = .false.
    error_text = ''
    keyword = ''

    read(message, *, iostat=ios) keyword
    if (ios.ne.0) then
        error_text = '无法解析 ZeroMQ 指令'
        return
    endif

    call uppercase_ascii(keyword)

    select case (trim(keyword))
    case ('RESET', 'START')
        if (len_trim(message) .eq. len_trim(keyword)) then
            row(1) = 0.0d0
        else
            read(message, *, iostat=ios) keyword, row(1)
            if (ios.ne.0) then
                error_text = 'RESET 指令缺少初始时间'
                return
            endif
        endif
        cmd = zmq_cmd_reset
        ok = .true.
    case ('STEP')
        read(message, *, iostat=ios) keyword, row(1:motion_stream_value_count)
        if (ios.ne.0) then
            error_text = 'STEP 指令缺少目标时刻或速度增量'
            return
        endif
        cmd = zmq_cmd_step
        ok = .true.
    case ('OUTPUT', 'SAVE')
        cmd = zmq_cmd_output
        ok = .true.
    case ('STOP', 'EXIT')
        cmd = zmq_cmd_stop
        ok = .true.
    case default
        error_text = '未知 ZeroMQ 指令: ' // trim(keyword)
    end select
end subroutine

subroutine zmq_recv_text(message, ok)
    implicit none

    character(len=*), intent(out) :: message
    logical, intent(out) :: ok

    character(kind=c_char) :: buffer(zmq_buffer_bytes)
    integer(c_int) :: rc
    integer :: i, recv_len

    ok = .false.
    message = ''

    rc = zmq_recv_c(zmq_reply_socket, buffer, int(zmq_buffer_bytes, c_size_t), 0_c_int)
    if (rc.lt.0_c_int) return

    recv_len = min(int(rc), len(message))
    do i = 1, recv_len
        message(i:i) = buffer(i)
    enddo
    if (recv_len.lt.len(message)) message(recv_len+1:) = ' '
    ok = .true.
end subroutine

subroutine zmq_send_text(message)
    implicit none

    character(len=*), intent(in) :: message

    character(kind=c_char), allocatable :: c_message(:)
    integer(c_int) :: rc
    integer :: message_len

    message_len = len_trim(message)
    if (message_len.le.0) return

    allocate(c_message(message_len))
    call copy_to_c_chars(message(1:message_len), c_message)
    rc = zmq_send_c(zmq_reply_socket, c_message, int(message_len, c_size_t), 0_c_int)
    deallocate(c_message)
end subroutine

subroutine to_c_string(text, c_text)
    implicit none

    character(len=*), intent(in) :: text
    character(kind=c_char), allocatable, intent(out) :: c_text(:)
    integer :: i, text_len

    text_len = len_trim(text)
    allocate(c_text(text_len + 1))
    do i = 1, text_len
        c_text(i) = text(i:i)
    enddo
    c_text(text_len + 1) = c_null_char
end subroutine

subroutine copy_to_c_chars(text, c_text)
    implicit none

    character(len=*), intent(in) :: text
    character(kind=c_char), intent(out) :: c_text(:)
    integer :: i

    do i = 1, len(text)
        c_text(i) = text(i:i)
    enddo
end subroutine

subroutine uppercase_ascii(text)
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

end module
