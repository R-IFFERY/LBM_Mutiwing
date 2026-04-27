subroutine move_mesh(dr,mesh)
    use mpi
    use com_date
    use meshtype
    implicit none
    
    type(mesh_dat) mesh(max_meshid)
    integer dr,i,l,j,k,m,i1,j1,i2,l1,l2,sp,ep,num,numt,status(MPI_STATUS_SIZE)
    Real (kind=ps), allocatable :: temp1(:,:),temp2(:,:)
    Real (kind=ps) tempf(Q),temp(2:mesh(1)%Ny-1,3),buf22(3*(mesh(1)%Ny-1))
    real(kind=ps), allocatable :: buf1(:),buf2(:)
    
    !$acc enter data create(temp)
    l=1
    if(dr.eq.1)then
        if(myid.eq.0)print*,'右移一格'
        do i=1,max_meshid
            if(mesh(l)%layer.eq.1)then
                if(numprocs.eq.1)then
                    if((l.eq.1).and.(periodical_bc.eq.1))then
                        !$acc kernels loop present(mesh,temp)
                        do j=2,mesh(l)%Ny-1
                            temp(j,1)=mesh(l)%u(2,j)
                            temp(j,2)=mesh(l)%v(2,j)
                            temp(j,3)=mesh(l)%rho(2,j)
                        enddo
                        !$acc end kernels loop
                    endif
                    !$acc kernels present(mesh)
                    !$acc loop independent 
                    do j=1,mesh(l)%Ny
                        !$acc loop independent 
                        do k=1,mesh(l)%Nx-1
                            mesh(l)%f0(k,j,:)=mesh(l)%f(k+1,j,:)
                            if(k.ge.1.and.k.le.mesh(l)%Nx-1)then
                                mesh(l)%u0(k,j)=mesh(l)%u(k+1,j)
                                mesh(l)%v0(k,j)=mesh(l)%v(k+1,j)
                                mesh(l)%pre(k,j)=mesh(l)%rho(k+1,j)
                            endif
                        enddo
                        !$acc end loop
                    enddo
                    !$acc end loop
                    !$acc loop independent 
                    do j=1,mesh(l)%Ny
                        !$acc loop independent 
                        do k=1,mesh(l)%Nx-1
                            mesh(l)%f(k,j,:)=mesh(l)%f0(k,j,:)
                            if(k.ge.1.and.k.le.mesh(l)%Nx-1)then
                                mesh(l)%u(k,j)=mesh(l)%u0(k,j)
                                mesh(l)%v(k,j)=mesh(l)%v0(k,j)
                                mesh(l)%rho(k,j)=mesh(l)%pre(k,j)
                            endif
                        enddo
                        !$acc end loop
                    enddo
                    !$acc end loop
                    !$acc end kernels
                else
                    if((periodical_bc.eq.1).and.(l.eq.1))then
                        numt=1
                        if(myid.eq.0)then
                            do j=2,mesh(l)%Ny-1
                                buf22(numt)=mesh(l)%u(2,j)
                                buf22(numt+1)=mesh(l)%v(2,j)
                                buf22(numt+2)=mesh(l)%rho(2,j)
                                numt=numt+3
                            enddo
                            numt=numt-1
                        endif
                    endif
                    if(myid.eq.0)then
                        do k=1,mesh(l)%Nx-1
                            mesh(l)%f(k,1,:)=mesh(l)%f(k+1,1,:)
                            mesh(l)%f(k,mesh(l)%Ny,:)=mesh(l)%f(k+1,mesh(l)%Ny,:)
                        enddo
                    endif
                    if(myid.eq.0) then
			            sp=1
			            ep=mesh(l)%Nx
		            else
			            sp=mesh(l)%start_row
			            ep=mesh(l)%end_row
                    endif
                    call transfanduandvandrho(1,dr,mesh(l)%f,mesh(l)%u,mesh(l)%v,mesh(l)%rho,sp,ep,mesh(l)%start_row,mesh(l)%end_row,mesh(l)%Nx,mesh(l)%Ny)
                    if(myid.eq.numprocs-1)then
                        do j=2,mesh(l)%Ny
                            mesh(l)%f(mesh(l)%Nx-1,j,:)=mesh(l)%f(mesh(l)%Nx,j,:)
                        enddo
                    endif
                    

                    if(myid.eq.0) then
			            sp=1
			            ep=mesh(l)%Nx
		            else
			            sp=mesh(l)%start_row
			            ep=mesh(l)%end_row
                    endif
                    
                    call scatter(mesh(l)%f,sp,ep,mesh(l)%Nx,mesh(l)%Ny)
                endif
            endif
        enddo
        if(periodical_bc.eq.1)then
            if(numprocs.eq.1)then
                !$acc kernels present(mesh,temp)
                !$acc loop
                do j=2,mesh(l)%Ny-1
                    mesh(1)%u(mesh(1)%Nx-1,j)=temp(j,1)
                    mesh(1)%v(mesh(1)%Nx-1,j)=temp(j,2)
                    mesh(1)%rho(mesh(1)%Nx-1,j)=temp(j,3)
                enddo
                !$acc end loop
                !$acc loop
                do k=1,Q
                    do j=1,mesh(1)%Ny
                        mesh(1)%f(mesh(1)%Nx,j,k)=mesh(1)%f(2,j,k)
                        mesh(1)%f(1,j,k)=mesh(1)%f(mesh(1)%Nx-1,j,k)
                    enddo
                enddo
                !$acc end loop
                !$acc end kernels
            else
                if(allocated(buf2)) deallocate(buf2)
                allocate(buf2(3*(mesh(1)%Ny-2)))
                numt=1
                if(myid.eq.0)then
                    call MPI_SEND(buf22,3*(mesh(1)%Ny-2),MPI_DOUBLE_PRECISION,numprocs-1,100,MPI_COMM_WORLD,ierr)
                elseif(myid.eq.numprocs-1)then
                    call MPI_RECV(buf2,3*(mesh(1)%Ny-2),MPI_DOUBLE_PRECISION,0,100,MPI_COMM_WORLD,status,ierr)
                    numt=1
                    do j=2,mesh(1)%Ny-1
                        mesh(1)%u(mesh(1)%Nx-1,j)=buf2(numt)
                        mesh(1)%v(mesh(1)%Nx-1,j)=buf2(numt+1)
                        mesh(1)%rho(mesh(1)%Nx-1,j)=buf2(numt+2)
                        numt=numt+3
                    enddo
                endif
                if(allocated(buf1)) deallocate(buf1)
                if(allocated(buf2)) deallocate(buf2)
                allocate(buf1(Q*mesh(1)%Ny),buf2(Q*mesh(1)%Ny))
                numt=1
                if(myid.eq.0)then
                    do j=1,mesh(1)%Ny
                        buf1(numt:numt+Q-1)=mesh(1)%f(2,j,:)
                        numt=numt+Q
                    enddo
                    numt=numt-1
                    call MPI_SENDRECV(buf1,numt,MPI_DOUBLE_PRECISION,numprocs-1,3,buf2,Q*mesh(1)%Ny,MPI_DOUBLE_PRECISION,numprocs-1,3,MPI_COMM_WORLD,status,ierr)
                    numt=1
                    do j=1,mesh(1)%Ny
                        mesh(1)%f(1,j,:)=buf2(numt:numt+Q-1)
                        numt=numt+Q
                    enddo                        
                elseif(myid.eq.numprocs-1)then
                    do j=1,mesh(1)%Ny
                        buf1(numt:numt+Q-1)=mesh(1)%f(mesh(1)%Nx-1,j,:)
                        numt=numt+Q
                    enddo
                    numt=numt-1
                    call MPI_SENDRECV(buf1,numt,MPI_DOUBLE_PRECISION,0,3,buf2,Q*mesh(1)%Ny,MPI_DOUBLE_PRECISION,0,3,MPI_COMM_WORLD,status,ierr)
                    numt=1
                    do j=1,mesh(1)%Ny
                        mesh(1)%f(mesh(1)%Nx,j,:)=buf2(numt:numt+Q-1)
                        numt=numt+Q
                    enddo    
                endif
            endif
        endif
    elseif(dr.eq.-1)then
        if(myid.eq.0)print*,'左移一格'
        do i=1,max_meshid
            if(mesh(l)%layer.eq.1)then
                if(numprocs.eq.1)then
                    if((l.eq.1).and.(periodical_bc.eq.1))then
                        !$acc kernels loop present(temp,mesh)
                        do j=2,mesh(l)%Ny-1
                            temp(j,1)=mesh(l)%u(mesh(l)%Nx-1,j)
                            temp(j,2)=mesh(l)%v(mesh(l)%Nx-1,j)
                            temp(j,3)=mesh(l)%rho(mesh(l)%Nx-1,j)
                        enddo
                        !$acc end kernels loop
                    endif
                    !$acc kernels present(mesh)
                    !$acc loop independent
                    do j=1,mesh(l)%Ny
                        !$acc loop independent
                        do k=2,mesh(l)%Nx
                            mesh(l)%f0(k,j,:)=mesh(l)%f(k-1,j,:)
                            if(k.ge.2.and.k.le.mesh(l)%Nx)then
                                mesh(l)%u0(k,j)=mesh(l)%u(k-1,j)
                                mesh(l)%v0(k,j)=mesh(l)%v(k-1,j)
                                mesh(l)%pre(k,j)=mesh(l)%rho(k-1,j)
                            endif
                        enddo
                        !$acc end loop
                    enddo
                    !$acc end loop
                    !$acc loop independent
                    do j=1,mesh(l)%Ny
                        !$acc loop independent
                        do k=2,mesh(l)%Nx
                            mesh(l)%f(k,j,:)=mesh(l)%f0(k,j,:)
                            if(k.ge.2.and.k.le.mesh(l)%Nx)then
                                mesh(l)%u(k,j)=mesh(l)%u0(k,j)
                                mesh(l)%v(k,j)=mesh(l)%v0(k,j)
                                mesh(l)%rho(k,j)=mesh(l)%pre(k,j)
                            endif
                        enddo
                        !$acc end loop
                    enddo
                    !$acc end loop
                    !$acc end kernels
                else
                    if((l.eq.1).and.(periodical_bc.eq.1))then
                        numt=1
                        if(myid.eq.numprocs-1)then
                            do j=2,mesh(l)%Ny-1
                                buf22(numt)=mesh(l)%u(mesh(l)%Nx-1,j)
                                buf22(numt+1)=mesh(l)%v(mesh(l)%Nx-1,j)
                                buf22(numt+2)=mesh(l)%rho(mesh(l)%Nx-1,j)
                                numt=numt+3
                            enddo
                            numt=numt-1
                        endif
                    endif
                    if(myid.eq.0)then
                        do k=mesh(l)%Nx,2,-1
                            mesh(l)%f(k,1,:)=mesh(l)%f(k-1,1,:)
                            mesh(l)%f(k,mesh(l)%Ny,:)=mesh(l)%f(k-1,mesh(l)%Ny,:)
                        enddo
                    endif
                    if(allocated(buf1)) deallocate(buf1)
                    if(allocated(buf2)) deallocate(buf2)
                    allocate(buf1(2*Q),buf2(2*Q))
                    numt=1

                    if(myid.eq.0) then
			            sp=1
			            ep=mesh(l)%Nx
		            else
			            sp=mesh(l)%start_row
			            ep=mesh(l)%end_row
                    endif
                    call transfanduandvandrho(1,dr,mesh(l)%f,mesh(l)%u,mesh(l)%v,mesh(l)%rho,sp,ep,mesh(l)%start_row,mesh(l)%end_row,mesh(l)%Nx,mesh(l)%Ny)
                    if(myid.eq.0)then
                        do j=2,mesh(l)%Ny
                            mesh(l)%f(2,j,:)=mesh(l)%f(1,j,:)
                        enddo
                    endif

                    if(allocated(buf1)) deallocate(buf1)
                    if(allocated(buf2)) deallocate(buf2)
                    allocate(buf1(Q*(mesh(l)%Ny-2)),buf2(Q*(mesh(l)%Ny-2)))                
                    if(myid.eq.0) then
			            sp=1
			            ep=mesh(l)%Nx
		            else
			            sp=mesh(l)%start_row
			            ep=mesh(l)%end_row
                    endif                    
                    call scatter(mesh(l)%f,sp,ep,mesh(l)%Nx,mesh(l)%Ny)
                endif
            endif
        enddo
        if(periodical_bc.eq.1)then
            if(numprocs.eq.1)then
                !$acc kernels present(temp,mesh)
                !$acc loop 
                do j=2,mesh(l)%Ny-1
                    mesh(1)%u(2,j)=temp(j,1)
                    mesh(1)%v(2,j)=temp(j,2)
                    mesh(1)%rho(2,j)=temp(j,3)
                enddo  
                !$acc end loop
                !$acc loop    
                do k=1,Q
                    do j=1,mesh(1)%Ny                 
                        mesh(1)%f(mesh(l)%Nx,j,k)=mesh(1)%f(2,j,k)
                        mesh(1)%f(1,j,k)=mesh(1)%f(mesh(1)%Nx-1,j,k)
                    enddo
                enddo
                !$acc end loop  
                !$acc end kernels
            else
                if(allocated(buf2)) deallocate(buf2)
                allocate(buf2(3*(mesh(1)%Ny-2)))
                numt=1
                if(myid.eq.numprocs-1)then
                    call MPI_SEND(buf22,3*(mesh(1)%Ny-2),MPI_DOUBLE_PRECISION,0,100,MPI_COMM_WORLD,ierr)
                elseif(myid.eq.0)then
                    call MPI_RECV(buf2,3*(mesh(1)%Ny-2),MPI_DOUBLE_PRECISION,numprocs-1,100,MPI_COMM_WORLD,status,ierr)
                    numt=1
                    do j=2,mesh(l)%Ny-1
                        mesh(l)%u(2,j)=buf2(numt)
                        mesh(l)%v(2,j)=buf2(numt+1)
                        mesh(l)%rho(2,j)=buf2(numt+2)
                        numt=numt+3
                    enddo
                endif
                    
                if(allocated(buf1)) deallocate(buf1)
                if(allocated(buf2)) deallocate(buf2)
                allocate(buf1(Q*mesh(l)%Ny),buf2(Q*mesh(l)%Ny))
                numt=1
                if(myid.eq.numprocs-1)then
                    do j=1,mesh(1)%Ny
                        buf1(numt:numt+Q-1)=mesh(1)%f(mesh(1)%Nx-1,j,:)
                        numt=numt+Q
                    enddo
                    numt=numt-1
                    call MPI_SENDRECV(buf1,numt,MPI_DOUBLE_PRECISION,numprocs-1,3,buf2,Q*mesh(1)%Ny,MPI_DOUBLE_PRECISION,numprocs-1,3,MPI_COMM_WORLD,status,ierr)
                    numt=1
                    do j=1,mesh(l)%Ny
                        mesh(1)%f(mesh(1)%Nx,j,:)=buf2(numt:numt+Q-1)
                        numt=numt+Q
                    enddo                        
                elseif(myid.eq.0)then
                    do j=1,mesh(1)%Ny
                        buf1(numt:numt+Q-1)=mesh(1)%f(2,j,:)
                        numt=numt+Q
                    enddo
                    numt=numt-1
                    call MPI_SENDRECV(buf1,numt,MPI_DOUBLE_PRECISION,0,3,buf2,Q*mesh(1)%Ny,MPI_DOUBLE_PRECISION,0,3,MPI_COMM_WORLD,status,ierr)
                    numt=1
                    do j=1,mesh(1)%Ny
                        mesh(1)%f(1,j,:)=buf2(numt:numt+Q-1)
                        numt=numt+Q
                    enddo    
                endif                
            endif
        endif
    endif
    !$acc exit data delete(temp)
end