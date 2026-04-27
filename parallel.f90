!subroutine space_distribution
!	use mpi
!	use com_date
!	implicit none

!	if (myid.eq.0) then
!		allocate (rho(Nx,Ny),u(Nx,Ny),v(Nx,Ny),u0(Nx,Ny),v0(Nx,Ny),pre(Nx,Ny),f(Nx,Ny,Q),f0(Nx,Ny,Q),feq(Nx,Ny,Q),sa(Nx,Ny,Q))		
!	else
!		allocate (rho(start_row:end_row,Ny),u(start_row:end_row,Ny),v(start_row:end_row,Ny),u0(start_row:end_row,Ny),v0(start_row:end_row,Ny),pre(start_row:end_row,Ny),f(start_row:end_row,Ny,Q),f0(start_row:end_row,Ny,Q),feq(start_row:end_row,Ny,Q),sa(start_row:end_row,Ny,Q))
!	endif
!end subroutine

subroutine gather(tem,sr,er,Nx,Ny,start_row,end_row)
	use mpi
	use com_date
	implicit none

	integer sr,er
	integer Nx,Ny,start_row,end_row
	real (kind=ps) tem(sr:er,Ny)
	integer sendcount,i,j,k1,k2,temsend,i2,i3,ii,i4,start_row1,end_row1
	integer,allocatable :: recvcounts(:),displs(:)
	real(kind=ps),allocatable :: sendbuf(:),recvbuf(:)

	allocate (recvcounts(numprocs),displs(numprocs))
	if (myid.eq.0) then
		k1=1
		k2=end_row-1
		sendcount=(end_row-1)*Ny
	elseif (myid.eq.numprocs-1) then
		k1=start_row+1
		k2=end_row
		sendcount=(end_row-start_row)*Ny
	else
		k1=start_row+1
		k2=end_row-1
		sendcount=((end_row-1)-(start_row+1)+1)*Ny
	endif
	
	allocate(sendbuf(sendcount),recvbuf(Nx*Ny))
	recvcounts=Nx*Ny
	displs(1)=0

	do i=1,numprocs
		if(i.ne.1) displs(i)=displs(i-1)+temsend
		start_row1=(i-1)*int(Nx/numprocs)
		end_row1=i*int(Nx/numprocs)-1
		if(i.eq.1) start_row1=1
		if(i.ge.numprocs) end_row1=Nx
        temsend=(end_row1-start_row1+1)*Ny
		recvcounts(i)=temsend
	enddo

	i2=1
	do j=1,Ny
		do i=k1,k2
			sendbuf(i2)=tem(i,j)
			i2=i2+1
		enddo
	enddo

	call MPI_BARRIER(MPI_COMM_WORLD,ierr)
	call mpi_gatherv(sendbuf,sendcount,mpi_real8,recvbuf,recvcounts,displs,mpi_real8,0,MPI_COMM_WORLD,ierr)
	call MPI_BARRIER(MPI_COMM_WORLD,ierr)

	if (myid.eq.0) then
		do i=1,numprocs
			start_row1=(i-1)*int(Nx/numprocs)
			end_row1=i*int(Nx/numprocs)-1
			if(i.eq.1) start_row1=1
			if(i.ge.numprocs) end_row1=Nx
!			if (i.eq.1) then
!				k1=1
!				k2=end_row1-1
!			elseif (i.eq.numprocs) then
!				k1=start_row1+1
!				k2=end_row1
!			else
!				k1=start_row1+1
!				k2=end_row1-1
!			endif
			ii=displs(i)
            k1=start_row1
            k2=end_row1
			do j=1,Ny
				do i4=k1,k2
					tem(i4,j)=recvbuf(ii+1)
					ii=ii+1
				enddo
			enddo

		enddo
	endif
	deallocate (recvcounts,displs,sendbuf,recvbuf)

end subroutine

subroutine transferdata(f,sp,ep,start_row,end_row,Nx,Ny,id,l1,l2)
	use mpi
	use com_date
	implicit none
	integer sp,ep,start_row,end_row,Nx,Ny,id
	integer spp,epp,l1,l2
	real (kind=ps) f(sp:ep,Ny,Q),cpu_ts1,cpu_es11
	integer cou,numt,i,j,k,l,status(MPI_STATUS_SIZE)
	real (kind=ps), allocatable :: buf1(:),buf2(:)
	integer, allocatable :: recvcounts(:),displs(:)
!	buf1(3*Q*(Ny-2)),buf2(3*Q*(Ny-2))

    if(myid.eq.0) call CPU_TIME(cpu_ts1)
	allocate (buf1(Q*(Ny-2)),buf2(Q*(Ny-2)))

	!从左向右平移数据	
	cou=0
	numt=Q*(Ny-2)	

    do k=1,Q
        do j=2,Ny-1
            buf1(cou+1)=f(end_row-1,j,k)
            cou=cou+1
        enddo
    enddo
	call MPI_BARRIER(MPI_COMM_WORLD,ierr) 
	call MPI_SENDRECV(buf1,numt,MPI_DOUBLE_PRECISION,rightpro,1,buf2,numt,MPI_DOUBLE_PRECISION,leftpro,1,MPI_COMM_WORLD,status,ierr)
	if(myid.ne.0)then
		cou=0
        do k=1,Q
            do j=2,Ny-1
                f(start_row,j,k)=buf2(cou+1)
                cou=cou+1
            enddo
        enddo
	endif
	!从右向左平移数据
	cou=0

    do k=1,Q
        do j=2,Ny-1
            buf1(cou+1)=f(start_row+1,j,k)
            cou=cou+1
        enddo
    enddo
	call MPI_BARRIER(MPI_COMM_WORLD,ierr)
	call MPI_SENDRECV(buf1,numt,MPI_DOUBLE_PRECISION,leftpro,2,buf2,numt,MPI_DOUBLE_PRECISION,rightpro,2,MPI_COMM_WORLD,status,ierr)
	if(myid.ne.numprocs-1)then
		cou=0
        do k=1,Q
            do j=2,Ny-1
                f(end_row,j,k)=buf2(cou+1)
                cou=cou+1
            enddo
        enddo
	endif
	call MPI_BARRIER(MPI_COMM_WORLD,ierr)
	deallocate (buf1,buf2)
    
    allocate (buf1(2*(Ny-2)*Q))

	if(myid.eq.(numprocs-1)) then
		cou=1
		do k=1,Q
			do j=2,Ny-1
				buf1(cou)=f(ep-1,j,k)
                buf1(cou+1)=f(ep-2,j,k)
				cou=cou+2
			enddo
		enddo
		cou=cou-2	
		call MPI_SEND(buf1,cou,MPI_DOUBLE_PRECISION,0,100,MPI_COMM_WORLD,ierr)
	elseif(myid.eq.0) then
		call MPI_RECV(buf1,2*(Ny-2)*Q,MPI_DOUBLE_PRECISION,numprocs-1,100,MPI_COMM_WORLD,status,ierr)
		cou=1
		do k=1,Q
			do j=2,Ny-1
				f(ep-1,j,k)=buf1(cou)
                f(ep-2,j,k)=buf1(cou+1)
				cou=cou+2
			enddo
		enddo
    endif
    deallocate (buf1)
	call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    
    spp=myid*int(Nx/numprocs)
    epp=(myid+1)*int(Nx/numprocs)-1
    if(myid.eq.0) spp=2
    if(myid.eq.(numprocs-1)) epp=Nx-1
    allocate (buf1(4*(Nx-2)*Q),buf2(4*(epp-spp+1)*Q))
    allocate (recvcounts(numprocs),displs(numprocs))
    cou=1
    numt=4*(epp-spp+1)*Q
    if(myid.eq.0) then
        displs(1)=0
        do i=0,numprocs-1
            recvcounts(i+1)=(i+1)*int(Nx/numprocs)-i*int(Nx/numprocs)
            if (i.eq.0) recvcounts(i+1)=int(Nx/numprocs)-2
            if (i.eq.numprocs-1) recvcounts(i+1)= Nx-i*int(Nx/numprocs)
            recvcounts(i+1)=4*recvcounts(i+1)*Q
            if(i.ne.0) displs(i+1)=displs(i)+recvcounts(i)
         enddo
    endif
    
    call MPI_BCAST(recvcounts,numprocs,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(displs,numprocs,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    
    do k=1,Q
        do i=spp,epp
            buf2(cou)=f(i,2,k)
            buf2(cou+1)=f(i,3,k)
            buf2(cou+2)=f(i,Ny-1,k)
            buf2(cou+3)=f(i,Ny-2,k)
            cou=cou+4
        enddo
    enddo
    call MPI_GATHERV(buf2,numt,MPI_DOUBLE_PRECISION,buf1,recvcounts,displs,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    cou=0
    
    if(myid.eq.0) then
        do l=0,numprocs-1
            spp=l*int(Nx/numprocs)
            epp=(l+1)*int(Nx/numprocs)-1
            if(l.eq.0) spp=2
            if(l.eq.(numprocs-1)) epp=Nx-1
            do k=1,Q
                do i=spp,epp
                    f(i,2,k)=buf1(cou+1)
                    f(i,3,k)=buf1(cou+2)
                    f(i,Ny-1,k)=buf1(cou+3)
                    f(i,Ny-2,k)=buf1(cou+4)
                    cou=cou+4
                enddo
            enddo
        enddo
    endif
    
    deallocate (buf1,buf2,recvcounts,displs)
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    
    if((mod(id,2).eq.1).and.(id.ne.max_meshid))then
        spp=myid*int(Nx/numprocs)
	    epp=(myid+1)*int(Nx/numprocs)-1
	    if(myid.eq.0) spp=2
	    if(myid.eq.(numprocs-1)) epp=Nx-1
	    allocate (buf1((Nx-2)*Q),buf2((epp-spp+1)*Q))
	    allocate (recvcounts(numprocs),displs(numprocs))
	    cou=1
	    numt=(epp-spp+1)*Q
    !	print*, numt,spp,epp,Ny,myid

	    if(myid.eq.0) then
		    displs(1)=0
		    do i=0,numprocs-1
			    recvcounts(i+1)=(i+1)*int(Nx/numprocs)-i*int(Nx/numprocs)
			    if (i.eq.0) recvcounts(i+1)=int(Nx/numprocs)-2
			    if (i.eq.numprocs-1) recvcounts(i+1)= Nx-i*int(Nx/numprocs)
			    recvcounts(i+1)=recvcounts(i+1)*Q
			    if(i.ne.0) displs(i+1)=displs(i)+recvcounts(i)
		    enddo
	    endif

	    call MPI_BCAST(recvcounts,numprocs,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	    call MPI_BCAST(displs,numprocs,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  
	    do k=1,Q
            do i=spp,epp
                buf2(cou)=f(i,l1,k)
                cou=cou+1
            enddo
        enddo
	    call MPI_GATHERV(buf2,numt,MPI_DOUBLE_PRECISION,buf1,recvcounts,displs,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
	    cou=0
    
	    if(myid.eq.0) then
		    do l=0,numprocs-1
			    spp=l*int(Nx/numprocs)
			    epp=(l+1)*int(Nx/numprocs)-1
			    if(l.eq.0) spp=2
			    if(l.eq.(numprocs-1)) epp=Nx-1
			    do k=1,Q
                     do i=spp,epp
                         f(i,l1,k)=buf1(cou+1)
                        cou=cou+1
                    enddo
			    enddo
		    enddo
        endif
        
        
        cou=1
        spp=myid*int(Nx/numprocs)
	    epp=(myid+1)*int(Nx/numprocs)-1
	    if(myid.eq.0) spp=2
	    if(myid.eq.(numprocs-1)) epp=Nx-1 
        
        do k=1,Q
            do i=spp,epp
                buf2(cou)=f(i,l2,k)
                cou=cou+1
            enddo
        enddo
        
	    call MPI_GATHERV(buf2,numt,MPI_DOUBLE_PRECISION,buf1,recvcounts,displs,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        
	    cou=0
        
	    if(myid.eq.0) then
		    do l=0,numprocs-1
			    spp=l*int(Nx/numprocs)
			    epp=(l+1)*int(Nx/numprocs)-1
			    if(l.eq.0) spp=2
			    if(l.eq.(numprocs-1)) epp=Nx-1
			    do k=1,Q
                     do i=spp,epp
                         f(i,l2,k)=buf1(cou+1)
                        cou=cou+1
                    enddo
			    enddo
		    enddo
        endif 
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        deallocate (buf1,buf2,recvcounts,displs)
    endif

    end subroutine
    
subroutine transferdatapv(f,rho11,u11,v11,p11,sp,ep,start_row,end_row,Nx,Ny,id,l1,l2)
    use mpi
	use com_date
	implicit none
    integer sp,ep,Nx,Ny,start_row,end_row,spp,epp,id,l1,l2
    real (kind=ps) f(sp:ep,Ny,Q),rho11(sp:ep,Ny),u11(sp:ep,Ny),v11(sp:ep,Ny),p11(sp:ep,Ny)
    integer cou,numt,i,j,k,l,status(MPI_STATUS_SIZE)
	real (kind=ps), allocatable :: buf1(:),buf2(:)
	integer, allocatable :: recvcounts(:),displs(:)
    
!    start_row=sp
!    end_row=ep
        
    allocate (buf1(Q*(Ny-2)),buf2(Q*(Ny-2)))
    
    
    !!传递分布函数
    !从左向右平移数据	
	cou=0
	numt=Q*(Ny-2)	

    do k=1,Q
        do j=2,Ny-1
            buf1(cou+1)=f(end_row-1,j,k)
            cou=cou+1
        enddo
    enddo
	call MPI_BARRIER(MPI_COMM_WORLD,ierr) 
	call MPI_SENDRECV(buf1,numt,MPI_DOUBLE_PRECISION,rightpro,1,buf2,numt,MPI_DOUBLE_PRECISION,leftpro,1,MPI_COMM_WORLD,status,ierr)
	if(myid.ne.0)then
		cou=0
        do k=1,Q
            do j=2,Ny-1
                f(start_row,j,k)=buf2(cou+1)
                cou=cou+1
            enddo
        enddo
    endif
    
	!从右向左平移数据
	cou=0

    do k=1,Q
        do j=2,Ny-1
            buf1(cou+1)=f(start_row+1,j,k)
            cou=cou+1
        enddo
    enddo
	call MPI_BARRIER(MPI_COMM_WORLD,ierr)
	call MPI_SENDRECV(buf1,numt,MPI_DOUBLE_PRECISION,leftpro,2,buf2,numt,MPI_DOUBLE_PRECISION,rightpro,2,MPI_COMM_WORLD,status,ierr)
	if(myid.ne.numprocs-1)then
		cou=0
        do k=1,Q
            do j=2,Ny-1
                f(end_row,j,k)=buf2(cou+1)
                cou=cou+1
            enddo
        enddo
	endif
	call MPI_BARRIER(MPI_COMM_WORLD,ierr)
	deallocate (buf1,buf2)
    
    
    !收集最右侧三列的数据
    allocate (buf1(2*(Ny-2)*Q))

	if(myid.eq.(numprocs-1)) then
		cou=1
		do k=1,Q
			do j=2,Ny-1
				buf1(cou)=f(ep-1,j,k)
                buf1(cou+1)=f(ep-2,j,k)
				cou=cou+2
			enddo
		enddo
		cou=cou-2	
		call MPI_SEND(buf1,cou,MPI_DOUBLE_PRECISION,0,100,MPI_COMM_WORLD,ierr)
	elseif(myid.eq.0) then
		call MPI_RECV(buf1,2*(Ny-2)*Q,MPI_DOUBLE_PRECISION,numprocs-1,100,MPI_COMM_WORLD,status,ierr)
		cou=1
		do k=1,Q
			do j=2,Ny-1
				f(ep-1,j,k)=buf1(cou)
                f(ep-2,j,k)=buf1(cou+1)
				cou=cou+2
			enddo
		enddo
    endif
    deallocate (buf1)
	call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    
    !收集上下最外侧三行的数据
    spp=myid*int(Nx/numprocs)
    epp=(myid+1)*int(Nx/numprocs)-1
    if(myid.eq.0) spp=2
    if(myid.eq.(numprocs-1)) epp=Nx-1
    allocate (buf1(4*(Nx-2)*Q),buf2(4*(epp-spp+1)*Q))
    allocate (recvcounts(numprocs),displs(numprocs))
    cou=1
    numt=4*(epp-spp+1)*Q
    if(myid.eq.0) then
        displs(1)=0
        do i=0,numprocs-1
            recvcounts(i+1)=(i+1)*int(Nx/numprocs)-i*int(Nx/numprocs)
            if (i.eq.0) recvcounts(i+1)=int(Nx/numprocs)-2
            if (i.eq.numprocs-1) recvcounts(i+1)= Nx-i*int(Nx/numprocs)
            recvcounts(i+1)=4*recvcounts(i+1)*Q
            if(i.ne.0) displs(i+1)=displs(i)+recvcounts(i)
         enddo
    endif
    
    call MPI_BCAST(recvcounts,numprocs,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(displs,numprocs,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    
    do k=1,Q
        do i=spp,epp
            buf2(cou)=f(i,2,k)
            buf2(cou+1)=f(i,3,k)
            buf2(cou+2)=f(i,Ny-1,k)
            buf2(cou+3)=f(i,Ny-2,k)
            cou=cou+4
        enddo
    enddo
    call MPI_GATHERV(buf2,numt,MPI_DOUBLE_PRECISION,buf1,recvcounts,displs,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    cou=0
    
    if(myid.eq.0) then
        do l=0,numprocs-1
            spp=l*int(Nx/numprocs)
            epp=(l+1)*int(Nx/numprocs)-1
            if(l.eq.0) spp=2
            if(l.eq.(numprocs-1)) epp=Nx-1
            do k=1,Q
                do i=spp,epp
                    f(i,2,k)=buf1(cou+1)
                    f(i,3,k)=buf1(cou+2)
                    f(i,Ny-1,k)=buf1(cou+3)
                    f(i,Ny-2,k)=buf1(cou+4)
                    cou=cou+4
                enddo
            enddo
        enddo
    endif
    
    deallocate (buf1,buf2,recvcounts,displs)
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    
    if((mod(id,2).eq.1).and.(id.ne.max_meshid))then
        spp=myid*int(Nx/numprocs)
	    epp=(myid+1)*int(Nx/numprocs)-1
	    if(myid.eq.0) spp=2
	    if(myid.eq.(numprocs-1)) epp=Nx-1
	    allocate (buf1((Nx-2)*Q),buf2((epp-spp+1)*Q))
	    allocate (recvcounts(numprocs),displs(numprocs))
	    cou=1
	    numt=(epp-spp+1)*Q
    !	print*, numt,spp,epp,Ny,myid

	    if(myid.eq.0) then
		    displs(1)=0
		    do i=0,numprocs-1
			    recvcounts(i+1)=(i+1)*int(Nx/numprocs)-i*int(Nx/numprocs)
			    if (i.eq.0) recvcounts(i+1)=int(Nx/numprocs)-2
			    if (i.eq.numprocs-1) recvcounts(i+1)= Nx-i*int(Nx/numprocs)
			    recvcounts(i+1)=recvcounts(i+1)*Q
			    if(i.ne.0) displs(i+1)=displs(i)+recvcounts(i)
		    enddo
	    endif

	    call MPI_BCAST(recvcounts,numprocs,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	    call MPI_BCAST(displs,numprocs,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  
	    do k=1,Q
            do i=spp,epp
                buf2(cou)=f(i,l1,k)
                cou=cou+1
            enddo
        enddo
	    call MPI_GATHERV(buf2,numt,MPI_DOUBLE_PRECISION,buf1,recvcounts,displs,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
	    cou=0
    
	    if(myid.eq.0) then
		    do l=0,numprocs-1
			    spp=l*int(Nx/numprocs)
			    epp=(l+1)*int(Nx/numprocs)-1
			    if(l.eq.0) spp=2
			    if(l.eq.(numprocs-1)) epp=Nx-1
			    do k=1,Q
                     do i=spp,epp
                         f(i,l1,k)=buf1(cou+1)
                        cou=cou+1
                    enddo
			    enddo
		    enddo
        endif
        
        
        cou=1
        spp=myid*int(Nx/numprocs)
	    epp=(myid+1)*int(Nx/numprocs)-1
	    if(myid.eq.0) spp=2
	    if(myid.eq.(numprocs-1)) epp=Nx-1 
        
        do k=1,Q
            do i=spp,epp
                buf2(cou)=f(i,l2,k)
                cou=cou+1
            enddo
        enddo
        
	    call MPI_GATHERV(buf2,numt,MPI_DOUBLE_PRECISION,buf1,recvcounts,displs,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        
	    cou=0
        
	    if(myid.eq.0) then
		    do l=0,numprocs-1
			    spp=l*int(Nx/numprocs)
			    epp=(l+1)*int(Nx/numprocs)-1
			    if(l.eq.0) spp=2
			    if(l.eq.(numprocs-1)) epp=Nx-1
			    do k=1,Q
                     do i=spp,epp
                         f(i,l2,k)=buf1(cou+1)
                        cou=cou+1
                    enddo
			    enddo
		    enddo
        endif 
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        deallocate (buf1,buf2,recvcounts,displs)
    endif
    
    !!传递u速度
    !从左向右平移数据
    allocate (buf1(Ny),buf2(Ny))
	cou=0
	numt=Ny	

    do j=1,Ny
        buf1(cou+1)=u11(end_row-1,j)
        cou=cou+1
    enddo
	call MPI_BARRIER(MPI_COMM_WORLD,ierr) 
	call MPI_SENDRECV(buf1,numt,MPI_DOUBLE_PRECISION,rightpro,1,buf2,numt,MPI_DOUBLE_PRECISION,leftpro,1,MPI_COMM_WORLD,status,ierr)
	if(myid.ne.0)then
		cou=0
        do j=1,Ny
            u11(start_row,j)=buf2(cou+1)
            cou=cou+1
        enddo
	endif
	!从右向左平移数据
	cou=0

    do j=1,Ny
        buf1(cou+1)=u11(start_row+1,j)
        cou=cou+1
    enddo
	call MPI_BARRIER(MPI_COMM_WORLD,ierr)
	call MPI_SENDRECV(buf1,numt,MPI_DOUBLE_PRECISION,leftpro,2,buf2,numt,MPI_DOUBLE_PRECISION,rightpro,2,MPI_COMM_WORLD,status,ierr)
	if(myid.ne.numprocs-1)then
		cou=0
        do j=1,Ny
            u11(end_row,j)=buf2(cou+1)
            cou=cou+1
        enddo
    endif
    deallocate (buf1,buf2)
	call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    
    !收集最右侧三列的数据
    allocate (buf1(3*(Ny-2)))

	if(myid.eq.(numprocs-1)) then
		cou=1
		do j=2,Ny-1
            buf1(cou)=u11(ep,j)
            buf1(cou+1)=u11(ep-1,j)
            buf1(cou+2)=u11(ep-2,j)
            cou=cou+3
		enddo
		cou=cou-3	
		call MPI_SEND(buf1,cou,MPI_DOUBLE_PRECISION,0,100,MPI_COMM_WORLD,ierr)
	elseif(myid.eq.0) then
		call MPI_RECV(buf1,3*(Ny-2),MPI_DOUBLE_PRECISION,numprocs-1,100,MPI_COMM_WORLD,status,ierr)
		cou=1
		do j=2,Ny-1
            u11(ep,j)=buf1(cou)
            u11(ep-1,j)=buf1(cou+1)
            u11(ep-2,j)=buf1(cou+2)
            cou=cou+3
		enddo
    endif
    deallocate (buf1)
	call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    
    !收集上下最外侧三行的数据
    spp=myid*int(Nx/numprocs)
    epp=(myid+1)*int(Nx/numprocs)-1
    if(myid.eq.0) spp=1
    if(myid.eq.(numprocs-1)) epp=Nx
    allocate (buf1(6*Nx),buf2(6*(epp-spp+1)))
    allocate (recvcounts(numprocs),displs(numprocs))
    cou=1
    numt=6*(epp-spp+1)
    if(myid.eq.0) then
        displs(1)=0
        do i=0,numprocs-1
            recvcounts(i+1)=(i+1)*int(Nx/numprocs)-i*int(Nx/numprocs)
            if (i.eq.0) recvcounts(i+1)=int(Nx/numprocs)-1
            if (i.eq.numprocs-1) recvcounts(i+1)= Nx-i*int(Nx/numprocs)+1
            recvcounts(i+1)=6*recvcounts(i+1)
            if(i.ne.0) displs(i+1)=displs(i)+recvcounts(i)
         enddo
    endif
    
    call MPI_BCAST(recvcounts,numprocs,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(displs,numprocs,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    
    do i=spp,epp
        buf2(cou)=u11(i,1)
        buf2(cou+1)=u11(i,2)
        buf2(cou+2)=u11(i,3)
        buf2(cou+3)=u11(i,Ny)
        buf2(cou+4)=u11(i,Ny-1)
        buf2(cou+5)=u11(i,Ny-2)
        cou=cou+6
    enddo
    call MPI_GATHERV(buf2,numt,MPI_DOUBLE_PRECISION,buf1,recvcounts,displs,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    cou=0
    
    if(myid.eq.0) then
        do l=0,numprocs-1
            spp=l*int(Nx/numprocs)
            epp=(l+1)*int(Nx/numprocs)-1
            if(l.eq.0) spp=1
            if(l.eq.(numprocs-1)) epp=Nx
            do i=spp,epp
                u11(i,1)=buf1(cou+1)
                u11(i,2)=buf1(cou+2)
                u11(i,3)=buf1(cou+3)
                u11(i,Ny)=buf1(cou+4)
                u11(i,Ny-1)=buf1(cou+5)
                u11(i,Ny-2)=buf1(cou+6)
                cou=cou+6
            enddo
        enddo
    endif
    
    deallocate (buf1,buf2,recvcounts,displs)
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    
    if((mod(id,2).eq.1).and.(id.ne.max_meshid))then
        spp=myid*int(Nx/numprocs)
	    epp=(myid+1)*int(Nx/numprocs)-1
	    if(myid.eq.0) spp=1
	    if(myid.eq.(numprocs-1)) epp=Nx
	    allocate (buf1(Nx),buf2((epp-spp+1)))
	    allocate (recvcounts(numprocs),displs(numprocs))
	    cou=1
	    numt=(epp-spp+1)
    !	print*, numt,spp,epp,Ny,myid

	    if(myid.eq.0) then
		    displs(1)=0
		    do i=0,numprocs-1
			    recvcounts(i+1)=(i+1)*int(Nx/numprocs)-i*int(Nx/numprocs)
			    if (i.eq.0) recvcounts(i+1)=int(Nx/numprocs)-1
			    if (i.eq.numprocs-1) recvcounts(i+1)= Nx-i*int(Nx/numprocs)+1
!			    recvcounts(i+1)=recvcounts(i+1)*Q
			    if(i.ne.0) displs(i+1)=displs(i)+recvcounts(i)
		    enddo
	    endif

	    call MPI_BCAST(recvcounts,numprocs,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	    call MPI_BCAST(displs,numprocs,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  
        do i=spp,epp
            buf2(cou)=u11(i,l1)
            cou=cou+1
        enddo
	    call MPI_GATHERV(buf2,numt,MPI_DOUBLE_PRECISION,buf1,recvcounts,displs,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
	    cou=0
    
	    if(myid.eq.0) then
		    do l=0,numprocs-1
			    spp=l*int(Nx/numprocs)
			    epp=(l+1)*int(Nx/numprocs)-1
			    if(l.eq.0) spp=1
			    if(l.eq.(numprocs-1)) epp=Nx
                do i=spp,epp
                    u11(i,l1)=buf1(cou+1)
                    cou=cou+1
                enddo
		    enddo
        endif
        
        
        cou=1
        spp=myid*int(Nx/numprocs)
	    epp=(myid+1)*int(Nx/numprocs)-1
	    if(myid.eq.0) spp=1
	    if(myid.eq.(numprocs-1)) epp=Nx
        
        do i=spp,epp
            buf2(cou)=u11(i,l2)
            cou=cou+1
        enddo
        
	    call MPI_GATHERV(buf2,numt,MPI_DOUBLE_PRECISION,buf1,recvcounts,displs,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        
	    cou=0
        
	    if(myid.eq.0) then
		    do l=0,numprocs-1
			    spp=l*int(Nx/numprocs)
			    epp=(l+1)*int(Nx/numprocs)-1
			    if(l.eq.0) spp=1
			    if(l.eq.(numprocs-1)) epp=Nx
                do i=spp,epp
                    u11(i,l2)=buf1(cou+1)
                    cou=cou+1
                enddo
		    enddo
        endif 
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        deallocate (buf1,buf2,recvcounts,displs)
    endif
    
    !!传递v速度
    !从左向右平移数据
    allocate (buf1(Ny),buf2(Ny))
	cou=0
	numt=Ny	

    do j=1,Ny
        buf1(cou+1)=v11(end_row-1,j)
        cou=cou+1
    enddo
	call MPI_BARRIER(MPI_COMM_WORLD,ierr) 
	call MPI_SENDRECV(buf1,numt,MPI_DOUBLE_PRECISION,rightpro,1,buf2,numt,MPI_DOUBLE_PRECISION,leftpro,1,MPI_COMM_WORLD,status,ierr)
	if(myid.ne.0)then
		cou=0
        do j=1,Ny
            v11(start_row,j)=buf2(cou+1)
            cou=cou+1
        enddo
	endif
	!从右向左平移数据
	cou=0

    do j=1,Ny
        buf1(cou+1)=v11(start_row+1,j)
        cou=cou+1
    enddo
	call MPI_BARRIER(MPI_COMM_WORLD,ierr)
	call MPI_SENDRECV(buf1,numt,MPI_DOUBLE_PRECISION,leftpro,2,buf2,numt,MPI_DOUBLE_PRECISION,rightpro,2,MPI_COMM_WORLD,status,ierr)
	if(myid.ne.numprocs-1)then
		cou=0
        do j=1,Ny
            v11(end_row,j)=buf2(cou+1)
            cou=cou+1
        enddo
    endif
    deallocate (buf1,buf2)
	call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    
    !收集最右侧三列的数据
    allocate (buf1(3*(Ny-2)))

	if(myid.eq.(numprocs-1)) then
		cou=1
		do j=2,Ny-1
            buf1(cou)=v11(ep,j)
            buf1(cou+1)=v11(ep-1,j)
            buf1(cou+2)=v11(ep-2,j)
            cou=cou+3
		enddo
		cou=cou-3	
		call MPI_SEND(buf1,cou,MPI_DOUBLE_PRECISION,0,100,MPI_COMM_WORLD,ierr)
	elseif(myid.eq.0) then
		call MPI_RECV(buf1,3*(Ny-2),MPI_DOUBLE_PRECISION,numprocs-1,100,MPI_COMM_WORLD,status,ierr)
		cou=1
		do j=2,Ny-1
            v11(ep,j)=buf1(cou)
            v11(ep-1,j)=buf1(cou+1)
            v11(ep-2,j)=buf1(cou+2)
            cou=cou+3
		enddo
    endif
    deallocate (buf1)
	call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    
    !收集上下最外侧三行的数据
    spp=myid*int(Nx/numprocs)
    epp=(myid+1)*int(Nx/numprocs)-1
    if(myid.eq.0) spp=1
    if(myid.eq.(numprocs-1)) epp=Nx
    allocate (buf1(6*Nx),buf2(6*(epp-spp+1)))
    allocate (recvcounts(numprocs),displs(numprocs))
    cou=1
    numt=6*(epp-spp+1)
    if(myid.eq.0) then
        displs(1)=0
        do i=0,numprocs-1
            recvcounts(i+1)=(i+1)*int(Nx/numprocs)-i*int(Nx/numprocs)
            if (i.eq.0) recvcounts(i+1)=int(Nx/numprocs)-1
            if (i.eq.numprocs-1) recvcounts(i+1)= Nx-i*int(Nx/numprocs)+1
            recvcounts(i+1)=6*recvcounts(i+1)
            if(i.ne.0) displs(i+1)=displs(i)+recvcounts(i)
         enddo
    endif
    
    call MPI_BCAST(recvcounts,numprocs,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(displs,numprocs,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    
    do i=spp,epp
        buf2(cou)=v11(i,1)
        buf2(cou+1)=v11(i,2)
        buf2(cou+2)=v11(i,3)
        buf2(cou+3)=v11(i,Ny)
        buf2(cou+4)=v11(i,Ny-1)
        buf2(cou+5)=v11(i,Ny-2)
        cou=cou+6
    enddo
    call MPI_GATHERV(buf2,numt,MPI_DOUBLE_PRECISION,buf1,recvcounts,displs,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    cou=0
    
    if(myid.eq.0) then
        do l=0,numprocs-1
            spp=l*int(Nx/numprocs)
            epp=(l+1)*int(Nx/numprocs)-1
            if(l.eq.0) spp=1
            if(l.eq.(numprocs-1)) epp=Nx
            do i=spp,epp
                v11(i,1)=buf1(cou+1)
                v11(i,2)=buf1(cou+2)
                v11(i,3)=buf1(cou+3)
                v11(i,Ny)=buf1(cou+4)
                v11(i,Ny-1)=buf1(cou+5)
                v11(i,Ny-2)=buf1(cou+6)
                cou=cou+6
            enddo
        enddo
    endif
    
    deallocate (buf1,buf2,recvcounts,displs)
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    
    if((mod(id,2).eq.1).and.(id.ne.max_meshid))then
        spp=myid*int(Nx/numprocs)
	    epp=(myid+1)*int(Nx/numprocs)-1
	    if(myid.eq.0) spp=1
	    if(myid.eq.(numprocs-1)) epp=Nx
	    allocate (buf1(Nx),buf2((epp-spp+1)))
	    allocate (recvcounts(numprocs),displs(numprocs))
	    cou=1
	    numt=(epp-spp+1)
    !	print*, numt,spp,epp,Ny,myid

	    if(myid.eq.0) then
		    displs(1)=0
		    do i=0,numprocs-1
			    recvcounts(i+1)=(i+1)*int(Nx/numprocs)-i*int(Nx/numprocs)
			    if (i.eq.0) recvcounts(i+1)=int(Nx/numprocs)-1
			    if (i.eq.numprocs-1) recvcounts(i+1)= Nx-i*int(Nx/numprocs)+1
!			    recvcounts(i+1)=recvcounts(i+1)*Q
			    if(i.ne.0) displs(i+1)=displs(i)+recvcounts(i)
		    enddo
	    endif

	    call MPI_BCAST(recvcounts,numprocs,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	    call MPI_BCAST(displs,numprocs,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  
        do i=spp,epp
            buf2(cou)=v11(i,l1)
            cou=cou+1
        enddo
	    call MPI_GATHERV(buf2,numt,MPI_DOUBLE_PRECISION,buf1,recvcounts,displs,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
	    cou=0
    
	    if(myid.eq.0) then
		    do l=0,numprocs-1
			    spp=l*int(Nx/numprocs)
			    epp=(l+1)*int(Nx/numprocs)-1
			    if(l.eq.0) spp=1
			    if(l.eq.(numprocs-1)) epp=Nx
                do i=spp,epp
                    v11(i,l1)=buf1(cou+1)
                    cou=cou+1
                enddo
		    enddo
        endif
        
        
        cou=1
        spp=myid*int(Nx/numprocs)
	    epp=(myid+1)*int(Nx/numprocs)-1
	    if(myid.eq.0) spp=1
	    if(myid.eq.(numprocs-1)) epp=Nx
        
        do i=spp,epp
            buf2(cou)=v11(i,l2)
            cou=cou+1
        enddo
        
	    call MPI_GATHERV(buf2,numt,MPI_DOUBLE_PRECISION,buf1,recvcounts,displs,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        
	    cou=0
        
	    if(myid.eq.0) then
		    do l=0,numprocs-1
			    spp=l*int(Nx/numprocs)
			    epp=(l+1)*int(Nx/numprocs)-1
			    if(l.eq.0) spp=1
			    if(l.eq.(numprocs-1)) epp=Nx
                do i=spp,epp
                    v11(i,l2)=buf1(cou+1)
                    cou=cou+1
                enddo
		    enddo
        endif 
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        deallocate (buf1,buf2,recvcounts,displs)
    endif
    
    !!传递密度
    !从左向右平移数据
    allocate (buf1(Ny),buf2(Ny))
	cou=0
	numt=Ny	

    do j=1,Ny
        buf1(cou+1)=rho11(end_row-1,j)
        cou=cou+1
    enddo
	call MPI_BARRIER(MPI_COMM_WORLD,ierr) 
	call MPI_SENDRECV(buf1,numt,MPI_DOUBLE_PRECISION,rightpro,1,buf2,numt,MPI_DOUBLE_PRECISION,leftpro,1,MPI_COMM_WORLD,status,ierr)
	if(myid.ne.0)then
		cou=0
        do j=1,Ny
            rho11(start_row,j)=buf2(cou+1)
            cou=cou+1
        enddo
	endif
	!从右向左平移数据
	cou=0

    do j=1,Ny
        buf1(cou+1)=rho11(start_row+1,j)
        cou=cou+1
    enddo
	call MPI_BARRIER(MPI_COMM_WORLD,ierr)
	call MPI_SENDRECV(buf1,numt,MPI_DOUBLE_PRECISION,leftpro,2,buf2,numt,MPI_DOUBLE_PRECISION,rightpro,2,MPI_COMM_WORLD,status,ierr)
	if(myid.ne.numprocs-1)then
		cou=0
        do j=1,Ny
            rho11(end_row,j)=buf2(cou+1)
            cou=cou+1
        enddo
    endif
    deallocate (buf1,buf2)
	call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    
    !收集最右侧三列的数据
    allocate (buf1(3*(Ny-2)))

	if(myid.eq.(numprocs-1)) then
		cou=1
		do j=2,Ny-1
            buf1(cou)=rho11(ep,j)
            buf1(cou+1)=rho11(ep-1,j)
            buf1(cou+2)=rho11(ep-2,j)
            cou=cou+3
		enddo
		cou=cou-3	
		call MPI_SEND(buf1,cou,MPI_DOUBLE_PRECISION,0,100,MPI_COMM_WORLD,ierr)
	elseif(myid.eq.0) then
		call MPI_RECV(buf1,3*(Ny-2),MPI_DOUBLE_PRECISION,numprocs-1,100,MPI_COMM_WORLD,status,ierr)
		cou=1
		do j=2,Ny-1
            rho11(ep,j)=buf1(cou)
            rho11(ep-1,j)=buf1(cou+1)
            rho11(ep-2,j)=buf1(cou+2)
            cou=cou+3
		enddo
    endif
    deallocate (buf1)
	call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    
    !收集上下最外侧三行的数据
    spp=myid*int(Nx/numprocs)
    epp=(myid+1)*int(Nx/numprocs)-1
    if(myid.eq.0) spp=1
    if(myid.eq.(numprocs-1)) epp=Nx
    allocate (buf1(6*Nx),buf2(6*(epp-spp+1)))
    allocate (recvcounts(numprocs),displs(numprocs))
    cou=1
    numt=6*(epp-spp+1)
    if(myid.eq.0) then
        displs(1)=0
        do i=0,numprocs-1
            recvcounts(i+1)=(i+1)*int(Nx/numprocs)-i*int(Nx/numprocs)
            if (i.eq.0) recvcounts(i+1)=int(Nx/numprocs)-1
            if (i.eq.numprocs-1) recvcounts(i+1)= Nx-i*int(Nx/numprocs)+1
            recvcounts(i+1)=6*recvcounts(i+1)
            if(i.ne.0) displs(i+1)=displs(i)+recvcounts(i)
         enddo
    endif
    
    call MPI_BCAST(recvcounts,numprocs,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(displs,numprocs,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    
    do i=spp,epp
        buf2(cou)=rho11(i,1)
        buf2(cou+1)=rho11(i,2)
        buf2(cou+2)=rho11(i,3)
        buf2(cou+3)=rho11(i,Ny)
        buf2(cou+4)=rho11(i,Ny-1)
        buf2(cou+5)=rho11(i,Ny-2)
        cou=cou+6
    enddo
    call MPI_GATHERV(buf2,numt,MPI_DOUBLE_PRECISION,buf1,recvcounts,displs,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    cou=0
    
    if(myid.eq.0) then
        do l=0,numprocs-1
            spp=l*int(Nx/numprocs)
            epp=(l+1)*int(Nx/numprocs)-1
            if(l.eq.0) spp=1
            if(l.eq.(numprocs-1)) epp=Nx
            do i=spp,epp
                rho11(i,1)=buf1(cou+1)
                rho11(i,2)=buf1(cou+2)
                rho11(i,3)=buf1(cou+3)
                rho11(i,Ny)=buf1(cou+4)
                rho11(i,Ny-1)=buf1(cou+5)
                rho11(i,Ny-2)=buf1(cou+6)
                cou=cou+6
            enddo
        enddo
    endif
    
    deallocate (buf1,buf2,recvcounts,displs)
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    
    if((mod(id,2).eq.1).and.(id.ne.max_meshid))then
        spp=myid*int(Nx/numprocs)
	    epp=(myid+1)*int(Nx/numprocs)-1
	    if(myid.eq.0) spp=1
	    if(myid.eq.(numprocs-1)) epp=Nx
	    allocate (buf1(Nx),buf2((epp-spp+1)))
	    allocate (recvcounts(numprocs),displs(numprocs))
	    cou=1
	    numt=(epp-spp+1)
    !	print*, numt,spp,epp,Ny,myid

	    if(myid.eq.0) then
		    displs(1)=0
		    do i=0,numprocs-1
			    recvcounts(i+1)=(i+1)*int(Nx/numprocs)-i*int(Nx/numprocs)
			    if (i.eq.0) recvcounts(i+1)=int(Nx/numprocs)-1
			    if (i.eq.numprocs-1) recvcounts(i+1)= Nx-i*int(Nx/numprocs)+1
!			    recvcounts(i+1)=recvcounts(i+1)*Q
			    if(i.ne.0) displs(i+1)=displs(i)+recvcounts(i)
		    enddo
	    endif

	    call MPI_BCAST(recvcounts,numprocs,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	    call MPI_BCAST(displs,numprocs,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  
        do i=spp,epp
            buf2(cou)=rho11(i,l1)
            cou=cou+1
        enddo
	    call MPI_GATHERV(buf2,numt,MPI_DOUBLE_PRECISION,buf1,recvcounts,displs,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
	    cou=0
    
	    if(myid.eq.0) then
		    do l=0,numprocs-1
			    spp=l*int(Nx/numprocs)
			    epp=(l+1)*int(Nx/numprocs)-1
			    if(l.eq.0) spp=1
			    if(l.eq.(numprocs-1)) epp=Nx
                do i=spp,epp
                    rho11(i,l1)=buf1(cou+1)
                    cou=cou+1
                enddo
		    enddo
        endif
        
        
        cou=1
        spp=myid*int(Nx/numprocs)
	    epp=(myid+1)*int(Nx/numprocs)-1
	    if(myid.eq.0) spp=1
	    if(myid.eq.(numprocs-1)) epp=Nx
        
        do i=spp,epp
            buf2(cou)=rho11(i,l2)
            cou=cou+1
        enddo
        
	    call MPI_GATHERV(buf2,numt,MPI_DOUBLE_PRECISION,buf1,recvcounts,displs,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        
	    cou=0
        
	    if(myid.eq.0) then
		    do l=0,numprocs-1
			    spp=l*int(Nx/numprocs)
			    epp=(l+1)*int(Nx/numprocs)-1
			    if(l.eq.0) spp=1
			    if(l.eq.(numprocs-1)) epp=Nx
                do i=spp,epp
                    rho11(i,l2)=buf1(cou+1)
                    cou=cou+1
                enddo
		    enddo
        endif 
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        deallocate (buf1,buf2,recvcounts,displs)
    endif
    
    !!传递压力
    !从左向右平移数据
    allocate (buf1(Ny),buf2(Ny))
	cou=0
	numt=Ny	

    do j=1,Ny
        buf1(cou+1)=p11(end_row-1,j)
        cou=cou+1
    enddo
	call MPI_BARRIER(MPI_COMM_WORLD,ierr) 
	call MPI_SENDRECV(buf1,numt,MPI_DOUBLE_PRECISION,rightpro,1,buf2,numt,MPI_DOUBLE_PRECISION,leftpro,1,MPI_COMM_WORLD,status,ierr)
	if(myid.ne.0)then
		cou=0
        do j=1,Ny
            p11(start_row,j)=buf2(cou+1)
            cou=cou+1
        enddo
	endif
	!从右向左平移数据
	cou=0

    do j=1,Ny
        buf1(cou+1)=p11(start_row+1,j)
        cou=cou+1
    enddo
	call MPI_BARRIER(MPI_COMM_WORLD,ierr)
	call MPI_SENDRECV(buf1,numt,MPI_DOUBLE_PRECISION,leftpro,2,buf2,numt,MPI_DOUBLE_PRECISION,rightpro,2,MPI_COMM_WORLD,status,ierr)
	if(myid.ne.numprocs-1)then
		cou=0
        do j=1,Ny
            p11(end_row,j)=buf2(cou+1)
            cou=cou+1
        enddo
	endif
	call MPI_BARRIER(MPI_COMM_WORLD,ierr)
	deallocate (buf1,buf2)
    
    !收集最右侧三列的数据
    allocate (buf1(3*(Ny-2)))

	if(myid.eq.(numprocs-1)) then
		cou=1
		do j=2,Ny-1
            buf1(cou)=p11(ep,j)
            buf1(cou+1)=p11(ep-1,j)
            buf1(cou+2)=p11(ep-2,j)
            cou=cou+3
		enddo
		cou=cou-3	
		call MPI_SEND(buf1,cou,MPI_DOUBLE_PRECISION,0,100,MPI_COMM_WORLD,ierr)
	elseif(myid.eq.0) then
		call MPI_RECV(buf1,3*(Ny-2),MPI_DOUBLE_PRECISION,numprocs-1,100,MPI_COMM_WORLD,status,ierr)
		cou=1
		do j=2,Ny-1
            p11(ep,j)=buf1(cou)
            p11(ep-1,j)=buf1(cou+1)
            p11(ep-2,j)=buf1(cou+2)
            cou=cou+3
		enddo
    endif
    deallocate (buf1)
	call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    
    !收集上下最外侧三行的数据
    spp=myid*int(Nx/numprocs)
    epp=(myid+1)*int(Nx/numprocs)-1
    if(myid.eq.0) spp=1
    if(myid.eq.(numprocs-1)) epp=Nx
    allocate (buf1(6*Nx),buf2(6*(epp-spp+1)))
    allocate (recvcounts(numprocs),displs(numprocs))
    cou=1
    numt=6*(epp-spp+1)
    if(myid.eq.0) then
        displs(1)=0
        do i=0,numprocs-1
            recvcounts(i+1)=(i+1)*int(Nx/numprocs)-i*int(Nx/numprocs)
            if (i.eq.0) recvcounts(i+1)=int(Nx/numprocs)-1
            if (i.eq.numprocs-1) recvcounts(i+1)= Nx-i*int(Nx/numprocs)+1
            recvcounts(i+1)=6*recvcounts(i+1)
            if(i.ne.0) displs(i+1)=displs(i)+recvcounts(i)
         enddo
    endif
    
    call MPI_BCAST(recvcounts,numprocs,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(displs,numprocs,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    
    do i=spp,epp
        buf2(cou)=p11(i,1)
        buf2(cou+1)=p11(i,2)
        buf2(cou+2)=p11(i,3)
        buf2(cou+3)=p11(i,Ny)
        buf2(cou+4)=p11(i,Ny-1)
        buf2(cou+5)=p11(i,Ny-2)
        cou=cou+6
    enddo
    call MPI_GATHERV(buf2,numt,MPI_DOUBLE_PRECISION,buf1,recvcounts,displs,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    cou=0
    
    if(myid.eq.0) then
        do l=0,numprocs-1
            spp=l*int(Nx/numprocs)
            epp=(l+1)*int(Nx/numprocs)-1
            if(l.eq.0) spp=1
            if(l.eq.(numprocs-1)) epp=Nx
            do i=spp,epp
                p11(i,1)=buf1(cou+1)
                p11(i,2)=buf1(cou+2)
                p11(i,3)=buf1(cou+3)
                p11(i,Ny)=buf1(cou+4)
                p11(i,Ny-1)=buf1(cou+5)
                p11(i,Ny-2)=buf1(cou+6)
                cou=cou+6
            enddo
        enddo
    endif
    
    deallocate (buf1,buf2,recvcounts,displs)
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    
    if((mod(id,2).eq.1).and.(id.ne.max_meshid))then
        spp=myid*int(Nx/numprocs)
	    epp=(myid+1)*int(Nx/numprocs)-1
	    if(myid.eq.0) spp=1
	    if(myid.eq.(numprocs-1)) epp=Nx
	    allocate (buf1(Nx),buf2((epp-spp+1)))
	    allocate (recvcounts(numprocs),displs(numprocs))
	    cou=1
	    numt=(epp-spp+1)
    !	print*, numt,spp,epp,Ny,myid

	    if(myid.eq.0) then
		    displs(1)=0
		    do i=0,numprocs-1
			    recvcounts(i+1)=(i+1)*int(Nx/numprocs)-i*int(Nx/numprocs)
			    if (i.eq.0) recvcounts(i+1)=int(Nx/numprocs)-1
			    if (i.eq.numprocs-1) recvcounts(i+1)= Nx-i*int(Nx/numprocs)+1
!			    recvcounts(i+1)=recvcounts(i+1)*Q
			    if(i.ne.0) displs(i+1)=displs(i)+recvcounts(i)
		    enddo
	    endif

	    call MPI_BCAST(recvcounts,numprocs,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	    call MPI_BCAST(displs,numprocs,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  
        do i=spp,epp
            buf2(cou)=p11(i,l1)
            cou=cou+1
        enddo
	    call MPI_GATHERV(buf2,numt,MPI_DOUBLE_PRECISION,buf1,recvcounts,displs,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
	    cou=0
    
	    if(myid.eq.0) then
		    do l=0,numprocs-1
			    spp=l*int(Nx/numprocs)
			    epp=(l+1)*int(Nx/numprocs)-1
			    if(l.eq.0) spp=1
			    if(l.eq.(numprocs-1)) epp=Nx
                do i=spp,epp
                    p11(i,l1)=buf1(cou+1)
                    cou=cou+1
                enddo
		    enddo
        endif
        
        
        cou=1
        spp=myid*int(Nx/numprocs)
	    epp=(myid+1)*int(Nx/numprocs)-1
	    if(myid.eq.0) spp=1
	    if(myid.eq.(numprocs-1)) epp=Nx
        
        do i=spp,epp
            buf2(cou)=p11(i,l2)
            cou=cou+1
        enddo
        
	    call MPI_GATHERV(buf2,numt,MPI_DOUBLE_PRECISION,buf1,recvcounts,displs,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        
	    cou=0
        
	    if(myid.eq.0) then
		    do l=0,numprocs-1
			    spp=l*int(Nx/numprocs)
			    epp=(l+1)*int(Nx/numprocs)-1
			    if(l.eq.0) spp=1
			    if(l.eq.(numprocs-1)) epp=Nx
                do i=spp,epp
                    p11(i,l2)=buf1(cou+1)
                    cou=cou+1
                enddo
		    enddo
        endif 
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        deallocate (buf1,buf2,recvcounts,displs)
    endif
    
end subroutine

subroutine scatter(f,sp,ep,Nx,Ny)
	use mpi
	use com_date
	implicit none
	integer sp,ep,Nx,Ny,n,status(MPI_STATUS_SIZE)
	real(kind=ps) f(sp:ep,Ny,Q)
	integer i,j,k
	integer spp,epp,numt,cou
	integer sendcounts(numprocs),displs(numprocs)
	real(kind=ps), allocatable :: buf1(:),buf2(:)

    allocate (buf1((Ny-2)*Q))

	if(myid.eq.0) then
		cou=1
		do k=1,Q
			do j=2,Ny-1
				buf1(cou)=f(ep,j,k)
				cou=cou+1
			enddo
		enddo
		cou=cou-1	
		call MPI_SEND(buf1,cou,MPI_DOUBLE_PRECISION,numprocs-1,100,MPI_COMM_WORLD,ierr)
	elseif(myid.eq.(numprocs-1)) then
		call MPI_RECV(buf1,(Ny-2)*Q,MPI_DOUBLE_PRECISION,0,100,MPI_COMM_WORLD,status,ierr)
		cou=1
		do k=1,Q
			do j=2,Ny-1
				f(ep,j,k)=buf1(cou)
				cou=cou+1
			enddo
		enddo
	endif
	call MPI_BARRIER(MPI_COMM_WORLD,ierr)

	do i=0,numprocs-1
		spp=i*int(Nx/numprocs)-1
		epp=(i+1)*int(Nx/numprocs)
		if(i.eq.0) spp=1
		if(i.ge.numprocs-1) epp=Nx
		if(i.eq.0) then
			sendcounts(i+1)=0
			displs(i+1)=0
		else
			displs(i+1)= displs(i)+sendcounts(i)
			sendcounts(i+1)= 2*(epp-spp+1)*Q
		endif
	enddo

	numt=0
	do i=1,numprocs
		numt=numt+sendcounts(i)
	enddo
	
	cou=1
	deallocate (buf1)
	allocate (buf1(numt),buf2(sendcounts(myid+1)))
	if(myid.eq.0) then
		do j=0,numprocs-1
			spp=j*int(Nx/numprocs)-1
			epp=(j+1)*int(Nx/numprocs)
			if(j.eq.0) spp=1
			if(j.ge.numprocs-1) epp=Nx
			if (j.ne.0) then
				do k=1,Q
					do	i=spp,epp
						buf1(cou)=f(i,1,k)
						buf1(cou+1)=f(i,Ny,k)
						cou=cou+2
					enddo
				enddo
			endif
		enddo
	endif

	call MPI_SCATTERV(buf1,sendcounts,displs,MPI_DOUBLE_PRECISION,buf2,sendcounts(myid+1),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

	if(myid.ne.0) then
		cou=1
		do k=1,Q
			do i=sp,ep
				f(i,1,k)=buf2(cou)
				f(i,Ny,k)=buf2(cou+1)
				cou=cou+2
			enddo
		enddo
    endif
end

subroutine garthervf (f,sp,ep,start_row,end_row,Nx,Ny)

	use mpi
	use com_date
	implicit none
	integer sp,ep,start_row,end_row,Nx,Ny
	integer spp,epp
	real (kind=ps) f(sp:ep,Ny,Q)
	integer cou,numt,i,j,k,l,status(MPI_STATUS_SIZE)
	real (kind=ps), allocatable :: buf1(:),buf2(:)
	integer, allocatable :: recvcounts(:),displs(:)

	!主进程收集各子进程的数据
	spp=myid*Nx/numprocs+1
	epp=(myid+1)*Nx/numprocs
	if(myid.eq.0) spp=1
	if(myid.eq.(numprocs-1)) epp=Nx
	allocate (buf1(Nx*Ny*Q),buf2((epp-spp+1)*Ny*Q))
	allocate (recvcounts(numprocs),displs(numprocs))
	cou=1
	numt=(epp-spp+1)*Ny*Q
!	print*, numt,spp,epp,Ny,myid
	if(myid.eq.0) then
		displs(1)=0
		do i=0,numprocs-1
			recvcounts(i+1)=(i+1)*Nx/numprocs-i*Nx/numprocs
			if (i.eq.0) recvcounts(i+1)=Nx/numprocs
			if (i.eq.numprocs-1) recvcounts(i+1)= Nx-i*Nx/numprocs
			recvcounts(i+1)=recvcounts(i+1)*Ny*Q
			if(i.ne.0) displs(i+1)=displs(i)+recvcounts(i)
		enddo
	endif

	call MPI_BCAST(recvcounts,numprocs,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(displs,numprocs,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!	if(myid.eq.0) print*, recvcounts
	do k=1,Q
		do j=1,Ny
			do i=spp,epp
				buf2(cou)=f(i,j,k)
				cou=cou+1
			enddo
		enddo
	enddo
	call MPI_GATHERV(buf2,numt,MPI_DOUBLE_PRECISION,buf1,recvcounts,displs,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
	cou=0
	if(myid.eq.0) then
		do l=0,numprocs-1
			spp=l*Nx/numprocs+1
			epp=(l+1)*Nx/numprocs
			if(l.eq.0) spp=1
			if(l.eq.(numprocs-1)) epp=Nx
			do k=1,Q
				do j=1,Ny
					do i=spp,epp
						f(i,j,k)=buf1(cou+1)
						cou=cou+1
					enddo
				enddo
			enddo
		enddo
	endif
	deallocate (buf1,buf2)
end
    
subroutine transfanduandvandrho(move_num,dr,f,u11,v11,rho11,sp,ep,start_row,end_row,Nx,Ny)
	use mpi
	use com_date
	implicit none
    
    integer move_num,dr,start_row,end_row,Nx,Ny,sp,ep,i,j,k,num,numt,status(MPI_STATUS_SIZE)
    real(kind=ps) f(sp:ep,Ny,Q),u11(sp:ep,Ny),v11(sp:ep,Ny),rho11(sp:ep,Ny)
    real(kind=ps), allocatable :: buf1(:),buf2(:)
    
    if(dr.eq.1)then
        allocate (buf1((3+Q)*(Ny-2)*(move_num+1)),buf2((3+Q)*(Ny-2)*(move_num+1)))
        num=1

        do j=2,Ny-1
            do i=start_row+1,start_row+move_num+1
                buf1(num)=u11(i,j)
                num=num+1
            enddo
        enddo
        do j=2,Ny-1
            do i=start_row+1,start_row+move_num+1
                buf1(num)=v11(i,j)
                num=num+1
            enddo
        enddo
        do j=2,Ny-1
            do i=start_row+1,start_row+move_num+1
                buf1(num)=rho11(i,j)
                num=num+1
            enddo
        enddo
        do k=1,Q
            do j=2,Ny-1
                do i=start_row+1,start_row+move_num+1
                    buf1(num)=f(i,j,k)
                    num=num+1
                enddo
            enddo
        enddo
        num=num-1           

        call MPI_SENDRECV(buf1,num,MPI_DOUBLE_PRECISION,leftpro,2,buf2,num,MPI_DOUBLE_PRECISION,rightpro,2,MPI_COMM_WORLD,status,ierr)
        
        do j=2,Ny-1
            do i=start_row,end_row-move_num-1
                u11(i,j)=u11(i+move_num,j)
                v11(i,j)=v11(i+move_num,j)
                rho11(i,j)=rho11(i+move_num,j)
                f(i,j,:)=f(i+move_num,j,:)            
            enddo
        enddo
        num=1
        if(myid.ne.numprocs-1)then
            do j=2,Ny-1
                do i=end_row-move_num,end_row
                    u11(i,j)=buf2(num)
                    num=num+1
                enddo
            enddo
            do j=2,Ny-1
                do i=end_row-move_num,end_row
                    v11(i,j)=buf2(num)
                    num=num+1
                enddo
            enddo
            do j=2,Ny-1
                do i=end_row-move_num,end_row
                    rho11(i,j)=buf2(num)
                    num=num+1
                enddo
            enddo
            do k=1,Q
                do j=2,Ny-1
                    do i=end_row-move_num,end_row
                        f(i,j,k)=buf2(num)
                        num=num+1
                    enddo
                enddo
            enddo
            num=num-1           
        endif
    elseif(dr.eq.-1) then
        allocate (buf1((3+Q)*(Ny-2)*(move_num+1)),buf2((3+Q)*(Ny-2)*(move_num+1)))
        num=1
        do j=2,Ny-1
            do i=end_row-move_num-1,end_row-1
                buf1(num)=u11(i,j)
                num=num+1
            enddo
        enddo
        do j=2,Ny-1
            do i=end_row-move_num-1,end_row-1
                buf1(num)=v11(i,j)
                num=num+1
            enddo
        enddo
        do j=2,Ny-1
            do i=end_row-move_num-1,end_row-1
                buf1(num)=rho11(i,j)
                num=num+1
            enddo
        enddo
        do k=1,Q
            do j=2,Ny-1
                do i=end_row-move_num-1,end_row-1
                    buf1(num)=f(i,j,k)
                    num=num+1
                enddo
            enddo
        enddo
        num=num-1           
        call MPI_SENDRECV(buf1,num,MPI_DOUBLE_PRECISION,rightpro,3,buf2,num,MPI_DOUBLE_PRECISION,leftpro,3,MPI_COMM_WORLD,status,ierr)
        
        do j=2,Ny-1
            do i=end_row,start_row+move_num+1,-1
                u11(i,j)=u11(i-move_num,j)
                v11(i,j)=v11(i-move_num,j)
                rho11(i,j)=rho11(i-move_num,j)
                f(i,j,:)=f(i-move_num,j,:)            
            enddo
        enddo

        num=1
        if(myid.ne.0)then
            do j=2,Ny-1
                do i=start_row,start_row+move_num
                    u11(i,j)=buf2(num)
                    num=num+1
                enddo
            enddo
            do j=2,Ny-1
                do i=start_row,start_row+move_num
                    v11(i,j)=buf2(num)
                    num=num+1
                enddo
            enddo
            do j=2,Ny-1
                do i=start_row,start_row+move_num
                    rho11(i,j)=buf2(num)
                    num=num+1
                enddo
            enddo
            do k=1,Q
                do j=2,Ny-1
                    do i=start_row,start_row+move_num
                        f(i,j,k)=buf2(num)
                        num=num+1
                    enddo
                enddo
            enddo
            num=num-1           
        endif
!        if(sp.eq.127.and.ep.eq.162.and.myid.eq.numprocs-1) print*,f(162,4,:)
        if(allocated(buf1)) deallocate(buf1)
        if(allocated(buf2)) deallocate(buf2)
        allocate(buf1(Q*(Ny-2)),buf2(Q*(Ny-2)))
        numt=1
        if(myid.eq.numprocs-1)then
            do j=2,Ny-1
                buf1(numt:numt+Q-1)=f(Nx,j,:)
                numt=numt+Q
            enddo
            numt=numt-1
            call MPI_SEND(buf1,numt,MPI_DOUBLE_PRECISION,0,110,MPI_COMM_WORLD,ierr)
        elseif(myid.eq.0)then
            call MPI_RECV(buf2,Q*(Ny-2),MPI_DOUBLE_PRECISION,numprocs-1,110,MPI_COMM_WORLD,status,ierr)
            do j=2,Ny-1
                f(Nx,j,:)=buf2(numt:numt+Q-1)
                numt=numt+Q
            enddo
        endif
    endif
end
    
subroutine scattervariables2d(var,sp,ep,Ynum)
	use mpi
	use com_date
	implicit none
    
    integer sp,ep,Ynum
    real(kind=ps) var(sp:ep,Ynum)
    integer sendcounts(numprocs),displs(numprocs),spp,epp,i,j,k,numt,buf2num
    real(kind=ps), allocatable :: buf1(:),buf2(:)
    
    if(myid.eq.0)then
        epp=int(ep/numprocs)
        allocate(buf1((ep+2*(numprocs-1))*Ynum),buf2(epp*Ynum))
    else
        allocate(buf1((ep+2*(numprocs-1))*Ynum),buf2((ep-sp+1)*Ynum))
    endif
    numt=1
    if(myid.eq.0)then
        do i=0,numprocs-1
		    spp=i*int(ep/numprocs)-1
		    epp=(i+1)*int(ep/numprocs)
		    if(i.eq.0) spp=1
		    if(i.ge.numprocs-1) epp=ep
		    if(i.eq.0) then
			    sendcounts(i+1)=0
			    displs(i+1)=0
		    else
			    displs(i+1)= displs(i)+sendcounts(i)
			    sendcounts(i+1)= Ynum*(epp-spp+1)
            endif
            if(i.ne.0)then
                do j=1,Ynum
                    do k=spp,epp
                        buf1(numt)=var(k,j)
                        numt=numt+1
                    enddo
                enddo
            endif
        enddo

        buf2num=0
    else
        buf2num=Ynum*(ep-sp+1)
    endif
    
    call MPI_SCATTERV(buf1,sendcounts,displs,MPI_DOUBLE_PRECISION,buf2,buf2num,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    
    if(myid.ne.0)then
        numt=1
        do j=1,Ynum
            do k=sp,ep
                var(k,j)=buf2(numt)
                numt=numt+1
            enddo
        enddo       
    endif

end
    
subroutine scattervariables3d(var,sp,ep,Ynum)
	use mpi
	use com_date
	implicit none
    
    integer sp,ep,Ynum
    real(kind=ps) var(sp:ep,Ynum,Q)
    integer sendcounts(numprocs),displs(numprocs),spp,epp,i,j,k,l,numt,buf2num
    real(kind=ps), allocatable :: buf1(:),buf2(:)
    
    if(myid.eq.0)then
        epp=int(ep/numprocs)
        allocate(buf1((ep+2*(numprocs-1))*Ynum*Q),buf2(epp*Ynum*Q))
    else
        allocate(buf1((ep+2*(numprocs-1))*Ynum*Q),buf2((ep-sp+1)*Ynum*Q))
    endif
    numt=1
    if(myid.eq.0)then
        do i=0,numprocs-1
		    spp=i*int(ep/numprocs)-1
		    epp=(i+1)*int(ep/numprocs)
		    if(i.eq.0) spp=1
		    if(i.ge.numprocs-1) epp=ep
		    if(i.eq.0) then
			    sendcounts(i+1)=0
			    displs(i+1)=0
		    else
			    displs(i+1)= displs(i)+sendcounts(i)
			    sendcounts(i+1)= Q*Ynum*(epp-spp+1)
            endif
            if(i.ne.0)then
                do l=1,Q
                    do j=1,Ynum
                        do k=spp,epp
                            buf1(numt)=var(k,j,l)
                            numt=numt+1
                        enddo
                    enddo
                enddo
            endif
        enddo

        buf2num=0
    else
        buf2num=Q*Ynum*(ep-sp+1)
    endif
    
    call MPI_SCATTERV(buf1,sendcounts,displs,MPI_DOUBLE_PRECISION,buf2,buf2num,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    
    if(myid.ne.0)then
        numt=1
        do l=1,Q
            do j=1,Ynum
                do k=sp,ep
                    var(k,j,l)=buf2(numt)
                    numt=numt+1
                enddo
            enddo 
        enddo
    endif

end