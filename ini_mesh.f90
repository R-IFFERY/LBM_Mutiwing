subroutine ini_mesh(n,mesh)
    use mpi
    use com_date
    use meshtype
    implicit none
    integer :: n,k,total
    real(kind=ps) :: width,height
    type(mesh_dat) :: mesh(n)

    if (n.ne.max_meshid.or.numprocs.gt.((nint(40d0/mesh_max_scale)+1)/4)) then
        if(myid.eq.0) print*, 'Invalid mesh count or too many MPI ranks'
        call MPI_ABORT(MPI_COMM_WORLD,910,ierr)
    endif
    total=0
    do k=1,n
        width=grid_xmax(k)-grid_xmin(k)
        height=grid_ymax(k)-grid_ymin(k)
	!网格1初始化
	mesh(k)%meshid=k
	mesh(k)%layer=k
	mesh(k)%scale=mesh_max_scale/2**(k-1)
	mesh(k)%Nx=nint(width/mesh(k)%scale)+1
	mesh(k)%Ny=nint(height/mesh(k)%scale)+1
	mesh(k)%dx=1.0d0
	mesh(k)%dy=1.0d0
	mesh(k)%Lx=mesh(k)%dx*dble(mesh(k)%Nx)
	mesh(k)%Ly=mesh(k)%dy*dble(mesh(k)%Ny)
	mesh(k)%dt=mesh(k)%dx
	mesh(k)%c=mesh(k)%dx/mesh(k)%dt
	mesh(k)%niu=UU*(1.0d0/mesh(k)%scale)/Re
	mesh(k)%tau_f=3.0d0*mesh(k)%niu+0.5d0
    allocate(mesh(k)%s(Q))
	allocate(mesh(k)%body(mesh(k)%Nx,mesh(k)%Ny))
	allocate(mesh(k)%left(mesh(k)%Ny,4),mesh(k)%right(mesh(k)%Ny,4),mesh(k)%up(mesh(k)%Nx,4),mesh(k)%down(mesh(k)%Nx,4))	!(:,1)与几号网格邻接；=0，则为数值边界条件，(:,2)是否有对应的临界网格交换点，=0，没有交换点，=1有交换点，(:,3)邻接网格的x编号，(:,4)邻接网格的y编号
	allocate(mesh(k)%leftf2(2:mesh(k)%Ny-1,Q),mesh(k)%rightf2(2:mesh(k)%Ny-1,Q),mesh(k)%upf2(mesh(k)%Nx,Q),mesh(k)%downf2(mesh(k)%Nx,Q))
	mesh(k)%ox=grid_xmin(k)
	mesh(k)%oy=grid_ymin(k)
	!内存分配
	mesh(k)%start_row=myid*int(mesh(k)%Nx/numprocs)-1
	mesh(k)%end_row=(myid+1)*int(mesh(k)%Nx/numprocs)
	if(myid.eq.0) mesh(k)%start_row=1
	if(myid.ge.numprocs-1) mesh(k)%end_row=mesh(k)%Nx
	if(myid.eq.0) then
		allocate(mesh(k)%rho(mesh(k)%Nx,mesh(k)%Ny),mesh(k)%u(mesh(k)%Nx,mesh(k)%Ny),mesh(k)%v(mesh(k)%Nx,mesh(k)%Ny),mesh(k)%u0(mesh(k)%Nx,mesh(k)%Ny),mesh(k)%v0(mesh(k)%Nx,mesh(k)%Ny),mesh(k)%pre(mesh(k)%Nx,mesh(k)%Ny),mesh(k)%f(mesh(k)%Nx,mesh(k)%Ny,Q),mesh(k)%f0(mesh(k)%Nx,mesh(k)%Ny,Q),mesh(k)%vor(mesh(k)%Nx,mesh(k)%Ny),mesh(k)%vor1(mesh(k)%Nx,mesh(k)%Ny),mesh(k)%vor2(mesh(k)%Nx,mesh(k)%Ny),mesh(k)%dvordt(mesh(k)%Nx,mesh(k)%Ny))
	else
		allocate(mesh(k)%rho(mesh(k)%start_row:mesh(k)%end_row,mesh(k)%Ny),mesh(k)%u(mesh(k)%start_row:mesh(k)%end_row,mesh(k)%Ny),mesh(k)%v(mesh(k)%start_row:mesh(k)%end_row,mesh(k)%Ny),mesh(k)%u0(mesh(k)%start_row:mesh(k)%end_row,mesh(k)%Ny),mesh(k)%v0(mesh(k)%start_row:mesh(k)%end_row,mesh(k)%Ny),mesh(k)%pre(mesh(k)%start_row:mesh(k)%end_row,mesh(k)%Ny),mesh(k)%f(mesh(k)%start_row:mesh(k)%end_row,mesh(k)%Ny,Q),mesh(k)%f0(mesh(k)%start_row:mesh(k)%end_row,mesh(k)%Ny,Q),mesh(k)%vor(mesh(k)%start_row:mesh(k)%end_row,mesh(k)%Ny),mesh(k)%vor1(mesh(k)%start_row:mesh(k)%end_row,mesh(k)%Ny),mesh(k)%vor2(mesh(k)%start_row:mesh(k)%end_row,mesh(k)%Ny),mesh(k)%dvordt(mesh(k)%start_row:mesh(k)%end_row,mesh(k)%Ny))
    endif
    allocate(mesh(k)%f_old(lbound(mesh(k)%f,1):ubound(mesh(k)%f,1),mesh(k)%Ny,Q))
	mesh(k)%left=0
	mesh(k)%right=0
	mesh(k)%up=0
	mesh(k)%down=0

        mesh(k)%s=(/0d0,1.1d0,1.05d0,0d0,1.2d0,0d0,1.2d0, &
                    1d0/mesh(k)%tau_f,1d0/mesh(k)%tau_f/)
        mesh(k)%body=0
        mesh(k)%vor=0d0
        mesh(k)%vor1=0d0
        mesh(k)%vor2=0d0
        mesh(k)%dvordt=0d0
        mesh(k)%nodescounts=mesh(k)%Nx*mesh(k)%Ny
        mesh(k)%facecounts=(mesh(k)%Nx-1)*(mesh(k)%Ny-1)
        mesh(k)%start_num=total
        total=total+mesh(k)%nodescounts
        if (myid.eq.0) then
            print*, 'Level',k,' Nx,Ny=',mesh(k)%Nx,mesh(k)%Ny
            print*, ' x=',grid_xmin(k),grid_xmax(k),' y=',grid_ymin(k),grid_ymax(k)
            print*, ' spacing=',mesh(k)%scale,' physical dt=',UU*mesh(k)%scale,' tau=',mesh(k)%tau_f
        endif
    enddo
    if(myid.eq.0) then
        print*, 'Total stored mesh nodes:',total
        call show_mesh(n,mesh)
    endif
end subroutine ini_mesh

real*8 function stret(m,smin,smax,ss)
	integer m,ii
	real*8 smin,smax,ss(m),q0,a,b,d,eps,x,temp
	real*8 F

	F(x)=x**(m-1)-d*x+d-1
	if(m==2)then
		ss(1)=0.0
		ss(2)=smax
		stret=1.5
		return
	elseif (m==1) then
		ss(1)=smax
		stret=1.0
		return
	elseif(m==3) then
		ss(1)=0.0
		ss(2)=smin
		ss(3)=smax
		if(((smax-smin)/smin)<1) then
			stret=1.5
			return
		else
			stret=(smax-smin)/smin
			return
		endif
	else
		eps=0.000001
		ss(1)=0.0
		a=1.0d0
		b=((smax-smin)/(m-2.0))/smin
		d=smax/smin
		!在[a，b]范围内用二分法求解F(x)=0的解，解即等比数列公比q
		if((F(b)<0.0d0).or.(b<1.0d0))then
			print*,'没有存在满足条件的公比值，请重新输入！'
			stret=0.0
			return
		else
			do
				temp=F((a+b)/2)
				if(abs(temp)<eps) then
					stret=(a+b)/2
					do ii=2,(m-1)
						ss(ii)=ss(ii-1)+smin*(((a+b)/2)**(ii-2))
					enddo
					ss(m)=smax
					return
				elseif(temp>0) then
					b=(a+b)/2
				else
					a=(a+b)/2
				endif

				if(abs(b-a)<eps) then
					stret=(a+b)/2
					do ii=2,(m-1)
						ss(ii)=ss(ii-1)+smin*(((a+b)/2)**(ii-2))
					enddo
					ss(m)=smax
					return
				endif
			enddo
		endif
	endif
end function stret

subroutine show_mesh(n,mesh)
	use mpi
	use com_date
    use meshtype
	implicit none


	integer n,m,i,j
	type(mesh_dat) mesh(n)

	open(290,file='mesh.dat',status='unknown')

	do m=1,n
		write(290,*) 'zone,i=',mesh(m)%Nx,'j=',mesh(m)%Ny
		do J=1,mesh(m)%Ny
			do I=1,mesh(m)%Nx
				WRITE(290,*) mesh(m)%ox+dble(i-1)*mesh(m)%scale,mesh(m)%oy+dble(j-1)*mesh(m)%scale
            enddo
		enddo
	enddo

	close(290)

endsubroutine

subroutine allocate_nodes_and_faces_numbers(n,mesh)

	use mpi
	use com_date
    use meshtype
	implicit none


	integer n,m,i,j,is,ie,js,je
	type(mesh_dat) mesh(n)

	mesh(1)%start_num=0
	mesh(1)%nodescounts=(mesh(1)%Nx-2)*(mesh(1)%Ny-2)
	mesh(1)%facecounts=(mesh(1)%Nx-3)*(mesh(1)%Ny-3)+mesh(1)%Ny-3

	do i=2,n
		mesh(i)%start_num=mesh(i-1)%start_num+mesh(i-1)%nodescounts
		is=2
		ie=mesh(i)%Nx-1
		js=2
		je=mesh(i)%Ny-1
		if(mesh(i)%layer.ne.1) then
			if (mesh(mesh(i)%left(2,1))%layer.lt.mesh(i)%layer) is=3
			if (mesh(mesh(i)%right(2,1))%layer.lt.mesh(i)%layer) ie=mesh(i)%Nx-2
			if (mesh(mesh(i)%down(2,1))%layer.lt.mesh(i)%layer) js=3
			if (mesh(mesh(i)%up(2,1))%layer.lt.mesh(i)%layer) je=mesh(i)%Ny-2
		endif
		mesh(i)%nodescounts=(ie-is+1)*(je-js+1)
		mesh(i)%facecounts=(ie-is)*(je-js)

		if(mesh(i)%layer.ne.4) then
			if((mod(mesh(i)%meshid,4).eq.1).or.(mod(mesh(i)%meshid,4).eq.3)) mesh(i)%facecounts=mesh(i)%facecounts+je-js
			if((mod(mesh(i)%meshid,4).eq.2).or.(mod(mesh(i)%meshid,4).eq.0)) mesh(i)%facecounts=mesh(i)%facecounts+ie-is
		endif
	enddo

endsubroutine
