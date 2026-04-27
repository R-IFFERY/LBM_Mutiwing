Subroutine ini_mesh(n,mesh)	
	use mpi
	use com_date
    use meshtype
	implicit none


	integer n,i,j,k,m,sum
	type(mesh_dat) mesh(n)
	real(kind=ps) spacel(5),spacer(5),spaceu(5),spaced(5)
	real(kind=ps) width,height,qqq,minspace1,minspace2,xo,yo,s1,s2,temp
	Real (kind=ps), external :: stret

	!网格范围
	width=40.0d0
	height=20.0d0

	!网格1初始化
	mesh(1)%meshid=1
	mesh(1)%layer=1
	mesh(1)%scale=mesh_max_scale
	mesh(1)%Nx=nint(width/mesh(1)%scale)+1
	mesh(1)%Ny=nint(height/mesh(1)%scale)+1
	mesh(1)%dx=1.0d0
	mesh(1)%dy=1.0d0
	mesh(1)%Lx=mesh(1)%dx*dble(mesh(1)%Nx)
	mesh(1)%Ly=mesh(1)%dy*dble(mesh(1)%Ny)
	mesh(1)%dt=mesh(1)%dx
	mesh(1)%c=mesh(1)%dx/mesh(1)%dt
	mesh(1)%niu=UU*(1.0d0/mesh(1)%scale)/Re
	mesh(1)%tau_f=3.0d0*mesh(1)%niu+0.5d0
    allocate(mesh(1)%s(Q))
	allocate(mesh(1)%body(mesh(1)%Nx,mesh(1)%Ny))	
	allocate(mesh(1)%left(mesh(1)%Ny,4),mesh(1)%right(mesh(1)%Ny,4),mesh(1)%up(mesh(1)%Nx,4),mesh(1)%down(mesh(1)%Nx,4))	!(:,1)与几号网格邻接；=0，则为数值边界条件，(:,2)是否有对应的临界网格交换点，=0，没有交换点，=1有交换点，(:,3)邻接网格的x编号，(:,4)邻接网格的y编号
	allocate(mesh(1)%leftf2(2:mesh(1)%Ny-1,Q),mesh(1)%rightf2(2:mesh(1)%Ny-1,Q),mesh(1)%upf2(mesh(1)%Nx,Q),mesh(1)%downf2(mesh(1)%Nx,Q))
	mesh(1)%ox=-0.5d0*width
	mesh(1)%oy=-0.5d0*height
	!内存分配
	mesh(1)%start_row=myid*int(mesh(1)%Nx/numprocs)-1
	mesh(1)%end_row=(myid+1)*int(mesh(1)%Nx/numprocs)
	if(myid.eq.0) mesh(1)%start_row=1
	if(myid.ge.numprocs-1) mesh(1)%end_row=mesh(1)%Nx
	if(myid.eq.0) then
		allocate(mesh(1)%rho(mesh(1)%Nx,mesh(1)%Ny),mesh(1)%u(mesh(1)%Nx,mesh(1)%Ny),mesh(1)%v(mesh(1)%Nx,mesh(1)%Ny),mesh(1)%u0(mesh(1)%Nx,mesh(1)%Ny),mesh(1)%v0(mesh(1)%Nx,mesh(1)%Ny),mesh(1)%pre(mesh(1)%Nx,mesh(1)%Ny),mesh(1)%f(mesh(1)%Nx,mesh(1)%Ny,Q),mesh(1)%f0(mesh(1)%Nx,mesh(1)%Ny,Q),mesh(1)%vor(mesh(1)%Nx,mesh(1)%Ny),mesh(1)%vor1(mesh(1)%Nx,mesh(1)%Ny),mesh(1)%vor2(mesh(1)%Nx,mesh(1)%Ny),mesh(1)%dvordt(mesh(1)%Nx,mesh(1)%Ny))
	else
		allocate(mesh(1)%rho(mesh(1)%start_row:mesh(1)%end_row,mesh(1)%Ny),mesh(1)%u(mesh(1)%start_row:mesh(1)%end_row,mesh(1)%Ny),mesh(1)%v(mesh(1)%start_row:mesh(1)%end_row,mesh(1)%Ny),mesh(1)%u0(mesh(1)%start_row:mesh(1)%end_row,mesh(1)%Ny),mesh(1)%v0(mesh(1)%start_row:mesh(1)%end_row,mesh(1)%Ny),mesh(1)%pre(mesh(1)%start_row:mesh(1)%end_row,mesh(1)%Ny),mesh(1)%f(mesh(1)%start_row:mesh(1)%end_row,mesh(1)%Ny,Q),mesh(1)%f0(mesh(1)%start_row:mesh(1)%end_row,mesh(1)%Ny,Q),mesh(1)%vor(mesh(1)%start_row:mesh(1)%end_row,mesh(1)%Ny),mesh(1)%vor1(mesh(1)%start_row:mesh(1)%end_row,mesh(1)%Ny),mesh(1)%vor2(mesh(1)%start_row:mesh(1)%end_row,mesh(1)%Ny),mesh(1)%dvordt(mesh(1)%start_row:mesh(1)%end_row,mesh(1)%Ny))
    endif
	mesh(1)%left=0
	mesh(1)%right=0
	mesh(1)%up=0
	mesh(1)%down=0

!	print*, mesh(7).up(32,:)
	call allocate_nodes_and_faces_numbers(n,mesh)
	if(myid.eq.0)then
		call show_mesh(n,mesh)
    endif

    mesh(max_meshid)%s(1)=0.0d0
    mesh(max_meshid)%s(2)=1.1d0
    mesh(max_meshid)%s(3)=1.05d0
    mesh(max_meshid)%s(4)=0.0d0
    mesh(max_meshid)%s(5)=1.2d0
    mesh(max_meshid)%s(6)=0.0d0
    mesh(max_meshid)%s(7)=1.2d0
    mesh(max_meshid)%s(8)=1.0d0/mesh(max_meshid)%tau_f
    mesh(max_meshid)%s(9)=1.0d0/mesh(max_meshid)%tau_f         
    sum=0
    do i=1,max_meshid
        if(myid.eq.0) write(*,'(a,i3.2,a)') '网格',mesh(i)%meshid,'的信息:'
        if(myid.eq.0) write(*,'(a,f8.4,a,f8.4,a,f8.4,a,f8.4)') '网格无量纲范围，x：',mesh(i)%ox,'->',mesh(i)%ox+mesh(i)%scale*dble(mesh(i)%Nx-2),'，y：',mesh(i)%oy,'->',mesh(i)%oy+mesh(i)%scale*dble(mesh(i)%Ny-2)
        if(myid.eq.0) write(*,'(a,i5.4,a,i5.4,a,i8.7,a )') '节点数：',mesh(i)%Nx,'×',mesh(i)%Ny,'，共',mesh(i)%Nx*mesh(i)%Ny,'个节点'
        if(myid.eq.0) print*,'松弛因子为:'
        if(myid.eq.0) write(*,'(a,f6.4,a,f6.4,a,f6.4)') 's1=',mesh(i)%s(1),'   s2=',mesh(i)%s(2),'   s3=',mesh(i)%s(3)
        if(myid.eq.0) write(*,'(a,f6.4,a,f6.4,a,f6.4)') 's4=',mesh(i)%s(4),'   s5=',mesh(i)%s(5),'   s6=',mesh(i)%s(6)
        if(myid.eq.0) write(*,'(a,f6.4,a,f6.4,a,f6.4)') 's7=',mesh(i)%s(7),'   s8=',mesh(i)%s(8),'   s9=',mesh(i)%s(9)
        sum=sum+mesh(i)%Nx*mesh(i)%Ny   
    enddo
    if(myid.eq.0) print*,'本次计算的总网格节点数为',sum
endsubroutine

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
