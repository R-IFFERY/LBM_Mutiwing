!圆柱
subroutine body_point(mesh,Nx1,Ny1)
	use com_date
    use meshtype
	implicit none
	

	integer i,j,k,ip,jp,n,Nx1,Ny1
	real(kind=ps) x,y,x0,y0,x1,y1,x2,y2,res
	logical,allocatable :: ifbody(:,:)
    logical ifplate !是平板还是圆柱
	type(mesh_dat) mesh
	allocate (ifbody(Nx1,Ny1))
    
	ifplate=.true.
	ifbody=.false.
    if (ifplate.eq..true.) then
        do i=1,Nx1
		    do j=1,Ny1
			    x=dble(i-1)*mesh%scale+mesh%ox
			    y=dble(j-1)*mesh%scale+mesh%oy
			    x0=(x-xc)-sa_arr(1)
			    y0=(y-yc)-sb_arr(1) 
			    x1=cos(-angle_arr(1))*x0+sin(-angle_arr(1))*y0
			    y1=-sin(-angle_arr(1))*x0+cos(-angle_arr(1))*y0
                do k=1,Q
                    ip=i+e(k,1)
                    jp=j+e(k,2)
                    x=dble(ip-1)*mesh%scale+mesh%ox
			        y=dble(jp-1)*mesh%scale+mesh%oy
			        x0=(x-xc)-sa_arr(1)
			        y0=(y-yc)-sb_arr(1) 
			        x2=cos(-angle_arr(1))*x0+sin(-angle_arr(1))*y0
			        y2=-sin(-angle_arr(1))*x0+cos(-angle_arr(1))*y0
!			        res=1d-9
			        if((((x1.le.0.5d0).and.(x1.ge.-0.5d0)).or.((x2.le.0.5d0).and.(x2.ge.-0.5d0))).and.(y1*y2.le.0.0d0)) then
                        if(y1.le.0.0d0) then
				            mesh%body(i,j)=2 !平板下方的边界点
                            exit
                        else
                            mesh%body(i,j)=1 !平板上方的边界点
                            exit
                        endif
                    else
				        mesh%body(i,j)=3 !流体的内点
                    endif
                enddo
		    enddo
	    enddo
    else
	    do i=1,Nx1
		    do j=1,Ny1
			    x=dble(i-1)*mesh%scale+mesh%ox
			    y=dble(j-1)*mesh%scale+mesh%oy
			    x0=(x-xc)-sa_arr(1)
			    y0=(y-yc)-sb_arr(1) 
			    x1=cos(-angle_arr(1))*x0+sin(-angle_arr(1))*y0
			    y1=-sin(-angle_arr(1))*x0+cos(-angle_arr(1))*y0
			    res=1d-9
			    if(x1*x1+y1*y1.le.(rr*rr+res)) then
				    ifbody(i,j)=.true.
				    mesh%body(i,j)=1 !固体的内点
			    else
				    mesh%body(i,j)=3 !流体的内点
			    endif
		    enddo
	    enddo
	    do i=2,Nx1-1
		    do j=2,Ny1-1
			    do k=1,Q
				    ip=i+e(k,1)
				    jp=j+e(k,2)
				    if(ifbody(i,j)==.true..and.(ifbody(ip,jp)==.false.)) then
					    mesh%body(i,j)=2 !固体边界点
					    exit
				    elseif (ifbody(i,j)==.false..and.(ifbody(ip,jp)==.true.)) then
					    mesh%body(i,j)=4 !流体边界点
					    exit
				    endif
			    enddo
		    enddo
        enddo
    endif
	deallocate (ifbody)
end subroutine

Real (kind=ps) function distanceq(i,j,ip,jp,xb,yb,scale,ox,oy)
	use com_date
	implicit none

	integer i,j,ip,jp
	real (kind=ps) x,y,x1,y1,x2,y2,x0,y0,x3,y3,xb,yb,lamda1,lamda2,qr,qr1,qr2,scale,ox,oy
	real (kind=ps) a1,a2,r1,r2,r3,eps
	
	x=ox+dble(i-1)*scale
	y=oy+dble(j-1)*scale	
	x0=(x-xc)-sa_arr(1)
	y0=(y-xc)-sb_arr(1)
	x1=cos(-angle_arr(1))*x0+sin(-angle_arr(1))*y0
	y1=-sin(-angle_arr(1))*x0+cos(-angle_arr(1))*y0
	x=ox+dble(ip-1)*scale
	y=oy+dble(jp-1)*scale
	x0=(x-xc)-sa_arr(1)
	y0=(y-xc)-sb_arr(1)
	x2=cos(-angle_arr(1))*x0+sin(-angle_arr(1))*y0
	y2=-sin(-angle_arr(1))*x0+cos(-angle_arr(1))*y0

	qr=0.5d0
	qr1=0.0d0
	qr2=1.0d0
	eps=1.0d0
	do
		x0=x1+qr*(x2-x1)
		y0=y1+qr*(y2-y1)
		x=x1+qr1*(x2-x1)
		y=y1+qr1*(y2-y1)
		x3=x1+qr2*(x2-x1)
		y3=y1+qr2*(y2-y1)
		r1=sqrt(x0*x0+y0*y0)-rr
		r2=sqrt(x*x+y*y)-rr
		r3=sqrt(x3*x3+y3*y3)-rr
		if (abs(r1).lt.1d-9) then
			exit
		elseif(abs(r2).lt.1d-9) then
			qr=qr1
			exit
		elseif(abs(r3).lt.1d-9) then
			qr=qr2
			exit
		elseif (r1*r2.gt.0.0d0)then
			qr1=qr
			qr=0.5d0*(qr+qr2)
		else
			qr2=qr
			qr=0.5d0*(qr+qr1)
		endif
	enddo
	if(qr.gt.1.0d0) then
		qr=1.0d0
		print*, "qr>1!"
	endif
	if(qr.lt.0.0d0) then
		qr=0.0d0
		print*,"qr<0"
	endif
	xb=x1+qr*(x2-x1)
	yb=y1+qr*(y2-y1)
	distanceq=qr
!	if((i.eq.192).and.(j.eq.206).and.(ip.eq.193).and.(jp.eq.205)) print*, x1,y1, lamda1,lamda2

end function
    
Real (kind=ps) function distanceqline(i,j,ip,jp,xb,yb,scale,ox,oy)
	use com_date
	implicit none

	integer i,j,ip,jp
	real (kind=ps) x,y,x1,y1,x2,y2,x0,y0,x3,y3,xb,yb,lamda1,lamda2,qr,qr1,qr2,scale,ox,oy
	real (kind=ps) a1,a2,r1,r2,r3,eps
	
	x=ox+dble(i-1)*scale
	y=oy+dble(j-1)*scale	
	x0=(x-xc)-sa_arr(1)
	y0=(y-xc)-sb_arr(1)
	x1=cos(-angle_arr(1))*x0+sin(-angle_arr(1))*y0
	y1=-sin(-angle_arr(1))*x0+cos(-angle_arr(1))*y0
	x=ox+dble(ip-1)*scale
	y=oy+dble(jp-1)*scale
	x0=(x-xc)-sa_arr(1)
	y0=(y-xc)-sb_arr(1)
	x2=cos(-angle_arr(1))*x0+sin(-angle_arr(1))*y0
	y2=-sin(-angle_arr(1))*x0+cos(-angle_arr(1))*y0

    xb=-y1*(x2-x1)/(y2-y1)+x1
    yb=0.0d0
    qr=(xb-x1)/(x2-x1)
	if(qr.gt.1.0d0) then
        
		qr=1.0d0
!		print*, "qr>1!"
	endif
	if(qr.lt.0.0d0) then
		qr=0.0d0
!		print*,"qr<0"
	endif
	distanceqline=qr
!	if((i.eq.192).and.(j.eq.206).and.(ip.eq.193).and.(jp.eq.205)) print*, x1,y1, lamda1,lamda2

end function    

subroutine build_boundary_filename(wing_idx, filename)
	use com_date
	implicit none

	integer, intent(in) :: wing_idx
	character(len=*), intent(out) :: filename
	character(len=32) :: idx_text

	if (wing_idx.eq.1) then
		filename='boundary.dat'
	else
		write(idx_text,'(i0)') wing_idx-1
		filename='boundary'//trim(idx_text)//'.dat'
	endif
end subroutine

subroutine read_boundary(mesh)
	use com_date
    use meshtype
	implicit none
	
	type(mesh_dat) mesh
	integer i0,wing,point_start,seg_start,unit_id,ios
	character(len=256) :: filename

	xb0=0.0d0
	yb0=0.0d0
	arc_l=0.0d0
	do wing=1,wing_count
		point_start=wing_point_start(wing)
		seg_start=wing_segment_start(wing)
		call build_boundary_filename(wing,filename)
		unit_id=421+wing
		open(unit_id,file=trim(filename),status='old',iostat=ios)
		if (ios.ne.0) then
			write(*,*) '无法打开边界文件: ', trim(filename)
			stop 1
		endif
		read(unit_id,'(2f14.8)',iostat=ios) ((xb0(i0),yb0(i0)),i0=point_start,point_start+boundary_points_per_wing-1)
		close(unit_id)
		if (ios.ne.0) then
			write(*,*) '边界文件格式错误: ', trim(filename)
			stop 1
		endif
			do i0=0,boundary_segments_per_wing-1
				arc_l(seg_start+i0)=sqrt((xb0(point_start+i0+1)-xb0(point_start+i0))*(xb0(point_start+i0+1)-xb0(point_start+i0)) &
				    +(yb0(point_start+i0+1)-yb0(point_start+i0))*(yb0(point_start+i0+1)-yb0(point_start+i0)))/mesh%scale
			enddo
		enddo
    xb0(1:b_point)=xb0(1:b_point)+(xc-xp)
    yb0(1:b_point)=yb0(1:b_point)+(yc-yp)
end subroutine

subroutine move_boundary(mesh)
    use mpi
	use com_date
    use meshtype
	implicit none
	
	type(mesh_dat) mesh
	integer i0,i1,i2,j1,j2,i,j,wing,point_start,seg_start
	real(kind=ps) temp,temp1
    real(kind=ps) :: check_x,check_y,margin
	
    ! Keep the complete five-point IBM kernel inside the finest interior.
    margin=3.5d0*mesh%scale
    do wing=1,wing_count
        point_start=wing_point_start(wing)
        do i0=point_start,point_start+boundary_points_per_wing-1
            check_x=cos(angle_arr(wing))*xb0(i0)+sin(angle_arr(wing))*yb0(i0) &
                +wing_initial_x_origin+(wing-1)*wing_initial_x_spacing+sa_arr(wing)+xc+body_extra_s(wing)*mesh%scale
            check_y=-sin(angle_arr(wing))*xb0(i0)+cos(angle_arr(wing))*yb0(i0)+wing_initial_y_origin+sb_arr(wing)+yc
            if(check_x.lt.mesh%ox+margin.or.check_x.gt.mesh%ox+(mesh%Nx-1)*mesh%scale-margin.or. &
               check_y.lt.mesh%oy+margin.or.check_y.gt.mesh%oy+(mesh%Ny-1)*mesh%scale-margin) then
                if(myid.eq.0) print*, 'Wing',wing,'leaves finest mesh IBM support at time',time,check_x,check_y
                call MPI_ABORT(MPI_COMM_WORLD,913,ierr)
            endif
        enddo
    enddo
	!$acc update device (angle_arr,sa_arr,sb_arr,body_extra_s,movenum)
	!$acc kernels present(xb0,yb0,xb1,yb1,xb2,yb2,angle_arr,sa_arr,sb_arr,body_extra_s,mesh,movenum,adjoint_p,adjoint_wing_id,gforce,adjoint_n)
	!$acc loop seq
	do wing=1,wing_count
		point_start=wing_point_start(wing)
		seg_start=wing_segment_start(wing)
		!$acc loop
		do i0=point_start,point_start+boundary_points_per_wing-1
				xb1(i0)=cos(angle_arr(wing))*xb0(i0)+sin(angle_arr(wing))*yb0(i0)+wing_initial_x_origin+dble(wing-1)*wing_initial_x_spacing+sa_arr(wing)+xc &
				    +body_extra_s(wing)/(1.0d0/mesh%scale)-dble(movenum)*mesh_max_scale
				yb1(i0)=-sin(angle_arr(wing))*xb0(i0)+cos(angle_arr(wing))*yb0(i0)+wing_initial_y_origin+sb_arr(wing)+yc
			enddo
		!$acc end loop
		!$acc loop
		do i0=seg_start,seg_start+boundary_segments_per_wing-1
			xb2(i0)=0.5d0*(xb1(point_start+i0-seg_start+1)+xb1(point_start+i0-seg_start))
			yb2(i0)=0.5d0*(yb1(point_start+i0-seg_start+1)+yb1(point_start+i0-seg_start))
		enddo
		!$acc end loop
	enddo
	!$acc end loop

	!$acc loop
	do i=1,adjoint_n
		i1=adjoint_p(i,1)
		j1=adjoint_p(i,2)
        mesh%body(i1,j1)=0
        gforce(i1,j1,1)=0.0d0
		gforce(i1,j1,2)=0.0d0
		adjoint_wing_id(i)=0
    enddo
	!$acc end loop

	adjoint_n=1
	!$acc loop seq
	do i0=1,b_segment_count
		wing=(i0-1)/boundary_segments_per_wing+1
		i1=int((xb2(i0)-mesh%ox)/mesh%scale+1.0d0)
		do i=i1,1,-1
			temp=mesh%ox+dble(i-1)*mesh%scale
			if(abs((temp-xb2(i0))/mesh%scale).gt.2.5d0) then
				i1=i+1
				exit
			endif
		enddo
		i2=int((xb2(i0)-mesh%ox)/mesh%scale+1.0d0)
		do i=i2,mesh%Nx
			temp=mesh%ox+dble(i-1)*mesh%scale
			if(abs((temp-xb2(i0))/mesh%scale).gt.2.5d0) then
				i2=i-1
				exit
			endif
		enddo
		j2=int((yb2(i0)-mesh%oy)/mesh%scale+1.0d0)
		do i=j2,mesh%Ny
			temp=mesh%oy+dble(i-1)*mesh%scale
			if(abs((temp-yb2(i0))/mesh%scale).gt.2.5d0) then
				j2=i-1
				exit
			endif
		enddo
		j1=int((yb2(i0)-mesh%oy)/mesh%scale+1.0d0)
		do i=j1,1,-1
			temp=mesh%oy+dble(i-1)*mesh%scale
			if(abs((temp-yb2(i0))/mesh%scale).gt.2.5d0) then
				j1=i+1
				exit
			endif
        enddo
		do i=i1,i2
			do j=j1,j2
				if(mesh%body(i,j).eq.0)then
					mesh%body(i,j)=1
	                    adjoint_p(adjoint_n,1)=i
	                    adjoint_p(adjoint_n,2)=j
	                    adjoint_wing_id(adjoint_n)=wing
	                    adjoint_n=adjoint_n+1
	                endif
			enddo
        enddo
	enddo
	!$acc end loop
    adjoint_n=adjoint_n-1
	!$acc end kernels
	!-------------------------------------------------------------------------------
end subroutine
