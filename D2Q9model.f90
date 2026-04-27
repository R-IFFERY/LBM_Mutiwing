Subroutine D2Q9model
	use com_date
	implicit none

	e(1,1)=0
	e(1,2)=0
	e(2,1)=1
	e(2,2)=0
	e(3,1)=0
	e(3,2)=1
	e(4,1)=-1
	e(4,2)=0
	e(5,1)=0
	e(5,2)=-1
	e(6,1)=1
	e(6,2)=1
	e(7,1)=-1
	e(7,2)=1
	e(8,1)=-1
	e(8,2)=-1
	e(9,1)=1
	e(9,2)=-1

	w(1)=4.0/9.0
	w(2)=1.0/9.0
	w(3)=1.0/9.0
	w(4)=1.0/9.0
	w(5)=1.0/9.0
	w(6)=1.0/36.0
	w(7)=1.0/36.0
	w(8)=1.0/36.0
	w(9)=1.0/36.0
end Subroutine

Real*8 function D2Q9feq(rh,u1,v1,i,j,k)
	!$acc routine seq
	use com_date
	implicit none
	integer i,j,k
	real (kind=ps) eu,utwo,rh,u1,v1

	eu=e(k,1)*u1+e(k,2)*v1
	utwo=u1*u1+v1*v1
	D2Q9feq= w(k)*rh*(1.0+3.0*eu+4.5*eu*eu-1.5*utwo)

end function D2Q9feq

subroutine boundary1D2Q9
	use com_date
    use meshtype
	implicit none
	
	type(mesh_dat) mesh(max_meshid)
	Real (kind=ps), external :: D2Q9feq,distanceq
	integer i,j,k,k1,ip,jp,sp,ep
	real (kind=ps) qr,rhob,ub,vb,xb,yb

	if(myid.eq.0) then
		sp=mesh(1)%start_row+1
		ep=mesh(1)%end_row-3
	elseif(myid.eq.numprocs-1) then
		sp=mesh(1)%start_row+3
		ep=mesh(1)%end_row-1
	else
		sp=mesh(1)%start_row+3
		ep=mesh(1)%end_row-3
	endif
	if (numprocs.eq.1) then
		sp=mesh(1)%start_row+1
		ep=mesh(1)%end_row-1
	endif

	!边界条件
	do j=2,mesh(1)%Ny-1 !左右边界
		do k=1,Q
			if(myid.eq.numprocs-1) then
				mesh(1)%rho(mesh(1)%Nx,j)=1.0d0
				mesh(1)%u(mesh(1)%Nx,j)=mesh(1)%u(mesh(1)%Nx-1,j)
				mesh(1)%v(mesh(1)%Nx,j)=mesh(1)%v(mesh(1)%Nx-1,j)
				mesh(1)%f(mesh(1)%Nx,j,k)=D2Q9feq(mesh(1)%rho(mesh(1)%Nx,j),mesh(1)%u(mesh(1)%Nx,j),mesh(1)%v(mesh(1)%Nx,j),mesh(1)%Nx,j,k)+mesh(1)%f(mesh(1)%Nx-1,j,k)-D2Q9feq(mesh(1)%rho(mesh(1)%Nx-1,j),mesh(1)%u(mesh(1)%Nx-1,j),mesh(1)%v(mesh(1)%Nx-1,j),mesh(1)%Nx-1,j,k)
			endif
			if(myid.eq.0) then
				mesh(1)%rho(1,j)=1.0d0
				mesh(1)%u(1,j)=UU
				mesh(1)%v(1,j)=0.0d0
				mesh(1)%f(1,j,k)=D2Q9feq(mesh(1)%rho(1,j),mesh(1)%u(1,j),mesh(1)%v(1,j),1,j,k)+mesh(1)%f(2,j,k)-D2Q9feq(mesh(1)%rho(2,j),mesh(1)%u(2,j),mesh(1)%v(2,j),2,j,k)
			endif
		enddo
	enddo

	do i=mesh(1)%start_row,mesh(1)%end_row !上下边界
		do k=1,Q
			mesh(1)%rho(i,1)=1.0d0
			mesh(1)%u(i,1)=mesh(1)%u(i,2)
			mesh(1)%v(i,1)=0.0d0
			mesh(1)%f(i,1,k)=D2Q9feq(mesh(1)%rho(i,1),mesh(1)%u(i,1),mesh(1)%v(i,1),i,1,k)+mesh(1)%f(i,2,k)-D2Q9feq(mesh(1)%rho(i,2),mesh(1)%u(i,2),mesh(1)%v(i,2),i,2,k)
			mesh(1)%rho(i,mesh(1)%Ny)=1.0d0
			mesh(1)%u(i,mesh(1)%Ny)=mesh(1)%u(i,mesh(1)%Ny-1)
			mesh(1)%v(i,mesh(1)%Ny)=0.0d0
			mesh(1)%f(i,mesh(1)%Ny,k)=D2Q9feq(mesh(1)%rho(i,mesh(1)%Ny),mesh(1)%u(i,mesh(1)%Ny),mesh(1)%v(i,mesh(1)%Ny),i,mesh(1)%Ny,K)+mesh(1)%f(i,mesh(1)%Ny-1,k)-D2Q9feq(mesh(1)%rho(i,mesh(1)%Ny-1),mesh(1)%u(i,mesh(1)%Ny-1),mesh(1)%v(i,mesh(1)%Ny-1),i,mesh(1)%Ny-1,k)
		enddo
	enddo

	!圆柱边界
	myforce=0.0d0
	do i=mesh(1)%start_row,mesh(1)%end_row
		do j=2,mesh(1)%Ny-1
			do k=1,Q			
				ip=i-e(k,1)
				jp=j-e(k,2)			

				if ((mesh(1)%body(i,j).eq.4).and. (mesh(1)%body(ip,jp).eq.2)) then
					qr=distanceq(i,j,ip,jp,xb,yb,mesh(1)%scale,mesh(1)%ox,mesh(1)%oy)
					call reversedirection(k,k1)
					if(i.ge.sp.and.i.le.ep) then
						myforce(1)=myforce(1)+e(k1,1)*(mesh(1)%f(i,j,k1)+(D2Q9feq(mesh(1)%rho(i,j),mesh(1)%u(i,j),mesh(1)%v(i,j),i,j,k1)-mesh(1)%f(i,j,k1))/mesh(1)%tau_f)
						myforce(2)=myforce(2)+e(k1,2)*(mesh(1)%f(i,j,k1)+(D2Q9feq(mesh(1)%rho(i,j),mesh(1)%u(i,j),mesh(1)%v(i,j),i,j,k1)-mesh(1)%f(i,j,k1))/mesh(1)%tau_f)
					endif
!					myforce(1)=myforce(1)-e(k,1)*f(i,j,k)
!					myforce(2)=myforce(2)-e(k,2)*f(i,j,k)
!					qr=distanceq(i,j,ip,jp)
					if(qr.ge.0.75d0) then
						mesh(1)%rho(ip,jp)=mesh(1)%rho(i,j)
						mesh(1)%u(ip,jp)=(qr-1.0d0)*mesh(1)%u(i,j)/qr
						mesh(1)%v(ip,jp)=(sbv_arr(1)*UU+(qr-1.0d0)*mesh(1)%v(i,j))/qr
						mesh(1)%f(ip,jp,k)=D2Q9feq(mesh(1)%rho(ip,jp),mesh(1)%u(ip,jp),mesh(1)%v(ip,jp),ip,jp,k)+(1.0d0-1.0/mesh(1)%tau_f)*(mesh(1)%f(i,j,k)-D2Q9feq(mesh(1)%rho(i,j),mesh(1)%u(i,j),mesh(1)%v(i,j),i,j,k))
					else
						mesh(1)%rho(ip,jp)=mesh(1)%rho(i,j)
						mesh(1)%u(ip,jp)=(qr-1.0d0)*mesh(1)%u(i,j)+(1.0d0-qr)*(qr-1.0d0)/(1.0d0+qr)*mesh(1)%u(i+e(k,1),j+e(k,2))
						mesh(1)%v(ip,jp)=sbv_arr(1)*UU+(qr-1.0d0)*mesh(1)%v(i,j)+(1.0d0-qr)*(2.0*sbv_arr(1)*UU+(qr-1.0d0)*mesh(1)%v(i+e(k,1),j+e(k,2)))/(1.0d0+qr)
						mesh(1)%f(ip,jp,k)=D2Q9feq(mesh(1)%rho(ip,jp),mesh(1)%u(ip,jp),mesh(1)%v(ip,jp),ip,jp,k)+(1.0d0-1.0/mesh(1)%tau_f)*(qr*(mesh(1)%f(i,j,k)-D2Q9feq(mesh(1)%rho(i,j),mesh(1)%u(i,j),mesh(1)%v(i,j),i,j,k))+(1.0d0-qr)*(mesh(1)%f(i+e(k,1),j+e(k,2),k)-D2Q9feq(mesh(1)%rho(i+e(k,1),j+e(k,2)),mesh(1)%u(i+e(k,1),j+e(k,2)),mesh(1)%v(i+e(k,1),j+e(k,2)),i+e(k,1),j+e(k,2),k)))
					endif
					if(ip.ge.sp.and.ip.le.ep) then
						myforce(1)=myforce(1)-e(k,1)*mesh(1)%f(ip,jp,k)
						myforce(2)=myforce(2)-e(k,2)*mesh(1)%f(ip,jp,k)
					endif
!					if(myid.eq.0.and.ip.eq.321.and.jp.eq.444.and.k.eq.7) print*,f(ip,jp,7),i,j,"1"
				endif
			enddo
		enddo
	enddo
end

subroutine boundary2D2Q9(mesh)
	use com_date
    use meshtype
	implicit none
	
	type(mesh_dat) mesh(max_meshid)
	Real (kind=ps), external :: D2Q9feq,distanceq
	integer i,j,k,k1,ip,jp,sp,ep
	real (kind=ps) qr,rhob,ub,vb,wa,xb,yb

	if(myid.eq.0) then
		sp=mesh(1)%start_row+1
		ep=mesh(1)%end_row-3
	elseif(myid.eq.numprocs-1) then
		sp=mesh(1)%start_row+3
		ep=mesh(1)%end_row-1
	else
		sp=mesh(1)%start_row+3
		ep=mesh(1)%end_row-3
	endif
	if (numprocs.eq.1) then
		sp=mesh(1)%start_row+1
		ep=mesh(1)%end_row-1
	endif

	!边界条件
	do j=2,mesh(1)%Ny-1 !左右边界
		do k=1,Q
			if(myid.eq.numprocs-1) then
				mesh(1)%rho(mesh(1)%Nx,j)=1.0d0
				mesh(1)%u(mesh(1)%Nx,j)=mesh(1)%u(mesh(1)%Nx-1,j)
				mesh(1)%v(mesh(1)%Nx,j)=mesh(1)%v(mesh(1)%Nx-1,j)
				mesh(1)%f0(mesh(1)%Nx,j,k)=D2Q9feq(mesh(1)%rho(mesh(1)%Nx,j),mesh(1)%u(mesh(1)%Nx,j),mesh(1)%v(mesh(1)%Nx,j),mesh(1)%Nx,j,k)+mesh(1)%f(mesh(1)%Nx-1,j,k)-D2Q9feq(mesh(1)%rho(mesh(1)%Nx-1,j),mesh(1)%u(mesh(1)%Nx-1,j),mesh(1)%v(mesh(1)%Nx-1,j),mesh(1)%Nx-1,j,k)
			endif
			if(myid.eq.0) then
				mesh(1)%rho(1,j)=1.0d0
				mesh(1)%u(1,j)=UU
				mesh(1)%v(1,j)=0.0d0
				mesh(1)%f0(1,j,k)=D2Q9feq(mesh(1)%rho(1,j),mesh(1)%u(1,j),mesh(1)%v(1,j),1,j,k)+mesh(1)%f(2,j,k)-D2Q9feq(mesh(1)%rho(2,j),mesh(1)%u(2,j),mesh(1)%v(2,j),2,j,k)
			endif
		enddo
	enddo

	do i=mesh(1)%start_row,mesh(1)%end_row !上下边界
		do k=1,Q
			mesh(1)%rho(i,1)=1.0d0
			mesh(1)%u(i,1)=mesh(1)%u(i,2)
			mesh(1)%v(i,1)=0.0d0
			mesh(1)%f0(i,1,k)=D2Q9feq(mesh(1)%rho(i,1),mesh(1)%u(i,1),mesh(1)%v(i,1),i,1,k)+mesh(1)%f(i,2,k)-D2Q9feq(mesh(1)%rho(i,2),mesh(1)%u(i,2),mesh(1)%v(i,2),i,2,k)
			mesh(1)%rho(i,mesh(1)%Ny)=1.0d0
			mesh(1)%u(i,mesh(1)%Ny)=mesh(1)%u(i,mesh(1)%Ny-1)
			mesh(1)%v(i,mesh(1)%Ny)=0.0d0
			mesh(1)%f0(i,mesh(1)%Ny,k)=D2Q9feq(mesh(1)%rho(i,mesh(1)%Ny),mesh(1)%u(i,mesh(1)%Ny),mesh(1)%v(i,mesh(1)%Ny),i,mesh(1)%Ny,K)+mesh(1)%f(i,mesh(1)%Ny-1,k)-D2Q9feq(mesh(1)%rho(i,mesh(1)%Ny-1),mesh(1)%u(i,mesh(1)%Ny-1),mesh(1)%v(i,mesh(1)%Ny-1),i,mesh(1)%Ny-1,k)
		enddo
	enddo

	!圆柱边界
	myforce=0.0d0
	do i=sp,ep
		do j=2,mesh(1)%Ny-1
			do k=1,Q			
				ip=i-e(k,1)
				jp=j-e(k,2)			

				if ((mesh(1)%body(i,j).eq.4).and. (mesh(1)%body(ip,jp).eq.2)) then
					qr=distanceq(i,j,ip,jp,xb,yb,mesh(1)%scale,mesh(1)%ox,mesh(1)%oy)
					ub=UU*(-omega_arr(1)*sin(angle_arr(1))*xb+omega_arr(1)*cos(angle_arr(1))*yb-sau_arr(1))
					vb=UU*(-omega_arr(1)*cos(angle_arr(1))*xb-omega_arr(1)*sin(angle_arr(1))*yb+sbv_arr(1))
					call reversedirection(k,k1)
					myforce(1)=myforce(1)+e(k1,1)*mesh(1)%f(ip,jp,k1)
					myforce(2)=myforce(2)+e(k1,2)*mesh(1)%f(ip,jp,k1)
					if((k1.ge.2).and.(k1.le.5)) then
						wa=2.0/9.0
					elseif((k1.ge.6).and.(k1.le.9)) then
						wa=2.0/36.0
					endif
					if(qr.lt.0.5d0) then
						mesh(1)%f0(i,j,k)=qr*(1.0+2.0*qr)*mesh(1)%f(i-e(k,1),j-e(k,2),k1)+(1.0d0-4.0d0*qr*qr)*mesh(1)%f(i,j,k1)-qr*(1.0d0-2.0*qr)*mesh(1)%f(i+e(k,1),j+e(k,2),k1)+3.0*wa*mesh(1)%rho(i,j)*(dble(e(k1,1))*ub+dble(e(k1,2))*vb)
					else
						mesh(1)%f0(i,j,k)=mesh(1)%f(i-e(k,1),j-e(k,2),k1)/(qr*(2.0*qr+1.0))+(2.0*qr-1.0)/qr*mesh(1)%f(i+e(k,1),j+e(k,2),k)-(2.0*qr-1.0)/(2.0*qr+1.0)*mesh(1)%f(i+2*e(k,1),j+2*e(k,2),k)+3.0*wa*mesh(1)%rho(i,j)*(dble(e(k1,1))*ub+dble(e(k1,2))*vb)/(qr*(2.0*qr+1.0))
					endif
					myforce(1)=myforce(1)-e(k,1)*mesh(1)%f0(i,j,k)
					myforce(2)=myforce(2)-e(k,2)*mesh(1)%f0(i,j,k)
				endif
			enddo

		enddo
	enddo
end


subroutine boundary3D2Q9
	use com_date
    use meshtype
	implicit none
	
	type(mesh_dat) mesh(max_meshid)
	Real (kind=ps), external :: D2Q9feq,distanceq
	integer i,j,k,k1,ip,jp,sp,ep
	real (kind=ps) qr,rhob,ub,vb,wa,temp1,xb,yb

	if(myid.eq.0) then
		sp=mesh(1)%start_row+1
		ep=mesh(1)%end_row-3
	elseif(myid.eq.numprocs-1) then
		sp=mesh(1)%start_row+3
		ep=mesh(1)%end_row-1
	else
		sp=mesh(1)%start_row+3
		ep=mesh(1)%end_row-3
	endif
	if (numprocs.eq.1) then
		sp=mesh(1)%start_row+1
		ep=mesh(1)%end_row-1
	endif

	!边界条件
	do j=2,mesh(1)%Ny-1 !左右边界
		do k=1,Q
			if(myid.eq.numprocs-1) then
				mesh(1)%rho(mesh(1)%Nx,j)=1.0d0
				mesh(1)%u(mesh(1)%Nx,j)=mesh(1)%u(mesh(1)%Nx-1,j)
				mesh(1)%v(mesh(1)%Nx,j)=mesh(1)%v(mesh(1)%Nx-1,j)
				mesh(1)%f0(mesh(1)%Nx,j,k)=D2Q9feq(mesh(1)%rho(mesh(1)%Nx,j),mesh(1)%u(mesh(1)%Nx,j),mesh(1)%v(mesh(1)%Nx,j),mesh(1)%Nx,j,k)+mesh(1)%f(mesh(1)%Nx-1,j,k)-D2Q9feq(mesh(1)%rho(mesh(1)%Nx-1,j),mesh(1)%u(mesh(1)%Nx-1,j),mesh(1)%v(mesh(1)%Nx-1,j),mesh(1)%Nx-1,j,k)
			endif
			if(myid.eq.0) then
				mesh(1)%rho(1,j)=1.0d0
				mesh(1)%u(1,j)=uu
				mesh(1)%v(1,j)=0.0
				mesh(1)%f0(1,j,k)=D2Q9feq(mesh(1)%rho(1,j),mesh(1)%u(1,j),mesh(1)%v(1,j),1,j,k)+mesh(1)%f(2,j,k)-D2Q9feq(mesh(1)%rho(2,j),mesh(1)%u(2,j),mesh(1)%v(2,j),2,j,k)
			endif
		enddo
	enddo

	do i=mesh(1)%start_row,mesh(1)%end_row !上下边界
		do k=1,Q
			mesh(1)%rho(i,1)=1.0d0
			mesh(1)%u(i,1)=mesh(1)%u(i,2)
			mesh(1)%v(i,1)=0.0d0
			mesh(1)%f0(i,1,k)=D2Q9feq(mesh(1)%rho(i,1),mesh(1)%u(i,1),mesh(1)%v(i,1),i,1,k)+mesh(1)%f(i,2,k)-D2Q9feq(mesh(1)%rho(i,2),mesh(1)%u(i,2),mesh(1)%v(i,2),i,2,k)
			mesh(1)%rho(i,mesh(1)%Ny)=1.0d0
			mesh(1)%u(i,mesh(1)%Ny)=mesh(1)%u(i,mesh(1)%Ny-1)
			mesh(1)%v(i,mesh(1)%Ny)=0.0d0
			mesh(1)%f0(i,mesh(1)%Ny,k)=D2Q9feq(mesh(1)%rho(i,mesh(1)%Ny),mesh(1)%u(i,mesh(1)%Ny),mesh(1)%v(i,mesh(1)%Ny),i,mesh(1)%Ny,K)+mesh(1)%f(i,mesh(1)%Ny-1,k)-D2Q9feq(mesh(1)%rho(i,mesh(1)%Ny-1),mesh(1)%u(i,mesh(1)%Ny-1),mesh(1)%v(i,mesh(1)%Ny-1),i,mesh(1)%Ny-1,k)
		enddo
	enddo

	!圆柱边界
	myforce=0.0d0
	do i=sp,ep
		do j=2,mesh(1)%Ny-1
			do k=1,Q			
				ip=i-e(k,1)
				jp=j-e(k,2)			

				if ((mesh(1)%body(i,j).eq.4).and. (mesh(1)%body(ip,jp).eq.2)) then
					qr=distanceq(i,j,ip,jp,xb,yb,mesh(1)%scale,mesh(1)%ox,mesh(1)%oy)
					call reversedirection(k,k1)
					myforce(1)=myforce(1)+e(k1,1)*mesh(1)%f(ip,jp,k1)
					myforce(2)=myforce(2)+e(k1,2)*mesh(1)%f(ip,jp,k1)
					if((k1.ge.2).and.(k1.le.5)) then
						wa=2.0/9.0
					elseif((k1.ge.6).and.(k1.le.9)) then
						wa=2.0/36.0
					endif
					temp1=mesh(1)%f(i,j,k1)+qr*(mesh(1)%f(ip,jp,k1)-mesh(1)%f(i,j,k1))+0.5d0*qr*(qr-1.0d0)*(mesh(1)%f(ip,jp,k1)-2.0d0*mesh(1)%f(i,j,k1)+mesh(1)%f(i+e(k,1),j+e(k,2),k1))+3.0*wa*mesh(1)%rho(i,j)*(dble(e(k1,1))*(-sau_arr(1)*UU)+dble(e(k1,2))*(sbv_arr(1)*UU))
					mesh(1)%f0(i,j,k)=temp1+qr/(1.0d0+qr)*(mesh(1)%f(i+e(k,1),j+e(k,2),k)-temp1)-(mesh(1)%f(i+2*e(k,1),j+2*e(k,2),k)/((2.0d0+qr)/qr)-mesh(1)%f(i+e(k,1),j+e(k,2),k)/((1.0d0+qr)/qr)+qr*temp1/((2.0d0+qr)*(1.0d0+qr)))
					myforce(1)=myforce(1)-e(k,1)*mesh(1)%f0(i,j,k)
					myforce(2)=myforce(2)-e(k,2)*mesh(1)%f0(i,j,k)
				endif
			enddo
		enddo
	enddo
end
