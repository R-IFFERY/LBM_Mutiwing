Subroutine D2G9model
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

!Real (kind=ps) function D2G9feq(p,u1,v1,i,j,k)
!	use com_date
!	implicit none
!	real (kind=ps) p,u1,v1
!	integer i,j,k
!	real (kind=ps) eu,utwo,d0,d1,d2,ftemp 

!	d0=5.0/12.0
!	d1=1.0/3.0
!	d2=1.0/12.0

!	eu=e(k,1)*u1+e(k,2)*v1
!	utwo=u1*u1+v1*v1
!	sa(i,j,k)=w(k)*(3.0*eu+4.5*eu*eu-1.5*utwo)
!	if(k.ne.1) ftemp=ftemp+f(i,j,k)
!	if(k.eq.1) then
!		D2G9feq=rho0-4.0*d0*p+sa(i,j,k)
!	elseif ((k.le.5).and.(k.ge.2)) then
!		D2G9feq=d1*p+sa(i,j,k)
!	else
!		D2G9feq=d2*p+sa(i,j,k)
!	endif
!	return
!end function

!subroutine boundaryD2G9
!	use com_date
!	implicit none
!	Real (kind=ps), external :: D2G9feq,distanceq
!	integer i,j,k,ip,jp
!	real (kind=ps) qr,preb,ub,vb

	!边界条件
!	do j=2,Ny-1 !左右边界
!		do k=1,Q
!			pre(Nx,j)=pre(Nx-1,j)
!			u(Nx,j)=u(Nx-1,j)
!			v(Nx,j)=v(Nx-1,j)
!			f(Nx,j,k)=D2G9feq(Pre(Nx,j),u(Nx,j),v(Nx,j),Nx,j,k)+f(Nx-1,j,k)-D2G9feq(Pre(Nx-1,j),u(Nx-1,j),v(Nx-1,j),Nx-1,j,k)
!			pre(1,j)=pre(2,j)
!			u(1,j)=uu
!			v(1,j)=0.0d0
!			f(1,j,k)=D2G9feq(Pre(1,j),u(1,j),v(1,j),1,j,k)+f(2,j,k)-D2G9feq(Pre(2,j),u(2,j),v(2,j),2,j,k)
!		enddo
!	enddo

!	do i=1,Nx !上下边界
!		do k=1,Q
!			pre(i,1)=pre(i,2)
!			u(i,1)=u(i,2)
!			v(i,1)=0.0d0
!			f(i,1,k)=D2G9feq(Pre(i,1),u(i,1),v(i,1),i,1,k)+f(i,2,k)-D2G9feq(Pre(i,2),u(i,2),v(i,2),i,2,k)

!			u(i,Ny)=uu
!			pre(i,Ny)=pre(i,Ny-1)
!			u(i,Ny)=u(i,Ny-1)
!			v(i,Ny)=0.0d0
!			f(i,Ny,k)=D2G9feq(Pre(i,Ny),u(i,Ny),v(i,Ny),i,Ny,K)+f(i,Ny-1,k)-D2G9feq(Pre(i,Ny-1),u(i,Ny-1),v(i,Ny-1),i,Ny-1,k)
!		enddo
!	enddo

	!圆柱边界
!	do i=2,Nx-1
!		do j=2,Ny-1
!			do k=1,Q			
!				ip=i-e(k,1)
!				jp=j-e(k,2)			
!
!				if ((body(i,j).eq..false.).and. (body(ip,jp).eq..true.)) then
!					qr=distanceq(i,j,ip,jp)
!					if(qr.ge.0.75d0) then
!						preb=pre(i,j)
!						ub=(qr-1.0d0)/qr*u(i,j)
!						vb=(qr-1.0d0)/qr*v(i,j)
!						f(ip,jp,k)=D2G9feq(preb,ub,vb,ip,jp,k)+(1.0d0-1.0d0/tau_f)*(f(i,j,k)-D2G9feq(pre(i,j),u(i,j),v(i,j),i,j,k))
!					else
!						preb=pre(i,j)
!						ub=(qr-1.0d0)*u(i,j)+(1.0d0-qr)*(qr-1.0d0)/(1.0d0+qr)*u(i+e(k,1),j+e(k,2))
!						vb=(qr-1.0d0)*v(i,j)+(1.0d0-qr)*(qr-1.0d0)/(1.0d0+qr)*v(i+e(k,1),j+e(k,2))
!						f(ip,jp,k)=D2G9feq(preb,ub,vb,ip,jp,k)+(1.0d0-1.0d0/tau_f)*(qr*(f(i,j,k)-D2G9feq(pre(i,j),u(i,j),v(i,j),i,j,k))+(1.0d0-qr)*(f(i+e(k,1),j+e(k,2),k)-D2G9feq(pre(i+e(k,1),j+e(k,2)),u(i+e(k,1),j+e(k,2)),v(i+e(k,1),j+e(k,2)),i+e(k,1),j+e(k,2),k)))
!					endif
!				endif
!			enddo
!		enddo
!	enddo
!end