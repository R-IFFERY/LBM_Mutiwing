subroutine boundary_condition (mesh,shuzu,n) !远场边界条件
	use com_date
    use meshtype
	implicit none
	
	type(mesh_dat) mesh(max_meshid)
	integer n,i,j,k
    integer n1,n2,n3
	integer shuzu(n)
    real(kind=ps) rh,jx,jy,temp(Q),temp1(Q),temp2(Q),temp3(Q)
	Real (kind=ps), external :: D2Q9feq,distanceq

	i=1
    if(myid.eq.0) then
        if(ref_freestream.le.1)then
            if(mesh(1)%left(2,1).eq.0) then
                if(periodical_bc.eq.1)then
                    !$acc kernels loop independent present(mesh)
                    do k=1,Q
                        do j=2,mesh(1)%Ny-1
                            mesh(1)%f(1,j,k)=mesh(1)%f(mesh(1)%Nx-1,j,k)
                        enddo
                    enddo
                    !$acc end kernels loop
                else
                    do j=2,mesh(1)%Ny-1
                        mesh(1)%rho(1,j)=1.0d0
                        mesh(1)%u(1,j)=uu
                        mesh(1)%v(1,j)=0.0d0
                        rh=mesh(1)%rho(1,j)
                        jx=rh*mesh(1)%u(1,j)
                        jy=rh*mesh(1)%v(1,j)
                        temp1(1)=rh
                        temp1(2)=-2.0d0*rh+3.0d0*(jx*jx+jy*jy)/rh
                        temp1(3)=rh-3.0d0*(jx*jx+jy*jy)/rh
                        temp1(4)=jx
                        temp1(5)=-jx
                        temp1(6)=jy
                        temp1(7)=-jy
                        temp1(8)=(jx*jx-jy*jy)/rh
                        temp1(9)=(jx*jy)/rh
                        rh=mesh(1)%rho(2,j)
                        jx=rh*mesh(1)%u(2,j)
                        jy=rh*mesh(1)%v(2,j)
                        temp2(1)=rh
                        temp2(2)=-2.0d0*rh+3.0d0*(jx*jx+jy*jy)/rh
                        temp2(3)=rh-3.0d0*(jx*jx+jy*jy)/rh
                        temp2(4)=jx
                        temp2(5)=-jx
                        temp2(6)=jy
                        temp2(7)=-jy
                        temp2(8)=(jx*jx-jy*jy)/rh
                        temp2(9)=(jx*jy)/rh
                        temp=mesh(1)%f0(2,j,:)
                        call matmulvoc(temp3,matm,temp,Q,Q)
                        do k=1,Q
                            temp3(k)=temp1(k)+(1.0d0-mesh(1)%s(k))*(temp3(k)-temp2(k))
                        enddo
                        call matmulvoc(temp,inversematm,temp3,Q,Q)
                        mesh(1)%f(1,j,:)=temp
    !					    do k=1,Q
    !						    mesh(1)%rho(1,j)=1.0d0
    !						    mesh(1)%u(1,j)=UU
    !						    mesh(1)%v(1,j)=0.0d0
    !						    mesh(1)%f(1,j,k)=D2Q9feq(mesh(1)%rho(1,j),mesh(1)%u(1,j),mesh(1)%v(1,j),1,j,k)+(1.0d0-1.0d0/mesh(1)%tau_f)*(mesh(1)%f0(2,j,k)-D2Q9feq(mesh(1)%rho(2,j),mesh(1)%u(2,j),mesh(1)%v(2,j),2,j,k))
    !					    enddo
                    enddo
                endif
            endif

            if(mesh(1)%right(2,1).eq.0)then
                if(periodical_bc.eq.1)then
                    !$acc kernels loop independent present(mesh)
                    do k=1,Q
                        do j=2,mesh(1)%Ny-1
                            mesh(1)%f(mesh(1)%Nx,j,k)=mesh(1)%f(2,j,k)
                        enddo
                    enddo
                    !$acc end kernels loop
                else
                    do j=2,mesh(1)%Ny-1
                        mesh(1)%rho(mesh(1)%Nx,j)=mesh(1)%rho(mesh(1)%Nx-1,j) !(4.0d0*mesh(1)%rho(mesh(1)%Nx-1,j)-mesh(1)%rho(mesh(1)%Nx-2,j))/3.0d0
                        mesh(1)%u(mesh(1)%Nx,j)= mesh(1)%u(mesh(1)%Nx-1,j) !0.0d0 !mesh(1)%u(mesh(1)%Nx-1,j)!(4.0d0*mesh(1)%u(mesh(1)%Nx-1,j)-mesh(1)%u(mesh(1)%Nx-2,j))/3.0d0
                        mesh(1)%v(mesh(1)%Nx,j)=mesh(1)%v(mesh(1)%Nx-1,j)!(4.0d0*mesh(1)%v(mesh(1)%Nx-1,j)-mesh(1)%v(mesh(1)%Nx-2,j))/3.0d0
                        rh=mesh(1)%rho(mesh(1)%Nx,j)
                        jx=rh*mesh(1)%u(mesh(1)%Nx,j)
                        jy=rh*mesh(1)%v(mesh(1)%Nx,j)
                        temp1(1)=rh
                        temp1(2)=-2.0d0*rh+3.0d0*(jx*jx+jy*jy)/rh
                        temp1(3)=rh-3.0d0*(jx*jx+jy*jy)/rh
                        temp1(4)=jx
                        temp1(5)=-jx
                        temp1(6)=jy
                        temp1(7)=-jy
                        temp1(8)=(jx*jx-jy*jy)/rh
                        temp1(9)=(jx*jy)/rh
                        rh=mesh(1)%rho(mesh(1)%Nx-1,j)
                        jx=rh*mesh(1)%u(mesh(1)%Nx-1,j)
                        jy=rh*mesh(1)%v(mesh(1)%Nx-1,j)
                        temp2(1)=rh
                        temp2(2)=-2.0d0*rh+3.0d0*(jx*jx+jy*jy)/rh
                        temp2(3)=rh-3.0d0*(jx*jx+jy*jy)/rh
                        temp2(4)=jx
                        temp2(5)=-jx
                        temp2(6)=jy
                        temp2(7)=-jy
                        temp2(8)=(jx*jx-jy*jy)/rh
                        temp2(9)=(jx*jy)/rh
                        temp=mesh(1)%f0(mesh(1)%Nx-1,j,:)
                        call matmulvoc(temp3,matm,temp,Q,Q)
                        do k=1,Q
                            temp3(k)=temp1(k)+(1.0d0-mesh(1)%s(k))*(temp3(k)-temp2(k))
                        enddo
                        call matmulvoc(temp,inversematm,temp3,Q,Q)
                        mesh(1)%f(mesh(1)%Nx,j,:)=temp
    !					    do k=1,Q
    !						    mesh(1)%rho(mesh(1)%Nx,j)=1.0d0
    !						    mesh(1)%u(mesh(1)%Nx,j)=UU !(4.0d0*mesh(1)%u(mesh(1)%Nx-1,j)-mesh(1)%u(mesh(1)%Nx-2,j))/3.0d0
    !						    mesh(1)%v(mesh(1)%Nx,j)=0.0d0 !(4.0d0*mesh(1)%v(mesh(1)%Nx-1,j)-mesh(1)%v(mesh(1)%Nx-2,j))/3.0d0
    !						    mesh(1)%f(mesh(1)%Nx,j,k)=D2Q9feq(mesh(1)%rho(mesh(1)%Nx,j),mesh(1)%u(mesh(1)%Nx,j),mesh(1)%v(mesh(1)%Nx,j),mesh(1)%Nx,j,k)+(1.0d0-1.0d0/mesh(1)%tau_f)*(mesh(1)%f0(mesh(1)%Nx-1,j,k)-D2Q9feq(mesh(1)%rho(mesh(1)%Nx-1,j),mesh(1)%u(mesh(1)%Nx-1,j),mesh(1)%v(mesh(1)%Nx-1,j),mesh(1)%Nx-1,j,k))
    !					    enddo
                    enddo
                endif
            endif
        elseif(ref_freestream.eq.2) then 
            if(mesh(1)%left(2,1).eq.0) then
                do j=2,mesh(1)%Ny-1
                    do k=1,Q
                        mesh(1)%rho(1,j)=1.0d0
                        mesh(1)%u(1,j)=0.0d0
                        mesh(1)%v(1,j)=mesh(1)%v(2,j)
                        mesh(1)%f(1,j,k)=D2Q9feq(mesh(1)%rho(1,j),mesh(1)%u(1,j),mesh(1)%v(1,j),1,j,k)+mesh(1)%f(2,j,k)-D2Q9feq(mesh(1)%rho(2,j),mesh(1)%u(2,j),mesh(1)%v(2,j),2,j,k)
                    enddo
                enddo
            endif

            if(mesh(1)%right(2,1).eq.0)then
                do j=2,mesh(1)%Ny-1
                    do k=1,Q
                        mesh(1)%rho(mesh(1)%Nx,j)=1.0d0
                        mesh(1)%u(mesh(1)%Nx,j)=0.0d0
                        mesh(1)%v(mesh(1)%Nx,j)=mesh(1)%v(mesh(1)%Nx-1,j)
                        mesh(1)%f(mesh(1)%Nx,j,k)=D2Q9feq(mesh(1)%rho(mesh(1)%Nx,j),mesh(1)%u(mesh(1)%Nx,j),mesh(1)%v(mesh(1)%Nx,j),mesh(1)%Nx,j,k)+mesh(1)%f(mesh(1)%Nx-1,j,k)-D2Q9feq(mesh(1)%rho(mesh(1)%Nx-1,j),mesh(1)%u(mesh(1)%Nx-1,j),mesh(1)%v(mesh(1)%Nx-1,j),mesh(1)%Nx-1,j,k)
                    enddo
                enddo
            endif
        endif
            
        if(mesh(1)%up(1,1).eq.0) then
            !$acc kernels loop independent present(mesh) create(rh,jx,jy,temp1,temp2,temp3,temp)
            do j=1,mesh(1)%Nx
                mesh(1)%rho(j,mesh(1)%Ny)=mesh(1)%rho(j,mesh(1)%Ny-1) !1.0d0
                mesh(1)%u(j,mesh(1)%Ny)=0.0d0 !mesh(1)%u(j,mesh(1)%Ny-1)!(4.0d0*mesh(1)%u(j,mesh(1)%Ny-1)-mesh(1)%u(j,mesh(1)%Ny-2))/3.0d0
                mesh(1)%v(j,mesh(1)%Ny)=0.0d0 !mesh(1)%v(j,mesh(1)%Ny-1) !0.0d0 !(4.0d0*mesh(1)%v(j,mesh(1)%Ny-1)-mesh(1)%v(j,mesh(1)%Ny-2))/3.0d0
                if((j.eq.1).and.(periodical_bc.eq.1))then
                    if(i.eq.1)then
                        mesh(1)%rho(j,mesh(1)%Ny)=mesh(1)%rho(mesh(1)%Nx-1,mesh(1)%Ny-1) !1.0d0
                        mesh(1)%u(j,mesh(1)%Ny)=mesh(1)%u(mesh(1)%Nx-1,mesh(1)%Ny-1) !mesh(1)%u(j,mesh(1)%Ny-1)!(4.0d0*mesh(1)%u(j,mesh(1)%Ny-1)-mesh(1)%u(j,mesh(1)%Ny-2))/3.0d0
                        mesh(1)%v(j,mesh(1)%Ny)=mesh(1)%v(mesh(1)%Nx-1,mesh(1)%Ny-1) ! mesh(1)%v(j,mesh(1)%Ny-1) !0.0d0 !(4.0d0*mesh(1)%v(j,mesh(1)%Ny-1)-mesh(1)%v(j,mesh(1)%Ny-2))/3.0d0    
                        mesh(1)%rho(j,mesh(1)%Ny-1)=mesh(1)%rho(mesh(1)%Nx-1,mesh(1)%Ny-1) !1.0d0
                        mesh(1)%u(j,mesh(1)%Ny-1)=mesh(1)%u(mesh(1)%Nx-1,mesh(1)%Ny-1) !mesh(1)%u(j,mesh(1)%Ny-1)!(4.0d0*mesh(1)%u(j,mesh(1)%Ny-1)-mesh(1)%u(j,mesh(1)%Ny-2))/3.0d0
                        mesh(1)%v(j,mesh(1)%Ny-1)=mesh(1)%v(mesh(1)%Nx-1,mesh(1)%Ny-1) ! mesh(1)%v(j,mesh(1)%Ny-1) !0.0d0 !(4.0d0*mesh(1)%v(j,mesh(1)%Ny-1)-mesh(1)%v(j,mesh(1)%Ny-2))/3.0d0    
                    endif
                elseif((j.eq.mesh(1)%Nx).and.(periodical_bc.eq.1))then
                    if(i.eq.3)then
                        mesh(1)%rho(j,mesh(1)%Ny)=mesh(1)%rho(2,mesh(1)%Ny-1) !1.0d0
                        mesh(1)%u(j,mesh(1)%Ny)=mesh(1)%u(2,mesh(1)%Ny-1) !mesh(1)%u(j,mesh(1)%Ny-1)!(4.0d0*mesh(1)%u(j,mesh(1)%Ny-1)-mesh(1)%u(j,mesh(1)%Ny-2))/3.0d0
                        mesh(1)%v(j,mesh(1)%Ny)=mesh(1)%v(2,mesh(1)%Ny-1) ! mesh(1)%v(j,mesh(1)%Ny-1) !0.0d0 !(4.0d0*mesh(1)%v(j,mesh(1)%Ny-1)-mesh(1)%v(j,mesh(1)%Ny-2))/3.0d0    
                        mesh(1)%rho(j,mesh(1)%Ny-1)=mesh(1)%rho(2,mesh(1)%Ny-1) !1.0d0
                        mesh(1)%u(j,mesh(1)%Ny-1)=mesh(1)%u(2,mesh(1)%Ny-1) !mesh(1)%u(j,mesh(1)%Ny-1)!(4.0d0*mesh(1)%u(j,mesh(1)%Ny-1)-mesh(1)%u(j,mesh(1)%Ny-2))/3.0d0
                        mesh(1)%v(j,mesh(1)%Ny-1)=mesh(1)%v(2,mesh(1)%Ny-1) ! mesh(1)%v(j,mesh(1)%Ny-1) !0.0d0 !(4.0d0*mesh(1)%v(j,mesh(1)%Ny-1)-mesh(1)%v(j,mesh(1)%Ny-2))/3.0d0    
                    endif
                endif
                rh=mesh(1)%rho(j,mesh(1)%Ny)
                jx=rh*mesh(1)%u(j,mesh(1)%Ny)
                jy=rh*mesh(1)%v(j,mesh(1)%Ny)
                temp1(1)=rh
                temp1(2)=-2.0d0*rh+3.0d0*(jx*jx+jy*jy)/rh
                temp1(3)=rh-3.0d0*(jx*jx+jy*jy)/rh
                temp1(4)=jx
                temp1(5)=-jx
                temp1(6)=jy
                temp1(7)=-jy
                temp1(8)=(jx*jx-jy*jy)/rh
                temp1(9)=(jx*jy)/rh
                rh=mesh(1)%rho(j,mesh(1)%Ny-1)
                jx=rh*mesh(1)%u(j,mesh(1)%Ny-1)
                jy=rh*mesh(1)%v(j,mesh(1)%Ny-1)
                temp2(1)=rh
                temp2(2)=-2.0d0*rh+3.0d0*(jx*jx+jy*jy)/rh
                temp2(3)=rh-3.0d0*(jx*jx+jy*jy)/rh
                temp2(4)=jx
                temp2(5)=-jx
                temp2(6)=jy
                temp2(7)=-jy
                temp2(8)=(jx*jx-jy*jy)/rh
                temp2(9)=(jx*jy)/rh
                temp=mesh(1)%f0(j,mesh(1)%Ny-1,:)
                if((j.eq.1).and.(periodical_bc.eq.1))then
                    if(i.eq.1) temp= mesh(1)%f0(mesh(1)%Nx-1,mesh(1)%Ny-1,:)
                elseif((j.eq.mesh(1)%Nx).and.(periodical_bc.eq.1))then
                    if(i.eq.3)  temp= mesh(1)%f0(2,mesh(1)%Ny-1,:)
                endif
                call matmulvoc(temp3,matm,temp,Q,Q)
                do k=1,Q
                    temp3(k)=temp1(k)+(1.0d0-mesh(1)%s(k))*(temp3(k)-temp2(k))
                enddo
                call matmulvoc(temp,inversematm,temp3,Q,Q)
                mesh(1)%f(j,mesh(1)%Ny,:)=temp
!					do k=1,Q
!                        mesh(1)%rho(j,mesh(1)%Ny)=1.0d0
!                        mesh(1)%u(j,mesh(1)%Ny)=UU !(4.0d0*mesh(1)%u(j,mesh(1)%Ny-1)-mesh(1)%u(j,mesh(1)%Ny-2))/3.0d0
!                        mesh(1)%v(j,mesh(1)%Ny)=0.0d0 !(4.0d0*mesh(1)%v(j,mesh(1)%Ny-1)-mesh(1)%v(j,mesh(1)%Ny-2))/3.0d0
!                        mesh(1)%f(j,mesh(1)%Ny,k)=D2Q9feq(mesh(1)%rho(j,mesh(1)%Ny),mesh(1)%u(j,mesh(1)%Ny),mesh(1)%v(j,mesh(1)%Ny),j,mesh(1)%Ny,K)+(1.0d0-1.0d0/mesh(1)%tau_f)*(mesh(1)%f0(j,mesh(1)%Ny-1,k)-D2Q9feq(mesh(1)%rho(j,mesh(1)%Ny-1),mesh(1)%u(j,mesh(1)%Ny-1),mesh(1)%v(j,mesh(1)%Ny-1),j,mesh(1)%Ny-1,k))
!					enddo
            enddo
            !$acc end kernels loop
        endif

        if(mesh(1)%down(1,1).eq.0) then
            !$acc kernels loop independent present(mesh) create(rh,jx,jy,temp1,temp2,temp3,temp)
            do j=1,mesh(1)%Nx
                mesh(1)%rho(j,1)=mesh(1)%rho(j,2) !1.0d0
                mesh(1)%u(j,1)= 0.0d0 !mesh(1)%u(j,2) !(4.0d0*mesh(1)%u(j,2)-mesh(1)%u(j,3))/3.0d0
                mesh(1)%v(j,1)=0.0d0 !mesh(1)%v(j,2) !0.0d0 !(4.0d0*mesh(1)%v(j,2)-mesh(1)%v(j,3))/3.0d0
                if((j.eq.1).and.(periodical_bc.eq.1))then
                    if(i.eq.1)then
                        mesh(1)%rho(j,1)=mesh(1)%rho(mesh(1)%Nx-1,2) !1.0d0
                        mesh(1)%u(j,1)=mesh(1)%u(mesh(1)%Nx-1,2) !mesh(1)%u(j,mesh(1)%Ny-1)!(4.0d0*mesh(1)%u(j,mesh(1)%Ny-1)-mesh(1)%u(j,mesh(1)%Ny-2))/3.0d0
                        mesh(1)%v(j,1)=mesh(1)%v(mesh(1)%Nx-1,2) ! mesh(1)%v(j,mesh(1)%Ny-1) !0.0d0 !(4.0d0*mesh(1)%v(j,mesh(1)%Ny-1)-mesh(1)%v(j,mesh(1)%Ny-2))/3.0d0   
                        mesh(1)%rho(j,2)=mesh(1)%rho(mesh(1)%Nx-1,2) !1.0d0
                        mesh(1)%u(j,2)=mesh(1)%u(mesh(1)%Nx-1,2) !mesh(1)%u(j,mesh(1)%Ny-1)!(4.0d0*mesh(1)%u(j,mesh(1)%Ny-1)-mesh(1)%u(j,mesh(1)%Ny-2))/3.0d0
                        mesh(1)%v(j,2)=mesh(1)%v(mesh(1)%Nx-1,2) ! mesh(1)%v(j,mesh(1)%Ny-1) !0.0d0 !(4.0d0*mesh(1)%v(j,mesh(1)%Ny-1)-mesh(1)%v(j,mesh(1)%Ny-2))/3.0d0 
                    endif
                elseif((j.eq.mesh(1)%Nx).and.(periodical_bc.eq.1))then
                    if(i.eq.3)then
                        mesh(1)%rho(j,1)=mesh(1)%rho(2,2) !1.0d0
                        mesh(1)%u(j,1)=mesh(1)%u(2,2) !mesh(1)%u(j,mesh(1)%Ny-1)!(4.0d0*mesh(1)%u(j,mesh(1)%Ny-1)-mesh(1)%u(j,mesh(1)%Ny-2))/3.0d0
                        mesh(1)%v(j,1)=mesh(1)%v(2,2) ! mesh(1)%v(j,mesh(1)%Ny-1) !0.0d0 !(4.0d0*mesh(1)%v(j,mesh(1)%Ny-1)-mesh(1)%v(j,mesh(1)%Ny-2))/3.0d0  
                        mesh(1)%rho(j,2)=mesh(1)%rho(2,2) !1.0d0
                        mesh(1)%u(j,2)=mesh(1)%u(2,2) !mesh(1)%u(j,mesh(1)%Ny-1)!(4.0d0*mesh(1)%u(j,mesh(1)%Ny-1)-mesh(1)%u(j,mesh(1)%Ny-2))/3.0d0
                        mesh(1)%v(j,2)=mesh(1)%v(2,2) ! mesh(1)%v(j,mesh(1)%Ny-1) !0.0d0 !(4.0d0*mesh(1)%v(j,mesh(1)%Ny-1)-mesh(1)%v(j,mesh(1)%Ny-2))/3.0d0  
                    endif
                endif
                rh=mesh(1)%rho(j,1)
                jx=rh*mesh(1)%u(j,1)
                jy=rh*mesh(1)%v(j,1)
                temp1(1)=rh
                temp1(2)=-2.0d0*rh+3.0d0*(jx*jx+jy*jy)/rh
                temp1(3)=rh-3.0d0*(jx*jx+jy*jy)/rh
                temp1(4)=jx
                temp1(5)=-jx
                temp1(6)=jy
                temp1(7)=-jy
                temp1(8)=(jx*jx-jy*jy)/rh
                temp1(9)=(jx*jy)/rh
                rh=mesh(1)%rho(j,2)
                jx=rh*mesh(1)%u(j,2)
                jy=rh*mesh(1)%v(j,2)
                temp2(1)=rh
                temp2(2)=-2.0d0*rh+3.0d0*(jx*jx+jy*jy)/rh
                temp2(3)=rh-3.0d0*(jx*jx+jy*jy)/rh
                temp2(4)=jx
                temp2(5)=-jx
                temp2(6)=jy
                temp2(7)=-jy
                temp2(8)=(jx*jx-jy*jy)/rh
                temp2(9)=(jx*jy)/rh
                temp=mesh(1)%f0(j,2,:)
                if((j.eq.1).and.(periodical_bc.eq.1))then
                    if(i.eq.1) temp= mesh(1)%f0(mesh(1)%Nx-1,2,:)
                elseif((j.eq.mesh(1)%Nx).and.(periodical_bc.eq.1))then
                    if(i.eq.3)  temp= mesh(1)%f0(2,2,:)
                endif
                call matmulvoc(temp3,matm,temp,Q,Q)
                do k=1,Q
                    temp3(k)=temp1(k)+(1.0d0-mesh(1)%s(k))*(temp3(k)-temp2(k))
                enddo
                call matmulvoc(temp,inversematm,temp3,Q,Q)
                mesh(1)%f(j,1,:)=temp
            enddo
            !$acc end kernels loop
        endif
    endif
end


subroutine boundary_condition1 (mesh,shuzu,n) !远场边界条件,特征线边界条件
	use com_date
    use meshtype
	implicit none
	
	type(mesh_dat) mesh(max_meshid)
	integer n,i,j,k,i0,j0,k0
    integer n1,n2,n3
	integer shuzu(n)
    real(kind=ps) rh,jx,jy,temp(Q),temp1(Q),temp2(Q),temp3(Q),lx(3),ly(3),cs,temprho(mesh(1)%Ny),tempu(mesh(1)%Ny),tempv(mesh(1)%Ny),temprho1(mesh(1)%Nx),tempu1(mesh(1)%Nx),tempv1(mesh(1)%Nx),tempconr(4,3)
	Real (kind=ps), external :: D2Q9feq,distanceq
    real(kind=ps) rhox,rhoy,ux,uy,vx,vy

	!$acc enter data create(cs,lx,ly,rhox,rhoy,ux,uy,vx,vy,tempconr,temp1,temp,temp3,temp2,temprho1,tempu1,tempv1,temprho,tempu,tempv,rh,jx,jy,i,i0,j0,k0)
    if(myid.eq.0) then
        !$acc kernels present(cs,mesh,lx,ly,rhox,rhoy,ux,uy,vx,vy,tempconr) 
        cs=sqrt(3.0)/3.0d0      
        !左上角点处理
        rhox=(-3.0d0*mesh(1)%rho(1,mesh(1)%Ny)+4.0d0*mesh(1)%rho(2,mesh(1)%Ny)-mesh(1)%rho(3,mesh(1)%Ny))*0.5d0
        rhoy=(3.0d0*mesh(1)%rho(1,mesh(1)%Ny)-4.0d0*mesh(1)%rho(1,mesh(1)%Ny-1)+mesh(1)%rho(1,mesh(1)%Ny-2))*0.5d0  
        ux=(-3.0d0*mesh(1)%u(1,mesh(1)%Ny)+4.0d0*mesh(1)%u(2,mesh(1)%Ny)-mesh(1)%u(3,mesh(1)%Ny))*0.5d0 
        uy=(3.0d0*mesh(1)%u(1,mesh(1)%Ny)-4.0d0*mesh(1)%u(1,mesh(1)%Ny-1)+mesh(1)%u(1,mesh(1)%Ny-2))*0.5d0  
        vx=(-3.0d0*mesh(1)%v(1,mesh(1)%Ny)+4.0d0*mesh(1)%v(2,mesh(1)%Ny)-mesh(1)%v(3,mesh(1)%Ny))*0.5d0
        vy=(3.0d0*mesh(1)%v(1,mesh(1)%Ny)-4.0d0*mesh(1)%v(1,mesh(1)%Ny-1)+mesh(1)%v(1,mesh(1)%Ny-2))*0.5d0 
        lx(1)=(mesh(1)%u(1,mesh(1)%Ny)-cs)*(cs*cs*rhox-cs*mesh(1)%rho(1,mesh(1)%Ny)*ux)
        if(mesh(1)%u(1,mesh(1)%Ny).lt.0.0d0)then
            lx(2)=mesh(1)%u(1,mesh(1)%Ny)*vx
        else
            lx(2)=0.0d0
        endif
        lx(3)=0.0d0
        ly(1)=0.0d0 !(mesh(1)%u(1,j)-cs)*(cs*cs*rhox-cs*mesh(1)%rho(1,j)*ux)
        if(mesh(1)%v(1,mesh(1)%Ny).ge.0.0d0)then
            ly(2)=mesh(1)%v(1,mesh(1)%Ny)*uy
        else
            ly(2)=0.0d0
        endif
        ly(3)=(mesh(1)%v(1,mesh(1)%Ny)+cs)*(cs*cs*rhoy+cs*mesh(1)%rho(1,mesh(1)%Ny)*vy)
        tempconr(1,1)=mesh(1)%rho(1,mesh(1)%Ny)-(ly(1)+ly(3))*1.5d0-(lx(1)+lx(3))*1.5d0
        tempconr(1,2)=mesh(1)%u(1,mesh(1)%Ny)-ly(2)+(lx(1)-lx(3))/(2.0d0*cs*mesh(1)%rho(1,mesh(1)%Ny))
        tempconr(1,3)=mesh(1)%v(1,mesh(1)%Ny)+(ly(1)-ly(3))/(2.0d0*cs*mesh(1)%rho(1,mesh(1)%Ny))-lx(2)        
        !左下角点处理
        rhox=(-3.0d0*mesh(1)%rho(1,1)+4.0d0*mesh(1)%rho(2,1)-mesh(1)%rho(3,1))*0.5d0
        rhoy=(-3.0d0*mesh(1)%rho(1,1)+4.0d0*mesh(1)%rho(1,2)-mesh(1)%rho(1,3))*0.5d0  
        ux=(-3.0d0*mesh(1)%u(1,1)+4.0d0*mesh(1)%u(2,1)-mesh(1)%u(3,1))*0.5d0 
        uy=(-3.0d0*mesh(1)%u(1,1)+4.0d0*mesh(1)%u(1,2)-mesh(1)%u(1,3))*0.5d0  
        vx=(-3.0d0*mesh(1)%v(1,1)+4.0d0*mesh(1)%v(2,1)-mesh(1)%v(3,1))*0.5d0
        vy=(-3.0d0*mesh(1)%v(1,1)+4.0d0*mesh(1)%v(1,2)-mesh(1)%v(1,3))*0.5d0 
        lx(1)=(mesh(1)%u(1,1)-cs)*(cs*cs*rhox-cs*mesh(1)%rho(1,1)*ux)
        if(mesh(1)%u(1,1).lt.0.0d0)then
            lx(2)=mesh(1)%u(1,1)*vx
        else
            lx(2)=0.0d0
        endif
        lx(3)=0.0d0
        ly(1)=(mesh(1)%v(1,1)-cs)*(cs*cs*rhoy-cs*mesh(1)%rho(1,1)*vy)
        if(mesh(1)%v(1,1).lt.0.0d0)then
            ly(2)=mesh(1)%v(1,1)*uy
        else
            ly(2)=0.0d0
        endif
        ly(3)=0.0d0 !(mesh(1)%v(j,1)+cs)*(cs*cs*rhoy+cs*mesh(1)%rho(j,1)*vy)
        tempconr(3,1)=mesh(1)%rho(1,1)-(ly(1)+ly(3))*1.5d0-(lx(1)+lx(3))*1.5d0
        tempconr(3,2)=mesh(1)%u(1,1)-ly(2)+(lx(1)-lx(3))/(2.0d0*cs*mesh(1)%rho(1,1))
        tempconr(3,3)=mesh(1)%v(1,1)+(ly(1)-ly(3))/(2.0d0*cs*mesh(1)%rho(1,1))-lx(2)       
        !右上角点处理
        rhox=(3.0d0*mesh(1)%rho(mesh(1)%Nx,mesh(1)%Ny)-4.0d0*mesh(1)%rho(mesh(1)%Nx-1,mesh(1)%Ny)+mesh(1)%rho(mesh(1)%Nx-2,mesh(1)%Ny))*0.5d0
        rhoy=(3.0d0*mesh(1)%rho(mesh(1)%Nx,mesh(1)%Ny)-4.0d0*mesh(1)%rho(mesh(1)%Nx,mesh(1)%Ny-1)+mesh(1)%rho(mesh(1)%Nx,mesh(1)%Ny-2))*0.5d0  
        ux=(3.0d0*mesh(1)%u(mesh(1)%Nx,mesh(1)%Ny)-4.0d0*mesh(1)%u(mesh(1)%Nx-1,mesh(1)%Ny)+mesh(1)%u(mesh(1)%Nx-2,mesh(1)%Ny))*0.5d0 
        uy=(3.0d0*mesh(1)%u(mesh(1)%Nx,mesh(1)%Ny)-4.0d0*mesh(1)%u(mesh(1)%Nx,mesh(1)%Ny-1)+mesh(1)%u(mesh(1)%Nx,mesh(1)%Ny-2))*0.5d0  
        vx=(3.0d0*mesh(1)%v(mesh(1)%Nx,mesh(1)%Ny)-4.0d0*mesh(1)%v(mesh(1)%Nx-1,mesh(1)%Ny)+mesh(1)%v(mesh(1)%Nx-2,mesh(1)%Ny))*0.5d0
        vy=(3.0d0*mesh(1)%v(mesh(1)%Nx,mesh(1)%Ny)-4.0d0*mesh(1)%v(mesh(1)%Nx,mesh(1)%Ny-1)+mesh(1)%v(mesh(1)%Nx,mesh(1)%Ny-2))*0.5d0 
        lx(1)=0.0d0 !(mesh(1)%u(1,j)-cs)*(cs*cs*rhox-cs*mesh(1)%rho(1,j)*ux)
        if(mesh(1)%u(mesh(1)%Nx,mesh(1)%Ny).ge.0.0d0)then
            lx(2)=mesh(1)%u(mesh(1)%Nx,mesh(1)%Ny)*vx
        else
            lx(2)=0.0d0
        endif
        lx(3)=(mesh(1)%u(mesh(1)%Nx,mesh(1)%Ny)+cs)*(cs*cs*rhox+cs*mesh(1)%rho(mesh(1)%Nx,mesh(1)%Ny)*ux)
        ly(1)=0.0d0 !(mesh(1)%u(1,j)-cs)*(cs*cs*rhox-cs*mesh(1)%rho(1,j)*ux)
        if(mesh(1)%v(mesh(1)%Nx,mesh(1)%Ny).ge.0.0d0)then
            ly(2)=mesh(1)%v(mesh(1)%Nx,mesh(1)%Ny)*uy
        else
            ly(2)=0.0d0
        endif
        ly(3)=(mesh(1)%v(mesh(1)%Nx,mesh(1)%Ny)+cs)*(cs*cs*rhoy+cs*mesh(1)%rho(mesh(1)%Nx,mesh(1)%Ny)*vy)
        tempconr(2,1)=mesh(1)%rho(mesh(1)%Nx,mesh(1)%Ny)-(ly(1)+ly(3))*1.5d0-(lx(1)+lx(3))*1.5d0
        tempconr(2,2)=mesh(1)%u(mesh(1)%Nx,mesh(1)%Ny)-ly(2)+(lx(1)-lx(3))/(2.0d0*cs*mesh(1)%rho(mesh(1)%Nx,mesh(1)%Ny))
        tempconr(2,3)=mesh(1)%v(mesh(1)%Nx,mesh(1)%Ny)+(ly(1)-ly(3))/(2.0d0*cs*mesh(1)%rho(mesh(1)%Nx,mesh(1)%Ny))-lx(2)     
        !右下角点处理
        rhox=(3.0d0*mesh(1)%rho(mesh(1)%Nx,1)-4.0d0*mesh(1)%rho(mesh(1)%Nx-1,1)+mesh(1)%rho(mesh(1)%Nx-2,1))*0.5d0
        rhoy=(-3.0d0*mesh(1)%rho(mesh(1)%Nx,1)+4.0d0*mesh(1)%rho(mesh(1)%Nx,2)-mesh(1)%rho(mesh(1)%Nx,3))*0.5d0  
        ux=(3.0d0*mesh(1)%u(mesh(1)%Nx,1)-4.0d0*mesh(1)%u(mesh(1)%Nx-1,1)+mesh(1)%u(mesh(1)%Nx-2,1))*0.5d0 
        uy=(-3.0d0*mesh(1)%u(mesh(1)%Nx,1)+4.0d0*mesh(1)%u(mesh(1)%Nx,2)-mesh(1)%u(mesh(1)%Nx,3))*0.5d0  
        vx=(3.0d0*mesh(1)%v(mesh(1)%Nx,1)-4.0d0*mesh(1)%v(mesh(1)%Nx-1,1)+mesh(1)%v(mesh(1)%Nx-2,1))*0.5d0
        vy=(-3.0d0*mesh(1)%v(mesh(1)%Nx,1)+4.0d0*mesh(1)%v(mesh(1)%Nx,2)-mesh(1)%v(mesh(1)%Nx,3))*0.5d0 
        lx(1)=0.0d0 !(mesh(1)%u(1,j)-cs)*(cs*cs*rhox-cs*mesh(1)%rho(1,j)*ux)
        if(mesh(1)%u(mesh(1)%Nx,1).ge.0.0d0)then
            lx(2)=mesh(1)%u(mesh(1)%Nx,1)*vx
        else
            lx(2)=0.0d0
        endif
        lx(3)=(mesh(1)%u(mesh(1)%Nx,1)+cs)*(cs*cs*rhox+cs*mesh(1)%rho(mesh(1)%Nx,1)*ux)
        ly(1)=(mesh(1)%v(mesh(1)%Nx,1)-cs)*(cs*cs*rhoy-cs*mesh(1)%rho(mesh(1)%Nx,1)*vy)
        if(mesh(1)%v(mesh(1)%Nx,1).lt.0.0d0)then
            ly(2)=mesh(1)%v(mesh(1)%Nx,1)*uy
        else
            ly(2)=0.0d0
        endif
        ly(3)=0.0d0 !(mesh(1)%v(j,1)+cs)*(cs*cs*rhoy+cs*mesh(1)%rho(j,1)*vy)
        tempconr(4,1)=mesh(1)%rho(mesh(1)%Nx,1)-(ly(1)+ly(3))*1.5d0-(lx(1)+lx(3))*1.5d0
        tempconr(4,2)=mesh(1)%u(mesh(1)%Nx,1)-ly(2)+(lx(1)-lx(3))/(2.0d0*cs*mesh(1)%rho(mesh(1)%Nx,1))
        tempconr(4,3)=mesh(1)%v(mesh(1)%Nx,1)+(ly(1)-ly(3))/(2.0d0*cs*mesh(1)%rho(mesh(1)%Nx,1))-lx(2) 
        !$acc end kernels 
        
		do i=shuzu(1),shuzu(n)
            if(ref_freestream.le.1)then
			    if(mesh(1)%left(2,1).eq.0) then
                    if(periodical_bc.eq.1)then
                        !$acc kernels loop independent present(mesh,i)
                        do j=2,mesh(1)%Ny-1
                            mesh(1)%f(1,j,:)=mesh(1)%f(mesh(1)%Nx-1,j,:)
                        enddo
                        !$acc end kernels loop
                    else
                        !$acc kernels present(mesh,cs,lx,ly,rhox,rhoy,ux,uy,vx,vy,temprho1,tempu1,tempv1,temprho,tempu,tempv,temp,temp1,temp2,temp3,rh,jx,jy,i)
                        temprho=mesh(1)%rho(1,:)
                        tempu=mesh(1)%u(1,:)
                        tempv=mesh(1)%v(1,:)
                        !$acc loop independent
                        do j=2,mesh(1)%Ny-1
                            rhox=(-3.0d0*temprho(j)+4.0d0*mesh(1)%rho(2,j)-mesh(1)%rho(3,j))*0.5d0
                            rhoy=(temprho(j+1)-temprho(j-1))*0.5d0
                            ux=(-3.0d0*tempu(j)+4.0d0*mesh(1)%u(2,j)-mesh(1)%u(3,j))*0.5d0
                            uy=(tempu(j+1)-tempu(j-1))*0.5d0
                            vx=(-3.0d0*tempv(j)+4.0d0*mesh(1)%v(2,j)-mesh(1)%v(3,j))*0.5d0
                            vy=(tempv(j+1)-tempv(j-1))*0.5d0
                            lx(1)=(tempu(j)-cs)*(cs*cs*rhox-cs*temprho(j)*ux)
                            if(tempu(j).lt.0.0d0)then
                                lx(2)=tempu(j)*vx
                            else
                                lx(2)=0.0d0
                            endif
                            lx(3)=0.0d0
                            ly(1)=(tempv(j)-cs)*(cs*cs*rhoy-cs*temprho(j)*vy)
                            ly(2)=tempv(j)*uy
                            ly(3)=(tempv(j)+cs)*(cs*cs*rhoy+cs*temprho(j)*vy)
                            
                            mesh(1)%rho(1,j)=mesh(1)%rho(1,j)-(lx(1)+lx(3))*1.5d0 !-0.75d0*(mesh(1)%v(1,j)*rhoy+mesh(1)%rho(1,j)*vy) !-0.75d0*(ly(1)+ly(3))*1.5d0-(lx(1)+lx(3))*1.5d0
                            mesh(1)%u(1,j)=mesh(1)%u(1,j)+(lx(1)-lx(3))/(2.0d0*cs*mesh(1)%rho(1,j)) !-0.75d0*mesh(1)%v(1,j)*uy !-0.75d0*ly(2)+(lx(1)-lx(3))/(2.0d0*cs*mesh(1)%rho(1,j))
                            mesh(1)%v(1,j)=mesh(1)%v(1,j)-lx(2) !-0.75d0*(cs*cs*rhoy/mesh(1)%rho(1,j)+mesh(1)%v(1,j)*vy) !+0.75d0*(ly(1)-ly(3))/(2.0d0*cs*mesh(1)%rho(1,j))-lx(2)
!                            if(j.eq.425)print*,mesh(1)%rho(1,j),0.75d0*(ly(1)+ly(3))*1.5d0,-(lx(1)+lx(3))*1.5d0
                        enddo
                        !$acc end loop
                        !$acc loop independent     
				        do j=2,mesh(1)%Ny-1
                            rh=mesh(1)%rho(1,j)
                            jx=rh*mesh(1)%u(1,j)
                            jy=rh*mesh(1)%v(1,j)
                            temp1(1)=rh
                            temp1(2)=-2.0d0*rh+3.0d0*(jx*jx+jy*jy)/rh
                            temp1(3)=rh-3.0d0*(jx*jx+jy*jy)/rh
                            temp1(4)=jx
                            temp1(5)=-jx
                            temp1(6)=jy
                            temp1(7)=-jy
                            temp1(8)=(jx*jx-jy*jy)/rh
                            temp1(9)=(jx*jy)/rh
                            rh=mesh(1)%rho(2,j)
                            jx=rh*mesh(1)%u(2,j)
                            jy=rh*mesh(1)%v(2,j)
                            temp2(1)=rh
                            temp2(2)=-2.0d0*rh+3.0d0*(jx*jx+jy*jy)/rh
                            temp2(3)=rh-3.0d0*(jx*jx+jy*jy)/rh
                            temp2(4)=jx
                            temp2(5)=-jx
                            temp2(6)=jy
                            temp2(7)=-jy
                            temp2(8)=(jx*jx-jy*jy)/rh
                            temp2(9)=(jx*jy)/rh
                            temp=mesh(1)%f0(2,j,:)
                            call matmulvoc(temp3,matm,temp,Q,Q)
                            do k=1,Q
                                temp3(k)=temp1(k)+(1.0d0-mesh(1)%s(k))*(temp3(k)-temp2(k))
                            enddo
                            call matmulvoc(temp,inversematm,temp3,Q,Q)
                            mesh(1)%f(1,j,:)=temp
    !					    do k=1,Q
    !						    mesh(1)%rho(1,j)=1.0d0
    !						    mesh(1)%u(1,j)=UU
    !						    mesh(1)%v(1,j)=0.0d0
    !						    mesh(1)%f(1,j,k)=D2Q9feq(mesh(1)%rho(1,j),mesh(1)%u(1,j),mesh(1)%v(1,j),1,j,k)+(1.0d0-1.0d0/mesh(1)%tau_f)*(mesh(1)%f0(2,j,k)-D2Q9feq(mesh(1)%rho(2,j),mesh(1)%u(2,j),mesh(1)%v(2,j),2,j,k))
    !					    enddo
                        enddo
                        !$acc end loop
                        !$acc end kernels
                    endif
			    endif
                
			    if(mesh(1)%right(2,1).eq.0)then
                    if(periodical_bc.eq.1)then 
                        !$acc kernels loop independent present(mesh,i)  
                        do j=2,mesh(1)%Ny-1
                            mesh(1)%f(mesh(1)%Nx,j,:)=mesh(1)%f(2,j,:)
                        enddo
                        !$acc end kernels loop
                    else
                        !$acc kernels present(mesh,cs,lx,ly,rhox,rhoy,ux,uy,vx,vy,temprho1,tempu1,tempv1,temprho,tempu,tempv,temp,temp1,temp2,temp3,rh,jx,jy,i)
                        temprho=mesh(1)%rho(mesh(1)%Nx,:)
                        tempu=mesh(1)%u(mesh(1)%Nx,:)
                        tempv=mesh(1)%v(mesh(1)%Nx,:)
                        !$acc loop independent   
                        do j=2,mesh(1)%Ny-1
                            rhox=(3.0d0*temprho(j)-4.0d0*mesh(1)%rho(mesh(1)%Nx-1,j)+mesh(1)%rho(mesh(1)%Nx-2,j))*0.5d0
                            rhoy=(temprho(j+1)-temprho(j-1))*0.5d0
                            ux=(3.0d0*tempu(j)-4.0d0*mesh(1)%u(mesh(1)%Nx-1,j)+mesh(1)%u(mesh(1)%Nx-2,j))*0.5d0
                            uy=(tempu(j+1)-tempu(j-1))*0.5d0
                            vx=(3.0d0*tempv(j)-4.0d0*mesh(1)%v(mesh(1)%Nx-1,j)+mesh(1)%v(mesh(1)%Nx-2,j))*0.5d0
                            vy=(tempv(j+1)-tempv(j-1))*0.5d0
                            lx(1)=0.0d0 !(mesh(1)%u(1,j)-cs)*(cs*cs*rhox-cs*mesh(1)%rho(1,j)*ux)
                            if(tempu(j).ge.0.0d0)then
                                lx(2)=tempu(j)*vx
                            else
                                lx(2)=0.0d0
                            endif
                            lx(3)=(tempu(j)+cs)*(cs*cs*rhox+cs*temprho(j)*ux)
                            ly(1)=(tempv(j)-cs)*(cs*cs*rhoy-cs*temprho(j)*vy)
                            ly(2)=tempv(j)*uy
                            ly(3)=(tempv(j)+cs)*(cs*cs*rhoy+cs*temprho(j)*vy)
                            
                            mesh(1)%rho(mesh(1)%Nx,j)=mesh(1)%rho(mesh(1)%Nx,j)-(lx(1)+lx(3))*1.5d0 !-0.75d0*(ly(1)+ly(3))*1.5d0-(lx(1)+lx(3))*1.5d0
                            mesh(1)%u(mesh(1)%Nx,j)=mesh(1)%u(mesh(1)%Nx,j)+(lx(1)-lx(3))/(2.0d0*cs*mesh(1)%rho(mesh(1)%Nx,j)) !-0.75d0*ly(2)+(lx(1)-lx(3))/(2.0d0*cs*mesh(1)%rho(mesh(1)%Nx,j))
                            mesh(1)%v(mesh(1)%Nx,j)=mesh(1)%v(mesh(1)%Nx,j)-lx(2) !+0.75d0*(ly(1)-ly(3))/(2.0d0*cs*mesh(1)%rho(mesh(1)%Nx,j))-lx(2)
                        enddo
                        !$acc end loop
                        !$acc loop independent   
				        do j=2,mesh(1)%Ny-1
!                            mesh(1)%rho(mesh(1)%Nx,j)=1.0d0 !mesh(1)%rho(mesh(1)%Nx-1,j) !(4.0d0*mesh(1)%rho(mesh(1)%Nx-1,j)-mesh(1)%rho(mesh(1)%Nx-2,j))/3.0d0
!                            mesh(1)%u(mesh(1)%Nx,j)=0.0d0 ! mesh(1)%u(mesh(1)%Nx-1,j) !0.0d0 !mesh(1)%u(mesh(1)%Nx-1,j)!(4.0d0*mesh(1)%u(mesh(1)%Nx-1,j)-mesh(1)%u(mesh(1)%Nx-2,j))/3.0d0
!                            mesh(1)%v(mesh(1)%Nx,j)=0.0d0 !mesh(1)%v(mesh(1)%Nx-1,j)!(4.0d0*mesh(1)%v(mesh(1)%Nx-1,j)-mesh(1)%v(mesh(1)%Nx-2,j))/3.0d0
                            rh=mesh(1)%rho(mesh(1)%Nx,j)
                            jx=rh*mesh(1)%u(mesh(1)%Nx,j)
                            jy=rh*mesh(1)%v(mesh(1)%Nx,j)
                            temp1(1)=rh
                            temp1(2)=-2.0d0*rh+3.0d0*(jx*jx+jy*jy)/rh
                            temp1(3)=rh-3.0d0*(jx*jx+jy*jy)/rh
                            temp1(4)=jx
                            temp1(5)=-jx
                            temp1(6)=jy
                            temp1(7)=-jy
                            temp1(8)=(jx*jx-jy*jy)/rh
                            temp1(9)=(jx*jy)/rh
                            rh=mesh(1)%rho(mesh(1)%Nx-1,j)
                            jx=rh*mesh(1)%u(mesh(1)%Nx-1,j)
                            jy=rh*mesh(1)%v(mesh(1)%Nx-1,j)
                            temp2(1)=rh
                            temp2(2)=-2.0d0*rh+3.0d0*(jx*jx+jy*jy)/rh
                            temp2(3)=rh-3.0d0*(jx*jx+jy*jy)/rh
                            temp2(4)=jx
                            temp2(5)=-jx
                            temp2(6)=jy
                            temp2(7)=-jy
                            temp2(8)=(jx*jx-jy*jy)/rh
                            temp2(9)=(jx*jy)/rh
                            temp=mesh(1)%f0(mesh(1)%Nx-1,j,:)
                            call matmulvoc(temp3,matm,temp,Q,Q)
                            do k=1,Q
                                temp3(k)=temp1(k)+(1.0d0-mesh(1)%s(k))*(temp3(k)-temp2(k))
                            enddo
                            call matmulvoc(temp,inversematm,temp3,Q,Q)
                            mesh(1)%f(mesh(1)%Nx,j,:)=temp
    !					    do k=1,Q
    !						    mesh(1)%rho(mesh(1)%Nx,j)=1.0d0
    !						    mesh(1)%u(mesh(1)%Nx,j)=UU !(4.0d0*mesh(1)%u(mesh(1)%Nx-1,j)-mesh(1)%u(mesh(1)%Nx-2,j))/3.0d0
    !						    mesh(1)%v(mesh(1)%Nx,j)=0.0d0 !(4.0d0*mesh(1)%v(mesh(1)%Nx-1,j)-mesh(1)%v(mesh(1)%Nx-2,j))/3.0d0
    !						    mesh(1)%f(mesh(1)%Nx,j,k)=D2Q9feq(mesh(1)%rho(mesh(1)%Nx,j),mesh(1)%u(mesh(1)%Nx,j),mesh(1)%v(mesh(1)%Nx,j),mesh(1)%Nx,j,k)+(1.0d0-1.0d0/mesh(1)%tau_f)*(mesh(1)%f0(mesh(1)%Nx-1,j,k)-D2Q9feq(mesh(1)%rho(mesh(1)%Nx-1,j),mesh(1)%u(mesh(1)%Nx-1,j),mesh(1)%v(mesh(1)%Nx-1,j),mesh(1)%Nx-1,j,k))
    !					    enddo
                        enddo
                        !$acc end loop
                        !$acc end kernels
                    endif
                endif
            elseif(ref_freestream.eq.2) then 
                if(mesh(1)%left(2,1).eq.0) then 
                    !$acc kernels loop independent present(mesh)
				    do j=2,mesh(1)%Ny-1
					    do k=1,Q
						    mesh(1)%rho(1,j)=1.0d0
						    mesh(1)%u(1,j)=0.0d0
						    mesh(1)%v(1,j)=mesh(1)%v(2,j)
						    mesh(1)%f(1,j,k)=D2Q9feq(mesh(1)%rho(1,j),mesh(1)%u(1,j),mesh(1)%v(1,j),1,j,k)+mesh(1)%f(2,j,k)-D2Q9feq(mesh(1)%rho(2,j),mesh(1)%u(2,j),mesh(1)%v(2,j),2,j,k)
					    enddo
				    enddo
                    !$acc end kernels loop
			    endif

			    if(mesh(1)%right(2,1).eq.0)then 
                    !$acc kernels loop independent present(mesh)
				    do j=2,mesh(1)%Ny-1
					    do k=1,Q
						    mesh(1)%rho(mesh(1)%Nx,j)=1.0d0
						    mesh(1)%u(mesh(1)%Nx,j)=0.0d0
						    mesh(1)%v(mesh(1)%Nx,j)=mesh(1)%v(mesh(1)%Nx-1,j)
						    mesh(1)%f(mesh(1)%Nx,j,k)=D2Q9feq(mesh(1)%rho(mesh(1)%Nx,j),mesh(1)%u(mesh(1)%Nx,j),mesh(1)%v(mesh(1)%Nx,j),mesh(1)%Nx,j,k)+mesh(1)%f(mesh(1)%Nx-1,j,k)-D2Q9feq(mesh(1)%rho(mesh(1)%Nx-1,j),mesh(1)%u(mesh(1)%Nx-1,j),mesh(1)%v(mesh(1)%Nx-1,j),mesh(1)%Nx-1,j,k)
					    enddo
				    enddo
                    !$acc end kernels loop
                endif
            endif
            
			if(mesh(1)%up(1,1).eq.0) then 
                !$acc kernels default(present) create(j)
                !$acc loop independent 
                do j=1,mesh(1)%Nx
                    temprho1(j)=mesh(1)%rho(j,mesh(1)%Ny)
                    tempu1(j)=mesh(1)%u(j,mesh(1)%Ny)
                    tempv1(j)=mesh(1)%v(j,mesh(1)%Ny)
                enddo
                !$acc end loop
                !$acc loop independent 
                do j=1,mesh(1)%Nx
                    if((mesh(1)%left(mesh(1)%Ny-1,1).ne.0).and.(j.eq.1))then
                        k0=mesh(1)%left(mesh(1)%Ny-1,1)
                        i0=mesh(1)%left(mesh(1)%Ny-1,3)-1
                        j0=mesh(1)%left(mesh(1)%Ny-1,4)+1
                        rhox=(temprho1(j+1)-mesh(k0)%rho(i0,j0))*0.5d0
                        ux=(tempu1(j+1)-mesh(k0)%u(i0,j0))*0.5d0
                        vx=(tempv1(j+1)-mesh(k0)%v(i0,j0))*0.5d0
                    elseif((mesh(1)%right(mesh(1)%Ny-1,1).ne.0).and.(j.eq.mesh(1)%Nx))then
                        k0=mesh(1)%right(mesh(1)%Ny-1,1)
                        i0=mesh(1)%right(mesh(1)%Ny-1,3)+1
                        j0=mesh(1)%right(mesh(1)%Ny-1,4)+1  
                        rhox=(mesh(k0)%rho(i0,j0)-temprho1(j-1))*0.5d0
                        ux=(mesh(k0)%u(i0,j0)-tempu1(j-1))*0.5d0
                        vx=(mesh(k0)%v(i0,j0)-tempv1(j-1))*0.5d0
                    elseif((j.ge.2).and.(j.le.(mesh(1)%Nx-1)))then
                        rhox=(temprho1(j+1)-temprho1(j-1))*0.5d0
                        ux=(tempu1(j+1)-tempu1(j-1))*0.5d0
                        vx=(tempv1(j+1)-tempv1(j-1))*0.5d0
                    else
                        rhox=0.0d0
                        ux=0.0d0
                        vx=0.0d0
                    endif
                    rhoy=(3.0d0*temprho1(j)-4.0d0*mesh(1)%rho(j,mesh(1)%Ny-1)+mesh(1)%rho(j,mesh(1)%Ny-2))*0.5d0                    
                    uy=(3.0d0*tempu1(j)-4.0d0*mesh(1)%u(j,mesh(1)%Ny-1)+mesh(1)%u(j,mesh(1)%Ny-2))*0.5d0                    
                    vy=(3.0d0*tempv1(j)-4.0d0*mesh(1)%v(j,mesh(1)%Ny-1)+mesh(1)%v(j,mesh(1)%Ny-2))*0.5d0
                    ly(1)=0.0d0 !(mesh(1)%u(1,j)-cs)*(cs*cs*rhox-cs*mesh(1)%rho(1,j)*ux)
                    if(tempv1(j).ge.0.0d0)then
                        ly(2)=tempv1(j)*uy
                    else
                        ly(2)=0.0d0
                    endif
                    ly(3)=(tempv1(j)+cs)*(cs*cs*rhoy+cs*temprho1(j)*vy)
                    lx(1)=(tempu1(j)-cs)*(cs*cs*rhox-cs*temprho1(j)*ux)
                    lx(2)=tempu1(j)*vx
                    lx(3)=(tempu1(j)+cs)*(cs*cs*rhox+cs*temprho1(j)*ux)
                            
                    mesh(1)%rho(j,mesh(1)%Ny)=mesh(1)%rho(j,mesh(1)%Ny)-(ly(1)+ly(3))*1.5d0 !-0.75d0*(lx(1)+lx(3))*1.5d0-(ly(1)+ly(3))*1.5d0
                    mesh(1)%u(j,mesh(1)%Ny)=mesh(1)%u(j,mesh(1)%Ny)-ly(2) !+0.75d0*(lx(1)-lx(3))/(2.0d0*cs*mesh(1)%rho(j,mesh(1)%Ny))-ly(2)
                    mesh(1)%v(j,mesh(1)%Ny)=mesh(1)%v(j,mesh(1)%Ny)+(ly(1)-ly(3))/(2.0d0*cs*mesh(1)%rho(j,mesh(1)%Ny)) !-0.75d0*lx(2)+(ly(1)-ly(3))/(2.0d0*cs*mesh(1)%rho(j,mesh(1)%Ny))
                enddo
                !$acc end loop 
                mesh(1)%rho(1,mesh(1)%Ny)=tempconr(1,1)       
                mesh(1)%u(1,mesh(1)%Ny)= tempconr(1,2)
                mesh(1)%v(1,mesh(1)%Ny)=tempconr(1,3)   
                mesh(1)%rho(mesh(1)%Nx,mesh(1)%Ny)=tempconr(2,1)       
                mesh(1)%u(mesh(1)%Nx,mesh(1)%Ny)= tempconr(2,2)
                mesh(1)%v(mesh(1)%Nx,mesh(1)%Ny)=tempconr(2,3)                       
                !$acc loop independent 
				do j=1,mesh(1)%Nx
                    mesh(1)%rho(j,mesh(1)%Ny)=mesh(1)%rho(j,mesh(1)%Ny-1) !1.0d0
                    mesh(1)%u(j,mesh(1)%Ny)=0.0d0 !mesh(1)%u(j,mesh(1)%Ny-1)!(4.0d0*mesh(1)%u(j,mesh(1)%Ny-1)-mesh(1)%u(j,mesh(1)%Ny-2))/3.0d0
                    mesh(1)%v(j,mesh(1)%Ny)=0.0d0 !mesh(1)%v(j,mesh(1)%Ny-1) !0.0d0 !(4.0d0*mesh(1)%v(j,mesh(1)%Ny-1)-mesh(1)%v(j,mesh(1)%Ny-2))/3.0d0
                    if((j.eq.1).and.(periodical_bc.eq.1))then
                        if(i.eq.1)then
                            mesh(1)%rho(j,mesh(1)%Ny)=mesh(1)%rho(mesh(1)%Nx-1,mesh(1)%Ny-1) !1.0d0
                            mesh(1)%u(j,mesh(1)%Ny)=mesh(1)%u(mesh(1)%Nx-1,mesh(1)%Ny-1) !mesh(1)%u(j,mesh(1)%Ny-1)!(4.0d0*mesh(1)%u(j,mesh(1)%Ny-1)-mesh(1)%u(j,mesh(1)%Ny-2))/3.0d0
                            mesh(1)%v(j,mesh(1)%Ny)=mesh(1)%v(mesh(1)%Nx-1,mesh(1)%Ny-1) ! mesh(1)%v(j,mesh(1)%Ny-1) !0.0d0 !(4.0d0*mesh(1)%v(j,mesh(1)%Ny-1)-mesh(1)%v(j,mesh(1)%Ny-2))/3.0d0    
                            mesh(1)%rho(j,mesh(1)%Ny-1)=mesh(1)%rho(mesh(1)%Nx-1,mesh(1)%Ny-1) !1.0d0
                            mesh(1)%u(j,mesh(1)%Ny-1)=mesh(1)%u(mesh(1)%Nx-1,mesh(1)%Ny-1) !mesh(1)%u(j,mesh(1)%Ny-1)!(4.0d0*mesh(1)%u(j,mesh(1)%Ny-1)-mesh(1)%u(j,mesh(1)%Ny-2))/3.0d0
                            mesh(1)%v(j,mesh(1)%Ny-1)=mesh(1)%v(mesh(1)%Nx-1,mesh(1)%Ny-1) ! mesh(1)%v(j,mesh(1)%Ny-1) !0.0d0 !(4.0d0*mesh(1)%v(j,mesh(1)%Ny-1)-mesh(1)%v(j,mesh(1)%Ny-2))/3.0d0    
                        endif
                    elseif((j.eq.mesh(1)%Nx).and.(periodical_bc.eq.1))then
                        if(i.eq.3)then
                            mesh(1)%rho(j,mesh(1)%Ny)=mesh(1)%rho(2,mesh(1)%Ny-1) !1.0d0
                            mesh(1)%u(j,mesh(1)%Ny)=mesh(1)%u(2,mesh(1)%Ny-1) !mesh(1)%u(j,mesh(1)%Ny-1)!(4.0d0*mesh(1)%u(j,mesh(1)%Ny-1)-mesh(1)%u(j,mesh(1)%Ny-2))/3.0d0
                            mesh(1)%v(j,mesh(1)%Ny)=mesh(1)%v(2,mesh(1)%Ny-1) ! mesh(1)%v(j,mesh(1)%Ny-1) !0.0d0 !(4.0d0*mesh(1)%v(j,mesh(1)%Ny-1)-mesh(1)%v(j,mesh(1)%Ny-2))/3.0d0    
                            mesh(1)%rho(j,mesh(1)%Ny-1)=mesh(1)%rho(2,mesh(1)%Ny-1) !1.0d0
                            mesh(1)%u(j,mesh(1)%Ny-1)=mesh(1)%u(2,mesh(1)%Ny-1) !mesh(1)%u(j,mesh(1)%Ny-1)!(4.0d0*mesh(1)%u(j,mesh(1)%Ny-1)-mesh(1)%u(j,mesh(1)%Ny-2))/3.0d0
                            mesh(1)%v(j,mesh(1)%Ny-1)=mesh(1)%v(2,mesh(1)%Ny-1) ! mesh(1)%v(j,mesh(1)%Ny-1) !0.0d0 !(4.0d0*mesh(1)%v(j,mesh(1)%Ny-1)-mesh(1)%v(j,mesh(1)%Ny-2))/3.0d0    
                        endif
                    endif
                    rh=mesh(1)%rho(j,mesh(1)%Ny)
                    jx=rh*mesh(1)%u(j,mesh(1)%Ny)
                    jy=rh*mesh(1)%v(j,mesh(1)%Ny)
                    temp1(1)=rh
                    temp1(2)=-2.0d0*rh+3.0d0*(jx*jx+jy*jy)/rh
                    temp1(3)=rh-3.0d0*(jx*jx+jy*jy)/rh
                    temp1(4)=jx
                    temp1(5)=-jx
                    temp1(6)=jy
                    temp1(7)=-jy
                    temp1(8)=(jx*jx-jy*jy)/rh
                    temp1(9)=(jx*jy)/rh
                    rh=mesh(1)%rho(j,mesh(1)%Ny-1)
                    jx=rh*mesh(1)%u(j,mesh(1)%Ny-1)
                    jy=rh*mesh(1)%v(j,mesh(1)%Ny-1)
                    temp2(1)=rh
                    temp2(2)=-2.0d0*rh+3.0d0*(jx*jx+jy*jy)/rh
                    temp2(3)=rh-3.0d0*(jx*jx+jy*jy)/rh
                    temp2(4)=jx
                    temp2(5)=-jx
                    temp2(6)=jy
                    temp2(7)=-jy
                    temp2(8)=(jx*jx-jy*jy)/rh
                    temp2(9)=(jx*jy)/rh
                    temp=mesh(1)%f0(j,mesh(1)%Ny-1,:)
                    if((j.eq.1).and.(periodical_bc.eq.1))then
                        if(i.eq.1) temp= mesh(1)%f0(mesh(1)%Nx-1,mesh(1)%Ny-1,:)
                    elseif((j.eq.mesh(1)%Nx).and.(periodical_bc.eq.1))then
                        if(i.eq.3)  temp= mesh(1)%f0(2,mesh(1)%Ny-1,:)
                    endif
                    call matmulvoc(temp3,matm,temp,Q,Q)
                    do k=1,Q
                        temp3(k)=temp1(k)+(1.0d0-mesh(1)%s(k))*(temp3(k)-temp2(k))
                    enddo
                    call matmulvoc(temp,inversematm,temp3,Q,Q)
                    mesh(1)%f(j,mesh(1)%Ny,:)=temp
!					do k=1,Q
!                        mesh(1)%rho(j,mesh(1)%Ny)=1.0d0
!                        mesh(1)%u(j,mesh(1)%Ny)=UU !(4.0d0*mesh(1)%u(j,mesh(1)%Ny-1)-mesh(1)%u(j,mesh(1)%Ny-2))/3.0d0
!                        mesh(1)%v(j,mesh(1)%Ny)=0.0d0 !(4.0d0*mesh(1)%v(j,mesh(1)%Ny-1)-mesh(1)%v(j,mesh(1)%Ny-2))/3.0d0
!                        mesh(1)%f(j,mesh(1)%Ny,k)=D2Q9feq(mesh(1)%rho(j,mesh(1)%Ny),mesh(1)%u(j,mesh(1)%Ny),mesh(1)%v(j,mesh(1)%Ny),j,mesh(1)%Ny,K)+(1.0d0-1.0d0/mesh(1)%tau_f)*(mesh(1)%f0(j,mesh(1)%Ny-1,k)-D2Q9feq(mesh(1)%rho(j,mesh(1)%Ny-1),mesh(1)%u(j,mesh(1)%Ny-1),mesh(1)%v(j,mesh(1)%Ny-1),j,mesh(1)%Ny-1,k))
!					enddo
				enddo
                !$acc end loop
                !$acc end kernels
			endif
			if(mesh(1)%down(1,1).eq.0) then 
                !$acc kernels default(none) present(mesh,tempu1,tempv1,temprho1,ly,lx,tempconr,temp,temp1,temp2,temp3) create(j)
                !present(mesh,cs,lx,ly,rhox,rhoy,ux,uy,vx,vy,temprho1,tempu1,tempv1,temprho,tempu,tempv,temp,temp1,temp2,temp3,rh,jx,jy,tempconr,i0,j0,k0)
                temprho1=mesh(1)%rho(:,1)
                tempu1=mesh(1)%u(:,1)
                tempv1=mesh(1)%v(:,1)
                !$acc loop independent 
                do j=1,mesh(1)%Nx
                    if((mesh(1)%left(2,1).ne.0).and.(j.eq.1))then
                        k0=mesh(1)%left(2,1)
                        i0=mesh(1)%left(2,3)-1
                        j0=mesh(1)%left(2,4)-1
                        rhox=(temprho1(j+1)-mesh(k0)%rho(i0,j0))*0.5d0
                        ux=(tempu1(j+1)-mesh(k0)%u(i0,j0))*0.5d0
                        vx=(tempv1(j+1)-mesh(k0)%v(i0,j0))*0.5d0
                    elseif((mesh(1)%right(2,1).ne.0).and.(j.eq.mesh(1)%Nx))then
                        k0=mesh(1)%right(2,1)
                        i0=mesh(1)%right(2,3)+1
                        j0=mesh(1)%right(2,4)-1  
!                        print*,j0
                        rhox=(mesh(k0)%rho(i0,j0)-temprho1(j-1))*0.5d0
                        ux=(mesh(k0)%u(i0,j0)-tempu1(j-1))*0.5d0
                        vx=(mesh(k0)%v(i0,j0)-tempv1(j-1))*0.5d0
                    elseif((j.ge.2).and.(j.le.(mesh(1)%Nx-1)))then
                        rhox=(temprho1(j+1)-temprho1(j-1))*0.5d0
                        ux=(tempu1(j+1)-tempu1(j-1))*0.5d0
                        vx=(tempv1(j+1)-tempv1(j-1))*0.5d0
                    else
                        rhox=0.0d0
                        ux=0.0d0
                        vx=0.0d0
                    endif
                    rhoy=(-3.0d0*temprho1(j)+4.0d0*mesh(1)%rho(j,2)-mesh(1)%rho(j,3))*0.5d0                    
                    uy=(-3.0d0*tempu1(j)+4.0d0*mesh(1)%u(j,2)-mesh(1)%u(j,3))*0.5d0                    
                    vy=(-3.0d0*tempv1(j)+4.0d0*mesh(1)%v(j,2)-mesh(1)%v(j,3))*0.5d0
                    ly(1)=(tempv1(j)-cs)*(cs*cs*rhoy-cs*temprho1(j)*vy)
                    if(tempv1(j).lt.0.0d0)then
                        ly(2)=tempv1(j)*uy
                    else
                        ly(2)=0.0d0
                    endif
                    ly(3)=0.0d0 !(mesh(1)%v(j,1)+cs)*(cs*cs*rhoy+cs*mesh(1)%rho(j,1)*vy)
                    lx(1)=(tempu1(j)-cs)*(cs*cs*rhox-cs*temprho1(j)*ux)
                    lx(2)=tempu1(j)*vx
                    lx(3)=(tempu1(j)+cs)*(cs*cs*rhox+cs*temprho1(j)*ux)
!                            
                    mesh(1)%rho(j,1)=mesh(1)%rho(j,1)-(ly(1)+ly(3))*1.5d0 !-0.75d0*(lx(1)+lx(3))*1.5d0-(ly(1)+ly(3))*1.5d0
                    mesh(1)%u(j,1)=mesh(1)%u(j,1)-ly(2) !+0.75d0*(lx(1)-lx(3))/(2.0d0*cs*mesh(1)%rho(j,1))-ly(2)
                    mesh(1)%v(j,1)=mesh(1)%v(j,1)+(ly(1)-ly(3))/(2.0d0*cs*mesh(1)%rho(j,1)) !-0.75d0*lx(2)+(ly(1)-ly(3))/(2.0d0*cs*mesh(1)%rho(j,1))
                enddo
                !$acc end loop
                mesh(1)%rho(1,1)= tempconr(3,1)       
                mesh(1)%u(1,1)= tempconr(3,2)
                mesh(1)%v(1,1)=tempconr(3,3)   
                mesh(1)%rho(mesh(1)%Nx,1)= tempconr(4,1)       
                mesh(1)%u(mesh(1)%Nx,1)= tempconr(4,2)
                mesh(1)%v(mesh(1)%Nx,1)=tempconr(4,3)                       
                !$acc loop independent 
				do j=1,mesh(1)%Nx
                    mesh(1)%rho(j,1)=mesh(1)%rho(j,2) !1.0d0
                    mesh(1)%u(j,1)=0.0d0 ! mesh(1)%u(j,2) !(4.0d0*mesh(1)%u(j,2)-mesh(1)%u(j,3))/3.0d0
                    mesh(1)%v(j,1)=0.0d0 !mesh(1)%v(j,2) !0.0d0 !(4.0d0*mesh(1)%v(j,2)-mesh(1)%v(j,3))/3.0d0
                    if((j.eq.1).and.(periodical_bc.eq.1))then
                        if(i.eq.1)then
                            mesh(1)%rho(j,1)=mesh(1)%rho(mesh(1)%Nx-1,2) !1.0d0
                            mesh(1)%u(j,1)=mesh(1)%u(mesh(1)%Nx-1,2) !mesh(1)%u(j,mesh(1)%Ny-1)!(4.0d0*mesh(1)%u(j,mesh(1)%Ny-1)-mesh(1)%u(j,mesh(1)%Ny-2))/3.0d0
                            mesh(1)%v(j,1)=mesh(1)%v(mesh(1)%Nx-1,2) ! mesh(1)%v(j,mesh(1)%Ny-1) !0.0d0 !(4.0d0*mesh(1)%v(j,mesh(1)%Ny-1)-mesh(1)%v(j,mesh(1)%Ny-2))/3.0d0   
                            mesh(1)%rho(j,2)=mesh(1)%rho(mesh(1)%Nx-1,2) !1.0d0
                            mesh(1)%u(j,2)=mesh(1)%u(mesh(1)%Nx-1,2) !mesh(1)%u(j,mesh(1)%Ny-1)!(4.0d0*mesh(1)%u(j,mesh(1)%Ny-1)-mesh(1)%u(j,mesh(1)%Ny-2))/3.0d0
                            mesh(1)%v(j,2)=mesh(1)%v(mesh(1)%Nx-1,2) ! mesh(1)%v(j,mesh(1)%Ny-1) !0.0d0 !(4.0d0*mesh(1)%v(j,mesh(1)%Ny-1)-mesh(1)%v(j,mesh(1)%Ny-2))/3.0d0 
                        endif
                    elseif((j.eq.mesh(1)%Nx).and.(periodical_bc.eq.1))then
                        if(i.eq.3)then
                            mesh(1)%rho(j,1)=mesh(1)%rho(2,2) !1.0d0
                            mesh(1)%u(j,1)=mesh(1)%u(2,2) !mesh(1)%u(j,mesh(1)%Ny-1)!(4.0d0*mesh(1)%u(j,mesh(1)%Ny-1)-mesh(1)%u(j,mesh(1)%Ny-2))/3.0d0
                            mesh(1)%v(j,1)=mesh(1)%v(2,2) ! mesh(1)%v(j,mesh(1)%Ny-1) !0.0d0 !(4.0d0*mesh(1)%v(j,mesh(1)%Ny-1)-mesh(1)%v(j,mesh(1)%Ny-2))/3.0d0  
                            mesh(1)%rho(j,2)=mesh(1)%rho(2,2) !1.0d0
                            mesh(1)%u(j,2)=mesh(1)%u(2,2) !mesh(1)%u(j,mesh(1)%Ny-1)!(4.0d0*mesh(1)%u(j,mesh(1)%Ny-1)-mesh(1)%u(j,mesh(1)%Ny-2))/3.0d0
                            mesh(1)%v(j,2)=mesh(1)%v(2,2) ! mesh(1)%v(j,mesh(1)%Ny-1) !0.0d0 !(4.0d0*mesh(1)%v(j,mesh(1)%Ny-1)-mesh(1)%v(j,mesh(1)%Ny-2))/3.0d0  
                        endif
                    endif
                    rh=mesh(1)%rho(j,1)
                    jx=rh*mesh(1)%u(j,1)
                    jy=rh*mesh(1)%v(j,1)
                    temp1(1)=rh
                    temp1(2)=-2.0d0*rh+3.0d0*(jx*jx+jy*jy)/rh
                    temp1(3)=rh-3.0d0*(jx*jx+jy*jy)/rh
                    temp1(4)=jx
                    temp1(5)=-jx
                    temp1(6)=jy
                    temp1(7)=-jy
                    temp1(8)=(jx*jx-jy*jy)/rh
                    temp1(9)=(jx*jy)/rh
                    rh=mesh(1)%rho(j,2)
                    jx=rh*mesh(1)%u(j,2)
                    jy=rh*mesh(1)%v(j,2)
                    temp2(1)=rh
                    temp2(2)=-2.0d0*rh+3.0d0*(jx*jx+jy*jy)/rh
                    temp2(3)=rh-3.0d0*(jx*jx+jy*jy)/rh
                    temp2(4)=jx
                    temp2(5)=-jx
                    temp2(6)=jy
                    temp2(7)=-jy
                    temp2(8)=(jx*jx-jy*jy)/rh
                    temp2(9)=(jx*jy)/rh
                    temp=mesh(1)%f0(j,2,:)
                    if((j.eq.1).and.(periodical_bc.eq.1))then
                        if(i.eq.1) temp= mesh(1)%f0(mesh(1)%Nx-1,2,:)
                    elseif((j.eq.mesh(1)%Nx).and.(periodical_bc.eq.1))then
                        if(i.eq.3)  temp= mesh(1)%f0(2,2,:)
                    endif
                    call matmulvoc(temp3,matm,temp,Q,Q)
                    do k=1,Q
                        temp3(k)=temp1(k)+(1.0d0-mesh(1)%s(k))*(temp3(k)-temp2(k))
                    enddo
                    call matmulvoc(temp,inversematm,temp3,Q,Q)
                    mesh(1)%f(j,1,:)=temp
!					do k=1,Q
!                        mesh(1)%rho(j,1)=1.0d0
!                        mesh(1)%u(j,1)=UU !(4.0d0*mesh(1)%u(j,2)-mesh(1)%u(j,3))/3.0d0
!                        mesh(1)%v(j,1)=0.0d0 !(4.0d0*mesh(1)%v(j,2)-mesh(1)%v(j,3))/3.0d0
!                        mesh(1)%f(j,1,k)=D2Q9feq(mesh(1)%rho(j,1),mesh(1)%u(j,1),mesh(1)%v(j,1),j,1,K)+(1.0d0-1.0d0/mesh(1)%tau_f)*(mesh(1)%f0(j,2,k)-D2Q9feq(mesh(1)%rho(j,2),mesh(1)%u(j,2),mesh(1)%v(j,2),j,2,k))
!					enddo
				enddo
                !$acc end loop
                !$acc end kernels
			endif
!            if(i.eq.2)print*, mesh(1)%f0(mesh(1)%Nx,mesh(1)%Ny-1,5)
        enddo
	endif
    !$acc exit data delete(cs,lx,ly,rhox,rhoy,ux,uy,vx,vy,tempconr,temp1,temp,temp3,temp2,temprho1,tempu1,tempv1,temprho,tempu,tempv,rh,jx,jy,i,i0,j0,k0)
end
    
    
subroutine wall_boundary_condition(mesh)
	use com_date
    use meshtype
	implicit none
	
	type(mesh_dat) mesh
	integer sp,ep,ip,jp,i,j,k,k1
	real (kind=ps) qr,rhob,ub,vb,wa,xb,yb
	Real (kind=ps), external :: D2Q9feq,distanceq,distanceqline
	real (kind=ps), allocatable :: ff3(:,:,:)
    logical ifplate !是平板还是圆柱

	ifplate=.true.
    if(myid.eq.0) then
		sp=mesh%start_row+1
		ep=mesh%end_row-3
	elseif(myid.eq.numprocs-1) then
		sp=mesh%start_row+3
		ep=mesh%end_row-1
	else
		sp=mesh%start_row+3
		ep=mesh%end_row-3
	endif
	if (numprocs.eq.1) then
		sp=mesh%start_row+1
		ep=mesh%end_row-1
	endif

	allocate (ff3(sp:ep,mesh%Ny,Q))

	do i=sp,ep
		do j=2,mesh%Ny-1
			do k=1,Q			
				ip=i-e(k,1)
				jp=j-e(k,2)			

				if (((mesh%body(i,j).eq.1).and. (mesh%body(ip,jp).eq.2)).or.((mesh%body(i,j).eq.2).and. (mesh%body(ip,jp).eq.1))) then
!					qr=distanceq(i,j,ip,jp,xb,yb,mesh%scale,mesh%ox,mesh%oy)
                    qr=distanceqline(i,j,ip,jp,xb,yb,mesh%scale,mesh%ox,mesh%oy)
					ub=UU*(-omega_arr(1)*sin(angle_arr(1))*xb+omega_arr(1)*cos(angle_arr(1))*yb+sau_arr(1))
					vb=UU*(-omega_arr(1)*cos(angle_arr(1))*xb-omega_arr(1)*sin(angle_arr(1))*yb+sbv_arr(1))
					call reversedirection(k,k1)
					if((k1.ge.2).and.(k1.le.5)) then
						wa=2.0/9.0
					elseif((k1.ge.6).and.(k1.le.9)) then
						wa=2.0/36.0
					endif
					if(qr.lt.0.5d0) then
						ff3(i,j,k)=qr*(1.0+2.0*qr)*mesh%f0(i-e(k,1),j-e(k,2),k1)+(1.0d0-4.0d0*qr*qr)*mesh%f0(i,j,k1)-qr*(1.0d0-2.0*qr)*mesh%f0(i+e(k,1),j+e(k,2),k1)+3.0*wa*mesh%rho(i,j)*(dble(e(k1,1))*ub+dble(e(k1,2))*vb)
					else
						ff3(i,j,k)=mesh%f0(i-e(k,1),j-e(k,2),k1)/(qr*(2.0*qr+1.0))+(2.0*qr-1.0)/qr*mesh%f0(i+e(k,1),j+e(k,2),k)-(2.0*qr-1.0)/(2.0*qr+1.0)*mesh%f0(i+2*e(k,1),j+2*e(k,2),k)+3.0*wa*mesh%rho(i,j)*(dble(e(k1,1))*ub+dble(e(k1,2))*vb)/(qr*(2.0*qr+1.0))
					endif
					myforce(1)=myforce(1)-e(k,1)*ff3(i,j,k)
					myforce(2)=myforce(2)-e(k,2)*ff3(i,j,k)
				endif
			enddo
			
		enddo
		
	enddo
	
	do i=sp,ep
		do j=2,mesh%Ny-1
			do k=1,Q			
				ip=i-e(k,1)
				jp=j-e(k,2)			

				if (((mesh%body(i,j).eq.1).and. (mesh%body(ip,jp).eq.2)).or.((mesh%body(i,j).eq.2).and. (mesh%body(ip,jp).eq.1))) then
					mesh%f0(i,j,k)=ff3(i,j,k)
				endif
			enddo
		enddo
	enddo

	deallocate (ff3)

end

subroutine immersed_boundary_wang(mesh)
	use com_date
    use meshtype
	implicit none
	
	type(mesh_dat) mesh
	integer sp,ep,ip,jp,i,j,k,k1
	integer cout(mesh%Nx,mesh%Ny)
	real*8 constx,consty,sectionx,sectiony,ub,vb,temu,temv
	real*8 ustar(mesh%Nx,mesh%Ny),vstar(mesh%Nx,mesh%Ny)
	integer xt,yt,xt1,xt2,yt1,yt2
	real*8 x1,x2,x3,y1,y2,y3,x0,y0

	if(myid.eq.0) then
		sp=mesh%start_row+1
		ep=mesh%end_row-3
	elseif(myid.eq.numprocs-1) then
		sp=mesh%start_row+3
		ep=mesh%end_row-1
	else
		sp=mesh%start_row+3
		ep=mesh%end_row-3
	endif
	if (numprocs.eq.1) then
		sp=mesh%start_row+1
		ep=mesh%end_row-1
	endif

	if (myid.eq.0)then
		cout=0
		gforce=0.0d0
		ustar=mesh%u
		vstar=mesh%v
		do i=1,mesh%Nx
			constx=mesh%ox+mesh%scale*dble(i-1)
			do j=1,b_point-2
				if((constx.gt.min(xb1(j),xb1(j+1))).and.(constx.le.max(xb1(j),xb1(j+1)))) then
					sectionx=constx
					sectiony=yb1(j)+(sectionx-xb1(j))*(yb1(j+1)-yb1(j))/(xb1(j+1)-xb1(j))
					ub=UU*(-omega_arr(1)*sin(angle_arr(1))*(sectionx-xc)+omega_arr(1)*cos(angle_arr(1))*(sectiony-yc)+sau_arr(1))
					vb=UU*(-omega_arr(1)*cos(angle_arr(1))*(sectionx-xc)-omega_arr(1)*sin(angle_arr(1))*(sectiony-yc)+sbv_arr(1))
					yt1=int((sectiony-mesh%oy)/mesh%scale)+1
					yt2=int((sectiony-mesh%oy)/mesh%scale)+2
!					if((abs(mesh%oy+mesh%scale*dble(yt1-1)-sectiony).lt.1d-6).or.(abs(mesh%oy+mesh%scale*dble(yt2-1)-sectiony).lt.1d-6)) then
!						yt1=nint((sectiony-mesh%oy)/mesh%scale)
!						yt2=nint((sectiony-mesh%oy)/mesh%scale)+2
!						sectiony=mesh%oy+mesh%scale*dble(yt1)
!						ub=UU*(-omega_arr(1)*sin(angle_arr(1))*(sectionx-xc)+omega_arr(1)*cos(angle_arr(1))*(sectiony-yc)-sau_arr(1))
!						vb=UU*(-omega_arr(1)*cos(angle_arr(1))*(sectionx-xc)-omega_arr(1)*sin(angle_arr(1))*(sectiony-yc)+sbv_arr(1))
!						mesh%u(i,yt1+1)=ub
!						mesh%v(i,yt1+1)=vb
!						cout(i,yt1+1)=1
!						gforce(i,yt1+1,1)=2.0d0*mesh%rho(i,yt1+1)*(mesh%u(i,yt1+1)-ustar(i,yt1+1))
!						gforce(i,yt1+1,2)=2.0d0*mesh%rho(i,yt1+1)*(mesh%v(i,yt1+1)-vstar(i,yt1+1))
!						y0=mesh%oy+mesh%scale*dble(yt1-1)
!						y1=mesh%oy+mesh%scale*dble(yt1-2)
!						y2=mesh%oy+mesh%scale*dble(yt1)
!						y3=mesh%oy+mesh%scale*dble(yt2)
!						temu=(((y0-y2)*(y0-y3))/((y1-y2)*(y1-y3)))*ustar(i,yt1-1)+(((y0-y1)*(y0-y3))/((y2-y1)*(y2-y3)))*ub+(((y0-y1)*(y0-y2))/((y3-y1)*(y3-y2)))*ustar(i,yt2+1)
!						temv=(((y0-y2)*(y0-y3))/((y1-y2)*(y1-y3)))*vstar(i,yt1-1)+(((y0-y1)*(y0-y3))/((y2-y1)*(y2-y3)))*vb+(((y0-y1)*(y0-y2))/((y3-y1)*(y3-y2)))*vstar(i,yt2+1)
!						mesh%u(i,yt1)=(dble(cout(i,yt1))*mesh%u(i,yt1)+temu)/dble(cout(i,yt1)+1)
!						mesh%v(i,yt1)=(dble(cout(i,yt1))*mesh%v(i,yt1)+temv)/dble(cout(i,yt1)+1)
!						cout(i,yt1)=cout(i,yt1)+1
!						gforce(i,yt1,1)=2.0d0*mesh%rho(i,yt1)*(mesh%u(i,yt1)-ustar(i,yt1))
!						gforce(i,yt1,2)=2.0d0*mesh%rho(i,yt1)*(mesh%v(i,yt1)-vstar(i,yt1))
!						y0=mesh%oy+mesh%scale*dble(yt2-1)
!						temu=(((y0-y2)*(y0-y3))/((y1-y2)*(y1-y3)))*ustar(i,yt1-1)+(((y0-y1)*(y0-y3))/((y2-y1)*(y2-y3)))*ub+(((y0-y1)*(y0-y2))/((y3-y1)*(y3-y2)))*ustar(i,yt2+1)
!						temv=(((y0-y2)*(y0-y3))/((y1-y2)*(y1-y3)))*vstar(i,yt1-1)+(((y0-y1)*(y0-y3))/((y2-y1)*(y2-y3)))*vb+(((y0-y1)*(y0-y2))/((y3-y1)*(y3-y2)))*vstar(i,yt2+1)
!						mesh%u(i,yt2)=(dble(cout(i,yt2))*mesh%u(i,yt2)+temu)/dble(cout(i,yt2)+1)
!						mesh%v(i,yt2)=(dble(cout(i,yt2))*mesh%v(i,yt2)+temv)/dble(cout(i,yt2)+1)
!						cout(i,yt2)=cout(i,yt2)+1
!						gforce(i,yt2,1)=2.0d0*mesh%rho(i,yt2)*(mesh%u(i,yt2)-ustar(i,yt2))
!						gforce(i,yt2,2)=2.0d0*mesh%rho(i,yt2)*(mesh%v(i,yt2)-vstar(i,yt2))
!					else
						y0=mesh%oy+mesh%scale*dble(yt1-1)
						y1=mesh%oy+mesh%scale*dble(yt1-2)
						y2=sectiony
						y3=mesh%oy+mesh%scale*dble(yt2)
						temu=(((y0-y2)*(y0-y3))/((y1-y2)*(y1-y3)))*ustar(i,yt1-1)+(((y0-y1)*(y0-y3))/((y2-y1)*(y2-y3)))*ub+(((y0-y1)*(y0-y2))/((y3-y1)*(y3-y2)))*ustar(i,yt2+1)
						temv=(((y0-y2)*(y0-y3))/((y1-y2)*(y1-y3)))*vstar(i,yt1-1)+(((y0-y1)*(y0-y3))/((y2-y1)*(y2-y3)))*vb+(((y0-y1)*(y0-y2))/((y3-y1)*(y3-y2)))*vstar(i,yt2+1)
!						print*, y1,mesh%u(i,yt1-1),y2,mesh%u(i,yt1+1),y3,mesh%u(i,yt2+1),y0,temu
!						pause
!						if(cout(i,yt1).eq.2) print*,1,i,yt1,j,constx,min(xb1(j),xb1(j+1)),max(xb1(j),xb1(j+1))
						if(cout(i,yt1).eq.0) then
							mesh%u(i,yt1)=(dble(cout(i,yt1))*mesh%u(i,yt1)+temu)/dble(cout(i,yt1)+1)
							mesh%v(i,yt1)=(dble(cout(i,yt1))*mesh%v(i,yt1)+temv)/dble(cout(i,yt1)+1)
							cout(i,yt1)=cout(i,yt1)+1
							gforce(i,yt1,1)=2.0d0*mesh%rho(i,yt1)*(mesh%u(i,yt1)-ustar(i,yt1))
							gforce(i,yt1,2)=2.0d0*mesh%rho(i,yt1)*(mesh%v(i,yt1)-vstar(i,yt1))
						endif
						y0=mesh%oy+mesh%scale*dble(yt2-1)
						temu=(((y0-y2)*(y0-y3))/((y1-y2)*(y1-y3)))*ustar(i,yt1-1)+(((y0-y1)*(y0-y3))/((y2-y1)*(y2-y3)))*ub+(((y0-y1)*(y0-y2))/((y3-y1)*(y3-y2)))*ustar(i,yt2+1)
						temv=(((y0-y2)*(y0-y3))/((y1-y2)*(y1-y3)))*vstar(i,yt1-1)+(((y0-y1)*(y0-y3))/((y2-y1)*(y2-y3)))*vb+(((y0-y1)*(y0-y2))/((y3-y1)*(y3-y2)))*vstar(i,yt2+1)
!						if(cout(i,yt2).eq.2) print*,2,i,yt2,j,constx,min(xb1(j),xb1(j+1)),max(xb1(j),xb1(j+1))
						if(cout(i,yt2).eq.0) then
							mesh%u(i,yt2)=(dble(cout(i,yt2))*mesh%u(i,yt2)+temu)/dble(cout(i,yt2)+1)
							mesh%v(i,yt2)=(dble(cout(i,yt2))*mesh%v(i,yt2)+temv)/dble(cout(i,yt2)+1)
							cout(i,yt2)=cout(i,yt2)+1
							gforce(i,yt2,1)=2.0d0*mesh%rho(i,yt2)*(mesh%u(i,yt2)-ustar(i,yt2))
							gforce(i,yt2,2)=2.0d0*mesh%rho(i,yt2)*(mesh%v(i,yt2)-vstar(i,yt2))
						endif
!					endif
				endif
			enddo
		enddo
		do i=1,mesh%Ny
			consty=mesh%oy+mesh%scale*dble(i-1)
			do j=1,b_point-2
				if((consty.gt.min(yb1(j),yb1(j+1))).and.(consty.le.max(yb1(j),yb1(j+1)))) then
					sectiony=consty
					sectionx=xb1(j)+(sectiony-yb1(j))*(xb1(j+1)-xb1(j))/(yb1(j+1)-yb1(j))
					ub=UU*(-omega_arr(1)*sin(angle_arr(1))*(sectionx-xc)+omega_arr(1)*cos(angle_arr(1))*(sectiony-yc)+sau_arr(1))
					vb=UU*(-omega_arr(1)*cos(angle_arr(1))*(sectionx-xc)-omega_arr(1)*sin(angle_arr(1))*(sectiony-yc)+sbv_arr(1))
					xt1=int((sectionx-mesh%ox)/mesh%scale)+1
					xt2=int((sectionx-mesh%ox)/mesh%scale)+2
!					if((abs(mesh%ox+mesh%scale*dble(xt1-1)-sectionx).lt.1d-6).or.(abs(mesh%ox+mesh%scale*dble(xt2-1)-sectionx).lt.1d-6)) then
!						xt1=nint((sectionx-mesh%ox)/mesh%scale)
!						xt2=nint((sectionx-mesh%ox)/mesh%scale)+2
!						sectionx=mesh%ox+mesh%scale*dble(xt1)
!						ub=UU*(-omega_arr(1)*sin(angle_arr(1))*(sectionx-xc)+omega_arr(1)*cos(angle_arr(1))*(sectiony-yc)-sau_arr(1))
!						vb=UU*(-omega_arr(1)*cos(angle_arr(1))*(sectionx-xc)-omega_arr(1)*sin(angle_arr(1))*(sectiony-yc)+sbv_arr(1))
!						mesh%u(xt1+1,i)=ub
!						mesh%v(xt1+1,i)=vb
!						cout(xt1+1,i)=1
!						gforce(xt1+1,i,1)=2.0d0*mesh%rho(xt1+1,i)*(mesh%u(xt1+1,i)-ustar(xt1+1,i))
!						gforce(xt1+1,i,2)=2.0d0*mesh%rho(xt1+1,i)*(mesh%v(xt1+1,i)-vstar(xt1+1,i))
!						x0=mesh%ox+mesh%scale*dble(xt1-1)
!						x1=mesh%ox+mesh%scale*dble(xt1-2)
!						x2=mesh%ox+mesh%scale*dble(xt1)
!						x3=mesh%ox+mesh%scale*dble(xt2)
!						temu=(((x0-x2)*(x0-x3))/((x1-x2)*(x1-x3)))*ustar(xt1-1,i)+(((x0-x1)*(x0-x3))/((x2-x1)*(x2-x3)))*ub+(((x0-x1)*(x0-x2))/((x3-x1)*(x3-x2)))*ustar(xt2+1,i)
!						temv=(((x0-x2)*(x0-x3))/((x1-x2)*(x1-x3)))*vstar(xt1-1,i)+(((x0-y1)*(x0-x3))/((x2-x1)*(x2-x3)))*vb+(((x0-x1)*(x0-x2))/((x3-x1)*(x3-x2)))*vstar(xt2+1,i)
!						mesh%u(xt1,i)=(dble(cout(xt1,i))*mesh%u(xt1,i)+temu)/dble(cout(xt1,i)+1)
!						mesh%v(xt1,i)=(dble(cout(xt1,i))*mesh%v(xt1,i)+temv)/dble(cout(xt1,i)+1)
!						cout(xt1,i)=cout(xt1,i)+1
!						gforce(xt1,i,1)=2.0d0*mesh%rho(xt1,i)*(mesh%u(xt1,i)-ustar(xt1,i))
!						gforce(xt1,i,2)=2.0d0*mesh%rho(xt1,i)*(mesh%v(xt1,i)-vstar(xt1,i))
!						x0=mesh%ox+mesh%scale*dble(xt2-1)
!						temu=(((x0-x2)*(x0-x3))/((x1-x2)*(x1-x3)))*ustar(xt1-1,i)+(((x0-x1)*(x0-x3))/((x2-x1)*(x2-x3)))*ub+(((x0-x1)*(x0-x2))/((x3-x1)*(x3-x2)))*ustar(xt2+1,i)
!						temv=(((x0-x2)*(x0-x3))/((x1-x2)*(x1-x3)))*vstar(xt1-1,i)+(((x0-x1)*(x0-x3))/((x2-x1)*(x2-x3)))*vb+(((x0-x1)*(x0-x2))/((x3-x1)*(x3-x2)))*vstar(xt2+1,i)
!						mesh%u(xt2,i)=(dble(cout(xt2,i))*mesh%u(xt2,i)+temu)/dble(cout(xt2,i)+1)
!						mesh%v(xt2,i)=(dble(cout(xt2,i))*mesh%v(xt2,i)+temv)/dble(cout(xt2,i)+1)
!						cout(xt2,i)=cout(xt2,i)+1
!						gforce(xt2,i,1)=2.0d0*mesh%rho(xt2,i)*(mesh%u(xt2,i)-ustar(xt2,i))
!						gforce(xt2,i,2)=2.0d0*mesh%rho(xt2,i)*(mesh%v(xt2,i)-vstar(xt2,i))
!					else
						x0=mesh%ox+mesh%scale*dble(xt1-1)
						x1=mesh%ox+mesh%scale*dble(xt1-2)
						x2=sectionx
						x3=mesh%ox+mesh%scale*dble(xt2)
						temu=(((x0-x2)*(x0-x3))/((x1-x2)*(x1-x3)))*ustar(xt1-1,i)+(((x0-x1)*(x0-x3))/((x2-x1)*(x2-x3)))*ub+(((x0-x1)*(x0-x2))/((x3-x1)*(x3-x2)))*ustar(xt2+1,i)
						temv=(((x0-x2)*(x0-x3))/((x1-x2)*(x1-x3)))*vstar(xt1-1,i)+(((x0-x1)*(x0-x3))/((x2-x1)*(x2-x3)))*vb+(((x0-x1)*(x0-x2))/((x3-x1)*(x3-x2)))*vstar(xt2+1,i)
!						if(cout(xt1,i).eq.2) print*,3,xt1,i,j,constx,min(xb1(j),xb1(j+1)),max(xb1(j),xb1(j+1))
						if(cout(xt1,i).le.1) then
							mesh%u(xt1,i)=(dble(cout(xt1,i))*mesh%u(xt1,i)+temu)/dble(cout(xt1,i)+1)
							mesh%v(xt1,i)=(dble(cout(xt1,i))*mesh%v(xt1,i)+temv)/dble(cout(xt1,i)+1)
							cout(xt1,i)=cout(xt1,i)+1
							gforce(xt1,i,1)=2.0d0*mesh%rho(xt1,i)*(mesh%u(xt1,i)-ustar(xt1,i))
							gforce(xt1,i,2)=2.0d0*mesh%rho(xt1,i)*(mesh%v(xt1,i)-vstar(xt1,i))
						endif
						x0=mesh%ox+mesh%scale*dble(xt2-1)
						temu=(((x0-x2)*(x0-x3))/((x1-x2)*(x1-x3)))*ustar(xt1-1,i)+(((x0-x1)*(x0-x3))/((x2-x1)*(x2-x3)))*ub+(((x0-x1)*(x0-x2))/((x3-x1)*(x3-x2)))*ustar(xt2+1,i)
						temv=(((x0-x2)*(x0-x3))/((x1-x2)*(x1-x3)))*vstar(xt1-1,i)+(((x0-x1)*(x0-x3))/((x2-x1)*(x2-x3)))*vb+(((x0-x1)*(x0-x2))/((x3-x1)*(x3-x2)))*vstar(xt2+1,i)
!						if(cout(xt2,i).eq.2) print*,4,xt2,i,j,constx,min(xb1(j),xb1(j+1)),max(xb1(j),xb1(j+1))
						if(cout(xt2,i).le.1) then
							mesh%u(xt2,i)=(dble(cout(xt2,i))*mesh%u(xt2,i)+temu)/dble(cout(xt2,i)+1)
							mesh%v(xt2,i)=(dble(cout(xt2,i))*mesh%v(xt2,i)+temv)/dble(cout(xt2,i)+1)
							cout(xt2,i)=cout(xt2,i)+1
							gforce(xt2,i,1)=2.0d0*mesh%rho(xt2,i)*(mesh%u(xt2,i)-ustar(xt2,i))
							gforce(xt2,i,2)=2.0d0*mesh%rho(xt2,i)*(mesh%v(xt2,i)-vstar(xt2,i))
						endif
!					endif
				endif
			enddo
		enddo
	endif

!	print*, "Wang=",gforce(31,59,1),gforce(31,59,1)
	force(1)=-Sum(gforce(:,:,1))
	force(2)=-Sum(gforce(:,:,2))
end

subroutine immersed_boundary_suzuki_and_inamuro(mesh)
    use mpi
	use com_date
    use meshtype
	implicit none
	
	type(mesh_dat) mesh
	integer sp,ep,ip,jp,i,j,k,k1,num,temn,smooth,spp,epp
	real(kind=ps),  allocatable ::  ub(:),vb(:),ustar(:),vstar(:),rhos(:),ustar1(:),vstar1(:),rhos1(:),rhos2(:),buf1(:),buf2(:),as(:),ul(:),vl(:),glx(:),gly(:)
	real(kind=ps),  allocatable :: tem1(:,:),tem2(:,:),aa(:,:),aa1(:,:),bu(:),bv(:),tem11(:,:),tem22(:,:),tema(:,:),comi(:),x(:,:),xx1(:,:),xx2(:,:),xx3(:,:)
	real*8 disr,x0,y0,kernel1,kernel2,sum1,sum2
    integer nono,mn,nn,js,js1,ds,nnn
	integer,  allocatable :: non_zero(:),indexnum(:,:),recvcounts(:),displs(:),mrpos(:),nrpos(:),mrow(:),mcol(:),mrrow(:),mrcol(:),ia(:),ja(:)
    character(len=512) :: out_path

    smooth=1
    sp=myid*int((b_point-1)/numprocs)
	ep=(myid+1)*int((b_point-1)/numprocs)-1
    if(myid.eq.0) sp=1
    if(myid.eq.numprocs-1) ep=b_point-1
    
	num=adjoint_n
    nnn=1
    if(myid.eq.0) then
        call build_output_path('hhg.dat', out_path)
        open(2001,file=trim(out_path),status='unknown')
    endif
    
	allocate (ub(b_point-1),vb(b_point-1),ustar(num),vstar(num),rhos(num),ustar1(num),vstar1(num),rhos1(num),rhos2(num),tema(b_point-1,b_point-1),comi(b_point-1)) !,buf1(b_point-1))
	allocate (tem1(b_point-1,num),tem2(num,b_point-1),aa(b_point-1,b_point-1),aa1(b_point-1,b_point-1),bu(b_point-1),bv(b_point-1),x(b_point-1,2),xx1((b_point-1)*num,3),xx2((b_point-1)*num,3),xx3((b_point-1)*num,3),as((b_point-1)*(b_point-1)))
	allocate (non_zero(num),indexnum(num,2),recvcounts(numprocs),displs(numprocs),mrpos(b_point-1),nrpos(num),mcol(num),mrrow(b_point-1),mrow(b_point-1),mrcol(num),ia(b_point),ja((b_point-1)*(b_point-1)),ul(num),vl(num),glx(num),gly(num))
    
!    angle_arr(1)=-body_rel_theta(1)
!    omega_arr(1)=-body_rel_temp_theta(1)/(mesh%scale*UU)

    do i=1,b_point-1
        ub(i)=UU*(-omega_arr(1)*sin(angle_arr(1))*0.5d0*(xb0(i+1)+xb0(i))+omega_arr(1)*cos(angle_arr(1))*0.5d0*(yb0(i+1)+yb0(i))+sau_arr(1))
        vb(i)=UU*(-omega_arr(1)*cos(angle_arr(1))*0.5d0*(xb0(i+1)+xb0(i))-omega_arr(1)*sin(angle_arr(1))*0.5d0*(yb0(i+1)+yb0(i))+sbv_arr(1))
    enddo

    ub=ub+body_rel_v(1) !+body_rel_a(1)

    js=1
    mcol=0
    mrow=0 
    ustar1=0.0d0
    vstar1=0.0d0
    rhos1=0.0d0

    do k=sp,ep
        do i=1,adjoint_n
            nono=1
            x0=mesh%ox+mesh%scale*dble(adjoint_p(i,1)-1)
            y0=mesh%oy+mesh%scale*dble(adjoint_p(i,2)-1)
            if(smooth.eq.0) then
                disr=abs((x0-xb2(k))/mesh%scale)
                if(abs(disr).lt.1.0d0)then
                    kernel1=0.125d0*(3.0d0-2.0d0*abs(disr)+sqrt(1.0d0+4.0d0*abs(disr)-4.0d0*disr*disr))
                elseif ((abs(disr).lt.2.0d0).and.(abs(disr).ge.1.0d0))then
                    kernel1=0.125d0*(5.0d0-2.0d0*abs(disr)-sqrt(-7.0d0+12.0d0*abs(disr)-4.0d0*disr*disr))
                else
                    kernel1=0.0d0
                    nono=0
                endif
                disr=abs((y0-yb2(k))/mesh%scale)
                if(abs(disr).lt.1.0d0)then
                    kernel2=0.125d0*(3.0d0-2.0d0*abs(disr)+sqrt(1.0d0+4.0d0*abs(disr)-4.0d0*disr*disr))
                elseif ((abs(disr).lt.2.0d0).and.(abs(disr).ge.1.0d0))then
                    kernel2=0.125d0*(5.0d0-2.0d0*abs(disr)-sqrt(-7.0d0+12.0d0*abs(disr)-4.0d0*disr*disr))
                else
                    kernel2=0.0d0
                    nono=0
                endif
            else
                disr=abs((x0-xb2(k))/mesh%scale)
                if(abs(disr).le.0.5d0)then
                    kernel1=0.375d0+0.03125d0*pai-0.25d0*disr*disr
                elseif ((abs(disr).le.1.5d0).and.(abs(disr).ge.0.5d0)) then
                    kernel1=0.25d0+0.125d0*(1.0d0-disr)*sqrt(-2.0d0+8.0d0*disr-4.0d0*disr*disr)-0.125*asin(sqrt(2.0d0)*(disr-1.0d0))
                elseif ((abs(disr).le.2.5d0).and.(abs(disr).ge.1.5d0)) then
                    kernel1=1.0625d0-0.015625d0*pai-0.75d0*disr+0.125*disr*disr+0.0625d0*(disr-2.0d0)*sqrt(-14.0d0+16.0d0*disr-4.0d0*disr*disr)+0.0625d0*asin(sqrt(2.0d0)*(disr-2.0d0))
                else
                    kernel1=0.0d0
                    nono=0
                endif
                disr=abs((y0-yb2(k))/mesh%scale)
                if(abs(disr).le.0.5d0)then
                    kernel2=0.375d0+0.03125d0*pai-0.25d0*disr*disr
                elseif ((abs(disr).le.1.5d0).and.(abs(disr).ge.0.5d0)) then
                    kernel2=0.25d0+0.125d0*(1.0d0-disr)*sqrt(-2.0d0+8.0d0*disr-4.0d0*disr*disr)-0.125*asin(sqrt(2.0d0)*(disr-1.0d0))
                elseif ((abs(disr).le.2.5d0).and.(abs(disr).ge.1.5d0)) then
                    kernel2=1.0625d0-0.015625d0*pai-0.75d0*disr+0.125*disr*disr+0.0625d0*(disr-2.0d0)*sqrt(-14.0d0+16.0d0*disr-4.0d0*disr*disr)+0.0625d0*asin(sqrt(2.0d0)*(disr-2.0d0))
                else
                    kernel2=0.0d0
                    nono=0
                endif
            endif
            if(nono.eq.1) then
                xx3(js,1)=dble(k)
                xx3(js,2)=dble(i)
                xx3(js,3)=kernel1*kernel2*arc_l(k)
                xx1(js,1)=dble(k)
                xx1(js,2)=dble(i)
                xx1(js,3)=(kernel1*kernel2)
                mrow(k)=mrow(k)+1
                mcol(i)=mcol(i)+1
                js=js+1
            endif
        enddo
    enddo
    js=js-1
          
    if(numprocs.ne.1) then
        call MPI_ALLREDUCE(mrow,mrrow,b_point-1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(mcol,mrcol,adjoint_n,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)

        mrow=mrrow
        mcol=mrcol

        call MPI_ALLGATHER(js,1,MPI_INTEGER,recvcounts,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
        js1=sum(recvcounts)

        allocate(buf1(3*js),buf2(3*js1))

        recvcounts=3*recvcounts
        displs=0
        do j=1,numprocs-1
            displs(j+1)=displs(j)+recvcounts(j)
        enddo
        j=1
        do i=1,js
            buf1(j)=xx1(i,1)
            buf1(j+1)=xx1(i,2)
            buf1(j+2)=xx1(i,3)
            j=j+3
        enddo
 
        call MPI_ALLGATHERV(buf1,recvcounts(myid+1),MPI_DOUBLE_PRECISION,buf2,recvcounts,displs,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
        j=1
        do i=1,js1
            xx1(i,1)=buf2(j)
            xx1(i,2)=buf2(j+1)
            xx1(i,3)=buf2(j+2)
            j=j+3
        enddo

        j=1
        do i=1,js
            buf1(j)=xx3(i,1)
            buf1(j+1)=xx3(i,2)
            buf1(j+2)=xx3(i,3)
            j=j+3
        enddo
        call MPI_ALLGATHERV(buf1,recvcounts(myid+1),MPI_DOUBLE_PRECISION,buf2,recvcounts,displs,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
        j=1
        do i=1,js1
            xx3(i,1)=buf2(j)
            xx3(i,2)=buf2(j+1)
            xx3(i,3)=buf2(j+2)
            j=j+3
        enddo

        if(myid.eq.0) then
		    sp=mesh%start_row+1
		    ep=mesh%end_row-1
	    elseif(myid.eq.numprocs-1) then
		    sp=mesh%start_row+1
		    ep=mesh%end_row-1
	    else
		    sp=mesh%start_row+1
		    ep=mesh%end_row-1
	    endif
        do i=1,adjoint_n
            if((adjoint_p(i,1).ge.sp).and.(adjoint_p(i,1).le.ep))then
                ustar1(i)=mesh%u(adjoint_p(i,1),adjoint_p(i,2))
                vstar1(i)=mesh%v(adjoint_p(i,1),adjoint_p(i,2))
                rhos1(i)=mesh%rho(adjoint_p(i,1),adjoint_p(i,2))
            endif
        enddo

        call MPI_ALLREDUCE(ustar1,ustar,num,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(vstar1,vstar,num,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(rhos1,rhos,num,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

    else
        do i=1,adjoint_n       
            ustar(i)=mesh%u(adjoint_p(i,1),adjoint_p(i,2))
            vstar(i)=mesh%v(adjoint_p(i,1),adjoint_p(i,2))
        enddo
    endif
   
    mrpos(1)=1
    do k=2,b_point-1
        mrpos(k)=mrpos(k-1)+mrow(k-1)
    enddo
    nrpos(1)=1
    do k=2,num
         nrpos(k)=nrpos(k-1)+mcol(k-1)
    enddo

    if(numprocs.ne.1)then
        call fasttransposematrix(xx3,xx2,num,b_point-1,js1,nrpos)
    else
        call fasttransposematrix(xx3,xx2,num,b_point-1,js,nrpos)
    endif
    
!    if(numprocs.ne.1)then 
!        call multsmatrix(xx1,xx2,mrpos,nrpos,b_point-1,num,js1,as,ia,ja,ds)
!    else
!        call multsmatrix(xx1,xx2,mrpos,nrpos,b_point-1,num,js,as,ia,ja,ds)
!    endif

    if(numprocs.ne.1)then 
        call matmulvocxs(bu,xx1,mrpos,ustar,b_point-1,num,js1)
    else
        call matmulvocxs(bu,xx1,mrpos,ustar,b_point-1,num,js)
    endif
    bu=ub-bu  

    if(numprocs.ne.1)then 
        call matmulvocxs(bv,xx1,mrpos,vstar,b_point-1,num,js1)
    else
        call matmulvocxs(bv,xx1,mrpos,vstar,b_point-1,num,js)
    endif
   
    bv=vb-bv
  
    do k=1,10
        if(numprocs.ne.1)then 
            call matmulvocxs(glx,xx2,nrpos,bu,num,b_point-1,js1)
            call matmulvocxs(gly,xx2,nrpos,bv,num,b_point-1,js1)
        else
            call matmulvocxs(glx,xx2,nrpos,bu,num,b_point-1,js)
            call matmulvocxs(gly,xx2,nrpos,bv,num,b_point-1,js)
        endif
        
        ul=ustar+glx
        vl=vstar+gly
        
        if(numprocs.ne.1)then 
            call matmulvocxs(comj(:,1),xx1,mrpos,ul,b_point-1,num,js1)
            call matmulvocxs(comj(:,2),xx1,mrpos,vl,b_point-1,num,js1)
        else
            call matmulvocxs(comj(:,1),xx1,mrpos,ul,b_point-1,num,js)
            call matmulvocxs(comj(:,2),xx1,mrpos,vl,b_point-1,num,js)
        endif

        bu=bu+ub-comj(:,1)
        bv=bv+vb-comj(:,2)
    enddo
    
    !print*,maxval(comj(:,1))
    !pause
    
    spp=myid*int((adjoint_n)/numprocs)
	epp=(myid+1)*int((adjoint_n)/numprocs)-1
    if(myid.eq.0) spp=1
    if(myid.eq.numprocs-1) epp=adjoint_n

    do i=spp,epp
        sum1=0.0d0
        sum2=0.0d0
        if(i.lt.adjoint_n)then
            j=nrpos(i+1)
        else
            if(numprocs.ne.1)then
                j=js1+1
            else
                j=js+1
            endif
!            if(i.eq.1751)print*,j,nrpos(i),js1,myid
        endif

        if(numprocs.ne.1)then
            rhos1(i)=2.0d0*rhos(i)*glx(i)
            rhos2(i)=2.0d0*rhos(i)*gly(i)
        else
            gforce(adjoint_p(i,1),adjoint_p(i,2),1)=2.0d0*mesh%rho(adjoint_p(i,1),adjoint_p(i,2))*glx(i)
            gforce(adjoint_p(i,1),adjoint_p(i,2),2)=2.0d0*mesh%rho(adjoint_p(i,1),adjoint_p(i,2))*gly(i)
        endif
    enddo

    if(numprocs.ne.1)then
        do j=1,numprocs
            if(j.eq.1) then
                recvcounts(j)=int((adjoint_n)/numprocs)-1
                displs(j)=0
            elseif(j.eq.numprocs) then
                recvcounts(j)=adjoint_n-(numprocs-1)*int((adjoint_n)/numprocs)+1
                displs(j)=displs(j-1)+recvcounts(j-1)
            else
                recvcounts(j)=int((adjoint_n)/numprocs)
                displs(j)=displs(j-1)+recvcounts(j-1)
            endif
        enddo

        call MPI_ALLGATHERV(rhos1(spp:epp),recvcounts(myid+1),MPI_DOUBLE_PRECISION,ustar1,recvcounts,displs,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
        call MPI_ALLGATHERV(rhos2(spp:epp),recvcounts(myid+1),MPI_DOUBLE_PRECISION,vstar1,recvcounts,displs,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

        do i=1,adjoint_n 
            if((adjoint_p(i,1).ge.sp).and.(adjoint_p(i,1).le.ep)) then
                gforce(adjoint_p(i,1),adjoint_p(i,2),1)=ustar1(i)
                gforce(adjoint_p(i,1),adjoint_p(i,2),2)=vstar1(i)
            endif               
        enddo
    endif 
    
    if(myid.eq.0)then
		force=0.0d0
        mom=0.0d0
        if(numprocs.ne.1)then
            force(1)=-sum(ustar1)
            force(2)=-sum(vstar1)
            do k=1,adjoint_n
                x0=(mesh%ox+mesh%scale*dble(adjoint_p(k,1)-1)-(xc+wing_initial_x_origin))/mesh%scale-body_rel_pos(1)
                y0=(mesh%oy+mesh%scale*dble(adjoint_p(k,2)-1)-(yc+wing_initial_y_origin))/mesh%scale
                mom=mom-x0*vstar1(k)+y0*ustar1(k)
            enddo
        else
            force(1)=-sum(gforce(:,:,1))
            force(2)=-sum(gforce(:,:,2))
            do k=1,adjoint_n
                x0=(mesh%ox+mesh%scale*dble(adjoint_p(k,1)-1)-(xc+wing_initial_x_origin))/mesh%scale-body_rel_pos(1)
                y0=(mesh%oy+mesh%scale*dble(adjoint_p(k,2)-1)-(yc+wing_initial_y_origin))/mesh%scale
                mom=mom-x0*gforce(adjoint_p(k,1),adjoint_p(k,2),2)+y0*gforce(adjoint_p(k,1),adjoint_p(k,2),1)
            enddo
        endif

        body_rel_a_prev(1)=body_rel_a(1)
        if(ifdynamic.eq.1)then
            if(closed_boundary.eq..false.)then
                body_rel_a(1)=force(1)/(solid_den/(mesh%scale*mesh%scale))
            else
                body_rel_a(1)=force(1)/(solid_den*area_airfoil/(mesh%scale*mesh%scale))
            endif
        else
            body_rel_a(1)=0.0d0
        endif
        body_rel_v_prev(1)=body_rel_v(1)
        body_rel_v(1)=body_rel_a(1)+body_rel_v(1)
        body_rel_s_prev(1)=body_rel_s(1)
        body_rel_s(1)=body_rel_s(1)+body_rel_v_prev(1)
        body_rel_pos(1)=body_rel_pos(1)+body_rel_v_prev(1)
        body_rel_theta_prev(1)=body_rel_theta(1)
        body_rel_temp_theta_prev(1)=body_rel_temp_theta(1)
        body_rel_theta(1)=body_rel_theta_prev(1)+body_rel_temp_theta_prev(1)
!        if(closed_boundary.eq..false.)then
        if(ifdynamic.eq.1)then
            body_rel_temp_theta(1)=body_rel_temp_theta_prev(1)+par1*UU*UU*mesh%scale*mesh%scale*(body_rel_theta(1)+alphad)+par2*(mesh%scale*mesh%scale*mesh%scale*mesh%scale)*mom
        else
            body_rel_temp_theta(1)=0.0d0
        endif
!        else
!            body_rel_temp_theta(1)=body_rel_temp_theta_prev(1)
!        endif
!        print*, force(1),(solid_den*area_airfoil/(mesh%scale*mesh%scale))
!        body_rel_a(1)=0.0d0
    endif    
    call MPI_BCAST(body_rel_a,wing_count,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(body_rel_a_prev,wing_count,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(body_rel_v,wing_count,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(body_rel_v_prev,wing_count,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(body_rel_s,wing_count,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(body_rel_s_prev,wing_count,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(body_rel_pos,wing_count,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(body_rel_theta_prev,wing_count,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(body_rel_theta,wing_count,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(body_rel_temp_theta_prev,wing_count,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(body_rel_temp_theta,wing_count,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call update_body_extra_from_rel()
     
   
    do k=1,adjoint_n
        if(numprocs.ne.1)then
           if((adjoint_p(k,1).ge.sp).and.(adjoint_p(k,1).le.ep)) then
                i=adjoint_p(k,1)
                j=adjoint_p(k,2)
                mesh%u(i,j)=ul(k)
                mesh%v(i,j)=vl(k)
            endif                           
        else
            i=adjoint_p(k,1)
            j=adjoint_p(k,2)
    !        if(i.eq.127.and.j.eq.132) print*,k
            mesh%u(i,j)=ul(k)
            mesh%v(i,j)=vl(k)
        endif
    enddo

end

subroutine immersed_boundary_Wu(mesh,h)
    use mpi
    use cublas
    use cusolverDn
	use com_date
    use meshtype
	implicit none
	
	type(mesh_dat) mesh
	integer sp,ep,ip,jp,i,j,k,k1,num,temn,smooth,spp,epp,wing,seg_start,point_start
    real(kind=ps),  allocatable ::  ub(:),vb(:),ustar(:),vstar(:),glx(:),gly(:),bu(:),bv(:),xx1(:,:),xx2(:,:),as(:,:),rho_adj(:)
    real(kind=ps) disr,x0,y0,kernel1,kernel2,temp1,temp2
    real*8,  allocatable :: Workspace(:)
    integer :: devIpiv(max_boundary_segments),devInfo,lwork,stat
    type(cusolverDnHandle) :: h
    logical if_rk4 !是否采用龙格库塔方法求解动力学方程
    
    if(myid.eq.0) then
    if_rk4=.true.
    smooth=1
    
	!$acc update host(adjoint_n)
    num=adjoint_n
!	allocate (ub(b_segment_count),vb(b_segment_count),ustar(num),vstar(num),bu(b_segment_count),bv(b_segment_count),xx1(b_segment_count,num),xx2(num,b_segment_count),as((b_segment_count),(b_segment_count)),glx(num),gly(num))
    allocate (ub(b_segment_count),vb(b_segment_count),ustar(num),vstar(num),bu(b_segment_count),bv(b_segment_count),xx1(b_segment_count,num),xx2(num,b_segment_count),as(b_segment_count,b_segment_count),glx(num),gly(num),rho_adj(num))

    !$acc enter data create (ub,vb,ustar,vstar,bu,bv,xx1,xx2,as,glx,gly,rho_adj,devIpiv,devInfo)

    !$acc update device (omega_arr,angle_arr,sau_arr,sbv_arr,body_extra_v)
    !$acc kernels present(xb0,yb0,ub,vb,omega_arr,angle_arr,sau_arr,sbv_arr,body_extra_v,wing_segment_start,wing_point_start)
    !$acc loop seq
    do wing=1,wing_count
        seg_start=wing_segment_start(wing)
        point_start=wing_point_start(wing)
        !$acc loop
        do i=seg_start,seg_start+boundary_segments_per_wing-1
            ub(i)=UU*(-omega_arr(wing)*sin(angle_arr(wing))*0.5d0*(xb0(point_start+i-seg_start+1)+xb0(point_start+i-seg_start)) &
                +omega_arr(wing)*cos(angle_arr(wing))*0.5d0*(yb0(point_start+i-seg_start+1)+yb0(point_start+i-seg_start))+sau_arr(wing))+body_extra_v(wing)
            vb(i)=UU*(-omega_arr(wing)*cos(angle_arr(wing))*0.5d0*(xb0(point_start+i-seg_start+1)+xb0(point_start+i-seg_start)) &
                -omega_arr(wing)*sin(angle_arr(wing))*0.5d0*(yb0(point_start+i-seg_start+1)+yb0(point_start+i-seg_start))+sbv_arr(wing))
        enddo
        !$acc end loop
    enddo
    !$acc end loop
    !$acc end kernels

    call cpu_time(time_s5)
    !$acc kernels loop copyin(smooth) present(mesh,xx1,xx2,adjoint_p,arc_l,xb2,yb2) create(kernel1,kernel2,x0,y0)
    do k=1,b_segment_count
        do i=1,adjoint_n
            x0=mesh%ox+mesh%scale*dble(adjoint_p(i,1)-1)
            y0=mesh%oy+mesh%scale*dble(adjoint_p(i,2)-1)
            if(smooth.eq.0) then
                disr=abs((x0-xb2(k))/mesh%scale)
                if(abs(disr).lt.1.0d0)then
                    kernel1=0.125d0*(3.0d0-2.0d0*abs(disr)+sqrt(1.0d0+4.0d0*abs(disr)-4.0d0*disr*disr))
                elseif ((abs(disr).lt.2.0d0).and.(abs(disr).ge.1.0d0))then
                    kernel1=0.125d0*(5.0d0-2.0d0*abs(disr)-sqrt(-7.0d0+12.0d0*abs(disr)-4.0d0*disr*disr))
                else
                    kernel1=0.0d0
                endif
                disr=abs((y0-yb2(k))/mesh%scale)
                if(abs(disr).lt.1.0d0)then
                    kernel2=0.125d0*(3.0d0-2.0d0*abs(disr)+sqrt(1.0d0+4.0d0*abs(disr)-4.0d0*disr*disr))
                elseif ((abs(disr).lt.2.0d0).and.(abs(disr).ge.1.0d0))then
                    kernel2=0.125d0*(5.0d0-2.0d0*abs(disr)-sqrt(-7.0d0+12.0d0*abs(disr)-4.0d0*disr*disr))
                else
                    kernel2=0.0d0
                endif
            else
                disr=abs((x0-xb2(k))/mesh%scale)
                if(abs(disr).le.0.5d0)then
                    kernel1=0.375d0+0.03125d0*pai-0.25d0*disr*disr
                elseif ((abs(disr).le.1.5d0).and.(abs(disr).ge.0.5d0)) then
                    kernel1=0.25d0+0.125d0*(1.0d0-disr)*sqrt(-2.0d0+8.0d0*disr-4.0d0*disr*disr)-0.125*asin(sqrt(2.0d0)*(disr-1.0d0))
                elseif ((abs(disr).le.2.5d0).and.(abs(disr).ge.1.5d0)) then
                    kernel1=1.0625d0-0.015625d0*pai-0.75d0*disr+0.125*disr*disr+0.0625d0*(disr-2.0d0)*sqrt(-14.0d0+16.0d0*disr-4.0d0*disr*disr)+0.0625d0*asin(sqrt(2.0d0)*(disr-2.0d0))
                else
                    kernel1=0.0d0
                endif
                disr=abs((y0-yb2(k))/mesh%scale)
                if(abs(disr).le.0.5d0)then
                    kernel2=0.375d0+0.03125d0*pai-0.25d0*disr*disr
                elseif ((abs(disr).le.1.5d0).and.(abs(disr).ge.0.5d0)) then
                    kernel2=0.25d0+0.125d0*(1.0d0-disr)*sqrt(-2.0d0+8.0d0*disr-4.0d0*disr*disr)-0.125*asin(sqrt(2.0d0)*(disr-1.0d0))
                elseif ((abs(disr).le.2.5d0).and.(abs(disr).ge.1.5d0)) then
                    kernel2=1.0625d0-0.015625d0*pai-0.75d0*disr+0.125*disr*disr+0.0625d0*(disr-2.0d0)*sqrt(-14.0d0+16.0d0*disr-4.0d0*disr*disr)+0.0625d0*asin(sqrt(2.0d0)*(disr-2.0d0))
                else
                    kernel2=0.0d0
                endif
            endif
            xx1(k,i)=kernel1*kernel2
            xx2(i,k)=kernel1*kernel2*arc_l(k)
        enddo
    enddo
    !$acc end kernels loop
    call cpu_time(time_e5)

    !$acc kernels loop present(mesh,adjoint_p,ustar,vstar)
    do k=1,adjoint_n
        i=adjoint_p(k,1)
        j=adjoint_p(k,2)
        ustar(k)=mesh%u(i,j)
        vstar(k)=mesh%v(i,j)
    enddo
    !$acc end kernels loop

    !$acc host_data use_device (xx1,xx2,as,bu,bv,ustar,vstar)
    call cublasDgemm('N','N',b_segment_count,b_segment_count,adjoint_n,1.0d0,xx1,b_segment_count,xx2,adjoint_n,0.0d0,as,b_segment_count)
    call cublasDgemv('N',b_segment_count,adjoint_n,1.0d0,xx1,b_segment_count,ustar,1,0.0d0,bu,1)
    call cublasDgemv('N',b_segment_count,adjoint_n,1.0d0,xx1,b_segment_count,vstar,1,0.0d0,bv,1)
    !$acc end host_data

    !$acc kernels loop present(bu,ub,bv,vb,comj)
    do i=1,b_segment_count
        comj(i,1)=ub(i)-bu(i)   
        comj(i,2)=vb(i)-bv(i)
    enddo
    !$acc end kernels loop
	if (time.lt.0.5d0*mesh%scale*UU) then
		!$acc update self(as,comj(1:b_segment_count,1:2))
		write(*,*) 'IBM DEBUG A(min,max)=', minval(as), maxval(as)
		write(*,*) 'IBM DEBUG rhs-x(min,max)=', minval(comj(1:b_segment_count,1)), maxval(comj(1:b_segment_count,1))
		write(*,*) 'IBM DEBUG rhs-y(min,max)=', minval(comj(1:b_segment_count,2)), maxval(comj(1:b_segment_count,2))
	endif
    call cpu_time(time_s6)
    !$acc host_data use_device (as)
    stat=cusolverDnDgetrf_bufferSize(h, b_segment_count, b_segment_count, as, b_segment_count,lwork)
    !$acc end host_data
	if ((stat.ne.0).or.(lwork.le.0)) then
		write(*,*) 'cusolverDnDgetrf_bufferSize failed: status=',stat,' lwork=',lwork
		call MPI_ABORT(MPI_COMM_WORLD,904,ierr)
	endif
    allocate (Workspace(lwork))
    !$acc enter data create (Workspace)
    !$acc host_data use_device (as,Workspace,devIpiv,devInfo)
    stat=cusolverDnDgetrf(h,b_segment_count,b_segment_count,as,b_segment_count,Workspace,devIpiv,devInfo)
	!$acc end host_data
	!$acc update self(devInfo)
	if ((stat.ne.0).or.(devInfo.ne.0)) then
		write(*,*) 'cusolverDnDgetrf failed: status=',stat,' devInfo=',devInfo
		call MPI_ABORT(MPI_COMM_WORLD,901,ierr)
	endif
	!$acc host_data use_device (as,devIpiv,comj,devInfo)
    stat=cusolverDnDgetrs(h,CUBLAS_OP_N,b_segment_count,2,as,b_segment_count,devIpiv,comj,max_boundary_points,devInfo)
	!$acc end host_data
	!$acc update self(devInfo)
	if ((stat.ne.0).or.(devInfo.ne.0)) then
		write(*,*) 'cusolverDnDgetrs failed: status=',stat,' devInfo=',devInfo
		call MPI_ABORT(MPI_COMM_WORLD,902,ierr)
	endif
	if (time.lt.0.5d0*mesh%scale*UU) then
		!$acc update self(comj(1:b_segment_count,1:2))
		write(*,*) 'IBM DEBUG solution-x(min,max)=', minval(comj(1:b_segment_count,1)), maxval(comj(1:b_segment_count,1))
		write(*,*) 'IBM DEBUG solution-y(min,max)=', minval(comj(1:b_segment_count,2)), maxval(comj(1:b_segment_count,2))
	endif
	!$acc host_data use_device (xx2,comj,glx,gly)
    call cublasDgemv('N',adjoint_n,b_segment_count,1.0d0,xx2,adjoint_n,comj(:,1),1,0.0d0,glx,1)
    call cublasDgemv('N',adjoint_n,b_segment_count,1.0d0,xx2,adjoint_n,comj(:,2),1,0.0d0,gly,1)
    !$acc end host_data
    call cpu_time(time_e6)

    call cpu_time(time_s7)
    !$acc kernels present(mesh,ustar,vstar,adjoint_p,glx,gly,gforce,rho_adj) create(temp1,temp2,i,j)
    !$acc loop independent
    do k=1, adjoint_n
        i=adjoint_p(k,1)
        j=adjoint_p(k,2)
        mesh%u(i,j)=ustar(k)+glx(k)
        mesh%v(i,j)=vstar(k)+gly(k)
        temp1=2.0d0*mesh%rho(i,j)*glx(k)
        temp2=2.0d0*mesh%rho(i,j)*gly(k)
        rho_adj(k)=mesh%rho(i,j)
        gforce(i,j,1)=temp1
        gforce(i,j,2)=temp2
    enddo
    !$acc end loop
    !$acc end kernels
    call cpu_time(time_e7)
    !$acc update host(glx(1:adjoint_n),gly(1:adjoint_n),rho_adj(1:adjoint_n),adjoint_p(1:adjoint_n,:),adjoint_wing_id(1:adjoint_n))
    wing_force=0.0d0
    wing_moment=0.0d0
    do k=1,adjoint_n
        wing=adjoint_wing_id(k)
        temp1=2.0d0*rho_adj(k)*glx(k)
        temp2=2.0d0*rho_adj(k)*gly(k)
        wing_force(wing,1)=wing_force(wing,1)-temp1
        wing_force(wing,2)=wing_force(wing,2)-temp2
        i=adjoint_p(k,1)
        j=adjoint_p(k,2)
        x0=(mesh%ox+mesh%scale*dble(i-1)-(xc+wing_initial_x_origin+dble(wing-1)*wing_initial_x_spacing+sa_arr(wing)))/mesh%scale-body_extra_pos(wing)
        y0=(mesh%oy+mesh%scale*dble(j-1)-(yc+wing_initial_y_origin+sb_arr(wing)))/mesh%scale
        wing_moment(wing)=wing_moment(wing)-x0*temp2+y0*temp1
    enddo
    force=0.0d0
    do wing=1,wing_count
        force(1)=force(1)+wing_force(wing,1)
        force(2)=force(2)+wing_force(wing,2)
    enddo
    mom=sum(wing_moment(1:wing_count))
    force1=0.0d0
    force2=0.0d0
    mom2=0.0d0
    if (wing_count.ge.1) then
        force1=wing_force(1,:)
        mom=wing_moment(1)
    endif
    if (wing_count.ge.2) then
        force2=wing_force(2,:)
        mom2=wing_moment(2)
    endif

    if(myid.eq.0)then
        if(solve_dynamic_equations.eq..true.)then !动力学方程求解
            body_rel_a_prev(1)=body_rel_a(1)
            if(closed_boundary.eq..false.)then
                body_rel_a(1)=force(1)/(solid_den/(mesh%scale*mesh%scale))
            else
                body_rel_a(1)=force(1)/(solid_den*area_airfoil/(mesh%scale*mesh%scale))
            endif
            body_rel_v_prev(1)=body_rel_v(1)
            body_rel_s_prev(1)=body_rel_s(1)
            body_rel_theta_prev(1)=body_rel_theta(1)
            body_rel_temp_theta_prev(1)=body_rel_temp_theta(1)
            if(if_rk4.eq..false.)then
                body_rel_v(1)=body_rel_a(1)+body_rel_v(1)
                body_rel_s(1)=body_rel_s(1)+body_rel_v_prev(1)
                body_rel_pos(1)=body_rel_pos(1)+body_rel_v_prev(1)
                body_rel_theta(1)=body_rel_theta_prev(1)+body_rel_temp_theta_prev(1)
                body_rel_temp_theta(1)=body_rel_temp_theta_prev(1)+par1*UU*UU*mesh%scale*mesh%scale*(body_rel_theta(1)+alphad)+ &
                    par2*(mesh%scale*mesh%scale*mesh%scale*mesh%scale)*mom
            else
                call rk4(body_rel_s_prev(1),body_rel_s(1),body_rel_pos(1),body_rel_v_prev(1),body_rel_v(1),body_rel_a(1), &
                    body_rel_theta_prev(1),body_rel_theta(1),body_rel_temp_theta_prev(1),body_rel_temp_theta(1), &
                    par1*UU*UU*mesh%scale*mesh%scale,par2*(mesh%scale*mesh%scale*mesh%scale*mesh%scale)) !用四步四级龙格-库塔方法求解动力学方程
            endif
        endif
    endif 
    endif ! rank 0 owns the coupled IBM solve
    call MPI_BCAST(body_rel_a,wing_count,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(body_rel_a_prev,wing_count,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(body_rel_v,wing_count,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(body_rel_v_prev,wing_count,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(body_rel_s,wing_count,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(body_rel_s_prev,wing_count,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(body_rel_pos,wing_count,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(body_rel_theta_prev,wing_count,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(body_rel_theta,wing_count,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(body_rel_temp_theta_prev,wing_count,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(body_rel_temp_theta,wing_count,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call update_body_extra_from_rel()
    !$acc update device (body_rel_a,body_rel_a_prev,body_rel_v,body_rel_v_prev,body_rel_s,body_rel_s_prev,body_rel_pos)
    !$acc update device (body_rel_theta,body_rel_theta_prev,body_rel_temp_theta,body_rel_temp_theta_prev)
    !$acc update device (body_extra_v,body_extra_s,body_extra_pos)
    if(myid.eq.0) then
    !$acc exit data delete (ub,vb,ustar,vstar,bu,bv,xx1,xx2,as,glx,gly,rho_adj,devIpiv,devInfo,Workspace)
    endif
end

subroutine matmulvoc(c,a,b,num1,num2)
    !$acc routine seq
	implicit none
	integer num1,num2
	real*8 a(num1,num2),b(num2),c(num1)
	integer i,j,k
	real*8 sum

!    call dgemv('N',num1,num2,alpha,a,num1,b,1,bate,c,1)
	do i=1,num1
		sum=0.0d0
		do j=1,num2
			sum=sum+a(i,j)*b(j)
		enddo
		c(i)=sum
	enddo

end
    
subroutine fasttransposematrix(xx3,xx1,mu,nu,js,cpot)
	use com_date
    implicit none
    
    integer mu,nu,js,p,col,qq
    real (kind=ps) xx3(mu*nu,3),xx1(mu*nu,3) 
    integer cpot1(mu),cpot(mu)
    
    cpot1=cpot
    do p=1,js
        col=nint(xx3(p,2))
        qq=cpot1(col)
        xx1(qq,1)=xx3(p,2)
        xx1(qq,2)=xx3(p,1)
        xx1(qq,3)=xx3(p,3)
        cpot1(col)=cpot1(col)+1
    enddo
end
    
subroutine multsmatrix(xx1,xx2,mrpos,nrpos,mu,nu,js,as,ia,ja,ds)
    use mpi
	use com_date
    implicit none
    integer mu,nu,js,arow,brow,tp,p,t,qq,ccol,ds,ds1,i,j,sp,ep
    integer mrpos(mu),nrpos(nu),iia(mu),ia(mu+1),ja(mu*mu)
    real (kind=ps) xx1(mu*nu,3),xx2(mu*nu,3),as(mu*mu),ctemp(mu)
    integer,  allocatable :: buf1(:),buf2(:),recvcounts(:),displs(:)
    real (kind=ps),allocatable :: buf3(:) 
    
    sp=myid*int(mu/numprocs)
	ep=(myid+1)*int(mu/numprocs)-1
    if(myid.eq.0) sp=1
    if(myid.eq.numprocs-1) ep=mu
    
    allocate (buf1(sp:ep),buf2(mu*mu),recvcounts(numprocs),displs(numprocs)) 
    allocate (buf3(mu*mu))
    buf1=0

    ctemp=0.0d0
    iia=0
    ia=0
    ja=0
    ds1=1
    do arow=sp,ep
        if(arow.lt.mu)then
            tp=mrpos(arow+1)
        else
            tp=js+1
        endif
        do p=mrpos(arow),tp-1
            brow=nint(xx1(p,2))
            if(brow.lt.nu)then
                t=nrpos(brow+1)
            else
                t=js+1
            endif
            do qq=nrpos(brow),t-1
                ccol=nint(xx2(qq,2))
                ctemp(ccol)=ctemp(ccol)+xx1(p,3)*xx2(qq,3)
            enddo
        enddo
        do ccol=1,mu
            if(ctemp(ccol).ne.0.0d0)then
                if(numprocs.ne.1)then
                    buf1(arow)= buf1(arow)+1 
                    buf2(ds1)=ccol
                    buf3(ds1)=ctemp(ccol)
                    ctemp(ccol)=0.0d0
                    ds1=ds1+1
                else
                    iia(arow)= iia(arow)+1 
                    ja(ds1)=ccol
                    as(ds1)=ctemp(ccol)
                    ctemp(ccol)=0.0d0
                    ds1=ds1+1
                endif
            endif
        enddo
    enddo
    ds1=ds1-1

    if(numprocs.ne.1)then
        call MPI_ALLREDUCE(ds1,ds,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
        displs(1)=0
        do j=1,numprocs-1
            if(j.eq.1)then
                recvcounts(j)=int(mu/numprocs)-1
            else
                recvcounts(j)=int(mu/numprocs)
            endif
            displs(j+1)=displs(j)+recvcounts(j)
        enddo
        recvcounts(numprocs)=mu-(numprocs-1)*int(mu/numprocs)+1

        call MPI_ALLGATHERV(buf1(sp:ep),recvcounts(myid+1),MPI_INTEGER,iia,recvcounts,displs,MPI_INTEGER,MPI_COMM_WORLD,ierr)
        call MPI_ALLGATHER(ds1,1,MPI_INTEGER,recvcounts,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
        displs(1)=0
        do j=1,numprocs-1
            displs(j+1)=displs(j)+recvcounts(j)
        enddo
        call MPI_ALLGATHERV(buf2,recvcounts(myid+1),MPI_INTEGER,ja,recvcounts,displs,MPI_INTEGER,MPI_COMM_WORLD,ierr)
        call MPI_ALLGATHERV(buf3,recvcounts(myid+1),MPI_DOUBLE_PRECISION,as,recvcounts,displs,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
    else
        ds=ds1
    endif
    ia(1)=1
    do i=1,mu-1
        ia(i+1)=ia(i)+iia(i)
    enddo
    ia(mu+1)=ds+1
end
    
subroutine matmulvocxs(bb,xx1,mrpos,xx,mu,nu,js)
    use mpi
	use com_date
    implicit none
    integer mu,nu,js,arow,brow,tp,p,t,qq,ccol,j,sp,ep
    integer mrpos(mu)
    real (kind=ps) xx1(mu*nu,3),xx(nu),bb(mu)
    integer,  allocatable :: recvcounts(:),displs(:)
    real (kind=ps),allocatable :: buf(:) 
    
    sp=myid*int(mu/numprocs)
	ep=(myid+1)*int(mu/numprocs)-1
    if(myid.eq.0) sp=1
    if(myid.eq.numprocs-1) ep=mu
    allocate (recvcounts(numprocs),displs(numprocs)) 
    allocate (buf(sp:ep))
    bb=0.0d0
    buf=0.0d0
    do arow=sp,ep
        if(arow.lt.mu)then
            tp=mrpos(arow+1)
        else
            tp=js+1
        endif
        do p=mrpos(arow),tp-1
            brow=nint(xx1(p,2))
            if(numprocs.ne.1)then
                buf(arow)=buf(arow)+xx1(p,3)*xx(brow) 
            else
                bb(arow)=bb(arow)+xx1(p,3)*xx(brow)
            endif
        enddo
    enddo

    if(numprocs.ne.1)then
        displs(1)=0
        do j=1,numprocs-1
            if(j.eq.1)then
                recvcounts(j)=int(mu/numprocs)-1
            else
                recvcounts(j)=int(mu/numprocs)
            endif
            displs(j+1)=displs(j)+recvcounts(j)
        enddo
        recvcounts(numprocs)=mu-(numprocs-1)*int(mu/numprocs)+1
        call MPI_ALLGATHERV(buf(sp:ep),recvcounts(myid+1),MPI_DOUBLE_PRECISION,bb,recvcounts,displs,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
    endif
!    print*, buf(sp+5),bb(sp+5),myid
end

subroutine GMRES(x,a,ia,ja,ds,b,n)
    implicit none
    
    integer i,j,n,k,ds,m,l,s
    integer ia(n+1),ja(n*n)
    real*8 a(n*n,3),b(n),r(n),e1(n),beta(n),x(n)
    real*8 b_norm,r_norm,err,temp
    real*8, allocatable :: qq(:,:),h(:,:),cs(:),sn(:),y(:)
!    m=n-1
    m=20
    l=5
    allocate (qq(n,m),h(m+1,m+1),cs(m),sn(m),y(m)) 

!    x=1.0d0
 
    do s=1,l
        !call matmulvocxs(r,a,ia,x,n,n,ds)
        call matvol(a,x,r,n,ia,ja,ds) 
        r=b-r
        temp=maxval(abs(r))
        if(temp.le.1d-20) exit
        call norm(b,b_norm,n)
        call norm(r,r_norm,n)
        err=r_norm/b_norm
    
        e1=0.0d0
        e1(1)=1.0d0
        qq(:,1)=r/r_norm
        beta=r_norm*e1
        h=0.0d0
        sn=0.0d0
        cs=0.0d0

        do k=1,m
            call arnoldi(a,qq,k,n,h(:,k),ia,ja,ds,m)
            call apply_givens_rotation(h(:,k), cs, sn, k,n,m)
            beta(k+1) = -sn(k)*beta(k)
            beta(k)   = cs(k)*beta(k)
            err  = abs(beta(k+1)) / b_norm
 !           if(err.lt.1d-3) then
 !               print*, 'here!'
 !               exit
 !           endif
        enddo
!        open(2001,file='hhg.dat',status='unknown')    
!        do i=1,m+1
!            write(2001,'(21f14.8)') h(i,1:m+1)
!        enddo
!        pause

        y(k-1)=beta(k-1)/h(k-1,k-1)
        
        do i=k-2,1,-1
            temp=0.0d0
            do j=i+1,k-1
                temp=temp+h(i,j)*y(j)
            enddo
            y(i)=(beta(i)- temp)/h(i,i)
        enddo
        !y = h(1:k,1:k) \ beta(1:k)    
        !x = x + Q(:,1:k)*y
        do i=1,n
            temp=0.0d0
            do j=1,k-1
                temp=temp+qq(i,j)*y(j)
            enddo
            x(i)=x(i)+temp
        enddo
    enddo
!    call matmulvocxs(r,a,ia,x,n,n,ds)
!    print*,maxval(abs(b-r))
end
    
subroutine norm(a,a_norm,n)
    implicit none
    integer n,i,j
    real*8 a(n),tem,a_norm
    
    tem=0.0d0
    do i=1,n
        tem=tem+a(i)*a(i)
    enddo   
    a_norm=sqrt(tem)
end  
    
subroutine arnoldi(a, qq, k,n,h,ia,ja,ds,m)
    implicit none 
    integer n,k,i,j,ds,m
    integer ia(n+1),ja(n*n)
    real*8 a(n*n,3),qq(n,m),q1(n),h(m+1)
  
    call matvol(a,qq(:,k),q1,n,ia,ja,ds)
!    call matmulvocxs(q1,a,ia,qq(:,k),n,n,ds)
    do i = 1,k
        h(i)=0.0d0
        do j=1,n
            h(i)= h(i)+q1(j)*qq(j,i)
        enddo
        q1 = q1 - h(i)*qq(:,i)
    enddo
    call norm(q1,h(k+1),n)
    q1 = q1 / h(k+1)
    if (k.ne.m) qq(:,k+1)=q1
end
    
subroutine apply_givens_rotation(h, cs, sn, k,n,m)
    implicit none
    integer i,k,n,m
    real*8 cs(m),sn(m),temp,cs_k,sn_k,h(m+1)

    !apply for ith column
    do i = 1,k-1                              
        temp   =  cs(i)*h(i) + sn(i)*h(i+1)
        h(i+1) = -sn(i)*h(i) + cs(i)*h(i+1)
        h(i)   = temp
    enddo
  
    call givens_rotation(h(k), h(k+1),cs_k,sn_k)
    cs(k)=cs_k
    sn(k)=sn_k
  
    !eliminate H(i+1,i)
    h(k) = cs_k*h(k) + sn_k*h(k+1)
    h(k+1) = 0.0d0
end
    
subroutine givens_rotation (v1, v2,cs_k,sn_k)
    implicit none
    real*8 v1,v2,cs_k,sn_k,t
    
    if (v1==0)then
        cs_k = 0.0d0
        sn_k = 1.0d0
    else
        t=sqrt(v1*v1+v2*v2)
        cs_k = abs(v1) / t
        sn_k = cs_k * v2 / v1
    endif    
end

subroutine matvol(xx1,xx,bb,mu,ia,ja,js)
    use mpi
	use com_date
    implicit none
    integer mu,nu,js,arow,brow,tp,p,t,qq,ccol,j,sp,ep
    integer mrpos(mu)
    integer ia(mu+1),ja(mu*mu)
    real (kind=ps) xx1(mu*mu),xx(mu),bb(mu)
    integer,  allocatable :: recvcounts(:),displs(:)
    real (kind=ps),allocatable :: buf(:) 
    
    sp=myid*int(mu/numprocs)
	ep=(myid+1)*int(mu/numprocs)-1
    if(myid.eq.0) sp=1
    if(myid.eq.numprocs-1) ep=mu
    allocate (recvcounts(numprocs),displs(numprocs)) 
    allocate (buf(sp:ep))

    mrpos=ia(1:mu-1)
    bb=0.0d0
    buf=0.0d0
    do arow=sp,ep
        if(arow.lt.mu)then
            tp=mrpos(arow+1)
        else
            tp=js+1
        endif
        do p=mrpos(arow),tp-1
            brow=ja(p)
            if(numprocs.ne.1)then
                buf(arow)=buf(arow)+xx1(p)*xx(brow) 
            else
                bb(arow)=bb(arow)+xx1(p)*xx(brow)
            endif
        enddo
    enddo
    if(numprocs.ne.1)then
        displs(1)=0
        do j=1,numprocs-1
            if(j.eq.1)then
                recvcounts(j)=int(mu/numprocs)-1
            else
                recvcounts(j)=int(mu/numprocs)
            endif
            displs(j+1)=displs(j)+recvcounts(j)
        enddo
        recvcounts(numprocs)=mu-(numprocs-1)*int(mu/numprocs)+1
        call MPI_ALLGATHERV(buf(sp:ep),recvcounts(myid+1),MPI_DOUBLE_PRECISION,bb,recvcounts,displs,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
    endif
!    print*, buf(sp+5),bb(sp+5),myid
end

subroutine immersed_boundary_Wu1(mesh)
    use mpi
	use com_date
    use meshtype
	implicit none
	
	type(mesh_dat) mesh
	integer sp,ep,ip,jp,i,j,k,k1,num,temn,smooth,spp,epp
	real(kind=ps),  allocatable ::  ub(:),vb(:),ustar(:),vstar(:),rhos(:),ustar1(:),vstar1(:),rhos1(:),rhos2(:),buf1(:),buf2(:),as(:),ul(:),vl(:),glx(:),gly(:)
	real(kind=ps),  allocatable :: tem1(:,:),tem2(:,:),aa(:,:),aa1(:,:),bu(:),bv(:),tem11(:,:),tem22(:,:),tema(:,:),comi(:),x(:,:),xx1(:,:),xx2(:,:),xx3(:,:)
	real*8 disr,x0,y0,kernel1,kernel2,sum1,sum2
    integer nono,mn,nn,js,js1,ds,nnn
	integer,  allocatable :: non_zero(:),indexnum(:,:),recvcounts(:),displs(:),mrpos(:),nrpos(:),mrow(:),mcol(:),mrrow(:),mrcol(:),ia(:),ja(:)
    logical if_rk4 !是否采用龙格库塔方法求解动力学方程
    
    if_rk4=.true.
    smooth=1
    sp=myid*int((b_point-1)/numprocs)
	ep=(myid+1)*int((b_point-1)/numprocs)-1
    if(myid.eq.0) sp=1
    if(myid.eq.numprocs-1) ep=b_point-1
    
	!$acc update host(adjoint_n)
    num=adjoint_n
    nnn=1
	allocate (ub(b_point-1),vb(b_point-1),ustar(num),vstar(num),rhos(num),ustar1(num),vstar1(num),rhos1(num),rhos2(num),tema(b_point-1,b_point-1),comi(b_point-1)) !,buf1(b_point-1))
	allocate (tem1(b_point-1,num),tem2(num,b_point-1),aa(b_point-1,b_point-1),aa1(b_point-1,b_point-1),bu(b_point-1),bv(b_point-1),x(b_point-1,2),xx1((b_point-1)*num,3),xx2((b_point-1)*num,3),xx3((b_point-1)*num,3),as((b_point-1)*(b_point-1)))
	allocate (non_zero(num),indexnum(num,2),recvcounts(numprocs),displs(numprocs),mrpos(b_point-1),nrpos(num),mcol(num),mrrow(b_point-1),mrow(b_point-1),mrcol(num),ia(b_point),ja((b_point-1)*(b_point-1)),ul(num),vl(num),glx(num),gly(num))
    !$acc enter data create (js,js1,ds,ub,vb,ustar,vstar,rhos,ustar1,vstar1,rhos1,rhos2,tema,comi,tem1,tem2,aa,aa1,bu,bv,x,xx1,xx2,xx3,as,non_zero,indexnum,recvcounts,displs,mrpos,nrpos,mcol,mrrow,mrow,mrcol,ia,ja,ul,vl,glx,gly) copyin (num)
!    angle_arr(1)=-body_rel_theta(1)
!    omega_arr(1)=-body_rel_temp_theta(1)/(mesh%scale*UU)

    !$acc update device (omega_arr,angle_arr,sau_arr,sbv_arr,body_rel_v) 
    !$acc kernels loop present(xb0,yb0,ub,vb,omega_arr,angle_arr,sau_arr,sbv_arr,body_rel_v)
    do i=1,b_point-1
        ub(i)=UU*(-omega_arr(1)*sin(angle_arr(1))*0.5d0*(xb0(i+1)+xb0(i))+omega_arr(1)*cos(angle_arr(1))*0.5d0*(yb0(i+1)+yb0(i))+sau_arr(1))+body_rel_v(1)
        vb(i)=UU*(-omega_arr(1)*cos(angle_arr(1))*0.5d0*(xb0(i+1)+xb0(i))-omega_arr(1)*sin(angle_arr(1))*0.5d0*(yb0(i+1)+yb0(i))+sbv_arr(1))
    enddo
    !$acc end kernels loop

    !ub=ub+body_rel_v(1) !+body_rel_a(1)
    call cpu_time(time_s5)
    !$acc kernels copyin(smooth) present(mesh,xx1,xx3,js,mcol,mrow,ustar1,vstar1,rhos1)    
    js=0
    mcol=0
    mrow=0 
    ustar1=0.0d0
    vstar1=0.0d0
    rhos1=0.0d0

    !$acc loop
    do k=sp,ep
        do i=1,adjoint_n
            nono=1
            x0=mesh%ox+mesh%scale*dble(adjoint_p(i,1)-1)
            y0=mesh%oy+mesh%scale*dble(adjoint_p(i,2)-1)
            if(smooth.eq.0) then
                disr=abs((x0-xb2(k))/mesh%scale)
                if(abs(disr).lt.1.0d0)then
                    kernel1=0.125d0*(3.0d0-2.0d0*abs(disr)+sqrt(1.0d0+4.0d0*abs(disr)-4.0d0*disr*disr))
                elseif ((abs(disr).lt.2.0d0).and.(abs(disr).ge.1.0d0))then
                    kernel1=0.125d0*(5.0d0-2.0d0*abs(disr)-sqrt(-7.0d0+12.0d0*abs(disr)-4.0d0*disr*disr))
                else
                    kernel1=0.0d0
                    nono=0
                endif
                disr=abs((y0-yb2(k))/mesh%scale)
                if(abs(disr).lt.1.0d0)then
                    kernel2=0.125d0*(3.0d0-2.0d0*abs(disr)+sqrt(1.0d0+4.0d0*abs(disr)-4.0d0*disr*disr))
                elseif ((abs(disr).lt.2.0d0).and.(abs(disr).ge.1.0d0))then
                    kernel2=0.125d0*(5.0d0-2.0d0*abs(disr)-sqrt(-7.0d0+12.0d0*abs(disr)-4.0d0*disr*disr))
                else
                    kernel2=0.0d0
                    nono=0
                endif
            else
                disr=abs((x0-xb2(k))/mesh%scale)
                if(abs(disr).le.0.5d0)then
                    kernel1=0.375d0+0.03125d0*pai-0.25d0*disr*disr
                elseif ((abs(disr).le.1.5d0).and.(abs(disr).ge.0.5d0)) then
                    kernel1=0.25d0+0.125d0*(1.0d0-disr)*sqrt(-2.0d0+8.0d0*disr-4.0d0*disr*disr)-0.125*asin(sqrt(2.0d0)*(disr-1.0d0))
                elseif ((abs(disr).le.2.5d0).and.(abs(disr).ge.1.5d0)) then
                    kernel1=1.0625d0-0.015625d0*pai-0.75d0*disr+0.125*disr*disr+0.0625d0*(disr-2.0d0)*sqrt(-14.0d0+16.0d0*disr-4.0d0*disr*disr)+0.0625d0*asin(sqrt(2.0d0)*(disr-2.0d0))
                else
                    kernel1=0.0d0
                    nono=0
                endif
                disr=abs((y0-yb2(k))/mesh%scale)
                if(abs(disr).le.0.5d0)then
                    kernel2=0.375d0+0.03125d0*pai-0.25d0*disr*disr
                elseif ((abs(disr).le.1.5d0).and.(abs(disr).ge.0.5d0)) then
                    kernel2=0.25d0+0.125d0*(1.0d0-disr)*sqrt(-2.0d0+8.0d0*disr-4.0d0*disr*disr)-0.125*asin(sqrt(2.0d0)*(disr-1.0d0))
                elseif ((abs(disr).le.2.5d0).and.(abs(disr).ge.1.5d0)) then
                    kernel2=1.0625d0-0.015625d0*pai-0.75d0*disr+0.125*disr*disr+0.0625d0*(disr-2.0d0)*sqrt(-14.0d0+16.0d0*disr-4.0d0*disr*disr)+0.0625d0*asin(sqrt(2.0d0)*(disr-2.0d0))
                else
                    kernel2=0.0d0
                    nono=0
                endif
            endif
            if(nono.eq.1) then
                xx3(js+1,1)=dble(k)
                xx3(js+1,2)=dble(i)
                xx3(js+1,3)=kernel1*kernel2*arc_l(k)
                xx1(js+1,1)=dble(k)
                xx1(js+1,2)=dble(i)
                xx1(js+1,3)=(kernel1*kernel2)
                mrow(k)=mrow(k)+1
                mcol(i)=mcol(i)+1
                js=js+1
            endif
        enddo
    enddo
    !$acc end loop
    !js=js-1
    !$acc end kernels
    call cpu_time(time_e5)
    if(numprocs.ne.1) then
        call MPI_ALLREDUCE(mrow,mrrow,b_point-1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(mcol,mrcol,adjoint_n,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)

        mrow=mrrow
        mcol=mrcol

        call MPI_ALLGATHER(js,1,MPI_INTEGER,recvcounts,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
        js1=sum(recvcounts)

        allocate(buf1(3*js),buf2(3*js1))

        recvcounts=3*recvcounts
        displs=0
        do j=1,numprocs-1
            displs(j+1)=displs(j)+recvcounts(j)
        enddo
        j=1
        do i=1,js
            buf1(j)=xx1(i,1)
            buf1(j+1)=xx1(i,2)
            buf1(j+2)=xx1(i,3)
            j=j+3
        enddo
 
        call MPI_ALLGATHERV(buf1,recvcounts(myid+1),MPI_DOUBLE_PRECISION,buf2,recvcounts,displs,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
        j=1
        do i=1,js1
            xx1(i,1)=buf2(j)
            xx1(i,2)=buf2(j+1)
            xx1(i,3)=buf2(j+2)
            j=j+3
        enddo

        j=1
        do i=1,js
            buf1(j)=xx3(i,1)
            buf1(j+1)=xx3(i,2)
            buf1(j+2)=xx3(i,3)
            j=j+3
        enddo
        call MPI_ALLGATHERV(buf1,recvcounts(myid+1),MPI_DOUBLE_PRECISION,buf2,recvcounts,displs,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
        j=1
        do i=1,js1
            xx3(i,1)=buf2(j)
            xx3(i,2)=buf2(j+1)
            xx3(i,3)=buf2(j+2)
            j=j+3
        enddo

        if(myid.eq.0) then
		    sp=mesh%start_row+1
		    ep=mesh%end_row-1
	    elseif(myid.eq.numprocs-1) then
		    sp=mesh%start_row+1
		    ep=mesh%end_row-1
	    else
		    sp=mesh%start_row+1
		    ep=mesh%end_row-1
	    endif
        do i=1,adjoint_n
            if((adjoint_p(i,1).ge.sp).and.(adjoint_p(i,1).le.ep))then
                ustar1(i)=mesh%u(adjoint_p(i,1),adjoint_p(i,2))
                vstar1(i)=mesh%v(adjoint_p(i,1),adjoint_p(i,2))
                rhos1(i)=mesh%rho(adjoint_p(i,1),adjoint_p(i,2))
            endif
        enddo

        call MPI_ALLREDUCE(ustar1,ustar,num,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(vstar1,vstar,num,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(rhos1,rhos,num,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

    else
        !$acc kernels loop present(mesh,adjoint_p,ustar,vstar)
        do i=1,adjoint_n       
            ustar(i)=mesh%u(adjoint_p(i,1),adjoint_p(i,2))
            vstar(i)=mesh%v(adjoint_p(i,1),adjoint_p(i,2))
            !if(i.eq.1000)print*,ustar(i),adjoint_p(i,1),adjoint_p(i,2),mesh%u(adjoint_p(i,1),adjoint_p(i,2))
        enddo
        !$acc end kernels loop
    endif

    !$acc kernels present(mrpos,nrpos,mrow,mcol,js,js1,xx1,xx2,xx3)
    mrpos(1)=1
    do k=2,b_point-1
        mrpos(k)=mrpos(k-1)+mrow(k-1)
    enddo
    nrpos(1)=1
    
    do k=2,num
         nrpos(k)=nrpos(k-1)+mcol(k-1)
    enddo
    !$acc end kernels 
    !$acc update host(js1,js)
    if(numprocs.ne.1)then
        call fasttransposematrix(xx3,xx2,num,b_point-1,js1,nrpos)
    else
        call fasttransposematrix(xx3,xx2,num,b_point-1,js,nrpos)
    endif

    if(numprocs.ne.1)then 
        call multsmatrix(xx1,xx2,mrpos,nrpos,b_point-1,num,js1,as,ia,ja,ds)
    else
        call multsmatrix(xx1,xx2,mrpos,nrpos,b_point-1,num,js,as,ia,ja,ds)
    endif

    if(numprocs.ne.1)then 
        call matmulvocxs(bu,xx1,mrpos,ustar,b_point-1,num,js1)
    else
        call matmulvocxs(bu,xx1,mrpos,ustar,b_point-1,num,js)
    endif
    if(numprocs.ne.1)then 
        call matmulvocxs(bv,xx1,mrpos,vstar,b_point-1,num,js1)
    else
        call matmulvocxs(bv,xx1,mrpos,vstar,b_point-1,num,js)
    endif  

    !$acc kernels present(bu,ub,bv,vb)
    bu=ub-bu   
    bv=vb-bv
    !$acc end kernels 
    call cpu_time(time_s6)
!    call Jacobi_Iteration(comj,as,ia,ja,ds,bu,bv,b_point-1)
    call GMRES(comj(:,1),as,ia,ja,ds,bu,b_point-1)
    call GMRES(comj(:,2),as,ia,ja,ds,bv,b_point-1)
    call cpu_time(time_e6)
    call cpu_time(time_s7)
    if(numprocs.ne.1)then 
        call matmulvocxs(glx,xx2,nrpos,comj(:,1),num,b_point-1,js1)
        call matmulvocxs(gly,xx2,nrpos,comj(:,2),num,b_point-1,js1)
    else
        call matmulvocxs(glx,xx2,nrpos,comj(:,1),num,b_point-1,js)
        call matmulvocxs(gly,xx2,nrpos,comj(:,2),num,b_point-1,js)
    endif

!    if(numprocs.ne.1)then 
!        call matvol(as,comj(:,1),aa1(:,1),b_point-1,ia,ja,ds) 
!        call matvol(as,comj(:,2),aa1(:,2),b_point-1,ia,ja,ds) 
!    else
!        call matvol(as,comj(:,1),aa1(:,1),b_point-1,ia,ja,ds) 
!        call matvol(as,comj(:,2),aa1(:,2),b_point-1,ia,ja,ds) 
!    endif
!    print*, maxval(abs(bu-aa1(:,1))),maxval(abs(bv-aa1(:,2))),'1' 
    
    !$acc kernels present(ul,glx,ustar,vl,vstar,gly)
    ul=ustar+glx
    vl=vstar+gly
    !$acc end kernels 

!    if(numprocs.ne.1)then 
!        call matmulvocxs(bu,xx1,mrpos,ul,b_point-1,num,js1)
!    else
!        call matmulvocxs(bu,xx1,mrpos,ul,b_point-1,num,js)
!    endif   
!    if(numprocs.ne.1)then 
!        call matmulvocxs(bv,xx1,mrpos,vl,b_point-1,num,js1)
!    else
!        call matmulvocxs(bv,xx1,mrpos,vl,b_point-1,num,js)
!    endif   
!    print*, maxval(abs(bu-ub)),maxval(abs(bv-vb)),'2'

    spp=myid*int((adjoint_n)/numprocs)
	epp=(myid+1)*int((adjoint_n)/numprocs)-1
    if(myid.eq.0) spp=1
    if(myid.eq.numprocs-1) epp=adjoint_n

    !$acc kernels loop present(mesh) independent
    do i=spp,epp
        sum1=0.0d0
        sum2=0.0d0
        if(i.lt.adjoint_n)then
            j=nrpos(i+1)
        else
            if(numprocs.ne.1)then
                j=js1+1
            else
                j=js+1
            endif
!            if(i.eq.1751)print*,j,nrpos(i),js1,myid
        endif

        if(numprocs.ne.1)then
            rhos1(i)=2.0d0*rhos(i)*glx(i)
            rhos2(i)=2.0d0*rhos(i)*gly(i)
        else
            gforce(adjoint_p(i,1),adjoint_p(i,2),1)=2.0d0*mesh%rho(adjoint_p(i,1),adjoint_p(i,2))*glx(i)
            gforce(adjoint_p(i,1),adjoint_p(i,2),2)=2.0d0*mesh%rho(adjoint_p(i,1),adjoint_p(i,2))*gly(i)
        endif
    enddo
    !$acc end kernels loop

    if(numprocs.ne.1)then
        do j=1,numprocs
            if(j.eq.1) then
                recvcounts(j)=int((adjoint_n)/numprocs)-1
                displs(j)=0
            elseif(j.eq.numprocs) then
                recvcounts(j)=adjoint_n-(numprocs-1)*int((adjoint_n)/numprocs)+1
                displs(j)=displs(j-1)+recvcounts(j-1)
            else
                recvcounts(j)=int((adjoint_n)/numprocs)
                displs(j)=displs(j-1)+recvcounts(j-1)
            endif
        enddo

        call MPI_ALLGATHERV(rhos1(spp:epp),recvcounts(myid+1),MPI_DOUBLE_PRECISION,ustar1,recvcounts,displs,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
        call MPI_ALLGATHERV(rhos2(spp:epp),recvcounts(myid+1),MPI_DOUBLE_PRECISION,vstar1,recvcounts,displs,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

        do i=1,adjoint_n 
            if((adjoint_p(i,1).ge.sp).and.(adjoint_p(i,1).le.ep)) then
                gforce(adjoint_p(i,1),adjoint_p(i,2),1)=ustar1(i)
                gforce(adjoint_p(i,1),adjoint_p(i,2),2)=vstar1(i)
            endif               
        enddo
    endif 
    
    if(myid.eq.0)then
		force=0.0d0
        mom=0.0d0
        if(numprocs.ne.1)then
            force(1)=-sum(ustar1)
            force(2)=-sum(vstar1)
            do k=1,adjoint_n
                x0=(mesh%ox+mesh%scale*dble(adjoint_p(k,1)-1)-(xc+wing_initial_x_origin+sa_arr(1)))/mesh%scale-body_rel_pos(1)
                y0=(mesh%oy+mesh%scale*dble(adjoint_p(k,2)-1)-(yc+wing_initial_y_origin+sb_arr(1)))/mesh%scale
                mom=mom-x0*vstar1(k)+y0*ustar1(k)
            enddo
        else
            !$acc kernels present(mesh,adjoint_p,body_rel_pos,sa_arr,sb_arr,gforce,force) 
            force(1)=-sum(gforce(:,:,1))
            force(2)=-sum(gforce(:,:,2))
            !$acc loop private(x0,y0) reduction(+:mom)
            do k=1,adjoint_n
                x0=(mesh%ox+mesh%scale*dble(adjoint_p(k,1)-1)-(xc+wing_initial_x_origin+sa_arr(1)))/mesh%scale-body_rel_pos(1)
                y0=(mesh%oy+mesh%scale*dble(adjoint_p(k,2)-1)-(yc+wing_initial_y_origin+sb_arr(1)))/mesh%scale
                mom=mom-x0*gforce(adjoint_p(k,1),adjoint_p(k,2),2)+y0*gforce(adjoint_p(k,1),adjoint_p(k,2),1)
            enddo
            !$acc end loop
            !$acc end kernels
        endif

        !$acc update host(force,mom)
        if(solve_dynamic_equations.eq..true.)then !动力学方程求解
            body_rel_a_prev(1)=body_rel_a(1)
            if(closed_boundary.eq..false.)then
                body_rel_a(1)=force(1)/(solid_den/(mesh%scale*mesh%scale))
            else
                body_rel_a(1)=force(1)/(solid_den*area_airfoil/(mesh%scale*mesh%scale))
            endif
            body_rel_v_prev(1)=body_rel_v(1)
!            body_rel_v(1)=body_rel_a(1)+body_rel_v(1)
            body_rel_s_prev(1)=body_rel_s(1)
!            body_rel_s(1)=body_rel_s(1)+body_rel_v_prev(1)
!            body_rel_pos(1)=body_rel_pos(1)+body_rel_v_prev(1)
            body_rel_theta_prev(1)=body_rel_theta(1)
            body_rel_temp_theta_prev(1)=body_rel_temp_theta(1)
!            body_rel_theta(1)=body_rel_theta_prev(1)+body_rel_temp_theta_prev(1)
    !        if(closed_boundary.eq..false.)then
!            body_rel_temp_theta(1)=body_rel_temp_theta_prev(1)+par1*UU*UU*mesh%scale*mesh%scale*(body_rel_theta(1)+alphad)+par2*(mesh%scale*mesh%scale*mesh%scale*mesh%scale)*mom
    !        else
    !            body_rel_temp_theta(1)=body_rel_temp_theta_prev(1)
    !        endif
    !        print*, force(1),(solid_den*area_airfoil/(mesh%scale*mesh%scale))
    !        body_rel_a(1)=0.0d0
            if(if_rk4.eq..false.)then
                body_rel_v(1)=body_rel_a(1)+body_rel_v(1)
                body_rel_s(1)=body_rel_s(1)+body_rel_v_prev(1)
                body_rel_pos(1)=body_rel_pos(1)+body_rel_v_prev(1)
                body_rel_theta(1)=body_rel_theta_prev(1)+body_rel_temp_theta_prev(1)
                body_rel_temp_theta(1)=body_rel_temp_theta_prev(1)+par1*UU*UU*mesh%scale*mesh%scale*(body_rel_theta(1)+alphad)+par2*(mesh%scale*mesh%scale*mesh%scale*mesh%scale)*mom
            else
                call rk4(body_rel_s_prev(1),body_rel_s(1),body_rel_pos(1),body_rel_v_prev(1),body_rel_v(1),body_rel_a(1),body_rel_theta_prev(1),body_rel_theta(1),body_rel_temp_theta_prev(1),body_rel_temp_theta(1),par1*UU*UU*mesh%scale*mesh%scale,par2*(mesh%scale*mesh%scale*mesh%scale*mesh%scale)) !用四步四级龙格-库塔方法求解动力学方程
            endif
            call update_body_extra_from_rel()
        endif
    endif 
    if((solve_dynamic_equations.eq..true.).and.(numprocs.ne.1))then
        call MPI_BCAST(body_rel_a,wing_count,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(body_rel_a_prev,wing_count,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(body_rel_v,wing_count,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(body_rel_v_prev,wing_count,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(body_rel_s,wing_count,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(body_rel_s_prev,wing_count,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(body_rel_pos,wing_count,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(body_rel_theta_prev,wing_count,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(body_rel_theta,wing_count,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(body_rel_temp_theta_prev,wing_count,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(body_rel_temp_theta,wing_count,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call update_body_extra_from_rel()
    endif 
   
    !$acc kernels loop present(mesh) private(i,j)
    do k=1,adjoint_n
        if(numprocs.ne.1)then
           if((adjoint_p(k,1).ge.sp).and.(adjoint_p(k,1).le.ep)) then
                i=adjoint_p(k,1)
                j=adjoint_p(k,2)
                mesh%u(i,j)=ul(k)
                mesh%v(i,j)=vl(k)
            endif                           
        else
            i=adjoint_p(k,1)
            j=adjoint_p(k,2)
    !        if(i.eq.127.and.j.eq.132) print*,k
            mesh%u(i,j)=ul(k)
            mesh%v(i,j)=vl(k)
        endif
    enddo
    !$acc end kernels loop
    !$acc exit data delete (js,js1,ds,ub,vb,ustar,vstar,rhos,ustar1,vstar1,rhos1,rhos2,tema,comi,tem1,tem2,aa,aa1,bu,bv,x,xx1,xx2,xx3,as,non_zero,indexnum,recvcounts,displs,mrpos,nrpos,mcol,mrrow,mrow,mrcol,ia,ja,ul,vl,glx,gly,num)
    call cpu_time(time_e7)
end


subroutine LU_decomposition(A,b1,b2,Nm)
    implicit none
    !传入虚参
	integer::Nm,M(Nm)
    real(kind=8)::A(Nm,Nm),b1(Nm),b2(Nm),s(Nm),temp,tempy(Nm)
    !其余参数
    integer::i,j,k,t
    real(kind=8)::sum
	do k=1,Nm
		temp=maxval(A(k,:))
		A(k,:)=A(k,:)/temp
		b1(k)=b1(k)/temp
		b2(k)=b2(k)/temp
	enddo
    do k=1,Nm
        do i=k,Nm
			sum=0.d0
			do t=1,k-1
				sum=sum+A(i,t)*A(t,k)
			enddo
			s(i)=A(i,k)-sum
        enddo
		m(k)=k
		do i=k,Nm
			if(s(i).gt.s(m(k))) then
				m(k)=i
			endif
		enddo

		if(m(k).ne.k) then
			do t=1,Nm
				temp=a(k,t)
				a(k,t)=a(m(k),t)
				a(m(k),t)=temp
			enddo
			temp=s(k)
			s(k)=s(m(k))
			s(m(k))=temp
		endif

		a(k,k)=s(k)
		do j=k+1,Nm
			sum=0.0d0
			do t=1,k-1
				sum=sum+a(k,t)*a(t,j)
			enddo
			a(k,j)=a(k,j)-sum
		enddo
		do i=k+1,Nm
			a(i,k)=s(i)/a(k,k)
		enddo
    enddo

	do k=1,Nm-1
		t=m(k)
		temp=b1(k)
		b1(k)=b1(t)
		b1(t)=temp
		temp=b2(k)
		b2(k)=b2(t)
		b2(t)=temp
	enddo

	tempy(1)=b1(1)
	do i=2,Nm
		sum=0.0d0
		do t=1,i-1
			sum=sum+a(i,t)*tempy(t)
		enddo
		tempy(i)=b1(i)-sum
	enddo
	b1(Nm)=tempy(Nm)/a(Nm,Nm)
	do i=Nm-1,1,-1
		sum=0.0d0
		do t=i+1,Nm
			sum=sum+a(i,t)*b1(t)
		enddo
		b1(i)=(tempy(i)-sum)/a(i,i)
	enddo

	tempy(1)=b2(1)
	do i=2,Nm
		sum=0.0d0
		do t=1,i-1
			sum=sum+a(i,t)*tempy(t)
		enddo
		tempy(i)=b2(i)-sum
	enddo
	b2(Nm)=tempy(Nm)/a(Nm,Nm)
	do i=Nm-1,1,-1
		sum=0.0d0
		do t=i+1,Nm
			sum=sum+a(i,t)*b2(t)
		enddo
		b2(i)=(tempy(i)-sum)/a(i,i)
	enddo


    return
end

subroutine gauss_decomposition(A,b1,Nm)
    implicit none
	integer i,j,k,ik,Nm
	real*8 A(Nm,Nm),b1(Nm)
	real*8 temp,mik
	do k=1,Nm-1
		ik=k
		temp=abs(A(ik,k))
		do i =k,Nm
			if(abs(a(i,k)).gt.abs(temp)) then
				temp=abs(a(i,k))
				ik=i
			endif
		enddo
		do j=k,Nm
			temp=a(ik,j)
			a(ik,j)=a(k,j)
			a(k,j)=temp
		enddo
		temp=b1(k)
		b1(k)=b1(ik)
		b1(ik)=temp
		do i=k+1,Nm
			mik=a(i,k)/a(k,k)
			do j=k+1,Nm
				a(i,j)=a(i,j)-mik*a(k,j)
			enddo
			b1(i)=b1(i)-mik*b1(k)
		enddo
	enddo
	b1(Nm)=b1(Nm)/a(Nm,Nm)
	do k=Nm-1,1,-1
		temp=0.0d0
		do j=k+1,Nm
			temp=temp+a(k,j)*b1(j)
		enddo
		b1(k)=(b1(k)-temp)/a(k,k)
	enddo
end

subroutine matmulpe(c,a,b,num1,num2,num3)
	implicit none
	integer num1,num2,num3
	real*8 a(num1,num2),b(num2,num3),c(num1,num3) 
	integer i,j,k
	real*8 sum
	
	do i=1,num1
		do j=1,num3
			sum=0.0d0
			do k=1,num2
				sum=sum+a(i,k)*b(k,j)				
			enddo
			c(i,j)=sum
		enddo
	enddo
end

subroutine matmulvoc1(c,a,b,num1,num2)
    !$acc routine seq
	implicit none
	integer num1,num2
	real*8 a(num1,num2),b(num2),c(num1)
	integer i,j,k
	real*8 sum

!    call dgemv('N',num1,num2,alpha,a,num1,b,1,bate,c,1)
	do i=1,num1
		sum=0.0d0
		do j=1,num2
			sum=sum+a(i,j)*b(j)
		enddo
		c(i)=sum
	enddo

end
    
subroutine fasttransposematrix1(xx3,xx1,mu,nu,js,cpot)
	use com_date
    implicit none
    
    integer mu,nu,js,p,col,qq
    real (kind=ps) xx3(mu*nu,3),xx1(mu*nu,3) 
    integer cpot1(mu),cpot(mu)
    
    !$acc enter data create(p,col,qq,cpot1) 
    !$acc kernels present(p,xx3,js,cpot,qq,col,cpot1,xx1)
    cpot1=cpot
    do p=1,js
        col=nint(xx3(p,2))
        qq=cpot1(col)
        xx1(qq,1)=xx3(p,2)
        xx1(qq,2)=xx3(p,1)
        xx1(qq,3)=xx3(p,3)
        cpot1(col)=cpot1(col)+1
    enddo
    !$acc end kernels
    !$acc exit data delete(p,col,qq,cpot1) 
end
    
subroutine multsmatrix1(xx1,xx2,mrpos,nrpos,mu,nu,js,as,ia,ja,ds)
    use mpi
	use com_date
    implicit none
    integer mu,nu,js,arow,brow,tp,p,t,qq,ccol,ds,ds1,i,j,sp,ep
    integer mrpos(mu),nrpos(nu),iia(mu),ia(mu+1),ja(mu*mu)
    real (kind=ps) xx1(mu*nu,3),xx2(mu*nu,3),as(mu*mu),ctemp(mu)
    integer,  allocatable :: buf1(:),buf2(:),recvcounts(:),displs(:)
    real (kind=ps),allocatable :: buf3(:) 
    
    sp=myid*int(mu/numprocs)
	ep=(myid+1)*int(mu/numprocs)-1
    if(myid.eq.0) sp=1
    if(myid.eq.numprocs-1) ep=mu
    
    allocate (buf1(sp:ep),buf2(mu*mu),recvcounts(numprocs),displs(numprocs)) 
    allocate (buf3(mu*mu))
    buf1=0

    !$acc enter data create(ds1,iia,ctemp)
    !$acc kernels present (ja,xx1,xx2,mrpos,nrpos,as,ds,js,ds1,iia,ctemp) create(buf1,buf2,buf3)
    ctemp=0.0d0
    iia=0
    ia=0
    ja=0
    ds1=1
    do arow=sp,ep
        if(arow.lt.mu)then
            tp=mrpos(arow+1)
        else
            tp=js+1
        endif
        !$acc loop reduction(+:ctemp) private(t,brow)
        do p=mrpos(arow),tp-1
            brow=nint(xx1(p,2))
            if(brow.lt.nu)then
                t=nrpos(brow+1)
            else
                t=js+1
            endif
            do qq=nrpos(brow),t-1
                ccol=nint(xx2(qq,2))
                ctemp(ccol)=ctemp(ccol)+xx1(p,3)*xx2(qq,3)
            enddo
        enddo
        !$acc end loop
        !$acc loop 
        do ccol=1,mu
            if(ctemp(ccol).ne.0.0d0)then
                if(numprocs.ne.1)then
                    buf1(arow)= buf1(arow)+1 
                    buf2(ds1)=ccol
                    buf3(ds1)=ctemp(ccol)
                    ctemp(ccol)=0.0d0
                    ds1=ds1+1
                else
                    iia(arow)= iia(arow)+1 
                    ja(ds1)=ccol
                    as(ds1)=ctemp(ccol)
                    ctemp(ccol)=0.0d0
                    ds1=ds1+1
                endif
            endif
        enddo
        !$acc end loop
    enddo
    ds1=ds1-1
    !$acc end kernels

    if(numprocs.ne.1)then
        call MPI_ALLREDUCE(ds1,ds,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
        displs(1)=0
        do j=1,numprocs-1
            if(j.eq.1)then
                recvcounts(j)=int(mu/numprocs)-1
            else
                recvcounts(j)=int(mu/numprocs)
            endif
            displs(j+1)=displs(j)+recvcounts(j)
        enddo
        recvcounts(numprocs)=mu-(numprocs-1)*int(mu/numprocs)+1

        call MPI_ALLGATHERV(buf1(sp:ep),recvcounts(myid+1),MPI_INTEGER,iia,recvcounts,displs,MPI_INTEGER,MPI_COMM_WORLD,ierr)
        call MPI_ALLGATHER(ds1,1,MPI_INTEGER,recvcounts,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
        displs(1)=0
        do j=1,numprocs-1
            displs(j+1)=displs(j)+recvcounts(j)
        enddo
        call MPI_ALLGATHERV(buf2,recvcounts(myid+1),MPI_INTEGER,ja,recvcounts,displs,MPI_INTEGER,MPI_COMM_WORLD,ierr)
        call MPI_ALLGATHERV(buf3,recvcounts(myid+1),MPI_DOUBLE_PRECISION,as,recvcounts,displs,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
    else
        !$acc kernels present(ds,ds1)
        ds=ds1
        !$acc end kernels
    endif
    !$acc kernels present(ia,iia,ds)
    ia(1)=1
    do i=1,mu-1
        ia(i+1)=ia(i)+iia(i)
    enddo
    ia(mu+1)=ds+1
    !$acc end kernels
    !$acc exit data delete (ds1,iia,ctemp)
end
    
subroutine matmulvocxs1(bb,xx1,mrpos,xx,mu,nu,js)
    use mpi
	use com_date
    implicit none
    integer mu,nu,js,arow,brow,tp,p,t,qq,ccol,j,sp,ep
    integer mrpos(mu)
    real (kind=ps) xx1(mu*nu,3),xx(nu),bb(mu)
    integer,  allocatable :: recvcounts(:),displs(:)
    real (kind=ps),allocatable :: buf(:) 
    
    sp=myid*int(mu/numprocs)
	ep=(myid+1)*int(mu/numprocs)-1
    if(myid.eq.0) sp=1
    if(myid.eq.numprocs-1) ep=mu
    allocate (recvcounts(numprocs),displs(numprocs)) 
    allocate (buf(sp:ep))
    !$acc enter data create(buf)
    !$acc kernels present(buf,bb)
    bb=0.0d0
    buf=0.0d0
    do arow=sp,ep
        if(arow.lt.mu)then
            tp=mrpos(arow+1)
        else
            tp=js+1
        endif
        do p=mrpos(arow),tp-1
            brow=nint(xx1(p,2))
            if(numprocs.ne.1)then
                buf(arow)=buf(arow)+xx1(p,3)*xx(brow) 
            else
                bb(arow)=bb(arow)+xx1(p,3)*xx(brow)
            endif
        enddo
    enddo
    !$acc end kernels
    !$acc exit data delete (buf)
    if(numprocs.ne.1)then
        displs(1)=0
        do j=1,numprocs-1
            if(j.eq.1)then
                recvcounts(j)=int(mu/numprocs)-1
            else
                recvcounts(j)=int(mu/numprocs)
            endif
            displs(j+1)=displs(j)+recvcounts(j)
        enddo
        recvcounts(numprocs)=mu-(numprocs-1)*int(mu/numprocs)+1
        call MPI_ALLGATHERV(buf(sp:ep),recvcounts(myid+1),MPI_DOUBLE_PRECISION,bb,recvcounts,displs,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
    endif
!    print*, buf(sp+5),bb(sp+5),myid
end
    
subroutine yanzheng(as,bv,ia,ja,num1,buf1,xx1,mrpos,vstar,num,js,ds)
    implicit none
    integer num1,num,js,ds,k,j,tp,p
    real*8 as(num1*num1),bv(num1),buf1(num1),buf2(num1),vstar(num),xx1(num1*num,3)
    integer mrpos(num1),marpos,ia(num1+1),ja(num1*num1)
    real*8, allocatable :: ass(:,:)
    
    allocate (ass(num1*num1,3))
    
    call matmulvocxs(buf2,xx1,mrpos,vstar,num1,num,js)
    do k=1,num1
        if(k.lt.num1)then
            tp=ia(k+1)
        else
            tp=ds+1
        endif
        do p=ia(k),tp-1
            ass(p,1)=dble(k)
            ass(p,2)=dble(ja(p))
            ass(p,3)=as(p)
        enddo
    enddo
    call matmulvocxs(buf1,ass,ia,bv,num1,num1,ds)
    buf1=buf1+buf2
end

subroutine GMRES1(x,a,ia,ja,ds,b,n)
    implicit none
    
    integer i,j,n,k,ds,m,l,s
    integer ia(n+1),ja(n*n)
    real*8 a(n*n,3),b(n),r(n),e1(n),beta(n),x(n)
    real*8 b_norm,r_norm,err,temp
    real*8, allocatable :: qq(:,:),h(:,:),cs(:),sn(:),y(:)
    !$acc routine (norm) seq
!    m=n-1
    m=20
    l=5
    allocate (qq(n,m),h(m+1,m+1),cs(m),sn(m),y(m)) 
    !$acc enter data create(r,qq,h,cs,sn,y,b_norm,r_norm,err,e1,beta,y) copyin(n)

!    x=1.0d0
 
    do s=1,l
        !call matmulvocxs(r,a,ia,x,n,n,ds)
        call matvol(a,x,r,n,ia,ja,ds) 
        !$acc parallel present(b,r) copyout(temp)
        r=b-r
        temp=maxval(abs(r))
        !$acc end parallel
        if(temp.le.1d-20) exit
        !$acc kernels present(b,r,b_norm,r_norm,err,qq,e1,beta,n,sn,cs)
        call norm(b,b_norm,n)
        call norm(r,r_norm,n)
        err=r_norm/b_norm
    
        e1=0.0d0
        e1(1)=1.0d0
        qq(:,1)=r/r_norm
        beta=r_norm*e1
        h=0.0d0
        sn=0.0d0
        cs=0.0d0
        !$acc end kernels

        do k=1,m
            call arnoldi(a,qq,k,n,h(:,k),ia,ja,ds,m)
            call apply_givens_rotation(h(:,k), cs, sn, k,n,m)
            !update the residual vector
            !$acc kernels present (sn,cs,beta)
            beta(k+1) = -sn(k)*beta(k)
            beta(k)   = cs(k)*beta(k)
            err  = abs(beta(k+1)) / b_norm
            !$acc end kernels
 !           if(err.lt.1d-3) then
 !               print*, 'here!'
 !               exit
 !           endif
        enddo
!        open(2001,file='hhg.dat',status='unknown')    
!        do i=1,m+1
!            write(2001,'(21f14.8)') h(i,1:m+1)
!        enddo
!        pause

        !$acc kernels present(qq,h,x,y,beta)
        y(k-1)=beta(k-1)/h(k-1,k-1)
        !$acc loop
        do i=k-2,1,-1
            temp=0.0d0
            !$acc loop
            do j=i+1,k-1
                temp=temp+h(i,j)*y(j)
            enddo
            !$acc end loop
            y(i)=(beta(i)- temp)/h(i,i)
        enddo
        !$acc end loop
        !y = h(1:k,1:k) \ beta(1:k)    
        !x = x + Q(:,1:k)*y
        do i=1,n
            temp=0.0d0
            do j=1,k-1
                temp=temp+qq(i,j)*y(j)
            enddo
            x(i)=x(i)+temp
        enddo
        !$acc end kernels
    enddo
!    call matmulvocxs(r,a,ia,x,n,n,ds)
!    print*,maxval(abs(b-r))
    !$acc exit data delete(r,qq,h,cs,sn,y,b_norm,r_norm,err,e1,beta,n,y) 
end
    
subroutine norm2(a,a_norm,n)
    !$acc routine seq
    implicit none
    integer n,i,j
    real*8 a(n),tem,a_norm
    
    tem=0.0d0
    do i=1,n
        tem=tem+a(i)*a(i)
    enddo   
    a_norm=sqrt(tem)
end  

subroutine norm1(a,a_norm,n)
    implicit none
    integer n,i,j
    real*8 a(n),tem,a_norm
    
    !$acc kernels present(a,a_norm,n) create(tem)
    tem=0.0d0
    !$acc loop reduction(+:tem)
    do i=1,n
        tem=tem+a(i)*a(i)
    enddo  
    !$acc end loop  
    a_norm=sqrt(tem)
    !$acc end kernels
end  
    
subroutine arnoldi1(a, qq, k,n,h,ia,ja,ds,m)
    implicit none 
    integer n,k,i,j,ds,m
    integer ia(n+1),ja(n*n)
    real*8 a(n*n,3),qq(n,m),q1(n),h(m+1)
  
    !$acc routine (norm) seq
    !$acc enter data create(q1)
    call matvol(a,qq(:,k),q1,n,ia,ja,ds)
!    call matmulvocxs(q1,a,ia,qq(:,k),n,n,ds)
    !$acc kernels present (qq,n,h,q1)
    do i = 1,k
        h(i)=0.0d0
        do j=1,n
            h(i)= h(i)+q1(j)*qq(j,i)
        enddo
        q1 = q1 - h(i)*qq(:,i)
    enddo
    call norm(q1,h(k+1),n)
    q1 = q1 / h(k+1)
    if (k.ne.m) qq(:,k+1)=q1
    !$acc end kernels
    !$acc exit data delete(q1)
end
    
subroutine apply_givens_rotation1(h, cs, sn, k,n,m)
    implicit none
    integer i,k,n,m
    real*8 cs(m),sn(m),temp,cs_k,sn_k,h(m+1)

    !$acc routine (givens_rotation) seq
    !$acc kernels present(h,cs,sn) create(temp,cs_k,sn_k)
    !apply for ith column
    !$acc loop
    do i = 1,k-1                              
        temp   =  cs(i)*h(i) + sn(i)*h(i+1)
        h(i+1) = -sn(i)*h(i) + cs(i)*h(i+1)
        h(i)   = temp
    enddo
    !$acc end loop
  
    !update the next sin cos values for rotation
    call givens_rotation(h(k), h(k+1),cs_k,sn_k)
    cs(k)=cs_k
    sn(k)=sn_k
  
    !eliminate H(i+1,i)
    h(k) = cs_k*h(k) + sn_k*h(k+1)
    h(k+1) = 0.0d0
    !$acc end kernels
end
    
subroutine givens_rotation1 (v1, v2,cs_k,sn_k)
    !$acc routine seq
    implicit none
    real*8 v1,v2,cs_k,sn_k,t
    
    if (v1==0)then
        cs_k = 0.0d0
        sn_k = 1.0d0
    else
        t=sqrt(v1*v1+v2*v2)
        cs_k = abs(v1) / t
        sn_k = cs_k * v2 / v1
    endif    
end

subroutine matvol1(xx1,xx,bb,mu,ia,ja,js)
    use mpi
	use com_date
    implicit none
    integer mu,nu,js,arow,brow,tp,p,t,qq,ccol,j,sp,ep
    integer mrpos(mu)
    integer ia(mu+1),ja(mu*mu)
    real (kind=ps) xx1(mu*mu),xx(mu),bb(mu)
    integer,  allocatable :: recvcounts(:),displs(:)
    real (kind=ps),allocatable :: buf(:) 
    
    sp=myid*int(mu/numprocs)
	ep=(myid+1)*int(mu/numprocs)-1
    if(myid.eq.0) sp=1
    if(myid.eq.numprocs-1) ep=mu
    allocate (recvcounts(numprocs),displs(numprocs)) 
    allocate (buf(sp:ep))

    !$acc enter data create(buf,mrpos)
    !$acc kernels present(xx1,xx,bb,ia,ja,js,buf,mrpos)
    mrpos=ia(1:mu-1)
    bb=0.0d0
    buf=0.0d0
    !$acc loop private(tp,brow) 
    do arow=sp,ep
        if(arow.lt.mu)then
            tp=mrpos(arow+1)
        else
            tp=js+1
        endif
        do p=mrpos(arow),tp-1
            brow=ja(p)
            if(numprocs.ne.1)then
                buf(arow)=buf(arow)+xx1(p)*xx(brow) 
            else
                bb(arow)=bb(arow)+xx1(p)*xx(brow)
            endif
        enddo
    enddo
    !$acc end loop
    !$acc end kernels
    !$acc exit data delete(buf,mrpos)
    if(numprocs.ne.1)then
        displs(1)=0
        do j=1,numprocs-1
            if(j.eq.1)then
                recvcounts(j)=int(mu/numprocs)-1
            else
                recvcounts(j)=int(mu/numprocs)
            endif
            displs(j+1)=displs(j)+recvcounts(j)
        enddo
        recvcounts(numprocs)=mu-(numprocs-1)*int(mu/numprocs)+1
        call MPI_ALLGATHERV(buf(sp:ep),recvcounts(myid+1),MPI_DOUBLE_PRECISION,bb,recvcounts,displs,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
    endif
!    print*, buf(sp+5),bb(sp+5),myid
end

subroutine rk4(s1,s,s2,v1,v,a1,jd1,jd,jsd1,jsd,p1,p2)
    use mpi
	use com_date
    implicit none 

    real (kind=ps) s1,s,s2,v1,v,a1
    real (kind=ps) jd1,jd,jsd1,jsd,p1,p2
    real (kind=ps) k11,k12,k13,k14,k21,k22,k23,k24
    
    k11=v1
    k21=a1
    k12=v1+0.5d0*k21
    k22=a1
    k13=v1+0.5d0*k22
    k23=a1
    k14=v1+k23
    k24=a1
    
    s=s1+(k11+2.0d0*k12+2.0d0*k13+k14)/6.0d0
    v=v1+(k21+2.0d0*k22+2.0d0*k23+k24)/6.0d0
    s2=s2+(k11+2.0d0*k12+2.0d0*k13+k14)/6.0d0
    
    k11=jsd1
    k21=p1*(jd1+alphad)+p2*mom
    k12=jsd1+0.5d0*k21
    k22=p1*(jd1+alphad+0.5d0*k11)+p2*mom
    k13=jsd1+0.5d0*k22
    k23=p1*(jd1+alphad+0.5d0*k12)+p2*mom
    k14=jsd1+k23
    k24=p1*(jd1+alphad+k13)+p2*mom
    
    jd=jd1+(k11+2.0d0*k12+2.0d0*k13+k14)/6.0d0
    jsd=jsd1+(k21+2.0d0*k22+2.0d0*k23+k24)/6.0d0
    
end
