Subroutine init(mesh) 
	!$acc routine (D2Q9feq) seq
	use com_date
    use meshtype
	implicit none
	

	integer i,j,k,n,sp,ep
	type(mesh_dat) mesh(max_meshid)
	Real (kind=ps), external :: D2Q9feq

	time=0.0d0
	call mode
	call D2Q9model
	rho0=1.0d0

	matm(1,:)=(/1.0d0,1.0d0,1.0d0,1.0d0,1.0d0,1.0d0,1.0d0,1.0d0,1.0d0/)
    matm(2,:)=(/-4.0d0,-1.0d0,-1.0d0,-1.0d0,-1.0d0,2.0d0,2.0d0,2.0d0,2.0d0/)
    matm(3,:)=(/4.0d0,-2.0d0,-2.0d0,-2.0d0,-2.0d0,1.0d0,1.0d0,1.0d0,1.0d0/)
    matm(4,:)=(/0.0d0,1.0d0,0.0d0,-1.0d0,0.0d0,1.0d0,-1.0d0,-1.0d0,1.0d0/)
    matm(5,:)=(/0.0d0,-2.0d0,0.0d0,2.0d0,0.0d0,1.0d0,-1.0d0,-1.0d0,1.0d0/)
    matm(6,:)=(/0.0d0,0.0d0,1.0d0,0.0d0,-1.0d0,1.0d0,1.0d0,-1.0d0,-1.0d0/)
    matm(7,:)=(/0.0d0,0.0d0,-2.0d0,0.0d0,2.0d0,1.0d0,1.0d0,-1.0d0,-1.0d0/)
    matm(8,:)=(/0.0d0,1.0d0,-1.0d0,1.0d0,-1.0d0,0.0d0,0.0d0,0.0d0,0.0d0/)
    matm(9,:)=(/0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,1.0d0,-1.0d0,1.0d0,-1.0d0/)
    
    call nizhen(matm,inversematm,9)
!    print*,inversematm
!    print*, matmul(matm,inversematm)
!    rightmove=(/1,2,4,5,6,8,9,10,12,13,11,7,3/)
!    leftmove=(/3,2,4,7,6,8,11,10,12,13,9,5,1/)
    
    do n=1,max_meshid
		mesh(n)%u=0.0d0
		mesh(n)%v=0.0d0
		mesh(n)%u0=0.0d0
		mesh(n)%v0=0.0d0
		mesh(n)%rho=rho0
		mesh(n)%pre=1.0d0/3.0d0
        if(ref_freestream.eq.1) then
		    mesh(n)%u=uu
		    mesh(n)%u0=uu
        endif
		if(myid.eq.0) then
			sp=1
			ep=mesh(n)%Nx
		else
			sp=mesh(n)%start_row
			ep=mesh(n)%end_row
		endif

		Do i=sp,ep
			Do j=1,mesh(n)%Ny
				Do k=1,Q
					mesh(n)%f(i,j,k)=D2Q9feq(mesh(n)%rho(i,j),mesh(n)%u(i,j),mesh(n)%v(i,j),i,j,k)
					mesh(n)%f0(i,j,k)=mesh(n)%f(i,j,k)
				enddo
			enddo
        enddo
        if(myid.eq.0)mesh(n)%leftf2=mesh(n)%f(1,2:mesh(n)%Ny-1,:)
        if(myid.eq.0)mesh(n)%rightf2=mesh(n)%f(mesh(n)%Nx,2:mesh(n)%Ny-1,:)
        if(myid.eq.0)mesh(n)%upf2=mesh(n)%f(:,mesh(n)%Ny,:)
        if(myid.eq.0)mesh(n)%downf2=mesh(n)%f(:,1,:)

        mesh(n)%body=0
        mesh(n)%f_old=mesh(n)%f
        !$acc update device(mesh(n)%rho,mesh(n)%u,mesh(n)%v,mesh(n)%u0,mesh(n)%v0,mesh(n)%pre,mesh(n)%f,mesh(n)%f0,mesh(n)%f_old,mesh(n)%body)
	enddo
	    if (.not. allocated(gforce)) then
	        allocate(gforce(mesh(max_meshid)%Nx,mesh(max_meshid)%Ny,2))
	    endif
		gforce=0.0d0
        adjoint_n=0
	    mesh(max_meshid)%body=0
		!$acc update device (mesh)
		!$acc update device (gforce)
	    ! Boundary coordinates are loaded after init; update on the first solver tick.
end Subroutine

Subroutine erranalysis(mesh,bnh)
	use mpi
	use com_date
    use meshtype
	implicit none
	

	integer i,j,n,bnh
	integer sp,ep
	type(mesh_dat) mesh(max_meshid)
	real (kind=ps) mytemp1,mytemp2,temp1,temp2,tempt1,tempt2
	

	tempt1=0.0d0
	tempt2=0.0d0
	do n=1,max_meshid	
		if(myid.eq.0) then
			sp=mesh(n)%start_row+1
			ep=mesh(n)%end_row-1
		elseif(myid.eq.numprocs-1) then
			sp=mesh(n)%start_row+1
			ep=mesh(n)%end_row-1
		else
			sp=mesh(n)%start_row+1
			ep=mesh(n)%end_row-1
		endif
		mytemp1=0.0d0
		mytemp2=0.0d0
		!$acc kernels loop reduction(+:mytemp1,mytemp2)	
		do i= sp,ep
			do j=2,mesh(n)%Ny-1
                mytemp1=mytemp1+(mesh(n)%u(i,j)-mesh(n)%u0(i,j))*(mesh(n)%u(i,j)-mesh(n)%u0(i,j))+(mesh(n)%v(i,j)-mesh(n)%v0(i,j))*(mesh(n)%v(i,j)-mesh(n)%v0(i,j))
                mytemp2=mytemp2+mesh(n)%u(i,j)*mesh(n)%u(i,j)+mesh(n)%v(i,j)*mesh(n)%v(i,j)	
			enddo
        enddo
		!$acc end kernels loop	

		call MPI_REDUCE(mytemp1,temp1,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
		call MPI_REDUCE(mytemp2,temp2,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)

		if(myid.eq.0) then
			tempt1=tempt1+temp1
			tempt2=tempt2+temp2
		endif
	enddo	
	if (myid.eq.0) then
		tempt1=sqrt(tempt1)
		tempt2=sqrt(tempt2)
		err=tempt1/(tempt2+1d-30)
	endif

	call MPI_BCAST (err,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!	print*, err
end Subroutine

subroutine output(cycle1,mesh)
    use nested_io
    implicit none
    integer :: cycle1
    type(mesh_dat) :: mesh(max_meshid)
    call output_nested_flow(cycle1,mesh)
end subroutine

Subroutine    output_gforce(mm,mesh)
	use com_date
    use meshtype
	implicit none
    
    type(mesh_dat) mesh
    integer mm
    integer i,j,outfile
	integer i1,i2,i3,i4,nn,totalnodecounts,totalfacecounts
	integer is,ie,js,je,iis,iie,jjs,jje
    character(len=512) :: out_path

	parameter(outfile=1267)

	character(len=32) :: nm
    
    write(nm,'(i0)') mm
	    call build_flow_output_path("g"//trim(nm)//".dat", out_path)
    open(outfile,file=trim(out_path),status="unknown")
    write(outfile,'(a)') 'VARIABLES = "X" "Y" "gforcex" "gforcey" '
    write(outfile,*) 'zone,t="',mm,'",i=',mesh%Nx,',j=',mesh%Ny
    do i=1,mesh%Ny
        do j=1,mesh%Nx
            write(outfile,2000) mesh%ox+dble(i-1)*mesh%scale,mesh%oy+dble(j-1)*mesh%scale,gforce(i,j,1),gforce(i,j,2)
        enddo
    enddo
2000	format(25e25.17)    
 end

Subroutine evolution(mesh,cycle1,nforce,h)
	use mpi
	use cublas
    use cusolverDn
	use com_date
    use meshtype
    use nested_grid
!    USE DFPORT, ONLY: DTIME
	implicit none
	
	type(mesh_dat) mesh(max_meshid)
	Real (kind=ps), external :: D2Q9feq,distanceq

	integer i,j,k,k1,ip,jp,ip1,jp1,ip2,jp2,sp,ep,setnum,myhavebodyp,havebodyp,n
	integer m1,m2,m3,cycle1,nforce,l1,l2
	integer cou,numt,row_len,wing,col
	real (kind=ps) tempf,qr,rhob,ub,vb,total_moment
	real (kind=ps) force_row(4+12*max_wings_supported),motion_row(1+3*max_wings_supported)
	type(cusolverDnHandle) :: h

    l1=2
    l2=mesh(max_meshid)%Ny-1
    call prepare_nested_tick(mesh,cycle1)
	do n=1,size(layer1)
        time=dble(cycle1-1)*(mesh(layer1(1))%scale*UU)
        if(time.ge.0.0d0) then
            call mode
            !angle_arr(1)=-body_rel_theta(1)
            !omega_arr(1)=-body_rel_temp_theta(1)/(mesh(layer1(1))%scale*UU)
            !omgac_arr(1)=-(par1*body_rel_theta(1)+par2*mom/(UU*UU*(1.0d0/(mesh(layer1(1))%scale*mesh(layer1(1))%scale))))
            call cpu_time(time_s8)       
            call move_boundary(mesh(max_meshid))
			call cpu_time(time_e8)
        endif
		if(myid.eq.0) then
			sp=mesh(layer1(n))%start_row+1
			ep=mesh(layer1(n))%end_row-1
		elseif(myid.eq.numprocs-1) then
			sp=mesh(layer1(n))%start_row+1
			ep=mesh(layer1(n))%end_row-1
		else
			sp=mesh(layer1(n))%start_row+1
			ep=mesh(layer1(n))%end_row-1
		endif
		if (numprocs.eq.1) then
			sp=mesh(layer1(n))%start_row+1
			ep=mesh(layer1(n))%end_row-1
        endif

		call cpu_time(time_s1)
        call streaming(mesh(max_meshid),sp,ep,l1,l2)
		call cpu_time(time_e1)
		call cpu_time(time_s2)
		call fill_nested_boundary(mesh(max_meshid-1),mesh(max_meshid), &
            dble(mod(cycle1-1,2))/2d0+0.5d0)
		call physical_variables(mesh(max_meshid),l1,l2)
		call cpu_time(time_e2)
        !call immersed_boundary_suzuki_and_inamuro(mesh(max_meshid))
		call cpu_time(time_s3)
        call collect_ibm_velocity(mesh(max_meshid))
		call immersed_boundary_wu(mesh(max_meshid),h)
        call distribute_ibm_velocity(mesh(max_meshid))
		call cpu_time(time_e3)
		call cpu_time(time_s4)
		call colission(mesh(max_meshid),sp,ep,l1,l2)
		call cpu_time(time_e4)
	        do wing=1,wing_count
	            wing_force(wing,:)=2.0d0*wing_force(wing,:)/(UU*UU*(1.0d0/mesh(max_meshid)%scale))
	            wing_moment(wing)=-2.0d0*wing_moment(wing)/(UU*UU*(1.0d0/(mesh(max_meshid)%scale*mesh(max_meshid)%scale)))
	        enddo
	        force=0.0d0
	        do wing=1,wing_count
	            force(1)=force(1)+wing_force(wing,1)
	            force(2)=force(2)+wing_force(wing,2)
	        enddo
	        total_moment=sum(wing_moment(1:wing_count))
	        force1=0.0d0
	        force2=0.0d0
	        mom=0.0d0
	        mom2=0.0d0
	        if (wing_count.ge.1) then
	            force1=wing_force(1,:)
	            mom=wing_moment(1)
	        endif
	        if (wing_count.ge.2) then
	            force2=wing_force(2,:)
	            mom2=wing_moment(2)
	        endif
	        force_row=0.0d0
	        force_row(1)=time
	        col=2
	        do wing=1,wing_count
	            force_row(col)=wing_force(wing,1)
	            force_row(col+1)=wing_force(wing,2)
	            force_row(col+2)=wing_moment(wing)
	            col=col+3
	        enddo
	        force_row(col)=force(1)
	        force_row(col+1)=force(2)
	        force_row(col+2)=total_moment
	        col=col+3
		        do wing=1,wing_count
		            force_row(col)=wing_initial_x_origin + dble(wing-1)*wing_initial_x_spacing + sa_arr(wing) + &
		                body_extra_s(wing)/(1.0d0/mesh(layer1(n))%scale)
		            force_row(col+1)=sau_arr(wing)
		            force_row(col+2)=saac_arr(wing)
		            force_row(col+3)=wing_initial_y_origin + sb_arr(wing)
	            force_row(col+4)=sbv_arr(wing)
	            force_row(col+5)=sbac_arr(wing)
	            force_row(col+6)=angle_arr(wing)
	            force_row(col+7)=omega_arr(wing)
	            force_row(col+8)=omgac_arr(wing)
	            col=col+9
	        enddo
	        row_len=col-1
	        if(myid.eq.0) then
	            latest_force_row=0.0d0
	            latest_force_row(1:row_len)=force_row(1:row_len)
	            latest_force_sum_row(1:row_len)=latest_force_sum_row(1:row_len)+force_row(1:row_len)
	            latest_force_row_len=row_len
	            latest_force_sample_count=latest_force_sample_count+1
	        endif
	        if(myid.eq.0.and.(cycle1.gt.cycles)) then
	            write(nforce,'(1000f15.8)') force_row(1:row_len)
	        endif
	        nforce=102
		        if((myid.eq.0).and.(cycle1.gt.cycles))then
		            motion_row=0.0d0
		            motion_row(1)=time
		            col=2
		            do wing=1,wing_count
		                motion_row(col)=0.0d0
		                motion_row(col+1)=body_extra_v(wing)/UU
		                motion_row(col+2)=wing_initial_x_origin + dble(wing-1)*wing_initial_x_spacing + sa_arr(wing) + &
		                    body_extra_s(wing)/(1.0d0/mesh(layer1(n))%scale)
		                col=col+3
		            enddo
		            write(nforce,'(1000f15.8)') motion_row(1:col-1)
		        endif
	        nforce=101
	        if(numprocs.ne.1)then
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
	    enddo  
	call cpu_time(time_s9) 
	call finish_nested_tick(mesh,cycle1)
	call cpu_time(time_e9)
end Subroutine

Subroutine reversedirection(k,k1)
	implicit none
	integer k,k1

	if (k.eq.1) k1=1
	if (k.eq.2) k1=4
	if (k.eq.3) k1=5
	if (k.eq.4) k1=2
	if (k.eq.5) k1=3
	if (k.eq.6) k1=8
	if (k.eq.7) k1=9
	if (k.eq.8) k1=6
	if (k.eq.9) k1=7

end Subroutine

Subroutine streaming(mesh,sp,ep,l1,l2)
	use mpi
	use com_date
    use meshtype
	implicit none
	
	type(mesh_dat) mesh
	integer sp,ep
	integer i,j,k,ip,jp,l1,l2

	!$acc kernels present(mesh,mesh%f0,mesh%f)
	!$acc loop independent
	do k=1,Q	
		!$acc loop independent
		do j=2,mesh%Ny-1
			!$acc loop independent private(ip,jp)
			do i=sp,ep	
				ip=i-e(k,1)
				jp=j-e(k,2)
				mesh%f0(i,j,k)=mesh%f(ip,jp,k)  
            enddo
			!$acc end loop
		enddo
		!$acc end loop
    enddo
	!$acc end loop
	!$acc end kernels loop

	if(numprocs.ne.1) then
		if(myid.eq.0) then
			i=1
			j=mesh%Nx
		else
			i=mesh%start_row
			j=mesh%end_row
        endif
        
		call MPI_BARRIER(MPI_COMM_WORLD,ierr)
		!$acc update self(mesh%f0)
		call transferdata(mesh%f0,i,j,mesh%start_row,mesh%end_row,mesh%Nx,mesh%Ny,mesh%meshid,l1,l2)
        !$acc update device(mesh%f0)
		call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    endif
end Subroutine

Subroutine physical_variables(mesh,l1,l2)
	use mpi
	use com_date
    use meshtype
	implicit none
	
	type(mesh_dat) mesh
	integer sp,ep,l1,l2
	integer i,j,k,ip,jp
	real (kind=ps) temp,temprho,tempu,tempv
 	
	if(myid.eq.0) then
		sp=1
        if(numprocs.eq.1) then
            ep=mesh%Nx
        else
            ep=mesh%end_row-1
        endif
	elseif ((myid.eq.numprocs-1).and.(numprocs.ne.1)) then
		sp=mesh%start_row+1
		ep=mesh%Nx
	else
		sp=mesh%start_row+1
		ep=mesh%end_row-1
    endif
!    if((myid.eq.4).and.(mesh%meshid.eq.13)) print*, myid,mesh%f(403,129,:)

	!$acc kernels present(mesh)
	!$acc loop independent
	do j=1,mesh%Ny
!			if((mesh%body(i,j).eq.3).or.(mesh%body(i,j).eq.4)) then

	!				tempf=0.0
		!$acc loop independent
		do i=sp,ep
    		mesh%u0(i,j)=mesh%u(i,j)
    		mesh%v0(i,j)=mesh%v(i,j)
    		!mesh%rho(i,j)=0.0d0
    		!mesh%u(i,j)=0.0d0
    		!mesh%v(i,j)=0.0d0
			temprho=0.0d0
			tempu=0.0d0
			tempv=0.0d0			
    		!$acc loop independent
			do k=1,Q
				if((i.eq.1).or.(i.eq.mesh%Nx).or.(j.eq.1).or.(j.eq.mesh%Ny)) then
					temp=mesh%f0(i,j,k)
					temprho=temprho+temp
					tempu=tempu+e(k,1)*temp
					tempv=tempv+e(k,2)*temp
				else
					mesh%f(i,j,k)=mesh%f0(i,j,k)
					temp=mesh%f0(i,j,k)
					temprho=temprho+temp
					tempu=tempu+e(k,1)*temp
					tempv=tempv+e(k,2)*temp
		!			if(k.ne.1) tempf=tempf+f(i,j,k)
                endif
			enddo
			!$acc end loop		
	!				pre(i,j)=0.6*(tempf+sa(i,j,1))

!				endif
    		mesh%rho(i,j)=temprho
			mesh%pre(i,j)=temprho/3.0d0
    		mesh%u(i,j)=tempu/temprho
    		mesh%v(i,j)=tempv/temprho
        enddo
		!$acc end loop	
    enddo
	!$acc end loop	
	!$acc end kernels

    if (numprocs.ne.1) then
			if(myid.eq.0) then
				i=1
				j=mesh%Nx
			else
				i=mesh%start_row
				j=mesh%end_row
			endif
!			call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!			call garthervf(mesh%f,i,j,mesh%start_row,mesh%end_row,mesh%Nx,mesh%Ny)
!			call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!			if (myid.eq.0) then
!				call gather (mesh%u,1,mesh%Nx,mesh%Nx,mesh%Ny,mesh%start_row,mesh%end_row)
!				call gather (mesh%v,1,mesh%Nx,mesh%Nx,mesh%Ny,mesh%start_row,mesh%end_row)
!				call gather (mesh%rho,1,mesh%Nx,mesh%Nx,mesh%Ny,mesh%start_row,mesh%end_row)
!			else
!				call gather (mesh%u,mesh%start_row,mesh%end_row,mesh%Nx,mesh%Ny,mesh%start_row,mesh%end_row)
!				call gather (mesh%v,mesh%start_row,mesh%end_row,mesh%Nx,mesh%Ny,mesh%start_row,mesh%end_row)
!				call gather (mesh%rho,mesh%start_row,mesh%end_row,mesh%Nx,mesh%Ny,mesh%start_row,mesh%end_row)
!			endif
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        !$acc update self(mesh%f,mesh%rho,mesh%u,mesh%v,mesh%pre)
        call transferdatapv(mesh%f,mesh%rho,mesh%u,mesh%v,mesh%pre,i,j,mesh%start_row,mesh%end_row,mesh%Nx,mesh%Ny,mesh%meshid,l1,l2)
        !$acc update device(mesh%f,mesh%rho,mesh%u,mesh%v,mesh%pre)
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    endif
!     if(myid.eq.0.and.mesh%meshid.eq.13) print*,mesh%rho(59,99)
!    if((myid.eq.4).and.(mesh%meshid.eq.13)) print*, myid,mesh%f(403,129,:)
!    if(myid.eq.4.and.mesh%meshid.eq.11) print*, myid,mesh%f0(114,30,5)
end Subroutine

Subroutine colission(mesh,sp,ep,l1,l2)
	use mpi
	use com_date
    use meshtype
	implicit none
	
	type(mesh_dat) mesh
	integer sp,ep
	integer i,j,k,ip,jp,l1,l2
	real(kind=ps) gfx,gfy,temp(Q),temp1(Q),temmat(Q,Q),fp(Q),ut,vt,rt
	Real (kind=ps), external :: D2Q9feq,distanceq
	!$acc routine (distanceq) seq
	!$acc routine (matmulvoc) seq

	!$acc kernels present(mesh,gforce) create(gfx,gfy,ut,vt,rt)
	!$acc loop independent private (temp,temp1,fp,ut,vt,rt,gfx,gfy)
	do j=2,mesh%Ny-1
		!$acc loop independent
		do i=sp,ep
            !temp=mesh%f(i,j,:)
			!$acc loop seq
			do k=1,Q
				temp(k)=mesh%f(i,j,k)
			enddo
			!$acc end loop
            call matmulvoc(temp1,matm,temp,Q,Q)
			!mesh%m(i,j,:)=temp1

            ut=mesh%u(i,j)
            vt=mesh%v(i,j)
            rt=mesh%rho(i,j)
            temp(1)=0.0d0
            temp(2)=-2.0d0*rt+3.0d0*(rt*ut*rt*ut+rt*vt*rt*vt)
            temp(3)=rt-3.0d0*(rt*ut*rt*ut+rt*vt*rt*vt)
            temp(4)=0.0d0
            temp(5)=-rt*ut
            temp(6)=0.0d0
            temp(7)=-rt*vt
            temp(8)=(rt*ut*rt*ut-rt*vt*rt*vt)
            temp(9)=(rt*ut*rt*vt)

            fp=0.0d0
            if (mesh%layer.eq.max_meshid) then
				gfx=gforce(i,j,1)
				gfy=gforce(i,j,2)
                fp(1)=0.0d0
                fp(2)=6.0d0*(1.0d0-0.5d0*mesh%s(2))*(ut*gfx+vt*gfy)
                fp(3)=-6.0d0*(1.0d0-0.5d0*mesh%s(3))*(ut*gfx+vt*gfy)
                fp(4)=gfx
                fp(5)=-(1.0d0-0.5d0*mesh%s(5))*gfx
                fp(6)=gfy
                fp(7)=-(1.0d0-0.5d0*mesh%s(7))*gfy
                fp(8)=2.0d0*(1.0d0-0.5*mesh%s(8))*(ut*gfx-vt*gfy)
                fp(9)=(1.0d0-0.5*mesh%s(9))*(ut*gfy+vt*gfx)
            endif
         
			!$acc loop seq            
            !do k=1,Q			
            !    mesh%m(i,j,k)=mesh%m(i,j,k)+(temp(k)-mesh%m(i,j,k))*mesh%s(k)+fp(k)
            !enddo
			do k=1,Q			
                temp(k)=temp1(k)+(temp(k)-temp1(k))*mesh%s(k)+fp(k)
            enddo
			!$acc end loop
            !temp=mesh%m(i,j,:)
            call matmulvoc(temp1,inversematm,temp,Q,Q)
            mesh%f(i,j,:)=temp1
		enddo
		!$acc end loop
    enddo
	!$acc end loop
	!$acc end kernels

	if(numprocs.ne.1) then
		if(myid.eq.0) then
			i=1
			j=mesh%Nx
		else
			i=mesh%start_row
			j=mesh%end_row
		endif
		call MPI_BARRIER(MPI_COMM_WORLD,ierr)
		!$acc update self(mesh%f)
        call transferdata(mesh%f,i,j,mesh%start_row,mesh%end_row,mesh%Nx,mesh%Ny,mesh%meshid,l1,l2)
        !$acc update device(mesh%f)
		call MPI_BARRIER(MPI_COMM_WORLD,ierr)
	endif
end Subroutine


Subroutine transfer_coarse_to_fine1(mesh,n,cycle1)
	use mpi
	use com_date
    use meshtype
	implicit none
	
	type(mesh_dat) mesh(max_meshid)
	integer n,m,j,k,i1,j1,i2,sp,ep,spp,i,cycle1
	Real (kind=ps), external :: D2Q9feq,distanceq
	Real (kind=ps) temp
	Real (kind=ps), allocatable :: temp1(:,:),temp2(:,:)


	if(myid.eq.0) then
		m=mesh(n)%left(2,1) !\D7\F3\B2\E0\CD\F8\B8\F1\B1߽\E7\B5Ĵ\AB\B5\DD
		if(mesh(m)%layer.lt.mesh(n)%layer) then !\C8\E7\B9\FBmesh(n)\D7\F3\B1ߵ\C4\CD\F8\B8\F1\B5Ĳ\E3\CA\FDС\D3\DAmesh(n)\B5\C4\CD\F8\B8\F1\B5Ĳ\E3\CA\FD
			do j=2,mesh(n)%Ny-1
!				mesh(n)%leftf2(j,:)=mesh(n)%f(1,j,:)
				mesh(n)%f(1,j,:)=mesh(n)%f0(1,j,:)
			enddo
						
			i2=	mesh(m)%Ny-2
			if(allocated(temp1)) deallocate(temp1)
			if(allocated(temp2)) deallocate(temp2)
			allocate(temp1(i2,Q),temp2(mesh(n)%Ny-2,Q))
			do i=1,i2
				i1=mesh(m)%Nx-1
				j1=i+1
                
                call mrttransctf(temp1(i,:),mesh(m)%rho(i1,j1),mesh(m)%u(i1,j1),mesh(m)%v(i1,j1),mesh(m)%f(i1,j1,:),mesh(m)%s,mesh(n)%s,Q)
			enddo

			if(mesh(n)%left(2,2).ne.0) then
				sp=mesh(n)%left(2,4)-1
				ep=1
			else
				sp=mesh(n)%left(3,4)-2
				ep=0
			endif

			call sanwanju1(temp1,temp2,i2,mesh(n)%Ny-2,sp,ep,Q)
			do i=2,mesh(n)%Ny-1
				mesh(n)%f0(1,i,:)=temp2(i-1,:)
			enddo
		endif

		m=mesh(n)%right(2,1) !\D3Ҳ\E0\CD\F8\B8\F1\B1߽\E7\B5Ĵ\AB\B5\DD
		if(mesh(m)%layer.lt.mesh(n)%layer) then !\C8\E7\B9\FBmesh(n)\D3ұߵ\C4\CD\F8\B8\F1\B5Ĳ\E3\CA\FDС\D3\DAmesh(n)\B5\C4\CD\F8\B8\F1\B5Ĳ\E3\CA\FD
			do j=2,mesh(n)%Ny-1
!				mesh(n)%rightf2(j,:)=mesh(n)%f(mesh(n)%Nx,j,:)
				mesh(n)%f(mesh(n)%Nx,j,:)=mesh(n)%f0(mesh(n)%Nx,j,:)
			enddo
			
			i2=	mesh(m)%Ny-2

			if(allocated(temp1)) deallocate(temp1)
			if(allocated(temp2)) deallocate(temp2)
			allocate(temp1(i2,Q),temp2(mesh(n)%Ny-2,Q))

			do i=1,i2
				i1=2
				j1=i+1
                call mrttransctf(temp1(i,:),mesh(m)%rho(i1,j1),mesh(m)%u(i1,j1),mesh(m)%v(i1,j1),mesh(m)%f(i1,j1,:),mesh(m)%s,mesh(n)%s,Q)
			enddo

			if(mesh(n)%right(2,2).ne.0) then
				sp=mesh(n)%right(2,4)-1
				ep=1
			else
				sp=mesh(n)%right(3,4)-2
				ep=0
			endif

			call sanwanju1(temp1,temp2,i2,mesh(n)%Ny-2,sp,ep,Q)
			do i=2,mesh(n)%Ny-1
				mesh(n)%f0(mesh(n)%Nx,i,:)=temp2(i-1,:)
			enddo
		endif

		m=mesh(n)%up(1,1) !\C9ϲ\E0\CD\F8\B8\F1\B1߽\E7\B5Ĵ\AB\B5\DD
		if(mesh(m)%layer.lt.mesh(n)%layer) then !\C8\E7\B9\FBmesh(n)\C9ϱߵ\C4\CD\F8\B8\F1\B5Ĳ\E3\CA\FDС\D3\DAmesh(n)\B5\C4\CD\F8\B8\F1\B5Ĳ\E3\CA\FD
			do j=1,mesh(n)%Nx
!				mesh(n)%upf2(j,:)=mesh(n)%f(j,mesh(n)%Ny,:)
				mesh(n)%f(j,mesh(n)%Ny,:)=mesh(n)%f0(j,mesh(n)%Ny,:)
			enddo

			spp=4*(mesh(m)%layer-1)
			i2=	mesh(spp+1)%Nx+mesh(spp+2)%Nx+mesh(spp+3)%Nx-6			

			if(allocated(temp1)) deallocate(temp1)
			if(allocated(temp2)) deallocate(temp2)
			allocate(temp1(i2,Q),temp2(mesh(n)%Nx,Q))

			do i=1,i2
				if(i.le.mesh(spp+1)%Nx-2) then
					i1=i+1
					j1=mesh(spp+2)%left(2,4)
                    call mrttransctf(temp1(i,:),mesh(spp+1)%rho(i1,j1),mesh(spp+1)%u(i1,j1),mesh(spp+1)%v(i1,j1),mesh(spp+1)%f(i1,j1,:),mesh(spp+1)%s,mesh(n)%s,Q)
				elseif((i.gt.mesh(spp+1)%Nx-2).and.(i.le.(mesh(spp+1)%Nx+mesh(spp+2)%Nx-4)))then
					i1=i-(mesh(spp+1)%Nx-2)+1
					j1=2
                    call mrttransctf(temp1(i,:),mesh(spp+2)%rho(i1,j1),mesh(spp+2)%u(i1,j1),mesh(spp+2)%v(i1,j1),mesh(spp+2)%f(i1,j1,:),mesh(spp+2)%s,mesh(n)%s,Q)
				else
					i1=i-(mesh(spp+1)%Nx+mesh(spp+2)%Nx-4)+1
					j1=mesh(spp+2)%right(2,4)
                    call mrttransctf(temp1(i,:),mesh(spp+3)%rho(i1,j1),mesh(spp+3)%u(i1,j1),mesh(spp+3)%v(i1,j1),mesh(spp+3)%f(i1,j1,:),mesh(spp+3)%s,mesh(n)%s,Q)
				endif
			enddo

			if(mesh(n)%up(1,2).ne.0) then
				if(n.eq.spp+5) then
					sp=mesh(spp+1)%Nx-2
				elseif(n.eq.spp+6) then
					sp=mesh(n)%up(1,3)+mesh(spp+1)%Nx-3
				elseif(n.eq.spp+7) then
					sp=mesh(n)%up(1,3)+mesh(spp+1)%Nx-3
				endif
				ep=1
			else
				if(n.eq.spp+5) then
					sp=mesh(spp+1)%Nx-2
				elseif(n.eq.spp+6) then
					sp=mesh(n)%up(2,3)+mesh(spp+1)%Nx-4
				elseif(n.eq.spp+7) then
					sp=mesh(n)%up(2,3)+mesh(spp+1)%Nx-4
				endif
				ep=0
			endif

			call sanwanju1(temp1,temp2,i2,mesh(n)%Nx,sp,ep,Q)
			do i=1,mesh(n)%Nx
				mesh(n)%f0(i,mesh(n)%Ny,:)=temp2(i,:)
			enddo					
		endif

		m=mesh(n)%down(1,1) !\CF²\E0\CD\F8\B8\F1\B1߽\E7\B5Ĵ\AB\B5\DD
		if(mesh(m)%layer.lt.mesh(n)%layer) then !\C8\E7\B9\FBmesh(n)\CF±ߵ\C4\CD\F8\B8\F1\B5Ĳ\E3\CA\FDС\D3\DAmesh(n)\B5\C4\CD\F8\B8\F1\B5Ĳ\E3\CA\FD
			do j=1,mesh(n)%Nx
!				mesh(n)%downf2(j,:)=mesh(n)%f(j,1,:)
				mesh(n)%f(j,1,:)=mesh(n)%f0(j,1,:)
			enddo
			spp=4*(mesh(m)%layer-1)
			i2=	mesh(spp+1)%Nx+mesh(spp+3)%Nx+mesh(spp+4)%Nx-6
			if(allocated(temp1)) deallocate(temp1)
			if(allocated(temp2)) deallocate(temp2)
			allocate(temp1(i2,Q),temp2(mesh(n)%Nx,Q))

			do i=1,i2
				if(i.le.mesh(spp+1)%Nx-2) then
					i1=i+1
					j1=mesh(spp+4)%left(mesh(spp+4)%Ny-1,4)
                    call mrttransctf(temp1(i,:),mesh(spp+1)%rho(i1,j1),mesh(spp+1)%u(i1,j1),mesh(spp+1)%v(i1,j1),mesh(spp+1)%f(i1,j1,:),mesh(spp+1)%s,mesh(n)%s,Q)
				elseif((i.gt.mesh(spp+1)%Nx-2).and.(i.le.(mesh(spp+1)%Nx+mesh(spp+4)%Nx-4)))then
					i1=i-(mesh(spp+1)%Nx-2)+1
					j1=mesh(spp+4)%Ny-1
                    call mrttransctf(temp1(i,:),mesh(spp+4)%rho(i1,j1),mesh(spp+4)%u(i1,j1),mesh(spp+4)%v(i1,j1),mesh(spp+4)%f(i1,j1,:),mesh(spp+4)%s,mesh(n)%s,Q)
				else
					i1=i-(mesh(spp+1)%Nx+mesh(spp+4)%Nx-4)+1
					j1=mesh(spp+4)%right(mesh(spp+4)%Ny-1,4)
                    call mrttransctf(temp1(i,:),mesh(spp+3)%rho(i1,j1),mesh(spp+3)%u(i1,j1),mesh(spp+3)%v(i1,j1),mesh(spp+3)%f(i1,j1,:),mesh(spp+1)%s,mesh(n)%s,Q)
				endif
			enddo

			if(mesh(n)%down(1,2).ne.0) then
				if(n.eq.spp+5) then
					sp=mesh(spp+1)%Nx-2
				elseif(n.eq.spp+7) then
					sp=mesh(n)%down(1,3)+mesh(spp+1)%Nx-3
				elseif(n.eq.spp+8) then
					sp=mesh(n)%down(1,3)+mesh(spp+1)%Nx-3
				endif
				ep=1
			else
				if(n.eq.spp+5) then
					sp=mesh(spp+1)%Nx-2
				elseif(n.eq.spp+7) then
					sp=mesh(n)%down(2,3)+mesh(spp+1)%Nx-4
				elseif(n.eq.spp+8) then
					sp=mesh(n)%down(2,3)+mesh(spp+1)%Nx-4
				endif
				ep=0
			endif

			call sanwanju1(temp1,temp2,i2,mesh(n)%Nx,sp,ep,Q)
			do i=1,mesh(n)%Nx
				mesh(n)%f0(i,1,:)=temp2(i,:)
			enddo

!			if (cycle1.eq.100) then
!				call darw(temp1,0,i2,mesh(spp+1)%ox,mesh(spp+1)%scale)
!				call darw(temp2,1,mesh(n)%Nx,mesh(n)%ox,mesh(n)%scale)
!				pause
!			endif					
		endif

	endif

	if(numprocs.ne.1) then
		if(myid.eq.0) then
			sp=1
			ep=mesh(n)%Nx
		else
			sp=mesh(n)%start_row
			ep=mesh(n)%end_row
		endif
		call scatter(mesh(n)%f0,sp,ep,mesh(n)%Nx,mesh(n)%Ny)
        call scatter(mesh(n)%f,sp,ep,mesh(n)%Nx,mesh(n)%Ny)
    endif

!    if((cycle1.eq.1000).and.(myid.eq.0)) call darw

end Subroutine

Subroutine transfer_coarse_to_fine(mesh,n,cycle1)
	use mpi
	use com_date
    use meshtype
	implicit none
	
	type(mesh_dat) mesh(max_meshid)
	integer n,m,j,k,i1,j1,i2,sp,ep,i,cycle1
	Real (kind=ps), external :: D2Q9feq,distanceq
	Real (kind=ps) temp
	Real (kind=ps), allocatable :: temp1(:,:),temp2(:,:)


	if(myid.eq.0) then
		m=mesh(n)%left(2,1) !\D7\F3\B2\E0\CD\F8\B8\F1\B1߽\E7\B5Ĵ\AB\B5\DD
		if(mesh(m)%layer.lt.mesh(n)%layer) then !\C8\E7\B9\FBmesh(n)\D7\F3\B1ߵ\C4\CD\F8\B8\F1\B5Ĳ\E3\CA\FDС\D3\DAmesh(n)\B5\C4\CD\F8\B8\F1\B5Ĳ\E3\CA\FD
			do j=2,mesh(n)%Ny-1
				mesh(n)%leftf2(j,:)=mesh(n)%f(1,j,:)
				mesh(n)%f(1,j,:)=mesh(n)%f0(1,j,:)
			enddo
			do j=2,mesh(n)%Ny-1
				i1=mesh(n)%left(j,3)
				j1=mesh(n)%left(j,4)
				if(mesh(n)%left(j,2).eq.1) then
					do k=1,Q 
						temp=D2Q9feq(mesh(m)%rho(i1,j1),mesh(m)%u(i1,j1),mesh(m)%v(i1,j1),i1,j1,k)
						mesh(n)%f0(1,j,k)=temp+(mesh(n)%tau_f-1.0d0)/(2.0d0*(mesh(m)%tau_f-1.0d0))*(mesh(m)%f(i1,j1,k)-temp)
					enddo
				endif
			enddo
			
			if(mesh(n)%left(2,2).eq.0)then
				sp=1
			else
				sp=2
			endif
			if(mesh(n)%left(mesh(n)%Ny-1,2).eq.0)then
				ep=mesh(n)%Ny
			else
				ep=mesh(n)%Ny-1
			endif
			
			i2=	(ep-sp)/2+1
			if(allocated(temp1)) deallocate(temp1)
			if(allocated(temp2)) deallocate(temp2)
			allocate(temp1(i2,Q),temp2(2*i2-1,Q))
			do i=1,i2
				if ((i.eq.1).and.(sp.eq.1))then
					i1=mesh(n)%left(3,3)
					j1=mesh(n)%left(3,4)-1
					do k=1,Q
						temp=D2Q9feq(mesh(m)%rho(i1,j1),mesh(m)%u(i1,j1),mesh(m)%v(i1,j1),i1,j1,k)
						temp1(i,k)=temp+(mesh(n)%tau_f-1.0d0)/(2.0d0*(mesh(m)%tau_f-1.0d0))*(mesh(m)%f(i1,j1,k)-temp)
					enddo
				else if((i.eq.i2).and.(ep.eq.mesh(n)%Ny))then
					i1=mesh(n)%left(mesh(n)%Ny-2,3)
					j1=mesh(n)%left(mesh(n)%Ny-2,4)+1
					do k=1,Q
						temp=D2Q9feq(mesh(m)%rho(i1,j1),mesh(m)%u(i1,j1),mesh(m)%v(i1,j1),i1,j1,k)
						temp1(i,k)=temp+(mesh(n)%tau_f-1.0d0)/(2.0d0*(mesh(m)%tau_f-1.0d0))*(mesh(m)%f(i1,j1,k)-temp)
					enddo
				else
					temp1(i,:)=mesh(n)%f0(1,sp+(i-1)*2,:)
				endif
			enddo

			call sanwanju(temp1,temp2,i2,2*i2-1,Q)
			do i=2,mesh(n)%Ny-1
				mesh(n)%f0(1,i,:)=temp2(i+(1-sp),:)
			enddo
		endif

		m=mesh(n)%right(2,1) !\D3Ҳ\E0\CD\F8\B8\F1\B1߽\E7\B5Ĵ\AB\B5\DD
		if(mesh(m)%layer.lt.mesh(n)%layer) then !\C8\E7\B9\FBmesh(n)\D3ұߵ\C4\CD\F8\B8\F1\B5Ĳ\E3\CA\FDС\D3\DAmesh(n)\B5\C4\CD\F8\B8\F1\B5Ĳ\E3\CA\FD
			do j=2,mesh(n)%Ny-1
				mesh(n)%rightf2(j,:)=mesh(n)%f(mesh(n)%Nx,j,:)
				mesh(n)%f(mesh(n)%Nx,j,:)=mesh(n)%f0(mesh(n)%Nx,j,:)
			enddo
			do j=2,mesh(n)%Ny-1
				i1=mesh(n)%right(j,3)
				j1=mesh(n)%right(j,4)
				if(mesh(n)%right(j,2).eq.1) then
					do k=1,Q 
						temp=D2Q9feq(mesh(m)%rho(i1,j1),mesh(m)%u(i1,j1),mesh(m)%v(i1,j1),i1,j1,k)
						mesh(n)%f0(mesh(n)%Nx,j,k)=temp+(mesh(n)%tau_f-1.0d0)/(2.0d0*(mesh(m)%tau_f-1.0d0))*(mesh(m)%f(i1,j1,k)-temp)
					enddo
				endif
			enddo
			if(mesh(n)%right(2,2).eq.0)then
				sp=1
			else
				sp=2
			endif
			if(mesh(n)%right(mesh(n)%Ny-1,2).eq.0)then
				ep=mesh(n)%Ny
			else
				ep=mesh(n)%Ny-1
			endif
			
			i2=	(ep-sp)/2+1

			if(allocated(temp1)) deallocate(temp1)
			if(allocated(temp2)) deallocate(temp2)
			allocate(temp1(i2,Q),temp2(2*i2-1,Q))

			do i=1,i2
				if ((i.eq.1).and.(sp.eq.1))then
					i1=mesh(n)%right(3,3)
					j1=mesh(n)%right(3,4)-1
					do k=1,Q
						temp=D2Q9feq(mesh(m)%rho(i1,j1),mesh(m)%u(i1,j1),mesh(m)%v(i1,j1),i1,j1,k)
						temp1(i,k)=temp+(mesh(n)%tau_f-1.0d0)/(2.0d0*(mesh(m)%tau_f-1.0d0))*(mesh(m)%f(i1,j1,k)-temp)
					enddo
				else if((i.eq.i2).and.(ep.eq.mesh(n)%Ny))then
					i1=mesh(n)%right(mesh(n)%Ny-2,3)
					j1=mesh(n)%right(mesh(n)%Ny-2,4)+1
					do k=1,Q
						temp=D2Q9feq(mesh(m)%rho(i1,j1),mesh(m)%u(i1,j1),mesh(m)%v(i1,j1),i1,j1,k)
						temp1(i,k)=temp+(mesh(n)%tau_f-1.0d0)/(2.0d0*(mesh(m)%tau_f-1.0d0))*(mesh(m)%f(i1,j1,k)-temp)
					enddo
				else
					temp1(i,:)=mesh(n)%f0(mesh(n)%Nx,sp+(i-1)*2,:)
				endif
			enddo

			call sanwanju(temp1,temp2,i2,2*i2-1,Q)
			do i=2,mesh(n)%Ny-1
				mesh(n)%f0(mesh(n)%Nx,i,:)=temp2(i+(1-sp),:)
			enddo
		endif

		m=mesh(n)%up(1,1) !\C9ϲ\E0\CD\F8\B8\F1\B1߽\E7\B5Ĵ\AB\B5\DD
		if(mesh(m)%layer.lt.mesh(n)%layer) then !\C8\E7\B9\FBmesh(n)\C9ϱߵ\C4\CD\F8\B8\F1\B5Ĳ\E3\CA\FDС\D3\DAmesh(n)\B5\C4\CD\F8\B8\F1\B5Ĳ\E3\CA\FD
			do j=1,mesh(n)%Nx
				mesh(n)%upf2(j,:)=mesh(n)%f(j,mesh(n)%Ny,:)
				mesh(n)%f(j,mesh(n)%Ny,:)=mesh(n)%f0(j,mesh(n)%Ny,:)
			enddo
			do j=1,mesh(n)%Nx
				i1=mesh(n)%up(j,3)
				j1=mesh(n)%up(j,4)
				if(mesh(n)%up(j,2).eq.1) then
					do k=1,Q 
						temp=D2Q9feq(mesh(m)%rho(i1,j1),mesh(m)%u(i1,j1),mesh(m)%v(i1,j1),i1,j1,k)
						mesh(n)%f0(j,mesh(n)%Ny,k)=temp+(mesh(n)%tau_f-1.0d0)/(2.0d0*(mesh(m)%tau_f-1.0d0))*(mesh(m)%f(i1,j1,k)-temp)
					enddo
				endif
			enddo
			if(mesh(n)%up(1,2).eq.0)then
				sp=0
			else
				sp=1
			endif
			if(mesh(n)%up(mesh(n)%Nx,2).eq.0)then
				ep=mesh(n)%Nx+1
			else
				ep=mesh(n)%Nx
			endif
			
			i2=	(ep-sp)/2+1
			if(allocated(temp1)) deallocate(temp1)
			if(allocated(temp2)) deallocate(temp2)
			allocate(temp1(i2,Q),temp2(2*i2-1,Q))

			do i=1,i2
				if ((i.eq.1).and.(sp.eq.0))then
					i1=mesh(n)%up(2,3)-1
					j1=mesh(n)%up(2,4)
					do k=1,Q
						temp=D2Q9feq(mesh(m)%rho(i1,j1),mesh(m)%u(i1,j1),mesh(m)%v(i1,j1),i1,j1,k)
						temp1(i,k)=temp+(mesh(n)%tau_f-1.0d0)/(2.0d0*(mesh(m)%tau_f-1.0d0))*(mesh(m)%f(i1,j1,k)-temp)
					enddo
				else if((i.eq.i2).and.(ep.eq.mesh(n)%Nx+1))then
					i1=mesh(n)%up(mesh(n)%Nx-1,3)+1
					j1=mesh(n)%up(mesh(n)%Nx-1,4)
					do k=1,Q
						temp=D2Q9feq(mesh(m)%rho(i1,j1),mesh(m)%u(i1,j1),mesh(m)%v(i1,j1),i1,j1,k)
						temp1(i,k)=temp+(mesh(n)%tau_f-1.0d0)/(2.0d0*(mesh(m)%tau_f-1.0d0))*(mesh(m)%f(i1,j1,k)-temp)
					enddo
				else
					temp1(i,:)=mesh(n)%f0(sp+(i-1)*2,mesh(n)%Ny,:)
				endif
			enddo

			call sanwanju(temp1,temp2,i2,2*i2-1,Q)
			do i=1,mesh(n)%Nx
				mesh(n)%f0(i,mesh(n)%Ny,:)=temp2(i+(1-sp),:)
			enddo					
		endif

		m=mesh(n)%down(1,1) !\CF²\E0\CD\F8\B8\F1\B1߽\E7\B5Ĵ\AB\B5\DD
		if(mesh(m)%layer.lt.mesh(n)%layer) then !\C8\E7\B9\FBmesh(n)\CF±ߵ\C4\CD\F8\B8\F1\B5Ĳ\E3\CA\FDС\D3\DAmesh(n)\B5\C4\CD\F8\B8\F1\B5Ĳ\E3\CA\FD
			do j=1,mesh(n)%Nx
				mesh(n)%downf2(j,:)=mesh(n)%f(j,1,:)
				mesh(n)%f(j,1,:)=mesh(n)%f0(j,1,:)
			enddo
			do j=1,mesh(n)%Nx
				i1=mesh(n)%down(j,3)
				j1=mesh(n)%down(j,4)
				if(mesh(n)%down(j,2).eq.1) then
					do k=1,Q 
						temp=D2Q9feq(mesh(m)%rho(i1,j1),mesh(m)%u(i1,j1),mesh(m)%v(i1,j1),i1,j1,k)
						mesh(n)%f0(j,1,k)=temp+(mesh(n)%tau_f-1.0d0)/(2.0d0*(mesh(m)%tau_f-1.0d0))*(mesh(m)%f(i1,j1,k)-temp)
					enddo
				endif
			enddo
			if(mesh(n)%down(1,2).eq.0)then
				sp=0
			else
				sp=1
			endif
			if(mesh(n)%down(mesh(n)%Nx,2).eq.0)then
				ep=mesh(n)%Nx+1
			else
				ep=mesh(n)%Nx
			endif
			
			i2=	(ep-sp)/2+1
			if(allocated(temp1)) deallocate(temp1)
			if(allocated(temp2)) deallocate(temp2)
			allocate(temp1(i2,Q),temp2(2*i2-1,Q))

			do i=1,i2
				if ((i.eq.1).and.(sp.eq.0))then
					i1=mesh(n)%down(2,3)-1
					j1=mesh(n)%down(2,4)
					do k=1,Q
						temp=D2Q9feq(mesh(m)%rho(i1,j1),mesh(m)%u(i1,j1),mesh(m)%v(i1,j1),i1,j1,k)
						temp1(i,k)=temp+(mesh(n)%tau_f-1.0d0)/(2.0d0*(mesh(m)%tau_f-1.0d0))*(mesh(m)%f(i1,j1,k)-temp)
					enddo
				else if((i.eq.i2).and.(ep.eq.mesh(n)%Nx+1))then
					i1=mesh(n)%down(mesh(n)%Nx-1,3)+1
					j1=mesh(n)%down(mesh(n)%Nx-1,4)
					do k=1,Q
						temp=D2Q9feq(mesh(m)%rho(i1,j1),mesh(m)%u(i1,j1),mesh(m)%v(i1,j1),i1,j1,k)
						temp1(i,k)=temp+(mesh(n)%tau_f-1.0d0)/(2.0d0*(mesh(m)%tau_f-1.0d0))*(mesh(m)%f(i1,j1,k)-temp)
					enddo
				else
					temp1(i,:)=mesh(n)%f0(sp+(i-1)*2,1,:)
				endif
			enddo

			call sanwanju(temp1,temp2,i2,2*i2-1,Q)
			do i=1,mesh(n)%Nx
				mesh(n)%f0(i,1,:)=temp2(i+(1-sp),:)
			enddo					
		endif

	endif

	if(numprocs.ne.1) then
		if(myid.eq.0) then
			sp=1
			ep=mesh(n)%Nx
		else
			sp=mesh(n)%start_row
			ep=mesh(n)%end_row
		endif
		call scatter(mesh(n)%f0,sp,ep,mesh(n)%Nx,mesh(n)%Ny)
	endif



end Subroutine

subroutine sanwanju(f1,f3,n1,n3,qm)
	implicit none

	integer n1,n3,i,j,qm
	real*8 f1(n1,qm),f3(n3,qm)
	real*8 alpha(1:n1),bate(n1),gama(2:n1)
	real*8 p(n1),qqq(n1),m(n1)
	real*8 hi,x1,x2,x3
	integer nn

	hi=dble(n1-1)/dble(n3-1)	
	do j=1,9
		alpha=0.5d0
		gama=0.5d0

		do i=2,n1-1
			bate(i)=3.0d0*(f1(i+1,j)-2.0d0*f1(i,j)+f1(i-1,j))
		enddo
		alpha(1)=0.0d0
		bate(1)=0.0d0
		gama(n1)=0.0d0
		bate(n1)=0.0d0
		p(1)=2.0d0
		do i=1,n1-1
			qqq(i)=alpha(i)/p(i)
			p(i+1)=2.0d0-gama(i+1)*qqq(i)
		enddo
		bate(1)=bate(1)/p(1)
		do i=2,n1
			bate(i)=(bate(i)-gama(i)*bate(i-1))/p(i)
		enddo
		m(n1)=bate(n1)
		do i=n1-1,1,-1
			m(i)=bate(i)-qqq(i)*m(i+1)
		enddo		
		do i=1,n3-1
			x1=dble(i-1)*hi+1.0d0
			x2=dble(int(x1))
			x3=x2+1.0d0
			nn=nint(x3)
			f3(i,j)=m(nn-1)/6.0d0*((x3-x1)**3.0d0)+m(nn)/6.0d0*((x1-x2)**3.0d0)+(f1(nn-1,j)-m(nn-1)/6.0d0)*(x3-x1)+(f1(nn,j)-m(nn)/6.0d0)*(x1-x2)	
		enddo
		f3(n3,j)=f1(n1,j)
	enddo

end subroutine

subroutine sanwanju1(f1,f3,n1,n3,sp,ep,qm)
	implicit none

	integer n1,n3,i,j,qm,sp,ep
	real*8 f1(n1,qm),f3(n3,qm)
	real*8 alpha(1:n1),bate(n1),gama(2:n1)
	real*8 p(n1),qqq(n1),m(n1)
	real*8 hi,x1,x2,x3
	integer nn

	hi=0.5d0	
	do j=1,qm
		alpha=0.5d0
		gama=0.5d0

		do i=2,n1-1
			bate(i)=3.0d0*(f1(i+1,j)-2.0d0*f1(i,j)+f1(i-1,j))
		enddo
		alpha(1)=0.0d0
		bate(1)=0.0d0
		gama(n1)=0.0d0
		bate(n1)=0.0d0
		p(1)=2.0d0
		do i=1,n1-1
			qqq(i)=alpha(i)/p(i)
			p(i+1)=2.0d0-gama(i+1)*qqq(i)
		enddo
		bate(1)=bate(1)/p(1)
		do i=2,n1
			bate(i)=(bate(i)-gama(i)*bate(i-1))/p(i)
		enddo
		m(n1)=bate(n1)
		do i=n1-1,1,-1
			m(i)=bate(i)-qqq(i)*m(i+1)
		enddo		
		do i=1,n3
			if(ep.eq.0) then
				x1=dble(sp)+dble(i)*hi
				x2=dble(int(x1))
				x3=x2+1.0d0
				nn=nint(x3)
			else
				x1=dble(sp)+dble(i-1)*hi
				x2=dble(int(x1))
				x3=x2+1.0d0
				nn=nint(x3)
			endif
			f3(i,j)=m(nn-1)/6.0d0*((x3-x1)**3.0d0)+m(nn)/6.0d0*((x1-x2)**3.0d0)+(f1(nn-1,j)-m(nn-1)/6.0d0)*(x3-x1)+(f1(nn,j)-m(nn)/6.0d0)*(x1-x2)	
		enddo
	enddo

end subroutine

subroutine darw(ff,mm,nc,sp,dx)
	use com_date
	implicit none
	integer mm,nc,j,n
	real*8 ff(nc,9)
	real*8 sp,dx
    character(len=512) :: out_path


	if(mm.eq.0) then
        call build_output_path('frecordold.dat', out_path)
        open(1234,file=trim(out_path),status='unknown')
    endif
	if(mm.eq.1) then
        call build_output_path('frecordnew.dat', out_path)
        open(1235,file=trim(out_path),status='unknown')
    endif

	if(mm.eq.0) then
		do j=1,nc
			write(1234,3500) sp+j*dx,ff(j,1),ff(j,2),ff(j,3),ff(j,4),ff(j,5),ff(j,6),ff(j,7),ff(j,8),ff(j,9)
		enddo
3500 format(10e25.17)
	endif

	if(mm.eq.1) then
		do j=1,nc
			write(1235,3600) sp+(j-1)*dx,ff(j,1),ff(j,2),ff(j,3),ff(j,4),ff(j,5),ff(j,6),ff(j,7),ff(j,8),ff(j,9)
		enddo
3600 format(10e25.17)
	endif

	if(mm.eq.0) close(1234)
	if(mm.eq.1) close(1235)

end
    
subroutine drawu (n,id1,id2,mesh)
	use com_date
    use meshtype
	implicit none
	
	type(mesh_dat) mesh(max_meshid)
    integer n,id1,id2,i
    character(len=512) :: out_path
    
    call build_output_path('uprofile1.dat', out_path)
    open (2345,file=trim(out_path),status='unknown')
    call build_output_path('uprofile5.dat', out_path)
    open (2346,file=trim(out_path),status='unknown')
    
    do i=1,mesh(id1)%Ny
        write(2345,*) mesh(id1)%oy+dble(i-1)*mesh(id1)%scale, mesh(id1)%u(mesh(id1)%Nx-1,i), mesh(id1)%v(mesh(id1)%Nx-1,i)
    enddo
    
    do i=1,mesh(id2)%Ny
        write(2346,*) mesh(id2)%oy+dble(i-1)*mesh(id2)%scale, mesh(id2)%u(1,i),mesh(id2)%v(1,i)
    enddo
    
    close(2345)
    close(2346)
end
    
subroutine Lagrangian_difference_coarse_to_fine(mesh,n,cycle1)
	use mpi
	use com_date
    use meshtype
	implicit none
	
	type(mesh_dat) mesh(max_meshid)
	integer n,m,j,k,i1,j1,i2,sp,ep,i,cycle1
    Real (kind=ps) temp(Q)

	if(myid.eq.0) then

		m=mesh(n)%left(2,1) !\D7\F3\B2\E0\CD\F8\B8\F1\B1߽\E7\B5Ĵ\AB\B5\DD
		if(mesh(m)%layer.lt.mesh(n)%layer) then !\C8\E7\B9\FBmesh(n)\D7\F3\B1ߵ\C4\CD\F8\B8\F1\B5Ĳ\E3\CA\FDС\D3\DAmesh(n)\B5\C4\CD\F8\B8\F1\B5Ĳ\E3\CA\FD
			do j=2,mesh(n)%Ny-1
                temp=mesh(n)%f(1,j,:)
				mesh(n)%f(1,j,:)= -0.125d0*mesh(n)%leftf2(j,:)+0.75d0*mesh(n)%f(1,j,:)+0.375d0*mesh(n)%f0(1,j,:)
                mesh(n)%leftf2(j,:)=temp
			enddo			
		endif

		m=mesh(n)%right(2,1) !\D3Ҳ\E0\CD\F8\B8\F1\B1߽\E7\B5Ĵ\AB\B5\DD
		if(mesh(m)%layer.lt.mesh(n)%layer) then !\C8\E7\B9\FBmesh(n)\D3ұߵ\C4\CD\F8\B8\F1\B5Ĳ\E3\CA\FDС\D3\DAmesh(n)\B5\C4\CD\F8\B8\F1\B5Ĳ\E3\CA\FD
			do j=2,mesh(n)%Ny-1
                temp=mesh(n)%f(mesh(n)%Nx,j,:)
				mesh(n)%f(mesh(n)%Nx,j,:)= -0.125d0*mesh(n)%rightf2(j,:)+0.75d0*mesh(n)%f(mesh(n)%Nx,j,:)+0.375d0*mesh(n)%f0(mesh(n)%Nx,j,:)
                mesh(n)%rightf2(j,:)=temp
			enddo
		endif

		m=mesh(n)%up(1,1) !\C9ϲ\E0\CD\F8\B8\F1\B1߽\E7\B5Ĵ\AB\B5\DD
		if(mesh(m)%layer.lt.mesh(n)%layer) then !\C8\E7\B9\FBmesh(n)\C9ϱߵ\C4\CD\F8\B8\F1\B5Ĳ\E3\CA\FDС\D3\DAmesh(n)\B5\C4\CD\F8\B8\F1\B5Ĳ\E3\CA\FD
			do j=1,mesh(n)%Nx
                temp=mesh(n)%f(j,mesh(n)%Ny,:)
				mesh(n)%f(j,mesh(n)%Ny,:)= -0.125d0*mesh(n)%upf2(j,:)+0.75d0*mesh(n)%f(j,mesh(n)%Ny,:)+0.375d0*mesh(n)%f0(j,mesh(n)%Ny,:)
                mesh(n)%upf2(j,:)=temp
			enddo
		endif

		m=mesh(n)%down(1,1) !\CF²\E0\CD\F8\B8\F1\B1߽\E7\B5Ĵ\AB\B5\DD
		if(mesh(m)%layer.lt.mesh(n)%layer) then !\C8\E7\B9\FBmesh(n)\CF±ߵ\C4\CD\F8\B8\F1\B5Ĳ\E3\CA\FDС\D3\DAmesh(n)\B5\C4\CD\F8\B8\F1\B5Ĳ\E3\CA\FD
			do j=1,mesh(n)%Nx
                temp=mesh(n)%f(j,1,:)
				mesh(n)%f(j,1,:)= -0.125d0*mesh(n)%downf2(j,:)+0.75d0*mesh(n)%f(j,1,:)+0.375d0*mesh(n)%f0(j,1,:)
                mesh(n)%downf2(j,:)=temp
			enddo
		endif
		
	endif

	if(numprocs.ne.1) then
		if(myid.eq.0) then
			sp=1
			ep=mesh(n)%Nx
		else
			sp=mesh(n)%start_row
			ep=mesh(n)%end_row
		endif
		call scatter(mesh(n)%f,sp,ep,mesh(n)%Nx,mesh(n)%Ny)
	endif
end

Subroutine transfer_fine_to_coarse(mesh,n,cycle1)
	use mpi
	use com_date
    use meshtype
	implicit none
	
	type(mesh_dat) mesh(max_meshid)
	integer n,m,j,k,i1,j1,i2,sp,ep,i,cycle1
	Real (kind=ps), external :: D2Q9feq,distanceq
	Real (kind=ps) temp
	Real (kind=ps), allocatable :: temp1(:,:),temp2(:,:)

	if(myid.eq.0) then
		m=mesh(n)%left(2,1) !\D7\F3\B2\E0\CD\F8\B8\F1\B1߽\E7\B5Ĵ\AB\B5\DD
		if(m.ne.0) then
			if(mesh(m)%layer.ge.mesh(n)%layer) then !\C8\E7\B9\FBmesh(n)\D7\F3\B1ߵ\C4\CD\F8\B8\F1\B5Ĳ\E3\CA\FD\B4\F3\D3\DAmesh(n)\B5\C4\CD\F8\B8\F1\B5Ĳ\E3\CA\FD
				do j=2,mesh(n)%Ny-1
					m=mesh(n)%left(j,1)
					i1=mesh(n)%left(j,3)
					j1=mesh(n)%left(j,4)
					if(mesh(m)%layer.eq.mesh(n)%layer) then
						do k=1,Q 
							mesh(n)%f(1,j,k)=mesh(m)%f(i1,j1,k)
						enddo
                    else
                        call mrttransftc(mesh(n)%f(1,j,:),mesh(m)%rho(i1,j1),mesh(m)%u(i1,j1),mesh(m)%v(i1,j1),mesh(m)%f(i1,j1,:),mesh(m)%s,mesh(n)%s,Q)
					endif
				enddo
			endif
		endif
        
		m=mesh(n)%right(2,1) !\D3Ҳ\E0\CD\F8\B8\F1\B1߽\E7\B5Ĵ\AB\B5\DD
		if(m.ne.0) then
			if(mesh(m)%layer.ge.mesh(n)%layer) then !\C8\E7\B9\FBmesh(n)\D3ұߵ\C4\CD\F8\B8\F1\B5Ĳ\E3\CA\FDС\D3\DAmesh(n)\B5\C4\CD\F8\B8\F1\B5Ĳ\E3\CA\FD
				do j=2,mesh(n)%Ny-1
					m=mesh(n)%right(j,1)
					i1=mesh(n)%right(j,3)
					j1=mesh(n)%right(j,4)
					if(mesh(m)%layer.eq.mesh(n)%layer) then
						do k=1,Q 
							mesh(n)%f(mesh(n)%Nx,j,k)=mesh(m)%f(i1,j1,k)
						enddo
                    else
                        call mrttransftc(mesh(n)%f(mesh(n)%Nx,j,:),mesh(m)%rho(i1,j1),mesh(m)%u(i1,j1),mesh(m)%v(i1,j1),mesh(m)%f(i1,j1,:),mesh(m)%s,mesh(n)%s,Q)
					endif
				enddo
			endif
		endif

		m=mesh(n)%up(1,1) !\C9ϲ\E0\CD\F8\B8\F1\B1߽\E7\B5Ĵ\AB\B5\DD
		if(m.ne.0) then
			if(mesh(m)%layer.ge.mesh(n)%layer) then !\C8\E7\B9\FBmesh(n)\C9ϱߵ\C4\CD\F8\B8\F1\B5Ĳ\E3\CA\FD\B4\F3\D3ڵ\C8\D3\DAmesh(n)\B5\C4\CD\F8\B8\F1\B5Ĳ\E3\CA\FD
				do j=1,mesh(n)%Nx
					m=mesh(n)%up(j,1)
					i1=mesh(n)%up(j,3)
					j1=mesh(n)%up(j,4)
					if(mesh(m)%layer.eq.mesh(n)%layer) then
						do k=1,Q 
							mesh(n)%f(j,mesh(n)%Ny,k)=mesh(m)%f(i1,j1,k)
						enddo	
                    else
                        call mrttransftc(mesh(n)%f(j,mesh(n)%Ny,:),mesh(m)%rho(i1,j1),mesh(m)%u(i1,j1),mesh(m)%v(i1,j1),mesh(m)%f(i1,j1,:),mesh(m)%s,mesh(n)%s,Q)
					endif
				enddo				
			endif
		endif

		m=mesh(n)%down(1,1) !\CF²\E0\CD\F8\B8\F1\B1߽\E7\B5Ĵ\AB\B5\DD
		if(m.ne.0) then
			if(mesh(m)%layer.ge.mesh(n)%layer) then !\C8\E7\B9\FBmesh(n)\CF±ߵ\C4\CD\F8\B8\F1\B5Ĳ\E3\CA\FD\B4\F3\D3ڵ\C8\D3\DAmesh(n)\B5\C4\CD\F8\B8\F1\B5Ĳ\E3\CA\FD
				do j=1,mesh(n)%Nx
					m=mesh(n)%down(j,1)
					i1=mesh(n)%down(j,3)
					j1=mesh(n)%down(j,4)
					if(mesh(m)%layer.eq.mesh(n)%layer) then
						do k=1,Q 
							mesh(n)%f(j,1,k)=mesh(m)%f(i1,j1,k)
						enddo	
                    else
                        call mrttransftc(mesh(n)%f(j,1,:),mesh(m)%rho(i1,j1),mesh(m)%u(i1,j1),mesh(m)%v(i1,j1),mesh(m)%f(i1,j1,:),mesh(m)%s,mesh(n)%s,Q)
					endif
				enddo				
			endif
		endif
	endif

	
	
	if(numprocs.ne.1) then
		if(myid.eq.0) then
			sp=1
			ep=mesh(n)%Nx
		else
			sp=mesh(n)%start_row
			ep=mesh(n)%end_row
		endif
		call scatter(mesh(n)%f,sp,ep,mesh(n)%Nx,mesh(n)%Ny)
	endif
end
    
subroutine nizhen(AA,b1,Nm)
    use com_date
    implicit none
	integer i,j,k,ik,Nm
	real (kind=ps) A(Nm,Nm),b1(Nm,Nm),AA(Nm,Nm)
	real (kind=ps) temp,mik
    
    A=AA
    b1=0.0d0
    do k=1,Nm
        b1(k,k)=1.0d0
    enddo
    
	do k=1,Nm-1
		ik=k
		temp=abs(A(ik,k))
		do i =k,Nm
			if(abs(a(i,k)).gt.abs(temp)) then
				temp=abs(a(i,k))
				ik=i
			endif
		enddo
		do j=1,Nm
			temp=a(ik,j)
			a(ik,j)=a(k,j)
			a(k,j)=temp
            temp=b1(ik,j)
			b1(ik,j)=b1(k,j)
			b1(k,j)=temp
		enddo

		do i=1,Nm
			mik=a(i,k)/a(k,k)
            temp=a(k,k)
			do j=k,Nm
				if(i.ne.k) a(i,j)=a(i,j)-mik*a(k,j)
            enddo
            a(k,k:Nm)=a(k,k:Nm)/temp
            do j=1,Nm
                if(i.ne.k) b1(i,j)=b1(i,j)-mik*b1(k,j)
            enddo
            b1(k,:)=b1(k,:)/temp
		enddo
    enddo
    
    temp=a(Nm,Nm)
    do i=1,Nm-1
        mik=a(i,Nm)/a(Nm,Nm)
        a(i,Nm)=0.0d0
        do j=1,Nm
            b1(i,j)=b1(i,j)-mik*b1(Nm,j)
        enddo
    enddo
    a(Nm,Nm)=1.0d0
    b1(Nm,:)=b1(Nm,:)/temp
end
    
subroutine mrttransctf(res,rhoc,uc,vc,fc,sc,sf,n)
    use com_date
    implicit none
    integer n,j
    real (kind=ps) res(n),rhoc,uc,vc,fc(n),sc(n),sf(n)
    real (kind=ps) reqc(n),item1(n),item2(n),tempitem1(n),tempitem2(n)
    
    reqc(1)=rhoc
    reqc(2)=-2.0d0*rhoc+3.0d0*(uc*uc+vc*vc)*rhoc
    reqc(3)=rhoc-3.0d0*(uc*uc+vc*vc)*rhoc
    reqc(4)=uc*rhoc
    reqc(5)=-uc*rhoc
    reqc(6)=vc*rhoc
    reqc(7)=-vc*rhoc
    reqc(8)=(uc*uc-vc*vc)*rhoc
    reqc(9)=uc*vc*rhoc
    call matmulvoc(item1,inversematm,reqc,n,n)
    call matmulvoc(tempitem1,matm,fc,n,n)
    tempitem1=tempitem1-reqc

    do j=1,n
        if((j.eq.8).or.(j.eq.9))then
            tempitem2(j)=(1.0d0-sf(j))*(0.5d0*sc(j)/sf(j))/(1.0d0-sc(j))
        else
            tempitem2(j)=(1.0d0-sf(j))/(1.0d0-sc(j))
        endif
        tempitem1(j)=tempitem1(j)*tempitem2(j)
    enddo
    call matmulvoc(item2,inversematm,tempitem1,n,n)
    res=item1+item2

end 
    
subroutine mrttransftc(res,rhof,uf,vf,ff,sf,sc,n)
    use com_date
    implicit none
    integer n,j
    real (kind=ps) res(n),rhof,uf,vf,ff(n),sc(n),sf(n)
    real (kind=ps) reqf(n),item1(n),item2(n),tempitem1(n),tempitem2(n)
    
    reqf(1)=rhof
    reqf(2)=-2.0d0*rhof+3.0d0*(uf*uf+vf*vf)*rhof
    reqf(3)=rhof-3.0d0*(uf*uf+vf*vf)*rhof
    reqf(4)=uf*rhof
    reqf(5)=-uf*rhof
    reqf(6)=vf*rhof
    reqf(7)=-vf*rhof
    reqf(8)=(uf*uf-vf*vf)*rhof
    reqf(9)=uf*vf*rhof
    call matmulvoc(item1,inversematm,reqf,n,n)
    call matmulvoc(tempitem1,matm,ff,n,n)
    tempitem1=tempitem1-reqf
    
!    if(tempitem1(1).gt.1.0d-6) then
!        print*, tempitem1(1)
!        pause
!    endif
    
    do j=1,n
        if((j.eq.8).or.(j.eq.9))then
            tempitem2(j)=(1.0d0-sc(j))*(2.0d0*sf(j)/sc(j))/(1.0d0-sf(j))
        else
            tempitem2(j)=(1.0d0-sc(j))/(1.0d0-sf(j))
        endif
        tempitem1(j)=tempitem1(j)*tempitem2(j)
    enddo
    call matmulvoc(item2,inversematm,tempitem1,n,n)
    res=item1+item2

end 
    
subroutine integration_real_a (dt,xref)
    use com_date
    implicit none
    real (kind=ps) dt,xref,tem
    
    body_rel_v_prev(1)=body_rel_v(1)
    body_rel_v(1)=body_rel_v(1)+0.5*(body_rel_a_prev(1)+body_rel_a(1))
    body_rel_s_prev(1)=body_rel_s(1)
    body_rel_s(1)=body_rel_s(1)+0.5*(body_rel_v_prev(1)+body_rel_v(1))
    body_rel_pos(1)=body_rel_pos(1)+0.5*(body_rel_v_prev(1)+body_rel_v(1))
    call update_body_extra_from_rel()
    
end
    
subroutine saveflowfield(cycle1,mesh)
    use nested_io
    implicit none
    integer :: cycle1
    type(mesh_dat) :: mesh(max_meshid)
    call save_nested_checkpoint(cycle1,mesh)
end subroutine

subroutine readflowfield(cycle1,mesh)
    use nested_io
    implicit none
    integer :: cycle1
    type(mesh_dat) :: mesh(max_meshid)
    call read_nested_checkpoint(cycle1,mesh)
end subroutine

subroutine vor_mom_rate(cycle1,mesh)
	use com_date
    use meshtype
	implicit none
	
	type(mesh_dat) mesh(max_meshid)
    integer cycle1,n,i,j,k,m(4),is,ie,js,je
    real (kind=ps) dvx,duy,xr,yr,temp1(4),temp2(4),temp3(4),tempu,tempvor,temptwoder,tem,tem1,tem2

    if (numprocs.ne.1) then
		do n=1,max_meshid		
			if (myid.eq.0) then
				call gather (mesh(n)%u,1,mesh(n)%Nx,mesh(n)%Nx,mesh(n)%Ny,mesh(n)%start_row,mesh(n)%end_row)
				call gather (mesh(n)%v,1,mesh(n)%Nx,mesh(n)%Nx,mesh(n)%Ny,mesh(n)%start_row,mesh(n)%end_row)
				call gather (mesh(n)%pre,1,mesh(n)%Nx,mesh(n)%Nx,mesh(n)%Ny,mesh(n)%start_row,mesh(n)%end_row)
			else
				call gather (mesh(n)%u,mesh(n)%start_row,mesh(n)%end_row,mesh(n)%Nx,mesh(n)%Ny,mesh(n)%start_row,mesh(n)%end_row)
				call gather (mesh(n)%v,mesh(n)%start_row,mesh(n)%end_row,mesh(n)%Nx,mesh(n)%Ny,mesh(n)%start_row,mesh(n)%end_row)
				call gather (mesh(n)%pre,mesh(n)%start_row,mesh(n)%end_row,mesh(n)%Nx,mesh(n)%Ny,mesh(n)%start_row,mesh(n)%end_row)
			endif			
    	enddo
	endif
	
    if(myid.eq.0.and.periodical_bc.eq.1)then
		!$acc kernels present(mesh)
        mesh(1)%u(mesh(1)%Nx,2:mesh(1)%Ny-1) =mesh(1)%u(2,2:mesh(1)%Ny-1)
        mesh(1)%v(mesh(1)%Nx,2:mesh(1)%Ny-1) =mesh(1)%v(2,2:mesh(1)%Ny-1)
        mesh(1)%u(1,2:mesh(1)%Ny-1) =mesh(1)%u(mesh(1)%Nx-1,2:mesh(1)%Ny-1)
        mesh(1)%v(1,2:mesh(1)%Ny-1) =mesh(1)%v(mesh(1)%Nx-1,2:mesh(1)%Ny-1)
		!$acc end kernels
    endif
    !print*, "(1.5d0*(liy1(1)+liy2(1))-2.0d0*(liy1(2)+liy2(2))+0.5d0*(liy1(3)+liy2(3)))/(UU*mesh(layer1(n))%scale)=",(1.5d0*(liy1(1)+liy2(1))-2.0d0*(liy1(2)+liy2(2))+0.5d0*(liy1(3)+liy2(3)))/(UU*mesh(layer1(n))%scale)
    
!    if(myid.eq.0) print*, mesh(1)%v(mesh(1)%Nx,722)    
    if(myid.eq.0) then
		do k=1,max_meshid
            is=2
			ie=mesh(k)%Nx-1
			js=2
			je=mesh(k)%Ny-1
			!$acc kernels present(mesh) create(dvx,duy,tem,tem1,tem2)
			!$acc loop independent
			do j=js,je
				!$acc loop independent
				do i=is,ie
					tem1=mesh(k)%vor1(i,j)
					tem2=mesh(k)%vor(i,j)
					mesh(k)%vor2(i,j)=tem1
					mesh(k)%vor1(i,j)=tem2
                    dvx=(mesh(k)%v(i+1,j)-mesh(k)%v(i-1,j))/(UU*2.0*mesh(k)%scale)	
                    duy=(mesh(k)%u(i,j+1)-mesh(k)%u(i,j-1))/(UU*2.0*mesh(k)%scale)
                    tem=(dvx-duy) 
					mesh(k)%vor(i,j)=tem                  
                    mesh(k)%dvordt(i,j)=(1.5d0*tem-2.0d0*tem2+0.5d0*tem1)/(UU*mesh(1)%scale)
				enddo
				!$acc end loop
            enddo
			!$acc end loop
			!$acc end kernels
		enddo

2000	format(25e25.17)
2001	format(4i8.6)
	endif
!	close(outfile)
end

subroutine cal_sur_cp(mesh,cp_sur)
    use mpi
	use com_date
    use meshtype
	implicit none
	
	type(mesh_dat) mesh
	integer sp,ep,ip,jp,i,j,k,k1,num,temn,smooth,spp,epp
	real(kind=ps),  allocatable ::  ub(:),vb(:),ustar(:),vstar(:),rhos(:),ustar1(:),vstar1(:),rhos1(:),rhos2(:),buf1(:),buf2(:),as(:),ul(:),vl(:),glx(:),gly(:)
	real(kind=ps),  allocatable :: tem1(:,:),tem2(:,:),aa(:,:),aa1(:,:),bu(:),bv(:),tem11(:,:),tem22(:,:),tema(:,:),comi(:),x(:,:),xx1(:,:),xx2(:,:),xx3(:,:)
	real*8 disr,x0,y0,kernel1,kernel2,sum1,sum2,cp_sur(b_point-1)
    integer nono,mn,nn,js,js1,ds,nnn
	integer,  allocatable :: non_zero(:),indexnum(:,:),recvcounts(:),displs(:),mrpos(:),nrpos(:),mrow(:),mcol(:),mrrow(:),mrcol(:),ia(:),ja(:)
    logical if_rk4 !是否采用龙格库塔方法求解动力学方程
    character(len=512) :: out_path

    if_rk4=.true.

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

    !$acc update host(xb2,yb2,adjoint_p(1:adjoint_n,:))
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
                rhos1(i)=mesh%pre(adjoint_p(i,1),adjoint_p(i,2))
            endif
        enddo

        call MPI_ALLREDUCE(rhos1,rhos,num,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

    else
        do i=1,adjoint_n       
            rhos(i)=mesh%pre(adjoint_p(i,1),adjoint_p(i,2))
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

!    if(numprocs.ne.1)then 
!        call multsmatrix(xx1,xx2,mrpos,nrpos,b_point-1,num,js1,as,ia,ja,ds)
!    else
!        call multsmatrix(xx1,xx2,mrpos,nrpos,b_point-1,num,js,as,ia,ja,ds)
!    endif

    if(numprocs.ne.1)then 
        call matmulvocxs(cp_sur,xx1,mrpos,rhos,b_point-1,num,js1)
    else
        call matmulvocxs(cp_sur,xx1,mrpos,rhos,b_point-1,num,js)
    endif
end    
    
