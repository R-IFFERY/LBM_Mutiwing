!二维格子boltzmann程序
!作者：王逗
!2014-03-25
module com_date
	implicit none 
	integer Q,ps !Q：离散速度分量数目，依赖于LBM模型；Nx,Ny为x和y方向的格子数
	parameter (ps=8) !ps=4,单精度；ps=8,双精度
	real (kind=ps) pai
	parameter (pai=3.1415926535897932384626433832795d0)
	parameter (Q=9) 
	real (kind=ps) UU !最外层的流体格子速度
	parameter (UU=0.01d0)
	integer e(Q,2)
	real (kind=ps) w(Q)
		integer max_wings_supported,wing_count,motion_columns_per_wing,max_motion_columns !翼数量相关参数
		integer boundary_points_per_wing,boundary_segments_per_wing,max_boundary_points,max_boundary_segments
		parameter (max_wings_supported=8)
		parameter (wing_count=1)	!仿真翼数量
		parameter (motion_columns_per_wing=9)
		parameter (max_motion_columns=1+motion_columns_per_wing*max_wings_supported)
		parameter (boundary_points_per_wing=205)
		parameter (boundary_segments_per_wing=boundary_points_per_wing-1)
		parameter (max_boundary_points=boundary_points_per_wing*max_wings_supported)
		parameter (max_boundary_segments=boundary_segments_per_wing*max_wings_supported)
		real (kind=ps) myforce(2),force(2),force1(2),force2(2),mom,rho0,err,mom2,pm1,pm2
		real (kind=ps) wing_force(max_wings_supported,2),wing_moment(max_wings_supported)
	integer bc_num !选择物面边界条件的类型:1.非平衡外推(Guo-Zheng-Shi),2.插值反弹格式(Lallemand-Luo),3.插值反弹格式(Yu-Mei-shyy),4.浸润边界法(基于Wang et al.)
	parameter (bc_num=4)
	integer max_meshid
	parameter (max_meshid=4)
			real (kind=ps) mesh_max_scale !最外层网格尺寸
			parameter (mesh_max_scale=0.08d0)
    real(kind=ps), parameter :: mesh_min_scale=mesh_max_scale/2**(max_meshid-1)
    real(kind=ps), parameter :: grid_xmin(max_meshid)=(/-20d0,-16d0,-12d0,-8d0/)
    real(kind=ps), parameter :: grid_xmax(max_meshid)=(/ 20d0, 16d0, 12d0, 8d0/)
    real(kind=ps), parameter :: grid_ymin(max_meshid)=(/-10d0, -8d0, -6d0,-4d0/)
    real(kind=ps), parameter :: grid_ymax(max_meshid)=(/ 10d0,  8d0,  6d0, 4d0/)
	        real (kind=ps) wing_initial_x_origin,wing_initial_x_spacing,wing_initial_y_origin
	        parameter (wing_initial_x_origin=0.0d0)
	        parameter (wing_initial_x_spacing=1.5d0)
	        parameter (wing_initial_y_origin=0.0d0)
		integer output_stride_x,output_stride_y !流场输出抽样步长，仅影响输出文件密度，不影响仿真网格
		parameter (output_stride_x=2,output_stride_y=2)
	integer cycles !开始记录力与运动的时间步
	real (kind=ps) time_s1,time_e1,time_s2,time_e2,time_s3,time_e3,time_s4,time_e4,time_s5,time_e5,time_s6,time_e6,time_s7,time_e7,time_s8,time_e8,time_s9,time_e9,time_s10,time_e10,time_s11,time_e11
		!$acc declare create (e,w,myforce,force,force1,mom,pm1,force2,mom2,pm2,rho0,err,cycles,time_s1,time_e1,wing_force,wing_moment)
		!$acc declare copyin (Q,ps,pai,UU,bc_num,max_meshid,mesh_max_scale,wing_initial_x_origin,wing_initial_x_spacing,wing_initial_y_origin,output_stride_x,output_stride_y,max_wings_supported,wing_count,motion_columns_per_wing,max_motion_columns,boundary_points_per_wing,boundary_segments_per_wing,max_boundary_points,max_boundary_segments)
	!并行计算时的进程编号
	integer myid,numprocs,leftpro,rightpro,ierr,rc
	!$acc declare create (myid,numprocs,leftpro,rightpro,ierr,rc)

	!运动函数参数
	real (kind=ps) h0,h01,alphad,alphau,alphad1,alphau1,theta0,f1,phi,alpha,Re,time
		real (kind=ps) motion_end_time
	    real (kind=ps), allocatable :: motion_table(:,:)
	    real (kind=ps), allocatable :: motion_step_table(:,:)
	    character(len=256) motion_file
	    character(len=256) motion_source_name
	    character(len=256) zmq_endpoint
	    character(len=256) output_root
		    character(len=256) output_dir
		    character(len=256) episode_flow_output_dir
		    character(len=256) restart_input_dir
		    real (kind=ps) flow_save_k
			    integer motion_points
		    integer motion_step_count
		    integer motion_index
		    integer motion_column_count,motion_stream_value_count
		    integer motion_source_file,motion_source_zmq,motion_source_kind
		    integer zmq_cmd_invalid,zmq_cmd_reset,zmq_cmd_step,zmq_cmd_stop,zmq_cmd_output
		    integer zmq_flow_output_interval
		    real (kind=ps) motion_dt
		    logical motion_loaded
		    logical motion_stream_ready
		    integer latest_force_row_len,latest_force_sample_count,episode_index
		    real (kind=ps) motion_stream_start_row(max_motion_columns),motion_stream_target_row(max_motion_columns)
		    real (kind=ps) latest_force_row(4+12*max_wings_supported)
		    real (kind=ps) latest_force_sum_row(4+12*max_wings_supported)
		    integer ref_freestream !1：以自由来流做参考速度，0：不以自由来流做参考速度
		real (kind=ps) sb_arr(max_wings_supported),sbv_arr(max_wings_supported),sbac_arr(max_wings_supported)
		real (kind=ps) sa_arr(max_wings_supported),sau_arr(max_wings_supported),saac_arr(max_wings_supported)
		real (kind=ps) angle_arr(max_wings_supported),omega_arr(max_wings_supported),omgac_arr(max_wings_supported)
		real (kind=ps) body_rel_a(max_wings_supported),body_rel_v(max_wings_supported),body_rel_s(max_wings_supported),body_rel_pos(max_wings_supported)
		real (kind=ps) body_rel_theta(max_wings_supported),body_rel_temp_theta(max_wings_supported)
		real (kind=ps) body_rel_a_prev(max_wings_supported),body_rel_v_prev(max_wings_supported),body_rel_s_prev(max_wings_supported)
		real (kind=ps) body_rel_theta_prev(max_wings_supported),body_rel_temp_theta_prev(max_wings_supported)
		real (kind=ps) body_extra_v(max_wings_supported),body_extra_s(max_wings_supported),body_extra_pos(max_wings_supported)
	real (kind=ps) rr !圆半径
	real (kind=ps) xc,yc !流场中心坐标
    real (kind=ps) xp,yp !翼型旋转轴的坐标

	!参数的初始化
	    parameter (ref_freestream = 1)
	    parameter (motion_source_file = 1)
	    parameter (motion_source_zmq = 2)
		    parameter (zmq_cmd_invalid = 0)
		    parameter (zmq_cmd_reset = 1)
		    parameter (zmq_cmd_step = 2)
		    parameter (zmq_cmd_stop = 3)
		    parameter (zmq_cmd_output = 4)
	    parameter (h0=0.5d0)
	parameter (h01=0.5d0)
	parameter (rr=0.5d0)	
	parameter (alphad=10.00d0*pai/180.0d0)
	parameter (alphau=10.0d0*pai/180.0d0)
	parameter (alphad1=10.00d0*pai/180.0d0)
	parameter (alphau1=10.0d0*pai/180.0d0)
	parameter (theta0=45.0d0*pai/180.0d0)
	parameter (f1=0.5d0)
	parameter (phi=0.5d0*pai)
	parameter (alpha=0d0)
	parameter (Re=500.0d0)
	parameter (xc=0.0d0,yc=0.0d0)
	parameter (xp=0.0d0,yp=0.0d0)
!    parameter (xp1=15.250d0,yp1=7.75d0)
		!$acc declare create (sb_arr,sbv_arr,sbac_arr,sa_arr,sau_arr,saac_arr,angle_arr,omega_arr,omgac_arr)
		!$acc declare create (body_rel_a,body_rel_v,body_rel_s,body_rel_pos,body_rel_theta,body_rel_temp_theta)
		!$acc declare create (body_rel_a_prev,body_rel_v_prev,body_rel_s_prev,body_rel_theta_prev,body_rel_temp_theta_prev)
		!$acc declare create (body_extra_v,body_extra_s,body_extra_pos)
		!$acc declare copyin (ref_freestream,h0,h01,rr,alphad,alphau,alphad1,alphau1,Re,time,xc,yc,xp,yp)

	!记录多块网格每层的网格编号
	integer layer1(1)
	!浸润边界法参数
	real (kind=ps),allocatable :: gforce(:,:,:)  !gforce:体积力密度
		integer b_point1,b_point,b_segment_count !边界点点数
		integer wing_point_start(max_wings_supported),wing_segment_start(max_wings_supported)
	    logical closed_boundary
	    parameter (closed_boundary=.true.)
		logical solve_dynamic_equations
		parameter (solve_dynamic_equations=.false.)
		real (kind=ps) xb0(max_boundary_points),yb0(max_boundary_points)
		real (kind=ps) xb1(max_boundary_points),yb1(max_boundary_points)
	    real (kind=ps) xb2(max_boundary_segments),yb2(max_boundary_segments)
		real (kind=ps) arc_l(max_boundary_segments)
		real (kind=ps) comj(max_boundary_points,2)
	    integer adjoint_n,adjoint_p(36*max_boundary_points,2),adjoint_wing_id(36*max_boundary_points) !记录边界关联点的数目与序号
	    !$acc declare create (layer1,gforce,xb0,yb0,xb1,yb1,xb2,yb2,arc_l,comj,adjoint_n,adjoint_p,adjoint_wing_id,wing_point_start,wing_segment_start,b_point1,b_point,b_segment_count)
		!$acc declare copyin (closed_boundary,solve_dynamic_equations)

    !多松弛法引入的参数
    real (kind=ps) matm(9,9)
    real (kind=ps) inversematm(9,9)
	!$acc declare create (matm,inversematm)
    
    !动力学耦合求解研究引入的参数
    real (kind=ps) ifdynamic !是否引入动力学方程
    real (kind=ps) solid_den,ar_airfoil,area_airfoil,I_airfoil,fbi,par1,par2
    integer rightmove(max_meshid),leftmove(max_meshid),movenum,periodical_bc
    parameter (ifdynamic=0) !是否求解动力学方程
    parameter (solid_den= 10.0d0) !翼的密度
    parameter (ar_airfoil= 0.1d0) !椭圆翼型的长宽比,椭圆翼型时有效
    parameter (area_airfoil=  6.807946476488794d-2)!0.25d0*pai*ar_airfoil) !椭圆的面积,椭圆翼型时，用该注释前面的表达式
    parameter (I_airfoil= 0.01565)!0.25d0*pai*ar_airfoil) !椭圆的面积,椭圆翼型时，用该注释前面的表达式6.807946476488794d-2
    parameter (fbi=10.0d0) !无量纲弹簧刚度
    parameter (par1=0.0d0) !-fbi/(solid_den*I_airfoil))
    parameter (par2=0.0d0) !1.0d0/(solid_den*I_airfoil))
	parameter (periodical_bc=0) !出入流口是否使用周期边界条件
    !$acc declare create (rightmove,leftmove,movenum)
	!$acc declare copyin (ifdynamic,solid_den,ar_airfoil,area_airfoil,I_airfoil,fbi,par1,par2,periodical_bc)

    !是否续算
    logical goon
    !$acc declare create (goon)
end module
    
module meshtype
    implicit none
    type mesh_dat
	    integer meshid
	    integer layer	!网格的层数，最外层为1
	    integer Nx,Ny	!x向与y向网格的数目
        ! Post-collision state at the start of the current level step.
        real(kind=8), allocatable :: f_old(:,:,:)
	    real (kind=8),allocatable :: rho(:,:),u(:,:),v(:,:),u0(:,:),v0(:,:),pre(:,:),f(:,:,:),f0(:,:,:)	! rho：密度，u：x方向速度，v：y方向速度，f：速度分布函数
        real (kind=8),allocatable :: vor(:,:),vor1(:,:),vor2(:,:),dvordt(:,:) 	! 涡量时间变化率
	    real (kind=8) scale,c,dx,dy,Lx,Ly,dt,p0,tau_f,niu,err !本网格基本参数
	    real (kind=8) ox,oy	!网格的原点
	    integer ,allocatable :: body(:,:) !1.固体的内点，2.固体的边界点，3.流体的内点，4.流体的边界点
	    integer start_row,end_row	!并行计算所需参数
	    integer ,allocatable :: left(:,:),right(:,:),up(:,:),down(:,:)	!定义四条边线
	    real (kind=8) ,allocatable :: leftf2(:,:),rightf2(:,:),upf2(:,:),downf2(:,:) !储存四边界上一步的分布函数
	    integer nodescounts,facecounts,start_num !储存该节点拥有的节点数,面数以及开始的节点号
        real (kind=8) ,allocatable :: s(:) !松弛因子矩阵
    end type mesh_dat
    
end module

Program LBM
	use MPI
	use openacc
	use cublas
    use cusolverDn
	use com_date
!    USE DFPORT, ONLY: DTIME
    use meshtype
!    include 'mpif.h'
	implicit none
	

	type(mesh_dat) mesh(max_meshid)
		integer i,j,k,ip,jp,n,nmax,nforce,flag1,flag2,stat
		integer gpu_count,gpu_id
		integer save_second,nsave
	    real (kind=ps) cpu_ts,cpu_ts1,cpu_es1,cpu_es2,cpu_es3,ansstime
		real (kind=ps) streaming_time,ph_time,colission_time,bc_time,bc1_time,bc2_time,bc3_time,move_mesh_time,move_airfoil_time,outer_bc_time,vor_time
	    real (kind=ps) solver_dt
	    character(len=512) :: out_path
		type(cusolverDnHandle) :: h

	call MPI_INIT(ierr)
	call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
	call MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,ierr)

	gpu_count=acc_get_num_devices(acc_device_nvidia)
	if (gpu_count.le.0) then
		if (myid.eq.0) write(*,*) 'No NVIDIA GPU is available to OpenACC'
		call MPI_ABORT(MPI_COMM_WORLD,900,ierr)
	endif
	gpu_id=mod(myid,gpu_count)
	call acc_set_device_num(gpu_id,acc_device_nvidia)
	call acc_init(acc_device_nvidia)
	stat=cusolverDnCreate(h)
	if (stat.ne.0) then
		write(*,*) 'Rank',myid,'cusolverDnCreate failed: status=',stat
		call MPI_ABORT(MPI_COMM_WORLD,903,ierr)
	endif

	if(myid.gt.0) then
		leftpro=myid-1
	else
		leftpro=MPI_PROC_NULL
	endif
	if(myid.lt.numprocs-1) then
		rightpro=myid+1
	else
		rightpro=MPI_PROC_NULL
	endif
	print*,"Process",myid,"of",numprocs,"is alive"
	
    goon=.false.
    
	    motion_loaded=.false.
	    motion_points=0
	    motion_step_count=0
	    motion_index=1
	    motion_column_count=1+wing_count*motion_columns_per_wing
	    motion_stream_value_count=1+wing_count*3
	    motion_dt=0.0d0
	    motion_end_time=0.0d0
	    if ((wing_count.lt.1).or.(wing_count.gt.max_wings_supported)) then
	        if (myid.eq.0) write(*,*) 'wing_count 超出允许范围: ', wing_count
	        call MPI_FINALIZE(rc)
	        stop 1
	    endif
	    b_point1=boundary_points_per_wing
	    b_point=wing_count*boundary_points_per_wing
	    b_segment_count=wing_count*boundary_segments_per_wing
	    do i=1,max_wings_supported
	        wing_point_start(i)=0
	        wing_segment_start(i)=0
	    enddo
	    do i=1,wing_count
	        wing_point_start(i)=1+(i-1)*boundary_points_per_wing
	        wing_segment_start(i)=1+(i-1)*boundary_segments_per_wing
	    enddo
	    sb_arr=0.0d0
	    sbv_arr=0.0d0
	    sbac_arr=0.0d0
	    sa_arr=0.0d0
	    sau_arr=0.0d0
	    saac_arr=0.0d0
	    angle_arr=0.0d0
	    omega_arr=0.0d0
	    omgac_arr=0.0d0
	    body_rel_a=0.0d0
	    body_rel_v=0.0d0
	    body_rel_s=0.0d0
	    body_rel_pos=0.0d0
	    body_rel_theta=0.0d0
	    body_rel_temp_theta=0.0d0
	    body_rel_a_prev=0.0d0
	    body_rel_v_prev=0.0d0
	    body_rel_s_prev=0.0d0
	    body_rel_theta_prev=0.0d0
	    body_rel_temp_theta_prev=0.0d0
	    body_extra_v=0.0d0
	    body_extra_s=0.0d0
	    body_extra_pos=0.0d0
	    wing_force=0.0d0
	    wing_moment=0.0d0
	    call ini_mesh(max_meshid,mesh)
	print*, "myid=",myid,"start from",mesh(1)%start_row,"to",mesh(1)%end_row
!	call space_distribution 	
		!$acc enter data copyin (mesh)
	call init_motion_runtime()
    if(myid.eq.0) call init_output_dir()
    call MPI_BCAST(output_root,len(output_root),MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(output_dir,len(output_dir),MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(restart_input_dir,len(restart_input_dir),MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(motion_file,len(motion_file),MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
    if(myid.eq.0.and.motion_source_kind.eq.motion_source_file) print*,"Input motion file:",trim(motion_file)
    if(myid.eq.0) print*,"Output directory:",trim(output_dir)
	call MPI_BARRIER(MPI_COMM_WORLD,ierr)
	layer1(1)=max_meshid
	if (motion_source_kind.eq.motion_source_file) then
	    call load_motion_table()
	    call build_motion_step_table(mesh_min_scale*UU)
			nmax=10 !最大时间迭代步数
		cycles=0
			if(goon==.false.)then
				nforce=101
		        call build_output_path('force.dat', out_path)
				if(myid.eq.0) open(nforce,file=trim(out_path),status='replace')
				nforce=102
		        call build_output_path('airfoil_motion.dat', out_path)
				if(myid.eq.0) open(nforce,file=trim(out_path),status='replace')
				nforce=101	
			else
	            if(myid.eq.0) then
	                call build_restart_input_path('force.dat', out_path)
	                open(103,file=trim(out_path),status='old',iostat=ierr)
	            endif
				if(myid.eq.0.and.ierr.eq.0)then	
					do !while (not(eof(nforce)))
						read(103,*,iostat=ierr)
						if(ierr/=0) exit
						cycles=cycles+1
		            ENDDO
	                close(103)
					endif
					nforce=101
			        call build_output_path('force.dat', out_path)
					if(myid.eq.0) open(nforce,file=trim(out_path),status='old',position='append')
				    nforce=102
			        call build_output_path('airfoil_motion.dat', out_path)
					if(myid.eq.0) open(nforce,file=trim(out_path),status='old',position='append')
				    nforce=101
				endif	

	    n=1
		call init(mesh)
		call read_boundary(mesh(max_meshid))
		    if(goon.eq..false.)then
			    !call move_boundary(mesh(max_meshid))
		        body_rel_a=0.0d0
		        body_rel_a_prev=0.0d0
		        body_rel_v=0.0d0
		        body_rel_v_prev=0.0d0
		        body_rel_s=0.0d0
		        body_rel_s_prev=0.0d0
		        body_rel_pos=0.0d0
		        body_rel_theta=0.0d0
		        body_rel_theta_prev=0.0d0
		        body_rel_temp_theta=0.0d0
		        body_rel_temp_theta_prev=0.0d0
				mom=0.0d0
				mom2=0.0d0
				pm1=0.0d0
				pm2=0.0d0
		        body_extra_v=0.0d0
		        body_extra_s=0.0d0
		        body_extra_pos=0.0d0
		        wing_force=0.0d0
		        wing_moment=0.0d0
		        call update_body_extra_from_rel()

			    n=1
		        movenum=0
		    else
		        call readflowfield(n,mesh)
		        call update_body_extra_from_rel()
		        call output(n,mesh)
		        n=n+1
		    endif
	    flag1=0
	    flag2=1
		comj=0.0d0
		!$acc update device (mesh)
		!$acc update device (myid,numprocs,leftpro,rightpro,layer1,matm,inversematm,e,w,gforce,comj,wing_force,wing_moment) 
			!$acc update device (b_point1,b_point,b_segment_count,wing_point_start,wing_segment_start,adjoint_wing_id)
		!$acc update device (sb_arr,sbv_arr,sbac_arr,sa_arr,sau_arr,saac_arr,angle_arr,omega_arr,omgac_arr)
		!$acc update device (body_rel_a,body_rel_v,body_rel_s,body_rel_pos,body_rel_theta,body_rel_temp_theta)
		!$acc update device (body_rel_a_prev,body_rel_v_prev,body_rel_s_prev,body_rel_theta_prev,body_rel_temp_theta_prev)
		!$acc update device (body_extra_v,body_extra_s,body_extra_pos)
		!$acc update device (xb0,yb0,xb1,yb1,xb2,yb2,arc_l,adjoint_n,adjoint_p)
		!$acc update device (mom,mom2,movenum,pm1,pm2)
		ansstime=0.0d0
		streaming_time=0.0d0
		ph_time=0.0d0
		colission_time=0.0d0
		bc_time=0.0d0
		bc1_time=0.0d0
		bc2_time=0.0d0
		bc3_time=0.0d0
		move_airfoil_time=0.0d0
		move_mesh_time=0.0d0
		outer_bc_time=0.0d0
		vor_time=0.0d0
        solver_dt=mesh_min_scale*UU
        if(myid.eq.0) print*, 'Offline flow output times: 30,31,32,33,34,35'
			call cpu_time(cpu_ts)
			do !n =1,1000
				time=dble(n-1)*(mesh_min_scale*UU)
			if(time.gt.motion_end_time+1.0d-12) exit
			call MPI_BARRIER(MPI_COMM_WORLD,ierr)
			call cpu_time(cpu_ts1)
			if(n.ge.2) call erranalysis(mesh,n)
			call cpu_time(cpu_es1)
			ansstime=ansstime+cpu_es1-cpu_ts1
			call evolution(mesh,n,nforce,h)
	!		print*,n
			streaming_time=streaming_time+time_e1-time_s1
			ph_time=ph_time+time_e2-time_s2
			bc_time=bc_time+time_e3-time_s3
			bc1_time=bc1_time+time_e5-time_s5
			bc2_time=bc2_time+time_e6-time_s6
			bc3_time=bc3_time+time_e7-time_s7
			colission_time=colission_time+time_e4-time_s4
			move_airfoil_time=move_airfoil_time+time_e8-time_s8
			outer_bc_time=outer_bc_time+time_e9-time_s9
        ! The four refinement boxes remain fixed in physical coordinates.
	        if(flag2.gt.0) flag2=flag2-1
			call cpu_time(time_s10)
	!        call vor_mom_rate(n,mesh)
			call cpu_time(time_e10)
			vor_time=vor_time+time_e10-time_s10
			if(myid.eq.0) then
				if (mod(n,100).eq.0) print*,"已迭代",n,"步，误差为：",err
	        endif

        ! Offline fixed-time snapshots (completed flow time = n*solver_dt).
        flag1=0
        do save_second=30,35
            nsave=nint(dble(save_second)/solver_dt)
            if(n.eq.nsave) flag1=1
        enddo
        ! End offline fixed-time snapshots.



	        if((flag1.eq.1).and.(flag2.eq.0))then
			!if(n.eq.1000)then
	            call output(n,mesh)
	            flag1=0
	        endif
		!if ((n.eq.916400).or.(n.eq.924400).or.(n.eq.917200).or.(n.eq.926400).or.(n.eq.916800).or.(n.eq.925200)) call output(n,mesh)
				if(n.eq.motion_step_count) then ! Save restart at the final offline motion step.
					call saveflowfield(n,mesh)
				endif
			n=n+1
		enddo
		!!$acc end data
		call cpu_time(cpu_es1)
		if(myid.eq.0)print*,"Run time =",cpu_es1-cpu_ts
		if(myid.eq.0)print*,"ansstime =",ansstime
		if(myid.eq.0)print*,"streaming_time =",streaming_time
		if(myid.eq.0)print*,"ph_time =",ph_time
		if(myid.eq.0)print*,"bc_time =",bc_time
		if(myid.eq.0)print*,"bc_time1 =",bc1_time
		if(myid.eq.0)print*,"bc_time2 =",bc2_time
		if(myid.eq.0)print*,"bc_time3 =",bc3_time
		if(myid.eq.0)print*,"colission_time =",colission_time
		if(myid.eq.0)print*,"move_airfoil_time =",move_airfoil_time
		if(myid.eq.0)print*,"move_mesh_time =",move_mesh_time
		if(myid.eq.0)print*,"outer_bc_time =",outer_bc_time
		if(myid.eq.0)print*,"vor_time =",vor_time
	else
	    call run_zmq_driver(mesh,h)
	endif
	stat = cusolverDnDestroy(h)
!	call output(n)
	call MPI_FINALIZE(rc)	
end Program
