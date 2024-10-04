!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use geometry,          only: cfg
   use hypre_str_class,   only: hypre_str
   use ddadi_class,       only: ddadi
   use tpns_class,        only: tpns
   use vfs_class,         only: vfs
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use surfmesh_class,    only: surfmesh
   use event_class,       only: event
   use monitor_class,     only: monitor
   implicit none
   private
   
   !> Single two-phase flow solver and volume fraction solver and corresponding time tracker
   type(hypre_str),   public :: ps
   type(ddadi),       public :: vs
   type(tpns),        public :: fs
   type(vfs),         public :: vf
   type(timetracker), public :: time
   
   !> Ensight postprocessing
   type(surfmesh) :: smesh
   type(ensight)  :: ens_out
   type(event)    :: ens_evt
   
   !> Simulation monitor file
   type(monitor) :: mfile,cflfile
   
   public :: simulation_init,simulation_run,simulation_final
   
   !> Private work arrays
   real(WP), dimension(:,:,:), allocatable :: resU,resV,resW
   real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi
   
   !> Problem definition and post-processing
   real(WP), dimension(3) :: Cdrop
   real(WP) :: Rdrop, Beta_NS
   integer :: CLsolver
   real(WP) :: height,R_wet,CArad,CAdeg,CLvel,C,alpha
   reaL(WP), dimension(:), allocatable :: all_time,all_rwet
   type(monitor) :: ppfile
   
contains
   
   
   !> Function that defines a level set function for a spreading drop problem
   function levelset_contact_drop(xyz,t) result(G)
      implicit none
      real(WP), dimension(3),intent(in) :: xyz
      real(WP), intent(in) :: t
      real(WP) :: G
      ! Create the droplet
      G=Rdrop-sqrt(sum((xyz-Cdrop)**2))
   end function levelset_contact_drop
   
   
   !> Specialized subroutine that outputs wetting information
   subroutine postproc_data()
      use irl_fortran_interface
      use mathtools, only: Pi
      use string,    only: str_medium
      use mpi_f08,   only: MPI_ALLREDUCE,MPI_SUM
      use parallel,  only: MPI_REAL_WP
      implicit none
      integer :: ierr,i,j,k,my_size
      real(WP) :: my_height,myR1,R1,myR2,R2,R_wet_old
      real(WP), dimension(:), allocatable :: temp
      ! Post-process height of drop
      my_height=0.0_WP
      do k=vf%cfg%kmin_,vf%cfg%kmax_
         do i=vf%cfg%imin_,vf%cfg%imax_
            ! Find closest vertical column to center
            if (vf%cfg%x(i).le.0.0_WP.and.vf%cfg%x(i+1).gt.0.0_WP.and.vf%cfg%z(k).le.0.0_WP.and.vf%cfg%z(k+1).gt.0.0_WP) then
               ! Integrate height
               do j=vf%cfg%jmin_,vf%cfg%jmax_
                  my_height=my_height+vf%VF(i,j,k)*vf%cfg%dy(j)
               end do
            end if
         end do
      end do
      call MPI_ALLREDUCE(my_height,height,1,MPI_REAL_WP,MPI_SUM,vf%cfg%comm,ierr)
      ! Post-process wetted radius, contact line velocity, and effective contact angle from bottom cell PLIC
      myR1=0.0_WP; myR2=0.0_WP
      if (vf%cfg%jproc.eq.1) then
         do k=vf%cfg%kmin_,vf%cfg%kmax_
            do i=vf%cfg%imin_,vf%cfg%imax_
               ! Compute the wetted area from VOF
               myR1=myR1+vf%VF(i,vf%cfg%jmin,k)*vf%cfg%dx(i)*vf%cfg%dz(k)
               ! Compute the wetted area from top of PLIC
               call getMoments(vf%polyface(2,i,vf%cfg%jmin+1,k),vf%liquid_gas_interface(i,vf%cfg%jmin,k),R2); myR2=myR2+abs(R2)
            end do
         end do
      end if
      call MPI_ALLREDUCE(myR1,R1,1,MPI_REAL_WP,MPI_SUM,vf%cfg%comm,ierr); R1=sqrt(R1/Pi)
      call MPI_ALLREDUCE(myR2,R2,1,MPI_REAL_WP,MPI_SUM,vf%cfg%comm,ierr); R2=sqrt(R2/Pi)
      R_wet_old=R_wet
      R_wet=2.0_WP*R1-R2
      if (time%t.eq.0.0_WP) then
         CLvel=0.0_WP
      else
         CLvel=(R_wet-R_wet_old)/time%dt
      end if
      CArad=atan2((R2-R1),0.5_WP*vf%cfg%dy(vf%cfg%jmin))+0.5_WP*Pi
      CAdeg=CArad*180.0_WP/Pi
      ! Also attempt to get C and alpha on the fly
      alpha=time%t*CLvel/R_wet
      C=R_wet/(time%t**alpha)
      ! Store time and amplitude series
      if (.not.allocated(all_time)) then
         my_size=0
      else
         my_size=size(all_time,dim=1)
      end if
      allocate(temp(my_size+1)); temp(1:my_size)=all_time; temp(my_size+1)=time%t; call MOVE_ALLOC(temp,all_time)
      allocate(temp(my_size+1)); temp(1:my_size)=all_rwet; temp(my_size+1)=R_wet ; call MOVE_ALLOC(temp,all_rwet)
   end subroutine postproc_data
   
   
   !> Function that localizes the top (x+) of the domain
   function xp_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.eq.pg%imax+1.and.pg%ym(j).ge.0.0_WP) isIn=.true.
   end function xp_locator
   
   
   !> Function that localizes the bottom (x-) of the domain
   function xm_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.eq.pg%imin.and.pg%ym(j).ge.0.0_WP) isIn=.true.
   end function xm_locator
   
   
   !> Function that localizes the top (z+) of the domain
   function zp_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (k.eq.pg%kmax+1.and.pg%ym(j).ge.0.0_WP) isIn=.true.
   end function zp_locator
   
   
   !> Function that localizes the bottom (z-) of the domain
   function zm_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (k.eq.pg%kmin.and.pg%ym(j).ge.0.0_WP) isIn=.true.
   end function zm_locator
   
   
   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none
      

      ! Allocate work arrays
      allocate_work_arrays: block
         allocate(resU(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resV(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resW(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Ui  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Vi  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Wi  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      end block allocate_work_arrays
      
      
      ! Initialize time tracker with 2 subiterations
      initialize_timetracker: block
         time=timetracker(amRoot=cfg%amRoot)
         call param_read('Max timestep size',time%dtmax)
         call param_read('Max cfl number',time%cflmax)
         call param_read('Max time',time%tmax)
         time%dt=time%dtmax
         time%itmax=2
      end block initialize_timetracker
      
      
      ! Initialize our VOF solver and field
      create_and_initialize_vof: block
         use mms_geom,  only: cube_refine_vol
         use vfs_class, only: lvira,VFhi,VFlo,remap
         use mathtools, only: Pi
         use, intrinsic :: iso_fortran_env, only: output_unit
         use string,    only: str_long
         use messager,  only: log
         integer :: i,j,k,n,si,sj,sk
         real(WP), dimension(3,8) :: cube_vertex
         real(WP), dimension(3) :: v_cent,a_cent
         real(WP) :: vol,area,contact
         integer, parameter :: amr_ref_lvl=4
         character(len=str_long) :: message
         ! Create a VOF solver
         call vf%initialize(cfg=cfg,reconstruction_method=lvira,transport_method=remap,name='VOF')
         ! Prepare the analytical calculation of a sphere on a wall
         call param_read('Initial drop radius',Rdrop,default=1.0_WP)
         call param_read('Initial contact angle',contact,default=180.0_WP); contact=contact*Pi/180.0_WP
         call param_read('Beta',Beta_NS,default=1.0_WP)
         call param_read('CLsolver',CLsolver,default=0)
         if (vf%cfg%nz.eq.1) then ! 2D analytical drop shape
            Rdrop=Rdrop*sqrt(Pi/(2.0_WP*(contact-sin(contact)*cos(contact))))
         else ! 3D analytical drop shape
            Rdrop=Rdrop*(4.0_WP/(2.0_WP-3.0_WP*cos(contact)+(cos(contact))**3))**(1.0_WP/3.0_WP)
         end if
         Cdrop=[0.0_WP,-Rdrop*cos(contact),0.0_WP]
         if (vf%cfg%amRoot) then
            write(message,'("Droplet initial radius is ",es12.5)') Rdrop; call log(message)
         end if
         ! Initialize the VOF field
         do k=vf%cfg%kmino_,vf%cfg%kmaxo_
            do j=vf%cfg%jmino_,vf%cfg%jmaxo_
               do i=vf%cfg%imino_,vf%cfg%imaxo_
                  ! Handle wall cells or cells below the plate surface
                  if (vf%mask(i,j,k).eq.1.or.vf%cfg%ym(j).lt.0.0_WP) then
                     vf%VF(i,j,k)=0.0_WP
                     vf%Lbary(:,i,j,k)=[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]
                     vf%Gbary(:,i,j,k)=[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]
                     cycle
                  end if
                  ! Set cube vertices
                  n=0
                  do sk=0,1
                     do sj=0,1
                        do si=0,1
                           n=n+1; cube_vertex(:,n)=[vf%cfg%x(i+si),vf%cfg%y(j+sj),vf%cfg%z(k+sk)]
                        end do
                     end do
                  end do
                  ! Call adaptive refinement code to get volume and barycenters recursively
                  vol=0.0_WP; area=0.0_WP; v_cent=0.0_WP; a_cent=0.0_WP
                  call cube_refine_vol(cube_vertex,vol,area,v_cent,a_cent,levelset_contact_drop,0.0_WP,amr_ref_lvl)
                  vf%VF(i,j,k)=vol/vf%cfg%vol(i,j,k)
                  if (vf%VF(i,j,k).ge.VFlo.and.vf%VF(i,j,k).le.VFhi) then
                     vf%Lbary(:,i,j,k)=v_cent
                     vf%Gbary(:,i,j,k)=([vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]-vf%VF(i,j,k)*vf%Lbary(:,i,j,k))/(1.0_WP-vf%VF(i,j,k))
                  else
                     vf%Lbary(:,i,j,k)=[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]
                     vf%Gbary(:,i,j,k)=[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]
                  end if
               end do
            end do
         end do
         ! Update the band
         call vf%update_band()
         ! Perform interface reconstruction from VOF field
         call vf%build_interface()
         ! Set interface planes at the boundaries
         call vf%set_full_bcond()
         ! Create discontinuous polygon mesh from IRL interface
         call vf%polygonalize_interface()
         ! Calculate distance from polygons
         call vf%distance_from_polygon()
         ! Calculate subcell phasic volumes
         call vf%subcell_vol()
         ! Calculate curvature
         call vf%get_curvature()
         ! Reset moments to guarantee compatibility with interface reconstruction
         call vf%reset_volume_moments()
      end block create_and_initialize_vof
      
      
      ! Create a two-phase flow solver without bconds
      create_and_initialize_flow_solver: block
         use tpns_class,      only: clipped_neumann,dirichlet
         use hypre_str_class, only: pcg_pfmg2
         use mathtools,       only: Pi
         integer :: i,j,k,n
         ! Create flow solver
         fs=tpns(cfg=cfg,name='Two-phase NS')
         ! Assign constant density to each phase
         fs%rho_l=1.0_WP; call param_read('Density ratio',fs%rho_g); fs%rho_g=fs%rho_l/fs%rho_g
         ! Read in surface tension coefficient and contact angle
         fs%sigma=1.0_WP; call param_read('CA',fs%contact_angle); fs%contact_angle=fs%contact_angle*Pi/180.0_WP
         ! Assign constant viscosity to each phase
         call param_read('Oh',fs%visc_l); call param_read('Viscosity ratio',fs%visc_g); fs%visc_g=fs%visc_l/fs%visc_g
         ! Setup boundary conditions
         call fs%add_bcond(name='bc_xp',type=clipped_neumann,face='x',dir=+1,canCorrect=.true.,locator=xp_locator)
         call fs%add_bcond(name='bc_xm',type=clipped_neumann,face='x',dir=-1,canCorrect=.true.,locator=xm_locator)
         call fs%add_bcond(name='bc_zp',type=clipped_neumann,face='z',dir=+1,canCorrect=.true.,locator=zp_locator)
         call fs%add_bcond(name='bc_zm',type=clipped_neumann,face='z',dir=-1,canCorrect=.true.,locator=zm_locator)
         ! Assign gravity from input file (will need to make non-dimensional form)
         fs%gravity=1.0_WP
         call param_read('Gravity',fs%gravity)
         ! Configure pressure solver
         ps=hypre_str(cfg=cfg,name='Pressure',method=pcg_pfmg2,nst=7)
         ps%maxlevel=10
         call param_read('Pressure iteration',ps%maxit)
         call param_read('Pressure tolerance',ps%rcvg)
         ! Configure implicit velocity solver
         vs=ddadi(cfg=cfg,name='Velocity',nst=7)
         ! Setup the solver
         call fs%setup(pressure_solver=ps,implicit_solver=vs)
         ! Zero initial field
         fs%U=0.0_WP; fs%V=0.0_WP; fs%W=0.0_WP
         ! Calculate cell-centered velocities and divergence
         call fs%interp_vel(Ui,Vi,Wi)
         call fs%get_div()
      end block create_and_initialize_flow_solver
      
      
      ! Create surfmesh object for interface polygon output
      create_smesh: block
         use irl_fortran_interface
         integer :: i,j,k,np,nplane
         smesh=surfmesh(nvar=1,name='plic')
         smesh%varname(1)='curv'
         call vf%update_surfmesh(smesh)
         np=0
         do k=vf%cfg%kmin_,vf%cfg%kmax_
            do j=vf%cfg%jmin_,vf%cfg%jmax_
               do i=vf%cfg%imin_,vf%cfg%imax_
                  do nplane=1,getNumberOfPlanes(vf%liquid_gas_interface(i,j,k))
                     if (getNumberOfVertices(vf%interface_polygon(nplane,i,j,k)).gt.0) then
                        np=np+1
                        smesh%var(1,np)=vf%curv(i,j,k)
                     end if
                  end do
               end do
            end do
         end do
      end block create_smesh
      
      
      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='twophase')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_scalar('pressure',fs%P)
         call ens_out%add_vector('velocity',Ui,Vi,Wi)
         call ens_out%add_scalar('VOF',vf%VF)
         call ens_out%add_scalar('curvature',vf%curv)
         call ens_out%add_surface('vofplic',smesh)
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
      end block create_ensight
      
      
      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         call fs%get_cfl(time%dt,time%cfl)
         call fs%get_max()
         call vf%get_max()
         ! Create simulation monitor
         mfile=monitor(fs%cfg%amRoot,'simulation')
         call mfile%add_column(time%n,'Timestep number')
         call mfile%add_column(time%t,'Time')
         call mfile%add_column(time%dt,'Timestep size')
         call mfile%add_column(time%cfl,'Maximum CFL')
         call mfile%add_column(fs%Umax,'Umax')
         call mfile%add_column(fs%Vmax,'Vmax')
         call mfile%add_column(fs%Wmax,'Wmax')
         call mfile%add_column(fs%Pmax,'Pmax')
         call mfile%add_column(vf%VFmax,'VOF maximum')
         call mfile%add_column(vf%VFmin,'VOF minimum')
         call mfile%add_column(vf%VFint,'VOF integral')
         call mfile%add_column(fs%divmax,'Maximum divergence')
         call mfile%add_column(fs%psolv%it,'Pressure iteration')
         call mfile%add_column(fs%psolv%rerr,'Pressure error')
         call mfile%write()
         ! Create CFL monitor
         cflfile=monitor(fs%cfg%amRoot,'cfl')
         call cflfile%add_column(time%n,'Timestep number')
         call cflfile%add_column(time%t,'Time')
         call cflfile%add_column(fs%CFLst,'STension CFL')
         call cflfile%add_column(fs%CFLc_x,'Convective xCFL')
         call cflfile%add_column(fs%CFLc_y,'Convective yCFL')
         call cflfile%add_column(fs%CFLc_z,'Convective zCFL')
         call cflfile%add_column(fs%CFLv_x,'Viscous xCFL')
         call cflfile%add_column(fs%CFLv_y,'Viscous yCFL')
         call cflfile%add_column(fs%CFLv_z,'Viscous zCFL')
         call cflfile%write()
      end block create_monitor
      
      
      ! Create a specialized post-processing file
      create_postproc: block
         ! Post-proc wetting information
         call postproc_data()
         ! Create monitor for liquid data
         ppfile=monitor(fs%cfg%amRoot,'dropinfo')
         call ppfile%add_column(time%n,'Timestep number')
         call ppfile%add_column(time%t,'Time')
         call ppfile%add_column(vf%VFmax,'VOF maximum')
         call ppfile%add_column(vf%VFmin,'VOF minimum')
         call ppfile%add_column(vf%VFint,'Total volume')
         call ppfile%add_column(height,'Drop height')
         call ppfile%add_column(R_wet,'Wetted radius')
         call ppfile%add_column(CLvel,'CL vel')
         call ppfile%add_column(CArad,'CA rad')
         call ppfile%add_column(CAdeg,'CA deg')
         call ppfile%add_column(C,'C')
         call ppfile%add_column(alpha,'alpha')
         call ppfile%write()
      end block create_postproc
      
      
      ! TEST THE REMOVAL OF WALL TREATMENT IN VOF TRANSPORT
      vf%vmask=0
      ! TEST THE REMOVAL OF CONSERVATIVE CORRECTION IN SL VOF TRANSPORT
      vf%cons_correct=.false.
      
      
   end subroutine simulation_init
   
   
   !> Perform an NGA2 simulation - this mimicks NGA's old time integration for multiphase
   subroutine simulation_run
      use tpns_class, only: static_contact,arithmetic_visc
      implicit none
      
      ! Perform time integration
      do while (.not.time%done())
         
         ! Increment time
         call fs%get_cfl(time%dt,time%cfl)
         call time%adjust_dt()
         call time%increment()
         
         ! Implement slip condition for contact line modeling here
         if (CLsolver.ne.0) then!enable slip model for solver numbers >0
            call contact_slip()
         end if
         ! Remember old VOF
         vf%VFold=vf%VF
         
         ! Remember old velocity
         fs%Uold=fs%U
         fs%Vold=fs%V
         fs%Wold=fs%W
         
         ! Prepare old staggered density (at n)
         call fs%get_olddensity(vf=vf)
         
         ! VOF solver step
         call vf%advance(dt=time%dt,U=fs%U,V=fs%V,W=fs%W)
         
         ! Prepare new staggered viscosity (at n+1)
         call fs%get_viscosity(vf=vf,strat=arithmetic_visc)
         
         ! Perform sub-iterations
         do while (time%it.le.time%itmax)
            
            ! Build mid-time velocity
            fs%U=0.5_WP*(fs%U+fs%Uold)
            fs%V=0.5_WP*(fs%V+fs%Vold)
            fs%W=0.5_WP*(fs%W+fs%Wold)
            
            ! Preliminary mass and momentum transport step at the interface
            call fs%prepare_advection_upwind(dt=time%dt)
            
            ! Explicit calculation of drho*u/dt from NS
            call fs%get_dmomdt(resU,resV,resW)
            ! Momentum source terms
            call fs%addsrc_gravity(resU,resV,resW)
            
            ! Assemble explicit residual
            resU=-2.0_WP*fs%rho_U*fs%U+(fs%rho_Uold+fs%rho_U)*fs%Uold+time%dt*resU
            resV=-2.0_WP*fs%rho_V*fs%V+(fs%rho_Vold+fs%rho_V)*fs%Vold+time%dt*resV
            resW=-2.0_WP*fs%rho_W*fs%W+(fs%rho_Wold+fs%rho_W)*fs%Wold+time%dt*resW
            
            ! Form implicit residuals
            call fs%solve_implicit(time%dt,resU,resV,resW)
            
            ! Apply these residuals
            fs%U=2.0_WP*fs%U-fs%Uold+resU
            fs%V=2.0_WP*fs%V-fs%Vold+resV
            fs%W=2.0_WP*fs%W-fs%Wold+resW
            
            ! Apply other boundary conditions
            call fs%apply_bcond(time%t,time%dt)
            
            ! Solve Poisson equation
            call fs%update_laplacian()
            call fs%correct_mfr()
            call fs%get_div()
            if (CLsolver.eq.0) then
               call fs%add_surface_tension_jump(dt=time%dt,div=fs%div,vf=vf,contact_model=static_contact)
            else if (CLsolver.gt.0) then
               call fs%add_surface_tension_jump(dt=time%dt,div=fs%div,vf=vf)
            end if
            fs%psolv%rhs=-fs%cfg%vol*fs%div/time%dt
            fs%psolv%sol=0.0_WP
            call fs%psolv%solve()
            call fs%shift_p(fs%psolv%sol)
            
            ! Correct velocity
            call fs%get_pgrad(fs%psolv%sol,resU,resV,resW)
            fs%P=fs%P+fs%psolv%sol
            fs%U=fs%U-time%dt*resU/fs%rho_U
            fs%V=fs%V-time%dt*resV/fs%rho_V
            fs%W=fs%W-time%dt*resW/fs%rho_W
            
            ! Increment sub-iteration counter
            time%it=time%it+1
            
         end do
         
         ! Recompute interpolated velocity and divergence
         call fs%interp_vel(Ui,Vi,Wi)
         call fs%get_div()
         
         ! Output to ensight
         if (ens_evt%occurs()) then
            update_smesh: block
               use irl_fortran_interface
               integer :: i,j,k,np,nplane
               call vf%update_surfmesh(smesh)
               np=0
               do k=vf%cfg%kmin_,vf%cfg%kmax_
                  do j=vf%cfg%jmin_,vf%cfg%jmax_
                     do i=vf%cfg%imin_,vf%cfg%imax_
                        do nplane=1,getNumberOfPlanes(vf%liquid_gas_interface(i,j,k))
                           if (getNumberOfVertices(vf%interface_polygon(nplane,i,j,k)).gt.0) then
                              np=np+1
                              smesh%var(1,np)=vf%curv(i,j,k)
                           end if
                        end do
                     end do
                  end do
               end do
            end block update_smesh
            call ens_out%write_data(time%t)
         end if
         
         ! Perform and output monitoring
         call fs%get_max()
         call vf%get_max()
         call mfile%write()
         call cflfile%write()
         
         ! Specialized post-processing
         call postproc_data(); call ppfile%write()
         
      end do
   end subroutine simulation_run
   
   
   !> Subroutine that updates the slip velocity based on contact line model
   subroutine contact_slip()
      use irl_fortran_interface
      implicit none
      integer :: i,j,k
      real(WP), dimension(3) :: nw
      real(WP) :: mysurf,mycos,cos_contact_angle
      real(WP) :: log_res_slip
      ! Precalculate cos(contact angle)
      cos_contact_angle=cos(fs%contact_angle)
      ! Force use of new Beta Factor
      log_res_slip=log(1e3*fs%cfg%dx(1)*fs%visc_l**2)  !fs%cfg%dx(1,1,1)*1e9   ... is an array
      Beta_NS=fs%contact_angle**2/(sin(fs%contact_angle)*3*log_res_slip*fs%visc_l)
      if (CLsolver.eq.1) then
         Beta_NS=1.0
      end if

      ! Loop over domain and identify cells that require contact angle model
      do k=fs%cfg%kmin_,fs%cfg%kmax_+1
         do j=fs%cfg%jmin_,fs%cfg%jmax_+1
            do i=fs%cfg%imin_,fs%cfg%imax_+1
               
               ! Check if there is a wall in y-
               if (fs%mask(i,j-1,k).eq.1.and.fs%mask(i-1,j-1,k).eq.1) then
                  ! Define wall normal
                  nw=[0.0_WP,+1.0_WP,0.0_WP]
                  ! Handle U-slip
                  mysurf=abs(calculateVolume(vf%interface_polygon(1,i-1,j,k)))+abs(calculateVolume(vf%interface_polygon(1,i,j,k)))
                  if (mysurf.gt.0.0_WP.and.fs%umask(i,j,k).eq.0) then
                     ! Surface-averaged local cos(CA)
                     mycos=(abs(calculateVolume(vf%interface_polygon(1,i-1,j,k)))*dot_product(calculateNormal(vf%interface_polygon(1,i-1,j,k)),nw)+&
                     &      abs(calculateVolume(vf%interface_polygon(1,i  ,j,k)))*dot_product(calculateNormal(vf%interface_polygon(1,i  ,j,k)),nw))/mysurf
                     ! Apply slip velocity
                     fs%U(i,j-1,k)=Beta_NS*fs%sigma*(mycos-cos_contact_angle)*sum(fs%divu_x(:,i,j,k)*vf%VF(i-1:i,j,k)*fs%cfg%dx(i))
                  end if
                  ! Handle W-slip
                  mysurf=abs(calculateVolume(vf%interface_polygon(1,i,j,k-1)))+abs(calculateVolume(vf%interface_polygon(1,i,j,k)))
                  if (mysurf.gt.0.0_WP.and.fs%wmask(i,j,k).eq.0) then
                     ! Surface-averaged local cos(CA)
                     mycos=(abs(calculateVolume(vf%interface_polygon(1,i,j,k-1)))*dot_product(calculateNormal(vf%interface_polygon(1,i,j,k-1)),nw)+&
                     &      abs(calculateVolume(vf%interface_polygon(1,i,j,k  )))*dot_product(calculateNormal(vf%interface_polygon(1,i,j,k  )),nw))/mysurf
                     ! Apply slip velocity
                     fs%W(i,j-1,k)=Beta_NS*fs%sigma*(mycos-cos_contact_angle)*sum(fs%divw_z(:,i,j,k)*vf%VF(i,j,k-1:k)*fs%cfg%dz(k))
                  end if
               end if
               
            end do
         end do
      end do
      
   end subroutine contact_slip
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      
      ! Get rid of all objects - need destructors
      ! monitor
      ! ensight
      ! bcond
      ! timetracker
      
      ! Deallocate work arrays
      deallocate(resU,resV,resW,Ui,Vi,Wi)
      
   end subroutine simulation_final
   
   
end module simulation
