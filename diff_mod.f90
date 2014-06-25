!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                                                                         !!
!!  eVsOneD_ICF.f90                                                                                                        !!   
!!  Description: Test of 1D ICF problem 					 															   !!
!!               spherical geometry																						   !!
!!               Lagrange hydro																				  			   !!
!!               option for laser power at outer boundary 																   !!
!!               plasma diffusivity of heat, momentum and mass 											 				   !!
!!  Revisions: Erik Vold																			   					   !!
!!  created: Sept. 2013																				   					   !!
!!  last modified: Sept 22, 2013																							   !!
!!  this version built on framework, select coding conventions, and namelist inputs created by Ryan Moll                   !!
!!      for the LANL Comp.Phys. Summer Workshop: July 22, 2013                                                             !!
!!                                                                                                                         !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                                                                         !!
!!  Modules:                                                                                                               !!
!!    prec_param    --> variable for setting precision                                                                     !!
!!    exp_params    --> contains variables for all experment parameters that do not vary across zones                      !!
!!    phys_const    --> defines values for necessary physical constants                                                    !!
!!    data_arrays   --> defines all values that vary across zones. Values include: radius, dr, velocity, density,          !!
!!                      internal energy, pressure, artificial viscosity, volume, sound speed, and composition as           !!
!!                      mass density of species 1 in binary diffusion:                                                     !!
!!                      arrays integrate din time are indexed (ir,1) for old time value and (ir,2) for new time value      !!
!!  Subroutines:                                                                                                           !!
!!    read_params      --> reads experiment parameters from file sphere_param                                              !!
!!    initialize       --> initializes experiment parameters, allocates memory for data arrays, opens data file for write  !!
!!    write_data_matlab--> writes to file values of r, u, rho, en, and P for current timestep for matlab input             !!
!!    write_data_pyout1--> writes to file values of r, u, rho, en, and P for current timestep for python input v.1         !!
!!    solve_eqns       --> does numerical integration of experiment variables                                              !!
!!    reset_arrays     --> sets arrays of old values to be new values calculated in timestep                               !!
!!    visc_mat_coeffs  --> calculates matrix coefficients for implicit solve of momentum equation. Plasma viscosity is     !!
!!                         calculated here.                                                                                !! 
!!    temp_mat_coeffs  --> calculates matrix coefficients for implicit solve of temperature equation. Plasma thermal       !!
!!                         conduction is calculated here.                                                                  !!
!!    comp_mat_coeffs  --> calculates matrix coefficients for implicit solve of mass diffusion equation. Diffusion         !!
!!                         coefficient is calculated here.                                                                 !!
!!    tridag           --> performs inversion of tridiagonal matrix to complete implicit solve of momentum, temperature    !!
!!                         and mass diffusion equations. Subroutine is called seperately for each equation                 !!     
!!    finish           --> closes data files and deallocates memory                                                        !!
!!                                                                                                                         !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                                                                         !!
!!   SI Units: r --> meters, Vol --> m^3, rho --> kg/m^3, u --> m/s, kT --> Joules, P --> kg/(m*s^2)                       !!
!!                                                                                                                         !!   
!!  Note1: In the plasma routines viscosity, thermal conduction and mass diffusion coefficients are calculated in cgs then !!
!!        multiplied by conversion factors to convert them to SI.     MODIFY AS ????                                       !!
!!  Note2: T is in eV so kT = 1.602e-19*T [ J/# ] and particle densities are in [ moles / cc ] so p = nkT =>               !!
!!         p[J/m3] = n[moles/cc]* 10^6*Avog*bk*T[eV], 																	   !!
!!            w/ pco = 10^6*Avog*bk =9.648e10  then: p [J/m3]= pco* n[moles/cc]* T[eV]									   !!
!!  Note3: amu [g/mole] so rho[kg/m3] = 10^(3)*amu[g/mole]*n[moles/cc] or n = 10^(-3)*rho/amu							   !!
!!                                                                                                                         !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


MODULE prec_param

  INTEGER, PARAMETER  ::  prec = 8

END MODULE prec_param

MODULE exp_params

  use prec_param
  
  INTEGER             ::  i, it, itprt, it_special_prtOut
  INTEGER             ::  nr, nmax_tsteps
  REAL(prec)          ::  time, dt, tmax, timeold, t_out_prt, Lr  

  LOGICAL             ::  nanDetected
  REAL(prec)          ::  gamma, cl, cQ, du, inv_rho_bar, dr_bar, ion_frac
  INTEGER             ::  forced_zones, in_bound, in_forced, visc_sw, diff_sw, hf_sw
  REAL(prec)          ::  init_T_gas, init_press, rho_fac, vol_fac, dtheta, dphi, pi, r_bar, r_bar_sq
  REAL(prec)          ::  out_press, out_density, out_temp, cfl_fac, init_dt, init_watts_pulse
  REAL(prec)          ::  highestPress, highestTemp, highestDensity, highEnergyTime, time_of_highestTemp
  REAL(prec)          ::  m1, m2, z1, z2
  CHARACTER(len=100)  ::  outfile, infile
  CHARACTER(len=7)    ::  dir
  CHARACTER(len=28)   ::  path = '/Users/archis/Documents/lanl/'

END MODULE exp_params

MODULE phys_const

  use prec_param
    
  ! defining physical parameters: k_boltz  = Boltzmann constant in [eV/Kelvin]
  !                               N_av     = Avogadro's number [particles/mol]
  !                               mm_DD    = molar mass of deuterium in [kg/mol]
  !                               mm_CH    = molar mass of the shell material [kg/mol]
  !                               JeV_conv = conversion factor between joules and eV [J/eV]
  REAL(prec),PARAMETER  :: k_boltz = 8.6173324e-5, N_av = 6.02214129e23,     &
                         &    JeV_conv = 1.602e-19, a_sbc = 7.5657e-16,      &
                         &    mm_DD = 4.029e-3,  mm_CH = 13.018947e-3,       &
                         &    eps = 1.e-12

END MODULE phys_const


MODULE data_arrays

  use prec_param

  REAL(prec), ALLOCATABLE  ::  r(:,:), u(:,:), rho(:,:), Temp(:,:), P(:,:), Comp(:,:)
  REAL(prec), ALLOCATABLE  ::  dr(:,:), Vol(:,:), q(:), Cs(:), N_part(:), M(:), vel_fac(:)
  REAL(prec), ALLOCATABLE  ::  ionized(:), N_elec(:), subdiag(:), diag(:), supdiag(:)
  REAL(prec), ALLOCATABLE  ::  vel_RHS(:), TempIC(:), temp_RHS(:), forcing(:), r_sq(:)

END MODULE data_arrays


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


PROGRAM eVsOneD_ICF

  use exp_params
  use phys_const
  use data_arrays


  implicit none

  ! read parameter file
  call read_params

  ! initialize variables
  call initialize

  ! write initial data to file
  call write_data_matlab
  call write_data_pyout1


  ! enter time stepping loop
!..................................................................................
  do it = 1, nmax_tsteps
    
    ! solve equations
    call solve_eqns

    ! reset indicies of data arrays
    call reset_arrays

    ! set max core values for diagnostics
    if (P(1,1) .GT. highestPress) then
      highestPress = P(1,1)
    end if

    if (Temp(1,1) .GT. highestTemp) then
      highestTemp = Temp(1,1)
      time_of_highestTemp = time
    end if

    if (rho(1,1) .GT. highestDensity) then
      highestDensity = rho(1,1)
    end if
    

    ! write to file (uncomment if statement to start saving data after certain time)
    !if (time .GT. 3.227e-9) then
    if (mod(it,itprt) .EQ. 0) then
      call write_data_matlab
      call write_data_pyout1
      print *, 'it, time, dt =',it, time, dt
    end if
    !end if
    ! special printout as re-compiled
    !if ( ( it_special_prtOut .eq. 0 ) .and. ( time .gt. t_out_prt )  ) then
    !    do i = 2, nr
    !        print *, r(i,2), comp(i,2)
    !    enddo
    !        it_special_prtOut = 1
    !endif

    ! Terminate program after certain amount of time (uncomment to activate)
    if (time .GT. tmax)  exit

  end do     !    ends it = 1, nnmax_tsteps
!..................................................................................


  ! complete program (close read/write files, deallocate arrays)
  call finish

END PROGRAM eVsOneD_ICF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE read_params

  use prec_param  
  use exp_params

  INTEGER             ::  number_r_cells, nmax_time_steps, it_steps_per_print, forced_region_zones
  INTEGER             ::  viscosity_sw, heatflux_sw, diffusion_sw
  REAL(prec)          ::  domain_size, cfl_factor, initial_time_step, adiabatic_index, time_max, t_out_special_print
  REAL(prec)          ::  inner_boundary, initial_pressure, initial_T_gas, density_change
  REAL(prec)          ::  initial_watts_pulse, laser_absorbed_frac, laser_watts_pulse
  REAL(prec)          ::  ionization_fraction, outerBC_density, outerBC_temp
  REAL(prec)          ::  mass1, mass2, charge1, charge2
  CHARACTER(len=100)  ::  output_file_name

  ! create namelist of items in parameter file 
  
  NAMELIST /input_values/ number_r_cells
  NAMELIST /input_values/ domain_size
  NAMELIST /input_values/ nmax_time_steps
  NAMELIST /input_values/ time_max
  NAMELIST /input_values/ cfl_factor
  NAMELIST /input_values/ initial_time_step
  NAMELIST /input_values/ adiabatic_index
  NAMELIST /input_values/ inner_boundary         !   as a fraction of domain size, Lr
  NAMELIST /input_values/ forced_region_zones    !   as an integer number of zones
  NAMELIST /input_values/ mass1
  NAMELIST /input_values/ charge1
  NAMELIST /input_values/ mass2
  NAMELIST /input_values/ charge2
  NAMELIST /input_values/ initial_pressure
  NAMELIST /input_values/ initial_T_gas
  NAMELIST /input_values/ density_change
  NAMELIST /input_values/ initial_watts_pulse
  NAMELIST /input_values/ laser_absorbed_frac
  NAMELIST /input_values/ ionization_fraction
  NAMELIST /input_values/ outerBC_density
  NAMELIST /input_values/ outerBC_temp
  NAMELIST /input_values/ viscosity_sw
  NAMELIST /input_values/ heatflux_sw
  NAMELIST /input_values/ diffusion_sw


  NAMELIST /input_values/ it_steps_per_print
  
  NAMELIST /input_values/ output_file_name  
  NAMELIST /input_values/ t_out_special_print


  ! open parameter file for reading and read data into namelist 
  infile = 'input'
  open(30, file=infile, action='read')
  read(30, NML=input_values)

  ! set experiment parameters to values read from parameter file

  nr                = number_r_cells                              ! number of grid zones in r (or coordinate 1)
  Lr                = domain_size                                 ! physical size of domain in meters
  nmax_tsteps       = nmax_time_steps                             ! number of time steps to use
  tmax              = time_max                                    ! maximum time for simulation run
  cfl_fac           = cfl_factor                                  ! set max timestep as fraction of CFL limit 
  init_dt           = initial_time_step                           ! initial time step size
  gamma             = adiabatic_index                             ! self explanatory
  in_bound          = floor(nr*inner_boundary)                    ! location of shell/fuel boundary as a fraction of total radius
  forced_zones      = forced_region_zones                         ! number of zones to apply laser forcing to
  m1                = mass1
  z1                = charge1
  m2                = mass2
  z2                = charge2
  init_press        = initial_pressure                            ! initial pressure in capsule 
  init_T_gas        = initial_T_gas                               ! initial gas temperature 
  rho_fac           = density_change                              ! ratio of shell density to gas density
  laser_watts_pulse = initial_watts_pulse                         ! Power in Watts supplied by laser
  init_watts_pulse  = laser_watts_pulse*laser_absorbed_frac       ! absorbed Power in Watts from laser
  ion_frac          = ionization_fraction                         ! fraction of atoms ionized by laser forcing  -? absorption fraction??
  out_density       = outerBC_density                             ! initialize outer BC  values 1.225
  out_temp          = outerBC_temp                                ! initialize outer BC  values 300.
  visc_sw           = viscosity_sw
  hf_sw             = heatflux_sw
  diff_sw           = diffusion_sw

  itprt             = it_steps_per_print         ! # of time steps between data output writes
  
  outfile           = output_file_name           ! name of output file
  t_out_prt         = t_out_special_print        ! time to dump a special write file

END SUBROUTINE read_params

SUBROUTINE initialize

  use exp_params
  use phys_const
  use data_arrays

  ! allocate memory for data arrays
  allocate( r(nr+2,2)      )
  allocate( u(nr+2,2)      )
  allocate( rho(nr+1,2)    )
  allocate( Temp(nr+1,2)   )
  allocate( TempIC(nr+1)   )
  allocate( P(nr+1,2)      )
  allocate( Comp(nr+1,2)   )
  allocate( dr(nr+1,2)     )
  allocate( Vol(nr+1,2)    )
  allocate( q(nr+1)        )
  allocate( cs(nr+1)       )
  allocate( N_part(nr+1)   )
  allocate( M(nr+1)        )
  allocate( vel_fac(nr+1)  )
  allocate( ionized(nr+1)  )
  allocate( N_elec(nr+1)   )
  allocate( subdiag(nr+1)  )
  allocate( diag(nr+1)     )
  allocate( supdiag(nr+1)  )
  allocate( vel_RHS(nr+2)  )
  allocate( temp_RHS(nr+1) )
  allocate( forcing(nr+1)  )
  allocate( r_sq(nr+2)     )

   ! open output files
     open(20,file=outfile,status='unknown',action='write')
  

  !  print *, dir, 'is dir'
   !   open file 11-16, etc to write multiple time density arrays for python scripts 
     
     !OPEN(UNIT = 11,FILE = path//dir//'time_xtOut',STATUS = 'unknown')
     !OPEN(UNIT = 12,FILE = dir//'radius_xtOut',STATUS = 'new')
     !OPEN(UNIT = 13,FILE = dir//'density_xtOut',STATUS = 'new')
     !OPEN(UNIT = 14,FILE = dir//'ein_xtOut',STATUS = 'new')
     !OPEN(UNIT = 15,FILE = dir//'pres_xtOut',STATUS = 'new')
     !OPEN(UNIT = 16,FILE = dir//'tev_xtOut',STATUS = 'new')
     !OPEN(UNIT = 17,FILE = dir//'velr_xtOut',STATUS = 'new')
     !OPEN(UNIT = 18,FILE = dir//'comp_xtOut',STATUS = 'new')
     !OPEN(UNIT = 19,FILE = dir//'dr_xtOut',STATUS = 'new')
     
     !OPEN(UNIT = 21,FILE = dir//'t_r_d_xtOut',STATUS = 'new')

     OPEN(UNIT = 11,FILE = 'time_xtOut',STATUS = 'unknown')
     OPEN(UNIT = 12,FILE = 'radius_xtOut',STATUS = 'unknown')
     OPEN(UNIT = 13,FILE = 'density_xtOut',STATUS = 'unknown')
     OPEN(UNIT = 14,FILE = 'ein_xtOut',STATUS = 'unknown')
     OPEN(UNIT = 15,FILE = 'pres_xtOut',STATUS = 'unknown')
     OPEN(UNIT = 16,FILE = 'tev_xtOut',STATUS = 'unknown')
     OPEN(UNIT = 17,FILE = 'velr_xtOut',STATUS = 'unknown')
     OPEN(UNIT = 18,FILE = 'comp_xtOut',STATUS = 'unknown')
     OPEN(UNIT = 19,FILE = 'dr_xtOut',STATUS = 'unknown')
     
     OPEN(UNIT = 21,FILE = 't_r_d_xtOut',STATUS = 'unknown')
     
 !    OPEN(UNIT = 98,FILE = 'diff_flux', STATUS = 'UNKNOWN')	
 !    OPEN(UNIT = 99,FILE = 'HeatFlux_Experimental', STATUS = 'UNKNOWN')


  ! initialize constants
  pi        = 2. * acos(0.0)
  dtheta    = pi/64.
  dphi      = pi/64.
  cQ        = 0.25*(gamma + 1)
  cL        = 0.5
  vol_fac   = 4./3.*pi    !   now it's vol of whole shell  1.   !   (dphi/3.)*(1. - cos(dtheta))   !   ???? 
  in_forced = nr-forced_zones+1   ! forced region from in_forced to nr, with outer BC state fixed in nr+1

  ! initialize matrix coefficients
  subdiag(:) = 0.
  diag(:)    = 0.
  supdiag(:) = 0.

  ! initialize diagnostic variables
  nanDetected = .FALSE.

  ! initialize velocity and grid positions
  u(:,:)  = 0.
  vel_fac(:) = 0.
  dr(:,1) = Lr/nr   ! uniform grid
  r(1,1)  = 0.
  do i = 2,(nr+2)
    r(i,1) = dr(i-1,1) * (i - 1)
  end do

  ! initialize volumes
  do i = 1,(nr+1)
    Vol(i,1) = vol_fac*( r(i+1,1)**3 - r(i,1)**3 )
  end do
 
  ! initialize density, cell masses, particle numbers
  rho(1:in_bound,1)                = (init_press*mm_DD)/(N_av*k_boltz*JeV_conv*init_T_gas)
  rho((in_bound+1):nr,1)           = rho_fac*rho(1,1)                       !  densities, rho [kg/m3]
  rho(nr+1,1)                      = out_density
  M(:)                             = rho(:,1)*Vol(:,1)                      !  total mass [kg] in cell
  N_part(1:in_bound)               = (N_av/mm_DD) * M(1:in_bound)           !  total number of particles in cell
  N_part((in_bound+1):(nr+1))      = (N_av/mm_CH) * M((in_bound+1):(nr+1)) 
  N_elec(1:in_bound)               =  z1    !   odd name - used as zbar  
  N_elec((in_bound+1):nr)          =  z2    !     but zbar by material not by zone - need to revise
  
  ! initialize ionization state
  ionized(:)                       =  0.0
  ionized(in_forced:nr)            =  ion_frac
  N_part(in_forced:nr)             =  (1+ionized(in_forced:nr)*N_elec(in_forced:nr)) *        &
                                      N_part(in_forced:nr)       !   add electrons to particle pressure 
                                                                 !      assuming in_forced region is ionized

  ! initialize temperature with forcing
  forcing(:)                       =  0
  Temp(1:in_bound,1)               =  k_boltz*JeV_conv*init_T_gas   ! init_T_gas[K], Temp[J]
  Temp((in_bound+1):nr,1)          =  (init_press * Vol((in_bound+1):nr,1)) / N_part((in_bound+1):nr)
  Temp((nr+1),1)                   =  k_boltz*JeV_conv*out_temp
  TempIC(1:nr+1)                   =  Temp(1:nr+1,1)
  
  ! initialize Deuterium mass fraction
  Comp(:,:) = 0
  Comp(1:in_bound,:) = 1

  ! initialize pressure
  P(:,1)        =  (N_part(:) / Vol(:,1)) * Temp(:,1)
  out_press     =  (N_part(nr+1) / Vol(nr+1,1)) * k_boltz * JeV_conv * out_temp      
  P((nr+1),1)   =  out_press 

  ! find sound speed
  cs(:) = sqrt( gamma * P(:,1) / rho(:,1) )
 
  ! initialize time variables
  time     = 0.
  dt       = min( init_dt, sqrt( (Lr/nr) / maxval(cs(:)) )  )
  timeold  = 0.
  it_special_prtOut = 0
  
  !  diagnostics for inputs
  print *,'first capsule zone & forced zone =',(in_bound+1),(nr-forced_zones+1)
  
  ! initialize diagnostic variables as needed
  highestPress   = P(1,1)
  highestTemp    = Temp(1,1)
  highestDensity = rho(1,1)

     
END SUBROUTINE initialize


SUBROUTINE solve_eqns

  use exp_params
  use phys_const
  use data_arrays

  !real compute_mass_fraction(nr+1)
  real dC(nr+1)

  ! update time
  time = time + dt

  ! find artificial viscosity
  do i = 1,(nr+1)
    du = u(i+1,1) - u(i,1)
    if (du .LT. 0) then
      q(i) = rho(i,1)*abs(du)*( cQ*abs(du) + cL*cs(i) )
    else
      q(i) = 0
    end if
  end do

  ! find new velocity
  do i = 2,(nr+1)
    dr_bar = .5 * (dr(i,1) + dr(i-1,1))
    inv_rho_bar =  ( (Vol(i,1)+Vol(i-1,1))/(M(i)+M(i-1)) ) 
    vel_fac(i) = inv_rho_bar * (P(i,1) + q(i) - P(i-1,1) - q(i-1)) / dr_bar
    call visc_mat_coeffs(subdiag,diag,supdiag)
  end do 
  
  vel_RHS(:) = u(:,1) - dt*vel_fac(:)

  u(:,2) = vel_RHS(:)  !!   test u w/o tridag
  
     call tri_diag(subdiag,diag,supdiag,vel_RHS(2:nr+1),u(2:nr+1,2),nr)

  u(nr+2,2) = u(nr+1,2)

  ! find new grid locations and cell lengths
  r(:,2)  = r(:,1) + .5*dt*(u(:,2) + u(:,1))
  r_sq(:) = r(:,2)*r(:,2)    !   r and r_sq at zone lower  face
  do i = 1,(nr+1)
    dr(i,2) = r(i+1,2) - r(i,2)
  end do

  ! find new volumes
  do i = 1,(nr+1)
    Vol(i,2) = vol_fac*(r(i+1,2)**3 - r(i,2)**3)
  end do

  ! find new densities
  rho(:,2) = rho(:,1) * (Vol(:,1)/Vol(:,2))
  rho((nr+1),2) = out_density
  
  ! continue forcing for 1 ns
    forcing(:) = 0
  if (time .LT. 1e-9) then
    forcing((nr-forced_zones+1):nr)  = init_watts_pulse / dr((nr-forced_zones+1):nr,2) /     &
   &                                       (4.*pi*r((nr-forced_zones+1):nr,2)**2)           !!  ( r(nr,2) - r(nr-forced_zones) )
  end if

  ! find new temperatures
  do i = 1,nr
    r_bar = .5*(r(i+1,2)+r(i,2))
    r_bar_sq = r_bar*r_bar     !   r_bar and r_bar_sq  at zone center
    call temp_mat_coeffs(subdiag,diag,supdiag) 
    
    temp_RHS(i) =  Temp(i,1) - dt* (gamma-1)* (Vol(i,2)/ N_part(i))*                          &
   &                          (0.5* P(i,1) + q(i))* ( r_sq(i+1)*u(i+1,2) - r_sq(i)*u(i,2) )   &
   &                           / (r_bar_sq*dr(i,2))
   
    temp_RHS(i) =  temp_RHS(i) +  dt* (gamma-1)* (Vol(i,2)/ N_part(i))* forcing(i)
  end do
  
    Temp(1:nr,2) = temp_RHS(:)   !  explicit temp from RMs tempRHS

  
    call tri_diag(subdiag,diag,supdiag,temp_RHS(1:nr),Temp(1:nr,2),nr)  
    Temp((nr+1),2) = k_boltz*JeV_conv*out_temp   !! units conversion from rewrite that was unfinished ??
  
  
  ! find new species mass fractions and number of particles
  do i = 1,nr
    r_bar = .5*(r(i+1,2)+r(i,2))
    r_bar_sq = r_bar*r_bar
    call comp_mat_coeffs(subdiag,diag,supdiag)  
  end do
  
     call tri_diag(subdiag,diag,supdiag,Comp(1:nr,1),Comp(1:nr,2),nr)
      do i = 1, nr
        if ( Comp(i,2) .lt. eps ) Comp(i,2) = 0.   !   python plotter doesn't like comp => 1.e-102
     enddo
   Comp((nr+1),2) = 0

   call compute_mass_fraction(dC)

  !if (mod(it,200).eq.0) then
  !  do i=1,nr
  !    print *,'dComp = ',dC(i)
  !    print *,'dCi =' , (Comp(i,2)-Comp(i,1))/dt
  !  end do

  !  print *,'step'
  !end if

   Comp(:,2) = Comp(:,1) + dt*dC(:)
   Comp((nr+1),2) = 0

  ! determine if zone is ionized based on temperature change and change N_part due to new Comp
  do i = 1,nr
    if ( (ionized(i) .EQ. 0) .AND. (Temp(i,2) .GT. (1.3*TempIC(i))  )  ) then
      ionized(i) = ion_frac;
    end if
  end do

  N_elec(:) = (z1*Comp(:,2)/mm_DD + z2*(1-Comp(:,2))/mm_CH)*ionized(:)*M(:)*N_av    !  ne - but not needed here 
  N_part(:) = N_av*M(:)*(  Comp(:,2)/ mm_DD* ( 1.+ionized(:)*z1 ) +(1.-Comp(:,2))/ mm_CH* ( 1.+ionized(:)*z2 )  )    !   z1 = 1 and z2 = z2 hardwired in
  !   do i=1,nr
  !    print *,'Npart = ',N_part(i)
  !    print *,'Nelec = ',N_elec(i)
  !    print *,'NT = ',N_elec(i)+N_part(i)
  !    print *,'oC =' , N_av*M(i)*(  Comp(i,2)/ mm_DD* ( 1.+ionized(i)*z1 ) + (1.-Comp(i,2))/ mm_CH* ( 1.+ionized(i)*z2 )  )
  !  end do
  ! find new pressures
  
  P(:,2) = (N_part(:)/Vol(:,2)) * Temp(:,2) 
  P((nr+1),2) = out_press

  ! find sound speed
  cs(:) = sqrt( gamma * P(:,2) / rho(:,2) )
  
  ! calculate new time step
  dt = cfl_fac * 0.5 * minval( dr(:,2)/(cs(:) + abs(u(:,2))) )
  dt = min( dt , 0.5*minval(dr(:,2)**2.0/1e-10))
  dt = min( dt, init_dt )


END SUBROUTINE solve_eqns 

SUBROUTINE reset_arrays

  use exp_params
  use data_arrays

  u(:,1)    = u(:,2)
  r(:,1)    = r(:,2)
  dr(:,1)   = dr(:,2)
  Vol(:,1)  = Vol(:,2)
  rho(:,1)  = rho(:,2)
  Temp(:,1) = Temp(:,2)
  P(:,1)    = P(:,2)
  Comp(:,1) = Comp(:,2)

END SUBROUTINE reset_arrays

SUBROUTINE visc_mat_coeffs (a,b,c)

  use prec_param
  use exp_params
  use phys_const
  use data_arrays

  REAL(prec), DIMENSION(nr), INTENT(out)  ::  a,b,c
  REAL(prec)                              ::  comb_dr, alpha, beta1, mu_plus, mu_minus
  REAL(prec)                              ::  eta_plus, eta_minus, r_plus, r_minus
  REAL(prec)                              ::  Temp_eV, Temp_eV_minus
  REAL(prec), PARAMETER                   ::  kb_ergeV = 1.6e-12

  Temp_eV        = Temp(i,1)/JeV_conv    ! convert temperature to eV
  Temp_eV_minus  = Temp(i-1,1)/JeV_conv  ! convert temperature to eV 

  mu_plus  = 1000*(mm_DD*Comp(i,1)/(1+ionized(i)) + mm_CH*(1-Comp(i,1))/(1+ionized(i)))
  mu_minus = 1000*(mm_DD*Comp(i-1,1)/(1+ionized(i-1)) + mm_CH*(1-Comp(i-1,1))/(1+ionized(i-1)))
  
  comb_dr        = .5*(dr(i,1) + dr(i-1,1))
  alpha          = dt*inv_rho_bar
  r_plus         = .5*(r(i+1,1) + r(i,1))
  r_minus        = .5*(r(i,1) + r(i-1,1))

  if (visc_sw .eq. 1) then
    eta_plus       = (1e-1)*(2.0064e7)*sqrt(mu_plus)*(.2)*kb_ergeV*(Temp_eV**2.5)         ! plasma viscosity on right side of zone boundary
    eta_minus      = (1e-1)*(2.0064e7)*sqrt(mu_minus)*(.2)*kb_ergeV*(Temp_eV_minus**2.5)  ! plasma viscosity on left side of zoen boundary
  else 
    eta_plus = 0
    eta_minus = 0
  end if
  beta1          = 1/((r(i,1)**2)*comb_dr)

  a(i-1) = -alpha*beta1*(r_minus**2)*eta_minus/dr(i-1,1)
  c(i-1) = -alpha*beta1*(r_plus**2)*eta_plus/dr(i,1)
  b(i-1) = 1 - (a(i-1) + c(i-1))

END SUBROUTINE visc_mat_coeffs

SUBROUTINE temp_mat_coeffs (a,b,c)

  use prec_param
  use exp_params
  use phys_const
  use data_arrays

  REAL(prec), DIMENSION(nr+1), INTENT(out)  ::  a,b,c
  REAL(prec)                                ::  dr_minus, dr_plus, alpha, beta1
  REAL(prec)                                ::  kappa, Temp_eV_minus, Temp_eV, Temp_eV_plus, logLambda
  REAL(prec), PARAMETER                     ::  kb_ergeV = 1.6e-12, ErgJ_conv = 1e-7, ccM_conv = 1e-6, gKg_conv = 1e-3, &
m_elec = 9.10938e-28

  

  !logLambda = log(8.799e29*(Temp(i,2)/(Comp(i,1)*z1+(1-Comp(i,1))*z2))**1.5/(ni+nj)**0.5)
  logLambda = 5.0

  if (logLambda .LE. 3.0) then
    logLambda = 3.0  
  else if (logLambda .GE. 1e4) then
    print *, logLambda, 'is L'
  end if
  
  Temp_eV        = Temp(i,1)/JeV_conv  ! convert temperature to eV
  if (i .EQ. 1) then
    dr_minus   = dr(i,2)
  else
    dr_minus   = .5*(dr(i,2) + dr(i-1,2))
  end if
    dr_plus      = .5*(dr(i+1,2) + dr(i,2))
  alpha        = dt*(gamma - 1)*(Vol(i,2)/N_part(i))
  kappa        = 1.0*1.93335e23*(Temp_eV**(5./2))/logLambda  ! plasma thermal conduction coefficient
  beta1        = kappa/(dr(i,2)*r_bar_sq)

  if (i .EQ. 1) then
    a(i) = 0.
  else
    a(i) = -alpha*beta1*r_sq(i)/dr_minus
  end if

  if (i .EQ. nr) then
    c(i) = 0.
  else
    c(i) = -alpha*beta1*r_sq(i+1)/dr_plus
  end if
  b(i) = 1  + .5*dt*(gamma - 1.)*(r_sq(i+1)*u(i+1,2) - r_sq(i)*u(i,2))/(r_bar_sq*dr(i,2)) - (a(i) + c(i))

END SUBROUTINE temp_mat_coeffs

SUBROUTINE comp_mat_coeffs (a,b,c)

  use prec_param
  use exp_params
  use phys_const
  use data_arrays

  REAL(prec), DIMENSION(nr+1), INTENT(out)  ::  a,b,c
  REAL(prec)                                ::  dr_minus, dr_plus, alpha, beta1, Temp_keV, Temp_eV
  REAL(prec)                                ::  diffc, rho_plus, rho_minus, DD_mass, CH_mass
  REAL(prec)                                ::  n1, n2, epsm, chi, coll_freq, logLambda, diffconst


  Temp_keV = Temp(i,2)/(1000.*JeV_conv)  ! convert temperature to keV
  Temp_eV = Temp(i,2)/(JeV_conv)
  DD_mass = mm_DD/N_av       ! calculating mass per particle
  CH_mass = mm_CH/N_av       ! calculating mass per particle

  n1 = (rho(i,2)*Comp(i,1))/mm_DD   !mol/zone
  n2 = (rho(i,2)*(1-Comp(i,1)))/mm_CH   ! mol/zone



  epsm = mm_CH/mm_DD
  chi = Comp(i,1)/(Comp(i,1)+(1.0-Comp(i,1))/epsm)
  logLambda = log(8.799e29*(Temp(i,2)/(Comp(i,1)*z1+(1-Comp(i,1))*z2))**1.5/(n1+n2)**0.5)

  if (logLambda .LE. 3.0) then
    logLambda = 3.0  
  else if (logLambda .GE. 1e4) then
    print *, logLambda, 'is L'
  end if
  
  !print *, logLambda, 'is L'

  !coll_freq = 3.81181e-18*logLambda/mm_DD/Temp(i,2)/((mm_DD+mm_CH)/(mm_CH*mm_DD)*Temp(i,2))**0.5
  diffconst = 3.2137e18/(z1*z2)**2/logLambda*(mm_CH*mm_DD*(mm_DD+mm_CH))**0.5*Temp(i,2)**1.5
  !if (coll_freq .EQ. 0.0) then
  !  coll_freq = 1e21
  !end if

  if (i .EQ. 1) then
    dr_minus = dr(i,2)
    rho_minus = rho(i,2)
  else 
    dr_minus = .5*(dr(i,2) + dr(i-1,2))
    rho_minus = .5*(rho(i,2) + rho(i-1,2))
  end if

  dr_plus  = .5*(dr(i+1,2) + dr(i,2))
  rho_plus = .5*(rho(i+1,2) + rho(i,2))
  alpha    = dt

  if (diff_sw .EQ. 1) then
    !diffc = epsm/(epsm+Comp(i,1)-epsm*Comp(i,1))**2*P(i,1)/rho(i,2)*mm_CH/coll_freq
    diffc = epsm/(epsm+Comp(i,1)-epsm*Comp(i,1))**2*P(i,1)/rho(i,2)*diffconst
  else
  	diffc = 0.0;
  end if

  !print *, diffc, 'is D'
  beta1    = 1.0/(rho(i,2)*r_bar_sq*dr(i,2))

  ! code for handling boundaries properly (shuld be implemented in other plasma routines)
  if (i .EQ. 1) then
    a(i) = 0
  else
    a(i)   = -alpha*beta1*r_sq(i)*rho_minus*diffc/dr_minus 
  end if

  if (i .EQ. nr) then
    c(i) = 0
  else
    c(i) = -alpha*beta1*r_sq(i+1)*rho_plus*diffc/dr_plus
  end if
  b(i)   = 1. - (a(i) + c(i))

END SUBROUTINE comp_mat_coeffs

SUBROUTINE tri_diag (a,b,c,d,u,n)

  use prec_param

  INTEGER                                ::  j
  INTEGER, INTENT(in)                    ::  n
  REAL(prec), DIMENSION(n), INTENT(in)   ::  a, b, c, d
  REAL(prec), DIMENSION(n), INTENT(out)  ::  u
  REAL(prec)                             ::  bet, gam(n)

  if (b(1) .EQ. 0) then
    print *, 'tridiag failed w/ b(1) = 0'
    read(*,*)
  end if

  bet = b(1)
  u(1) = d(1) / bet

  do j = 2,n
    gam(j) = c(j-1)/bet
    bet = b(j) - a(j)*gam(j)
    if (bet .EQ. 0) then
      print *, 'tridiag failed w/ bet=0'
      read(*,*)
    end if
    u(j) = (d(j) - a(j)*u(j-1)) / bet
  end do

  do j = n-1, 1, -1
    u(j) = u(j) - gam(j+1)*u(j+1)
  end do

END SUBROUTINE tri_diag

SUBROUTINE write_data_matlab     !     matlab format by Ryan Moll, July, 2013

  use exp_params
  use phys_const
  use data_arrays

  ! write data to file
  do i = 1,(nr+1)
    if ( (r(i,1) .NE. r(i,1))       .OR. &
         (u(i,1) .NE. u(i,1))       .OR. &
         (Temp(i,1) .NE. Temp(i,1)) .OR. &
         (rho(i,1) .NE. rho(i,1))   .OR. &
         (P(i,1) .NE. P(i,1))            &
       ) then
         nanDetected = .TRUE.
    end if 
    write(20,'(7e14.6)') time, r(i,1), u(i,1), Temp(i,1)/JeV_conv, (rho(i,1)*N_part(i))/((1.e6)*M(i)), P(i,1), Comp(i,1)
  end do

END SUBROUTINE write_data_matlab

SUBROUTINE write_data_pyout1     !     python output, eV, Sept, 2013

  use exp_params
  use phys_const
  use data_arrays

!  REAL(prec), ALLOCATABLE  ::  dTdr(:)
  REAL(prec) :: dTdr(nr-2,1), dCdr(nr-2,1)
  REAL(prec) :: kappa, diff, epsm, Temp_eV, chi, coll_freq
  REAL(prec) :: Temp_Interface, Comp_Interface
  integer    :: ir
       
	write (11,'(200e14.6)')   ( time*1.e9, i = 1,nr )  			    !  11 => 'time_xtOut' in ns
	write (12,'(200e14.6)')     r(1:nr,1)!*1.e6            	            !  12 => 'radius_xtOut'  in  microns
	write (13,'(200e14.6)')     rho(1:nr,1)          		            !  13 => 'density_xtOut'
	write (14,'(200e14.6)')     P(1:nr,1)/rho(1:nr,1)/(gamma - 1.)             !  14 => 'ein_xtOut'
	write (15,'(200e14.6)')     P(1:nr,1)                                      !  15 => 'pres_xtOut'
	write (16,'(200e14.6)')     Temp(1:nr,1)/Jev_conv                          !  16 => 'tev_xtOut'
  write (17,'(200e14.6)')     u(1:nr,1)         				    !  17 => 'velr_xtOut'
  write (18,'(200e14.6)')     Comp(1:nr,1)         			    !  18 => 'comp_xtOut'
  write (19,'(200e14.6)')     dr(1:nr,1)         			    !  19 => 'dr_xtOut'

!kappa        = (1.e2)*(.2)*(1.933e21)*(Temp_eV**(5./2))  ! plasma thermal conduction coefficient

!  do i=2,nr-1
!    Temp_Interface=0.5*(Temp(i-1,1)+Temp(i,1))/JeV_conv
!                  print *, Temp_Interface, 'eV'
!    kappa=(1.e2)*(.2)*(1.933e21)*(Temp_Interface**(5./2))
!                  print *, kappa, 'Plasma Thermal Conductivities'
!    dTdr(i-1,1)=-kappa*((Temp(i-1,1)-Temp(i,1))/dr(i,1))

!    Temp_eV = Temp(i,1)/(JeV_conv)
!	  epsm = mm_CH/mm_DD
!	  chi = Comp(i,1)/(Comp(i,1)+(1.0-Comp(i,1))/epsm)

!	  coll_freq = 4e16 * (1.0*z2)**2*5/ (mm_DD*(mm_DD*mm_CH)*1e-6/(mm_DD+mm_CH)*&
!	                (Temp_eV/mm_DD*1e-3+Temp_eV/mm_CH*1e-3)**1.5)

!	  diff = epsm/(epsm+Comp(i,1)-epsm*Comp(i,1))**2*P(i,1)*(1-Comp(i,1))/coll_freq
!    Comp_Interface=0.5*(Comp(i-1,1)+Comp(i,1))
!                  print *, Temp_Interface, 'eV'
!                  print *, kappa, 'Plasma Thermal Conductivities'
!    dCdr(i-1,1)=-1*((Comp(i-1,1)-Comp(i,1))/dr(i,1))
!  end do
		
!	write (98,'(200e14.6)')     dCdr(1:nr-2,1)
!  write (99,'(200e14.6)')     dTdr(1:nr-2,1)
         
  do ir = 2, nr
    if ( r(ir,1) .lt. 500.e-6)&      
      write (21,'(3e14.6)') time*1.e9, r(ir,1)*1.e6, Temp(ir,1)/Jev_conv
  end do


END SUBROUTINE write_data_pyout1


SUBROUTINE finish
  
  use exp_params
  use data_arrays
  use phys_const


  ! print diagnostic results
  print *, '***PROGRAM COMPLETE***'
  if (nanDetected) then
  print *, 'With fatal errors              : YES'
  else 
  print *, 'With fatal errors              : NO'
  end if
  print *, 'Final time                     : ', time
  print *, 'Highest pressure at origin     : ', highestPress
  print *, 'Highest temperature at origin  : ', highestTemp / JeV_conv
  Print *, 'Highest temperature occurs at  : ', time_of_highestTemp
  print *, 'Highest density at origin      : ', highestDensity * (N_part(1)/(1+ionized(1))*(1+ionized(1)*N_elec(1))) & 
                                                / ((1e6)*M(1))  

  print *,'optional write out solutions at last time step'
    !  call write_data_matlab
    !  call write_data_pyout1

  
  ! close parameter file
  close(30)
 
  ! close data file
  close(20)
  
  ! close python plot data files	  
	  close (unit=11)
	  close (unit=12)
	  close (unit=13)
	  close (unit=14)
	  close (unit=15)
	  close (unit=16)
	  close (unit=17)
	  close (unit=18)
          close (unit=19)

	  close (unit=21)
          close (unit=99)
	  
  ! deallocate memory from arrays
  deallocate( r        )
  deallocate( u        )
  deallocate( rho      )
  deallocate( Temp     )
  deallocate( P        )
  deallocate( Comp     )
  deallocate( dr       )
  deallocate( Vol      )
  deallocate( q        )
  deallocate( cs       )
  deallocate( N_part   )
  deallocate( M        )
  deallocate( vel_fac  )
  deallocate( ionized  )
  deallocate( N_elec   )
  deallocate( subdiag  )
  deallocate( diag     )
  deallocate( supdiag  )
  deallocate( vel_RHS  )
  deallocate( temp_RHS )
  deallocate( forcing  )
  deallocate( r_sq     )

END SUBROUTINE finish

!    ...............................................................................

   function rho2_per_nu12(m1,m2,z1,z2,tev1,tev2,teve)
   !   computes (rho2 / nu12) to modify binary diffusion coefficient adjusting to mass averaged u
   !   nu12 is collision rate of test particle 1 on 2 w/ kinetic factor 0.29 for low z - high z diffusion coefficient
   !   use n2 = n1, etc for self collision rates in visc. and thermal conduction
   
   implicit none
   real*4   n1,n2,m1,m2,z1,z2,tev1,tev2,teve,m12
   real*4   rho2_per_nu12, coulog
   real*4  ::  bk = 1.602e-19, avogn = 6.022e23, eps0 = 8.854e-12
   real*4  ::  KperP, debyeL, debyeco =  1.355e-11   !   = sqrt [ 2*eps0/(10^6*Avogn*bk) ]
   real*4  ::  fnuco = 4.1e16    !  4.1e-16 => 4 pi eps0^2 / e^4 /sqrt[mp] / bk**(1.5)
   real*4  ::  two = 2., one = 1., pi = 3.1415925
   
      debyeL = debyeco/ sqrt( ( n1*z1 + n2*z2) * ( one/tev1 + one/teve ) )
      KperP  = 12.*pi*eps0/bk*tev1*debyeL/z1/z2   !   kinetic particle energy to potential energy at debye radius
      coulog = 0.5*log( one + KperP*KperP )
      m12 = m1*m2/(m1+m2)
      print *,'tev1,tev2,teve = ',tev1,tev2,teve
      print *,'debyeL,KperP,coulog = ',debyeL,KperP,coulog
      rho2_per_nu12 =    fnuco * z1*z1*z2*z2*coulog/ ( m1*m12*(tev1/m1 + tev2/m2)**1.5 ) *0.29  ! hard wired high z-low z kinetic co
   end function rho2_per_nu12

subroutine compute_mass_fraction(dC)

use phys_const
use data_arrays
use exp_params
implicit none

real :: dC(nr+1), chi(nr+1), pions(nr+1), pelecs(nr+1)              !! Arrays


real :: oneoverrhorbarsqdr
real :: rhoD_p, rhoD_m, rsq_p, rsq_m, dr_p, dr_m   				    !! General variables for the differencing

real :: dchidcomp_p, dchidcomp_m, dCompdr_p, dCompdr_m				!! Molar Gradient
real :: pions_p, pions_m, XmY_p, XmY_m, dPionsdr_p, dPionsdr_m     	!! Ion Pressure

real :: epsm = mm_CH/mm_DD
real :: rho2overnu12(nr+1)
real :: logLambda

logLambda = 5.0 !log(8.799e29*(Temp(:,2)/(Comp(:,1)*z1+(1-Comp(:,1))*z2))**1.5/(ni+nj)**0.5)
rho2overnu12(:) = 3.2137e18/(z1*z2)**2/logLambda*(mm_CH*mm_DD*(mm_DD+mm_CH))**0.5*Temp(:,2)**1.5
!oneoverdr(:) = 1.0/dr(:,2)
!dchidcomp(:) = epsm/(epsm+Comp(:,1)-epsm*Comp(:,1))**2.0
chi(:) = Comp(i,1)/(Comp(i,1)+(1.0-Comp(i,1))/epsm)
pions(:) = P(:,2)-N_elec(:)/Vol(:,2)*Temp(i,2)

do i = 1,nr
    oneoverrhorbarsqdr = 1.0/rho(i,2)/(0.5*(r(i+1,2)+r(i,2)))**2.0/dr(i,2)

    if (i .eq. 1) then
      rsq_m = r(i,2)*r(i,2)
      rsq_p = r(i+1,2)*r(i+1,2)

      rhoD_p = 2.0/(1.0/(P(i+1,1)/rho(i+1,2)*rho2overnu12(i+1))+1.0/(P(i,1)/rho(i,2)*rho2overnu12(i)))
      rhoD_m = rhoD_p


      !! Molar Gradient
      dchidcomp_p = 0.5*(epsm/(Comp(i+1,1)+epsm*(1.0-Comp(i+1,1)))**2.0+&
      	epsm/(Comp(i,1)+epsm*(1.0-Comp(i,1)))**2.0)

      dchidcomp_m = dchidcomp_p

      dCompdr_m = 0.0
      dCompdr_p = (Comp(i+1,1) - Comp(i,1))/(0.5*(dr(i+1,2)+dr(i,2)))

      !! Ion Pressure Gradient
	  pions_p = 2.0/(1.0/pions(i+1) + 1.0/pions(i))
	  pions_m = pions_p

	  XmY_p = 2.0/(1.0/(chi(i+1)-Comp(i+1,1))+1.0/(chi(i)-Comp(i,1)))
	  XmY_m = XmY_p

	  dPionsdr_p = (pions(i+1) - pions(i))/(0.5*(dr(i+1,2)+dr(i,2)))
      dPionsdr_m = 0.0
    
    else if (i .eq. nr) then
      rsq_m = r(i,2)*r(i,2)
      rsq_p = r(i+1,2)*r(i+1,2)

      rhoD_m = 2.0/(1.0/(P(i-1,1)/rho(i-1,2)*rho2overnu12(i-1))+1.0/(P(i,1)/rho(i,2)*rho2overnu12(i)))
      rhoD_p = rhoD_m 

      !! Molar Gradient 
      dchidcomp_m = 0.5*(epsm/(Comp(i,1)+epsm*(1.0-Comp(i,1)))**2.0+&
      	epsm/(Comp(i-1,1)+epsm*(1.0-Comp(i-1,1)))**2.0)

      dchidcomp_p = dchidcomp_m

      dCompdr_p = 0
      dCompdr_m = (Comp(i,1) - Comp(i-1,1))/(0.5*(dr(i,2)+dr(i-1,2)))
      
      !! Ion Pressure gradient
      pions_m = 2.0/(1.0/pions(i-1) + 1.0/pions(i))
      pions_p = pions_m

	  XmY_m = 2.0/(1.0/(chi(i-1)-Comp(i-1,1))+1.0/(chi(i)-Comp(i,1)))
	  XmY_p = XmY_m

      dPionsdr_m = (pions(i) - pions(i-1))/(0.5*(dr(i,2)+dr(i-1,2)))
      dPionsdr_p = 0            !! Check this boundary

    else     
      rsq_m = r(i,2)*r(i,2)
      rsq_p = r(i+1,2)*r(i+1,2)

      rhoD_p = 2.0/(1.0/(P(i+1,1)/rho(i+1,2)*rho2overnu12(i+1))+1.0/(P(i,1)/rho(i,2)*rho2overnu12(i)))
      rhoD_m = 2.0/(1.0/(P(i-1,1)/rho(i-1,2)*rho2overnu12(i-1))+1.0/(P(i,1)/rho(i,2)*rho2overnu12(i)))

	!! Molar Gradient Diffusion      
      dchidcomp_p = 0.5*(epsm/(Comp(i+1,1)+epsm*(1.0-Comp(i+1,1)))**2.0+&
      	epsm/(Comp(i,1)+epsm*(1.0-Comp(i,1)))**2.0)

      dchidcomp_m = 0.5*(epsm/(Comp(i,1)+epsm*(1.0-Comp(i,1)))**2.0+&
      	epsm/(Comp(i-1,1)+epsm*(1.0-Comp(i-1,1)))**2.0)

      dCompdr_p = (Comp(i+1,1) - Comp(i,1))/(0.5*(dr(i+1,2)+dr(i,2)))
      dCompdr_m = (Comp(i,1) - Comp(i-1,1))/(0.5*(dr(i,2)+dr(i-1,2)))

	!! Ion Pressure Gradient
	  pions_p = 2.0/(1.0/pions(i+1) + 1.0/pions(i))
	  pions_m = 2.0/(1.0/pions(i-1) + 1.0/pions(i))

	  XmY_p = 2.0/(1.0/(chi(i+1)-Comp(i+1,1))+1.0/(chi(i)-Comp(i,1)))
	  XmY_m = 2.0/(1.0/(chi(i-1)-Comp(i-1,1))+1.0/(chi(i)-Comp(i,1)))

	  dPionsdr_p = (pions(i+1) - pions(i))/(0.5*(dr(i+1,2)+dr(i,2)))
      dPionsdr_m = (pions(i) - pions(i-1))/(0.5*(dr(i,2)+dr(i-1,2)))

      !YmZ_p = 2.0/(1.0/(Comp(i+1)-Comp(i+1,1))+1.0/(chi(i)-Comp(i,1)))
	  !YmZ_m = 2.0/(1.0/(chi(i-1)-Comp(i-1,1))+1.0/(chi(i)-Comp(i,1)))

	  !dPionsdr_p = (pions(i+1) - pions(i))/(0.5*(dr(i+1,2)+dr(i,2)))
      !dPionsdr_m = (pions(i) - pions(i-1))/(0.5*(dr(i,2)+dr(i-1,2)))
      
    end if

    dC(i) = oneoverrhorbarsqdr*(rhoD_p*rsq_p*&
    	(dchidcomp_p*dCompdr_p+0.0*1.0/pions_p*XmY_p*dPionsdr_p)-&
    	(dchidcomp_m*dCompdr_m+0.0*1.0/pions_m*XmY_m*dPionsdr_m)&
    	*rhoD_m*rsq_m)
    !print *,'oneoverrhorbarsqdr,dComp,rho2overnu12 = ',oneoverrhorbarsqdr(i),dComp(i),rho2overnu12(i)

end do
 


!dC(:) = oneoverrhorbarsqdr(:)*(P(:,1)/rho(:,2)*rho2overnu12(:))*(dchidcomp(:)*oneoverdr(:)*dComp(:))




end subroutine compute_mass_fraction