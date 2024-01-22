!#################################################
!#################################################          

! the module defines the global variables for a 
! variational MC for hard spheres in an Harmonic
! potential in 2D.

! The additional subroutine, 'input', read the 
! input gile and allocate the variables accordingly

!#################################################
!#################################################          

module DMC_parameters
    implicit none
    
    !define pi
    real*8, parameter :: PI=dacos(-1.d0)

    !space dimensions 
    integer           :: DIM = 2

    !MCparams
    !{
    integer           :: NMCsteps     !number of MC steps with default value 
    integer           :: NThermSteps  !number of thermalization steps       
    integer           :: NStabSteps   !number of stability steps    
    real*8            :: dt     
    !}

    !system parameters
    !{
    integer           :: N0walkers = 1 !number of walkers(in VMC = 1)
    integer           :: Natoms       !number of atoms/particles
    real*8            :: a_osc        !harmonic oscillator lenght
    real*8            :: Rv           !bessel function juntion
    real*8            :: k0           !k0 in trial WF
    real*8            :: alpha        !exponential decay parameter
    !}

    !other parameters
    !{
    integer           :: NdensProfileSteps
    real*8            :: densProfileStep
    integer           :: TWFNPartitions
    logical           :: USE_TABLE    !True if the simulation should use the table approx
    logical           :: INIT_CONF_FROM_FILE =.FALSE.
    !}



    !dimensional parameter
    !{
    real*8, parameter :: h2over2m = 1 !set to one
    real*8            :: D = 1        !unity of energy 
    !}

    !verbosity parameters
    !{
    logical           :: PRINT_DENSITY_PROFILE       = .FALSE.
    logical           :: PRINT_ENERGY_EVOLUTION      = .FALSE.
    logical           :: PRINT_INITIAL_PARAMETERS    = .FALSE.
    logical           :: PRINT_TRIAL_WAVEFUNCTION    = .FALSE.
    logical           :: PRINT_FINAL_CONFIGURATION   = .FALSE.
    
    integer,parameter  :: MAX_FILENAME_LENGHT = 50
    character(MAX_FILENAME_LENGHT) :: outfile_path = "./data/"
    character(MAX_FILENAME_LENGHT) :: infile_path  = "./data/"
    character(10)                  :: outfile_ext  = ".dat"
    character(MAX_FILENAME_LENGHT) :: dens_profile_filename     = "density_profile.dat"
    character(MAX_FILENAME_LENGHT) :: energy_evolution_filename = "energy_evolution.dat"
    character(MAX_FILENAME_LENGHT) :: init_gas_conf_filename    = "init_gas_conf.dat"
    character(MAX_FILENAME_LENGHT) :: twf_filename              = "twf.dat"
    character(MAX_FILENAME_LENGHT) :: final_conf_filename       = "final_conf.dat"
    
    !}
    
end module 
!######################################################

!######################################################
!               Input Subroutine
!######################################################
subroutine input(filename)
    use DMC_parameters
    use iofile,only: io_open,read_data,is_present

    implicit none

    
    character(*), intent(in) :: filename
    integer :: N_lines
    call io_open(input_filename = filename,N_lines = N_lines )

    !read MCparams
    call read_data("NMCsteps",NMCsteps)
    call read_data("NTsteps",NThermSteps)
    call read_data("NStabSteps",NStabSteps)
    call read_data("TWFNPartitions", TWFNPartitions)
    call read_data("use_table",USE_TABLE)
    call read_data("use_init_conf_from_file",INIT_CONF_FROM_FILE)
    call read_data("dt",dt)

    !read system parameters
    call read_data("N0walkers",N0walkers) !VMC supposed to be one 
    call read_data("Natoms",Natoms)
    call read_data("NdensProfileSteps",NdensProfileSteps)
    call read_data("a_osc",a_osc)
    call read_data("k0",k0)
    call read_data("Rv",Rv)
    call read_data("alpha",alpha)
    
    !read verbosity flag
    call read_data("DENSITY_PROFILE"      ,PRINT_DENSITY_PROFILE      )
    call read_data("ENERGY_EVOLUTION"     ,PRINT_ENERGY_EVOLUTION     )
    call read_data("INITIAL_PARAMETERS"   ,PRINT_INITIAL_PARAMETERS   )
    call read_data("TRIAL_WF"             ,PRINT_TRIAL_WAVEFUNCTION   )
    call read_data("FINAL_CONFIGURATION"   ,PRINT_FINAL_CONFIGURATION  )
    
    if (is_present("DENS_PROFILE_FILENAME")) then 
        call read_data("DENS_PROFILE_FILENAME", dens_profile_filename)
    end if 
    if (is_present("ENERGY_EVOLUTION_FILENAME")) then 
        call read_data("ENERGY_EVOLUTION_FILENAME", energy_evolution_filename)
    end if 
    if (is_present("TWF_FILENAME")) then 
        call read_data("TWF_FILENAME", twf_filename)
    end if
    if (is_present("FINAL_CONF_FILENAME")) then 
        call read_data("FINAL_CONF_FILENAME", final_conf_filename)
    end if
    if (PRINT_INITIAL_PARAMETERS) then 
        call print_parameters()
    endif    
end subroutine

subroutine print_parameters()
    use DMC_parameters
    implicit none
    character(LEN=13) :: FD = "(x,a,4x,f8.2)" !double format
    
    print *, "#####################################################"
    print *, " Montecarlo Simulation Using the following parameters"
    print *, "   - Number MC Steps       : ", NMCsteps
    print *, "   - Number Therm Steps    : ", NThermSteps
    print *, "   - Number Stability steps: ", NStabSteps
    print FD,"   - dt                    : ", dt
    print *, "   - Atoms Number          : ", Natoms
    print *, "   - Number initial walkers: ", N0walkers
    print *, "   - Number step dens prof : ", NdensProfileSteps
    print *, "   - TWF discret. N ofsteps: ", TWFNPartitions
    print *, "   - Use table             : ", use_table
    print *, "   - Read initconf fromfile: ", INIT_CONF_FROM_FILE
    print FD,"   - Rv                    : ", Rv
    print FD,"   - k0                    : ", k0
    print FD,"   - alpha                 : ", alpha
    print FD,"   - Harmonic lenght       : ", a_osc
    print *, " FLAGS:                                              "
    print *, "   - DENSITY PROFILE       : ", PRINT_DENSITY_PROFILE
    print *, "   - ENERGY EVOLUTION      : ", PRINT_ENERGY_EVOLUTION
    print *, "   - INITIAL PARAMETERS    : ", PRINT_INITIAL_PARAMETERS 
    print *, "   - TRIAL WAVEFUNCTION    : ", PRINT_TRIAL_WAVEFUNCTION
    print *, "   - FINAL CONFIGURATION   : ", PRINT_FINAL_CONFIGURATION
    print *, "######################################################"

end subroutine

subroutine read_initial_configuration_fromfile(R)
    Use, intrinsic :: iso_fortran_env, Only : iostat_end
    use DMC_parameters
    implicit none
    real*8,dimension(Natoms,DIM),intent(out) :: R
    integer :: IERROR
    character(MAX_FILENAME_LENGHT) :: filename
    character(len=100) :: ioerrmsg,line
    real*8 :: x,y
    integer :: i_atom, FID = 66
    filename = trim(outfile_path)//trim(init_gas_conf_filename)

    open(FID,file=filename,status='old',action='read')
    IERROR = 0
    read(FID,"(A)",iostat= IERROR,iomsg=ioerrmsg) line
    do while(IERROR == 0 .and. i_atom <= Natoms)
        read(FID,*,iostat= IERROR,iomsg=ioerrmsg) i_atom,x,y

        Select Case(IERROR)
        Case(0)
            R(i_atom,1) = x
            R(i_atom,2) = y
        Case(iostat_end)
            if(i_atom < Natoms) then 
                print*, "ERROR: Reading the initial configuration, less than expected"
                print*, "       expected  : ", Natoms
                print*, "       encountred: ", i_atom 
            end if 
            exit
        Case Default 
            print *, "IERROR: ", IERROR 
            print *, ioerrmsg
        End Select
    end do 
    close(FID)
    

end subroutine