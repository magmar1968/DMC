
module DMC
    use DMC_parameters
    use DMC_print
    implicit none
    
    integer                             :: MC_step, step
    integer                             :: NEW,OLD,Nwalkers
    real*8                              :: MAX_RADIUS
    real*8                              :: E_acc,E2_acc, Ekin_acc, Ekinfor_acc, Epot_acc
    real*8                              :: E_avg,E2_avg
    real*8,allocatable,dimension(:,:,:,:) :: walker
    real*8,allocatable,dimension(:,:,:,:) :: DRIFT   !to store DRIFT
    real*8,allocatable,dimension(:,:)   :: R
    real*8,allocatable,dimension(:,:)   :: EL
    real*8,allocatable,dimension(:)     :: density_profile
    real*8,dimension(2)                 :: TWF 


    type :: DMC_results
        real*8 :: E
        real*8 :: error
    end type

    contains

    !####################################################
    !#           Run Variational MC Simulation          #
    !####################################################    
    subroutine run_DMC(results)
        use random,Only:init_random_seed
        use HS_puregas
        use clock
        implicit none 
        
        type(DMC_results),intent(out) :: results
        type(time_res) :: Time
        real*8  :: c1
        real*8  :: EL_new,E
        real*8  :: E_shift !to control the population growth
        real*8,dimension(Natoms,DIM)  :: R_TMP   !to store trial position 
        integer :: COUNTER
        integer :: i_walker, new_Nwalkers, N_sons,son 

        allocate(walker(Natoms,DIM,N0walkers + N0walkers/3,2)) !some extra space for population fluctuations
        allocate( DRIFT(Natoms,DIM,N0walkers + N0walkers/3,2)) !some extra space for population fluctuations
        allocate(               EL(N0walkers + N0walkers/3,2))
        allocate(density_profile(NdensProfileSteps))
        
        !init all simulation parameters
        call init_patameters()
        COUNTER = 0
        call init_random_seed()
        if(USE_TABLE) then 
            call gen_TWF_tables(Npartitions=TWFNPartitions) 
        end if 

        !gen walker position and starting energies
        E_shift = 0
        do i_walker=1,Nwalkers
            if( .not. INIT_CONF_FROM_FILE) then 
                call gen_initial_configuration(walker(:,:,i_walker,OLD))
            else
                call read_initial_configuration_fromfile(walker(:,:,i_walker,OLD))
            end if  
            EL(i_walker,OLD) = Elocal(walker(:,:,i_walker,OLD))
            E_shift = E_shift + EL(i_walker,OLD)
        end do
        E_shift = E_shift/real(Nwalkers,8)
 
        do MC_step= -NStabSteps, NMCsteps
            call start_clock()
            do step = 1, NThermSteps
                new_Nwalkers = 0
                E_acc        = 0
                
                do i_walker = 1,Nwalkers

                    !generate trial move
                    call gen_new_particle_position(walker(:,:,i_walker,OLD), &
                                                   R_TMP,dt=dt)
                    
                    ! DRIFT(:,:,i_walker,OLD) = F(R_TMP)
                    
                    !check if the new position has some hard core crossing
                    if( .not. check_hcore_crosses(R_TMP)) then 
                        EL_new = Elocal(R_TMP)!DRIFT(:,:,i_walker,OLD) 
                        !new generation
                        call random_number(c1)
                        N_sons = floor( dexp(-dt* ((EL_new + EL(i_walker,OLD))/2. - E_shift)) + c1)
                    else 
                        N_sons = 0
                    end if 

                    !transmit to next generation 
                    do son = 1,N_sons
                        !update walkers number
                        new_Nwalkers = new_Nwalkers + 1
                        !append new walker to the end 
                        E_acc = E_acc + EL_new
                        walker(:,:,new_Nwalkers,NEW) =  R_TMP
                        ! DRIFT(:,:,new_Nwalkers,NEW)  = DRIFT(:,:,i_walker,OLD)
                        EL(new_Nwalkers,NEW) = EL_new
                    end do 
                end do 
                NEW = 3 - NEW
                OLD = 3 - OLD
                Nwalkers = new_Nwalkers
                !adjust the energy to keep the pop +- constant
                E_shift = E_acc/real(Nwalkers,8) &
                          - 0.1/dt * log( real(Nwalkers,8)/N0walkers)
                
            end do !end thermalization loop

            call stop_clock(Time)
            print *, "MC_step: ", MC_step,"pop:",Nwalkers, "in:" ,Time%minutes,Time%seconds
            if(MC_step > 0) then !skip all the stability steps
                E = E_acc/real(Nwalkers,8)
                E_avg  = E_avg  * real(MC_step-1,kind=8)/real(MC_step,8) &
                                        + (E)   /real(MC_step,8)
                E2_avg = E2_avg * real(MC_step-1,kind=8)/real(MC_step,8) &
                                        + (E**2)/real(MC_step,8)
                if(PRINT_ENERGY_EVOLUTION) &
                    call update_energy_accumulators()

                if(PRINT_DENSITY_PROFILE) &
                    call update_density_profile()
            end if 
        end do 

        results%E     = E_avg
        results%error = sqrt( (E2_avg - E_avg**2))

        if(PRINT_DENSITY_PROFILE) &
            call print_density_profile_toFile(density_profile)

        deallocate(walker)
        deallocate(DRIFT)
        deallocate(density_profile)
        deallocate(EL)
    end subroutine

    subroutine init_patameters()
        use DMC_parameters
        use DMC_print
        use HS_puregas,Only:set_HS_parameters
        implicit none
        
        call set_HS_parameters(alpha=alpha,Rv=Rv,k0=k0,a_osc=a_osc)
        !density profile stuff
        MAX_RADIUS      = 3*a_osc
        densProfileStep = MAX_RADIUS/NdensProfileSteps
        density_profile = 0

        !init swappers
        NEW = 2; OLD = 1
        !set all energies acc to 0
        E_acc       = 0
        E2_acc      = 0
        Ekin_acc    = 0
        Ekinfor_acc = 0
        Epot_acc    = 0  

        Nwalkers = N0walkers
        
        if (PRINT_TRIAL_WAVEFUNCTION) then 
            call print_TWF_tofile()
        end if

    end subroutine

    subroutine update_energy_accumulators()
        use HS_puregas, Only:Energies,get_energies
        implicit none
        type(Energies) :: acc_E,res_E !resulting energies
        integer :: i_walker
        
        
        do i_walker = 1, Nwalkers
            call get_energies(walker(:,:,i_walker,OLD),res_E,DRIFT(:,:,i_walker,OLD))
            acc_E%Ekin    = acc_E%Ekin   + res_E%Ekin
            acc_E%Ekinfor = acc_E%Ekinfor + res_E%Ekinfor
            acc_E%Epot    = acc_E%Epot    + res_E%Ekinfor
            acc_E%EL      = acc_E%EL      + res_E%EL 
        end do 
        acc_E%Ekin    = acc_E%Ekin   /Nwalkers
        acc_E%Ekinfor = acc_E%Ekinfor/Nwalkers
        acc_E%Epot    = acc_E%Epot   /Nwalkers
        acc_E%EL      = acc_E%EL     /Nwalkers
        call print_E_evolution_toFile(MC_step,E_avg,acc_E)

    end subroutine

    subroutine update_density_profile()
        use DMC_parameters
        use HS_puregas
        use array_utility,Only:increase_size
        implicit none

        integer :: i_walker
        real*8     :: radius 
        integer    :: i_atom,i_step
        integer    :: Nenlarging,new_size

        
        do i_walker = 1, Nwalkers
            do i_atom = 1, Natoms
                radius = norm2(walker(i_atom,:,i_walker,OLD)) 

                if (radius > MAX_RADIUS) then 
                    Nenlarging = floor( (radius - MAX_RADIUS)/densProfileStep)
                    Nenlarging = Nenlarging + 2 !adding some extra space 
                    new_size   = size(density_profile) + Nenlarging
                    
                    call increase_size(array=density_profile, new_size=new_size)
                    
                    NdensProfileSteps = size(density_profile)
                    MAX_RADIUS        = densProfileStep*new_size
                end if 
                i_step = floor(radius/densProfileStep) + 1
                density_profile(i_step) = (density_profile(i_step) + 1)/real(Nwalkers,8)
            end do
        end do 
        
    end subroutine



end module DMC