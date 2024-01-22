
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
    real*8,allocatable,dimension(:,:)   :: F1,F2   !to store DRIFT
    real*8,allocatable,dimension(:,:)   :: R
    real*8,allocatable,dimension(:,:)   :: EL
    real*8,allocatable,dimension(:)     :: density_profile, tmp_density_profile
    real*8,dimension(2)                 :: TWF 


    type :: DMC_results
        real*8 :: E
        real*8 :: error
        real*8 :: cross_rate
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
        real*8  :: c1,sigma
        real*8  :: EL_new,E
        real*8  :: E_shift !to control the population growth
        real*8,dimension(Natoms,DIM)  :: R_TMP,R_DIFFUSED   !to store trial position 
        integer :: COUNT_REJECTED,COUNT_ACCEPTED
        integer :: i_walker, new_Nwalkers,Maxwalkers, N_sons,son 
        character(len=40) :: E_FMT

        E_FMT = "(A,I3.1,A,F5.3,A,I3.2,A,I2.1,A,I2.1)"
        sigma = sqrt(2*D*dt)
        Maxwalkers = N0walkers + N0walkers/2
        allocate(walker(Natoms,DIM,Maxwalkers,2)) !some extra space for population fluctuations
        allocate(               EL(Maxwalkers,2))
        allocate(density_profile(NdensProfileSteps))
        allocate(tmp_density_profile(NdensProfileSteps))
        allocate(F1(Natoms,DIM))
        allocate(F2(Natoms,DIM))

        
        !init all simulation parameters
        call init_patameters()
        COUNT_REJECTED = 0
        call init_random_seed()
        if(USE_TABLE) then 
            call gen_TWF_tables(Npartitions=TWFNPartitions) 
        end if 

        !gen walker position and starting energies
        E_shift = 0; COUNT_ACCEPTED = 0
        do i_walker=1,Nwalkers
            if( .not. INIT_CONF_FROM_FILE) then 
                call gen_initial_configuration(walker(:,:,i_walker,OLD))
            else
                call read_initial_configuration_fromfile(walker(:,:,i_walker,OLD))
            end if  
        end do

        !compute initial local energy and E_shift
        do i_walker=1, Nwalkers
            EL(i_walker,OLD) = Elocal(walker(:,:,i_walker,OLD))
            E_shift = E_shift + EL(i_walker,OLD)
        end do 
        E_shift = E_shift/real(Nwalkers,8)
        do MC_step= -NStabSteps, NMCsteps
            call start_clock()
            do step = 1, NThermSteps
                new_Nwalkers = 0
                E_acc        = 0
                tmp_density_profile = 0
                
                do i_walker = 1,Nwalkers
                    !first diffuse
                    call diffuse(R_IN  = walker(:,:,i_walker,OLD),&
                                 R_OUT = R_DIFFUSED, sigma = sigma)
                   
                    if( .not. check_hcore_crosses(R_DIFFUSED)) then 
                        COUNT_ACCEPTED = COUNT_ACCEPTED + 1
                        !first step R' = R + D*dt*F(R)/2 + gammma.
                        F1 = F(walker(:,:,i_walker,OLD)) !F(R)
                        R_TMP = R_DIFFUSED + D*dt*F1/2.       
                       
                        !second step R'' = R + D*dt*(F(R) + F(R'))/4.
                        R_TMP = R_DIFFUSED + D*dt*(F1+F(R_TMP))/4.   
                       
                        !compute energy E(R'') 
                        EL_new = Elocal(R_TMP)

                        !new generation
                        call random_number(c1)
                        N_sons = floor( dexp(-dt* ((EL_new + EL(i_walker,OLD))/2. - E_shift)) + c1)
                        

                        if(PRINT_DENSITY_PROFILE .and. MC_step > 0 .and. N_sons > 0) &
                            call update_density_profile(R_TMP,N_sons)

                        !third move R(R''') = R + D*dt*F(R'')/4.
                        R_TMP = R_DIFFUSED + D*dt*F(R_TMP)/2.
                    else
                        COUNT_REJECTED = COUNT_REJECTED + 1 
                        N_sons = 0
                    end if 
                    
                    !transmit to next generation 
                    do son = 1,N_sons
                        !update walkers number
                        new_Nwalkers = new_Nwalkers + 1
                        !append new walker to the end 
                        E_acc = E_acc + EL_new
                        walker(:,:,new_Nwalkers,NEW) = R_TMP
                        EL(new_Nwalkers,NEW)         = EL_new
                    end do 
                end do 
                
                NEW = 3 - NEW
                OLD = 3 - OLD
                Nwalkers = new_Nwalkers
                !adjust the energy to keep the pop +- constant
                E_shift = E_acc/real(Nwalkers,8) &
                          - 0.1/dt * log( real(Nwalkers,8)/N0walkers)
            end do !end thermalization loop
            E = E_acc/real(Nwalkers,8)
            call stop_clock(Time)
            print E_FMT, "MC_step: ",MC_step, &
                         "     EL: ",E,       &
                         "    pop: ",Nwalkers,& 
                         "     in: ",Time%minutes,":",Time%seconds

            if(MC_step > 0) then !skip all the stability steps
                E_avg  = E_avg  * real(MC_step-1,kind=8)/real(MC_step,8) &
                         + (E)   /real(MC_step,8)
                E2_avg = E2_avg * real(MC_step-1,kind=8)/real(MC_step,8) &
                         + (E**2)/real(MC_step,8)
                if(PRINT_ENERGY_EVOLUTION) &
                    call update_energy_accumulators()

                if(PRINT_DENSITY_PROFILE) &
                    density_profile = density_profile + tmp_density_profile/Nwalkers
                
            end if 
        end do 

        results%E     = E_avg
        results%error = sqrt( (E2_avg - E_avg**2))
        results%cross_rate = COUNT_REJECTED/COUNT_ACCEPTED


        if(PRINT_DENSITY_PROFILE) &
            call print_density_profile_toFile(density_profile)

        deallocate(walker)
        deallocate(density_profile)
        deallocate(tmp_density_profile)
        deallocate(EL)
        deallocate(F1);deallocate(F2)
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
            call get_energies(walker(:,:,i_walker,NEW),res_E)
            acc_E%Ekin    = acc_E%Ekin    + res_E%Ekin
            acc_E%Ekinfor = acc_E%Ekinfor + res_E%Ekinfor
            acc_E%Epot    = acc_E%Epot    + res_E%Epot
            acc_E%EL      = acc_E%EL      + res_E%EL 
        end do 
        acc_E%Ekin    = acc_E%Ekin   /Nwalkers
        acc_E%Ekinfor = acc_E%Ekinfor/Nwalkers
        acc_E%Epot    = acc_E%Epot   /Nwalkers
        acc_E%EL      = acc_E%EL     /Nwalkers
        call print_E_evolution_toFile(MC_step,E_avg,acc_E)

    end subroutine

    subroutine update_density_profile(POS,son_number)
        use DMC_parameters
        use HS_puregas
        use array_utility,Only:increase_size
        implicit none

        real*8,dimension(Natoms,DIM) :: POS
        integer :: son_number,i_atom,i_step
        integer :: Nenlarging,new_size
        real*8  :: radius



          do i_atom = 1, Natoms
            radius = norm2(POS(i_atom,:)) 
            if (radius > MAX_RADIUS) then 
                Nenlarging = floor( (radius - MAX_RADIUS)/densProfileStep)
                Nenlarging = Nenlarging + 2 !adding some extra space 
                new_size   = size(density_profile) + Nenlarging
                
                call increase_size(array=density_profile,     new_size=new_size)
                call increase_size(array=tmp_density_profile, new_size=new_size)
                
                NdensProfileSteps = size(density_profile)
                MAX_RADIUS        = densProfileStep*new_size
            end if 
            i_step = floor(radius/densProfileStep) + 1
            density_profile(i_step) = density_profile(i_step) + son_number
        end do

        
        ! do i_walker = 1, Nwalkers
        !     do i_atom = 1, Natoms
        !         radius = norm2(walker(i_atom,:,i_walker,OLD)) 
        !         if (radius > MAX_RADIUS) then 
        !             Nenlarging = floor( (radius - MAX_RADIUS)/densProfileStep)
        !             Nenlarging = Nenlarging + 2 !adding some extra space 
        !             new_size   = size(density_profile) + Nenlarging
                    
        !             call increase_size(array=density_profile, new_size=new_size)
                    
        !             NdensProfileSteps = size(density_profile)
        !             MAX_RADIUS        = densProfileStep*new_size
        !         end if 
        !         i_step = floor(radius/densProfileStep) + 1
        !         density_profile(i_step) = density_profile(i_step) + 1/real(Nwalkers,8)
        !     end do
        ! end do 

        
    end subroutine



end module DMC