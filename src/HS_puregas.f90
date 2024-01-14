!#########################################################
!#########################################################

! Hard sphere, 2D pure gas implementation for variational
! MC algorithm. This module uses the module "parameter" to 
! read some of the data and take as argoument Rv,alpha and 
! k0 since eventually this code is used to variatonaly 
! evaluate the values of the first two parameters. 


!    0                                                            r < R_{core}
!    J_0(k_0 r) - J_0(k_0)/Y_0(k_0) * Y_0(k_0 * r)     R_{core} < r <= R_v
!    1 - \Gamma * \exp{-r/\alpha}                                 r > R_v


! File content:

! Subroutine                           Label               Line    
! ----------                           ---                 ----
! 
! Initialize Module                    init
!
! End Module                           end 
!
! Generate Vectorize TWF                gen_vec_TWF
!
! Generate Initial COnfiguration       gen_initial_configuration
!
! Check Hard core Crosses              check_hcore_crosses   
!
! TRIAL-WF
!   - Compute Trial Wafefunction       trial_WF
!   - Two Body Function                twobody_corr
!   - Harmonic oscillator Ground State harmonic_GS
!
! DERIVATIVES
!   - Twobody WF First Derivative      twobody_corrprime
!   - Twobody WF Second Derivative     twobody_corrdoubleprime
!   - Harmonic Osc. GS First Deriv.    harmonic_GSprime
!   - Harmonic Osc. GS Second Deriv.   harmonic_GSdoubleprime
!
!
! ELOCAL
!   - Compute Local Energy             Elocal
!   - Compute Potential Energy         Epot
!   - Compute Kinetic Energy           Ekin
!   
! DRIFT
!   - Compute Drift Force              F
!#################################################################
! Last Updated 29/11/23 L.Magnoni lorenzo.magnoni@studenti.unimi.it
module HS_puregas
    implicit none
       

    !trial Wavefunction parameters    
   
    real*8,private,save :: C          
    real*8,private,save :: Gamma
    real*8,private,save :: D = 1. !energy scale(hbar^2/2m Rcore^2)

    integer,private,save :: NTWFpartitions
    real*8,private,save  :: drTable
    
    type :: HS_parameters
        real*8 Rcore
        real*8 alpha
        real*8 Rv
        real*8 k0
        real*8 a_osc
    end type
    type(HS_parameters),private,save :: params

    type :: Energies
        real*8 :: EL      = 0.
        real*8 :: Ekin    = 0.
        real*8 :: Ekinfor = 0.
        real*8 :: Epot    = 0.
    end type



    real*8,allocatable,dimension(:) :: table_coor_func
    logical,private :: use_table = .FALSE.

    !public functions
    public :: set_HS_parameters
    public :: gen_TWF_tables
    public :: free_TWF_tables
    public :: gen_initial_configuration
    public :: trial_WF
    public :: gen_new_particle_position
    public :: diffuse
    public :: Elocal
    public :: get_energies
    public :: F

    contains

    !####################################################
    !           Set HS Parameters
    !####################################################
    subroutine set_HS_parameters(alpha,Rv,k0,a_osc,Rcore)
        implicit none
        real*8, intent(in) :: alpha,Rv,k0,a_osc
        real*8, intent(in), optional :: Rcore

        if( .not. present(Rcore)) then 
            params = HS_parameters(alpha=alpha,Rv=Rv,k0=k0,a_osc=a_osc,Rcore=1) 
        else 
            params = HS_parameters(alpha=alpha,Rv=Rv,k0=k0,a_osc=a_osc,Rcore=Rcore)
        end if 

        C = - bessel_j0(params%k0*params%Rcore)/&
              bessel_y0(params%k0*params%Rcore)

        Gamma =  -(params%alpha*params%k0) /exp(- params%Rv/ params%alpha) * &
                  (bessel_j1(params%k0*params%Rv) - bessel_j0(params%k0*params%Rcore)/&
                   bessel_y0(k0*params%Rcore) * bessel_y1(k0*params%Rv))
    end subroutine
    !####################################################

    subroutine free_TWF_tables()
        deallocate(table_coor_func)
    end subroutine free_TWF_tables


    !####################################################
    !#              Set energy scale                    #
    !####################################################
    subroutine set_energy_scale(newD)
        real*8,intent(in) :: newD
        D = newD 
    end subroutine


    !####################################################
    !#           get Trial Wave Function Terms          #
    !####################################################
    subroutine gen_TWF_tables(Npartitions)
        implicit none
        integer,intent(in) :: Npartitions
        real*8  :: r
        integer :: i_step

        use_table = .TRUE.
        NTWFpartitions = Npartitions

        drTable = 8*params%a_osc/Npartitions
        ALLOCATE(table_coor_func(Npartitions))

        do i_step = 1, Npartitions
            !harmonic part 
            r = (i_step-1)*drTable
            
            if( r <= params%Rcore) then 
                table_coor_func(i_step) = 0.
            else 
                table_coor_func(i_step) = twobody_corr( r )

            end if 
        end do   
    end subroutine


    !####################################################
    !#           Generate Initial Configuration         #
    !####################################################
    subroutine gen_initial_configuration(R)
        use random,Only:gauss
        implicit none
        
        real*8,intent(out) :: R(:,:)
        integer :: Natoms, DIM
        integer :: i_atom, j_dim
        logical :: regen

        Natoms = size(R,1);DIM = size(R,2)
        !generate randomly initial positions
        regen = .TRUE.
        do while (regen .eqv. .TRUE.)
            regen = .TRUE.
            do i_atom = 1, Natoms !for each atom and dimension    
                ! THERE MIGHT BE DIFFERENT WAY TO DO THAT
                do j_dim = 1, DIM !gen position
                    R(i_atom,j_dim) = gauss(sigma=(params%a_osc))
                end do
            end do 
            !check there's not hard core crossing with the previously generated atoms
            regen = check_hcore_crosses(R)
        end do
    end subroutine 
    !####################################################

    !####################################################
    !#            Check Hard Core Crosses               #
    !####################################################
    !check in an array of positions if at least two coordinates are closer than 
    !Rcore and return .TRUE. if they do
    logical function check_hcore_crosses(R)
        implicit none
        real*8, intent(in), dimension(:,:) :: R
        real*8 :: dist
        integer :: i_atom,j_atom        
        integer :: Natoms,DIM
        check_hcore_crosses = .FALSE.

#ifndef TURNOFF_INTERACTION
        Natoms = size(R,1);DIM = size(R,2)
        do i_atom=1, Natoms -1 
            do j_atom = i_atom + 1, Natoms
                dist = norm2(R(i_atom,:) - R(j_atom,:))
                !if two atoms are closer than Rcore return True
                if (dist <= params%Rcore) then 
                    check_hcore_crosses = .TRUE.
                    return 
                end if 
            end do 
        end do 
#endif
    end function check_hcore_crosses
    !####################################################

    
    !####################################################
    !#           Compute Trial Wavefunction             #
    !####################################################
    ! return the value of the trial WF for a specific configuration R
    ! summing both the interaction between particle and potential and 
    ! the interaction inbetween particles. Note that to we are exploiting
    ! the symmetry of the interaction to reduce the computational cost and 
    ! therefore we have to include a factore 2 in the eq. for the twobody 
    ! correlation. We use logarithm to trasform a production into a summatory
    ! which is much faster to compute.
    real*8 function trial_WF(R)
        implicit none
        
        real*8, dimension(:,:), intent(in) :: R
        integer :: i_atom,j_atom
        integer :: Natoms, DIM
        integer :: index
        real*8  :: frac,remainder
        real*8  :: dist, u,x

        Natoms = size(R,1); DIM = size(R,2)
        u = 0.

#ifndef TURNOFF_INTERACTION
        !particle interaction
        do i_atom = 1, Natoms -1 
            do j_atom = i_atom + 1, Natoms
                dist = norm2(R(i_atom,:)-R(j_atom,:))
                !exact method 
                if (use_table .eqv. .FALSE.) then 
                    u = u + 2*dlog(twobody_corr(dist))
                else 
                !tabular method 
                    frac  = dist/drTable
                    index = floor(frac)
                    remainder = frac - index
                    u = u + 2 * dlog( (1.-remainder)*table_coor_func(index)   + &
                                      remainder *table_coor_func(index+1))    
                end if 
            end do
        end do
#endif
        !potential 
        do i_atom = 1, Natoms
            dist = norm2(R(i_atom,:))
            x    = dist/params%a_osc
            u    = u + (-x**2/2.)
        end do 

        trial_WF = u        
        return 
    end function trial_WF


    !####################################################
    !#       Compute Twobody Correlation Function       #
    !####################################################
    !return the value of the two body wavefunction taking as input 
    !the distance between the two particles
    real*8 function twobody_corr(r)
        real*8, intent(in) :: r
        twobody_corr = 0

        if (r > params%Rv) then 
            twobody_corr = 1 - Gamma* exp(-r / params%alpha)
            return
        else if ( r > params%Rcore .AND. r <= params%Rv ) then
            twobody_corr = bessel_j0(params%k0*r) + C *&
                           bessel_y0(params%k0*r)
            return
        else
            twobody_corr = -1000
            return  
        end if 
    end function twobody_corr

    !return the value of the harmonic oscilator GS taking as input 
    !the distance from the potential well
    real*8 function harmonic_GS(r)
        real*8, intent(in) :: r
        real*8 N,x    
        if (r < 0.) then
            harmonic_GS = -1000
            return 
        end if 
        !GS expressed in natural lenght
        N = 1
        x   = r/params%a_osc    !reduce quantity
        harmonic_GS = N * exp(-(x**2)/2)
        return
    end function 
    !#####################################################


    !#####################################################
    !#                 Derivatives                       #
    !#####################################################  
    real*8 function twobody_corrprime(r)
        real*8, intent(in) :: r
        twobody_corrprime = 0
    
        if (r > params%Rv) then 
            twobody_corrprime =  (Gamma * exp(-r / params%alpha))/params%alpha
            return 
        else if ( r > params%Rcore .AND. r <= params%Rv ) then
            twobody_corrprime = -params%k0* ( bessel_j1(params%k0 *r) + &
                                          C * bessel_y1(params%k0*r))
            return
        else 
            twobody_corrprime = -1000
            return
        end if
    end function
    
    real*8 function twobody_corrdoubleprime(r)
        real*8, intent(in) :: r
        twobody_corrdoubleprime = 0.

        if( r > params%Rv) then 
            twobody_corrdoubleprime =  &
                 - (Gamma * exp(-r / params%alpha))/(params%alpha**2)
            return 
        else if (r > params%Rcore .AND. r <= params%Rv ) then 
            twobody_corrdoubleprime =  -params%k0**2/2. * &
                (&
                        ( bessel_j0(params%k0*r) - bessel_jn(2,params%k0*r))&
                   +  C*( bessel_y0(params%k0*r) - bessel_yn(2,params%k0*r))&
                )
            return 
        else 
            twobody_corrdoubleprime = -1000
            return 
        end if 
    end function 

    real*8 function harmonic_GSprime(r)
        real*8, intent(in) :: r
        real*8  :: x,J
        if (r < 0.) then
            harmonic_GSprime = -1000
            return 
        end if
        x = r/params%a_osc    !reduce quantity
        J = 1./params%a_osc   !jacobian
        harmonic_GSprime = J*( -x * exp(-(x**2)/2))
        return
    end function 

    real*8 function harmonic_GSdoubleprime(r)
        real*8, intent(in) :: r
        real*8  :: x,J      
        if (r < 0.) then
            harmonic_GSdoubleprime = -1000
            return 
        end if
        x = r/params%a_osc    !reduce quantity
        J = 1./params%a_osc   !Jacobian
        harmonic_GSdoubleprime = J**2*( - exp(-(x**2)/2) + (x**2) *exp(-(x**2)/2))
        return 
    end function
    !######################################################

    !####################################################
    !#           Get Trial Wave Function Terms          #
    !####################################################
    subroutine get_TWF_terms(TWF_corr,TWF_harm,rstep,Npartitions)
        implicit none
        integer,intent(in) :: Npartitions
        real*8,intent(out),dimension(3,Npartitions) :: TWF_corr,TWF_harm
        real*8,intent(out) :: rstep
        real*8  :: r
        integer :: i_step

        rstep = 4.*params%a_osc/Npartitions

        do i_step = 1, Npartitions
            !harmonic part 
            r = (i_step-1)*rstep
            
            if( r <= params%Rcore) then 
                TWF_corr(1,i_step) = 0.
                TWF_corr(2,i_step) = 0.
                TWF_corr(3,i_step) = 0.
            else 
                TWF_corr(1,i_step) = twobody_corr( r )
                TWF_corr(2,i_step) = twobody_corrprime( r )
                TWF_corr(3,i_step) = twobody_corrdoubleprime(r )
            end if 
            TWF_harm(1,i_step) = harmonic_GS(r)
            TWF_harm(2,i_step) = harmonic_GSprime(r)
            TWF_harm(3,i_step) = harmonic_GSdoubleprime(r)
        end do   
    end subroutine

    !####################################################
    !#            Gen new particle position             #
    !####################################################
    ! Diffuse and Drift a bunch of particles, optionally
    ! optionally the user can decide how many. The new
    ! state is registred in R_OUT while R_in remain 
    ! unchanged. 

    subroutine gen_new_particle_position(R_IN,R_OUT,DRIFT_IN,dt)
        implicit none
        real*8,  intent(in),  dimension(:,:) :: R_IN 
        real*8,  dimension(:,:),optional     :: DRIFT_IN
        real*8,  intent(out), dimension(:,:) :: R_OUT
        ! real*8, dimension(size(R_in,dim=1),size(R_in,dim=2)) :: R_TEMP
        real*8, dimension(size(R_in,dim=1),size(R_in,dim=2)) :: DRIFT
        real*8,  intent(in) ::  dt
        real*8 :: sigma 
        sigma = sqrt(2*D*dt)
        if( .not. present(DRIFT_IN)) then
            DRIFT = F(R_IN)
        else 
            DRIFT = DRIFT_IN
        end if 

        R_OUT = R_IN + D*dt*DRIFT/ 2.
        call diffuse(R_OUT,R_OUT,sigma)
        R_OUT = R_OUT + D*dt* ( F(R_OUT) + DRIFT )/4.
        call diffuse( R_OUT, R_OUT,sigma )
    end subroutine


    !####################################################
    !#             Diffuse the particles                #
    !####################################################
    ! Diffuse a bunch or particle, optionally the uses can 
    ! decide how many, trying a gaussian step using as 
    ! variance given in input. The new state is 
    ! registered in R_OUT while R_in remain untouched.
    subroutine diffuse(R_IN,R_OUT,sigma, NatomsToDiffuse)
        use random
        implicit none
        
        integer, intent(in),  optional :: NatomsToDiffuse
        real*8,  intent(inout),  dimension(:,:) :: R_IN 
        real*8,  intent(out), dimension(:,:) :: R_OUT
        real*8,  intent(in) :: sigma
        integer :: i_atom, j_dim, n_atom
        integer :: Natoms, DIM
        real*8  :: u

        Natoms = size(R_IN,dim=1); DIM = size(R_IN,dim=2) 

        if(.not. present(NatomsToDiffuse)) then 
            !move each atoms 
            do i_atom = 1, Natoms !for each atom and dimension    
                do j_dim = 1, DIM !gen position
                    R_OUT(i_atom,j_dim) = R_IN(i_atom,j_dim) + gauss(sigma)  
                end do
            end do    
        else
            R_OUT = R_IN !copy initial position to old one
            do n_atom = 1, NatomsToDiffuse !for each atom and dimension    
                call random_number(u)
                i_atom = floor(Natoms*u) + 1 !choose randomly one atom to diffuse
                do j_dim = 1, DIM !gen position
                    R_OUT(i_atom,j_dim) = R_OUT(i_atom,j_dim) + gauss(sigma)
                end do
            end do 
        end if 
    end subroutine 
    !###########################################################

    !####################################################
    !#              Compute Local Energy                #
    !#################################################### 
    ! Compute the local energy starting from configuration
    ! R and energy scale D=hsquare/2m
    real*8 function Elocal(R,DRIFT)
        implicit none
        real*8, intent(in), dimension(:,:) :: R
        real*8, intent(in), dimension(:,:),optional :: DRIFT
        if(.not. present(DRIFT)) then
            Elocal = Ekin(R) + Epot(R)
        else 
            Elocal = Ekin(R,DRIFT)
        end if 
         
        return 
    end function

    !####################################################
    !#           Compute Potential Energy               #
    !#################################################### 
    real*8 function Epot(R)
        implicit none 
        real*8,intent(in),dimension(:,:) :: R
        real*8  :: radius
        integer :: i_atom
        Epot = 0

        do i_atom=1,size(R,dim=1) 
            radius = norm2(R(i_atom,:)) 
            !potential energy in normalized quantities
            Epot = Epot + (params%Rcore**2 / params%a_osc**4) * radius**2
        end do 
        return 
    end function

    !####################################################
    !#             Compute Kinetic Energy               #
    !####################################################          
    real*8 function Ekin(R,DRIFT_IN)
        real*8,intent(in),dimension(:,:) :: R
        real*8, dimension(size(R,dim=1),size(R,2)),optional :: DRIFT_IN
        real*8, dimension(size(R,dim=1),size(R,2)) :: DRIFT
        real*8  :: radius,r1m,up,us,appo!,remainder,num,den
        integer :: i_atom,j_atom!, i_step
        integer :: Natoms
        ekin = 0
        Natoms = size(R,dim=1)

#ifndef TURNOFF_INTERACTION
        !pair interaction term 
        do i_atom=1,Natoms-1
            do j_atom=i_atom+1,Natoms
                radius  = norm2(R(i_atom,:) - R(j_atom,:))
                r1m     = 1./radius  
                up      = twobody_corrprime(radius)/twobody_corr(radius)
                us      = twobody_corrdoubleprime(radius)/twobody_corr(radius)
                
                !multiply by two to consider symmetric therms
                appo = - 2*( us + up * r1m +  - (up)**2)
                ekin = ekin + appo
            end do 
        end do 
#endif

        !potential interaction term 
        ekin = ekin - (-2/(params%a_osc**2))*Natoms
        
        !drift term
        if(.not. present(DRIFT_IN)) then
            DRIFT = F(R)
        else
            DRIFT = DRIFT_IN
        end if 
        DRIFT = F(R)
        do i_atom=1,Natoms
            ekin = ekin - 0.25d0*dot_product(DRIFT(i_atom,:),&
                                             DRIFT(i_atom,:))
        end do 
        ekin = D*Ekin
        return
    end function

    real*8 function Ekinfor(R,DRIFT)
        real*8,intent(in),dimension(:,:) :: R
        real*8, dimension(size(R,dim=1),size(R,dim=2)),optional :: DRIFT
        integer :: i_atom, Natoms
        ekinfor = 0
        Natoms = size(R,dim=1)

        if(.not. present(DRIFT)) DRIFT = F(R)
        do i_atom=1,Natoms
            ekinfor = ekinfor + 0.25*dot_product(DRIFT(i_atom,:),&
                                                 DRIFT(i_atom,:))
        end do 
        ekinfor = D*ekinfor
        return
    end function

    !####################################################
    !#                  Get Energies                    #
    !####################################################
    subroutine get_energies(R,OUT_Energies,DRIFT)
        implicit none
        
        real*8,intent(in),dimension(:,:) :: R
        type(Energies),intent(out) :: OUT_Energies
        real*8,optional,dimension(:,:) :: DRIFT !to speed up energy computing

        !not having to compute the drift speed up the process
        if( present(DRIFT)) then 
            OUT_Energies%Ekin    = Ekin(R,DRIFT)
            OUT_Energies%Ekinfor = Ekinfor(R,DRIFT)
        else 
            OUT_Energies%Ekin    = Ekin(R)
            OUT_Energies%Ekinfor = Ekinfor(R)
        end if 
        OUT_Energies%Epot = Epot(R)
        OUT_Energies%EL   = OUT_Energies%Epot +  OUT_Energies%Ekin
        

    end subroutine
    !###########################################################


    !####################################################
    !#             Compute Total Force                  #
    !#################################################### 
    function F(R) 
        implicit none
        real*8,intent(in) ,dimension(:,:) :: R
        real*8,dimension(size(R,dim=1),size(R,dim=2)) :: F
        real*8,dimension(size(R,dim=2))   :: r_hat !versor
        real*8  :: radius, up
        integer :: i_atom,j_atom,k_atom
        integer :: Natoms
        Natoms = size(R,dim=1)
        
        F = 0

#ifndef TURNOFF_INTERACTION
        !interaction force
        do i_atom = 1,Natoms-1
            do j_atom = i_atom+1,Natoms
                radius  = norm2(R(i_atom,:) - R(j_atom,:))       
                !normalized versor
                r_hat = (R(i_atom,:) - R(j_atom,:))/radius  
                up = twobody_corrprime(radius)/twobody_corr(radius)
                
                F(i_atom,:) = F(i_atom,:) + up*r_hat 
                F(j_atom,:) = F(j_atom,:) - up*r_hat 
            end do 
        end do 
#endif
        !field force
        do k_atom=1,Natoms
            radius = norm2(R(k_atom,:))
            !normalize versor
            r_hat = R(k_atom,:)/radius
            up    = - radius/(params%a_osc**2)

            F(k_atom,:) = F(k_atom,:) + up*r_hat    
        end do

        F = 2*F !from force defintion
        return 
    end function

end module
 
