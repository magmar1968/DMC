!author: Lorenzo Magnoni
!date:   28/09/23
!code: This module is intended to contain usefull functions for the 
!      generation of random number for all possible distribution
!
!#############################################################
! File contents:
!
! Subroutine                           Label               Line    
! ----------                           ---                 ----
! Initial Random Seed                  init_random_seed
!
! Uniform Random Number Generator      uniform
!
! Gaussian Random Number Generator     gauss         
!#######################################################################
!Last Update 07/11/2023 L.Magnoni lorenzo.magnoni@studenti.unimi.it

module random
    contains

    !###################################
    !       Initialize Random Seed
    !###################################
    subroutine init_random_seed()
        use iso_fortran_env, only: int64
        implicit none
        integer, allocatable :: seed(:)
        integer :: i, n, un, istat, dt(8), pid
        integer(int64) :: t
      
        call random_seed(size = n)
        allocate(seed(n))
        ! First try if the OS provides a random number generator
        open(newunit=un, file="/dev/urandom", access="stream", &
             form="unformatted", action="read", status="old", iostat=istat)
        if (istat == 0) then
           read(un) seed
           close(un)
        else
           ! Fallback to XOR:ing the current time and pid. The PID is
           ! useful in case one launches multiple instances of the same
           ! program in parallel.
           call system_clock(t)
           if (t == 0) then
              call date_and_time(values=dt)
              t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
                   + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
                   + dt(3) * 24_int64 * 60 * 60 * 1000 &
                   + dt(5) * 60 * 60 * 1000 &
                   + dt(6) * 60 * 1000 + dt(7) * 1000 &
                   + dt(8)
           end if
           pid = getpid()
           t = ieor(t, int(pid, kind(t)))
           do i = 1, n
              seed(i) = lcg(t)
           end do
        end if
        call random_seed(put=seed)
      contains
        ! This simple PRNG might not be good enough for real work, but is
        ! sufficient for seeding a better PRNG.
        function lcg(s)
          integer :: lcg
          integer(int64) :: s
          if (s == 0) then
             s = 104729
          else
             s = mod(s, 4294967296_int64)
          end if
          s = mod(s * 279470273_int64, 4294967291_int64)
          lcg = int(mod(s, int(huge(0), int64)), kind(0))
        end function lcg
      end subroutine init_random_seed
        

    !##############################################
    !     Uniform random number generator
    !##############################################
    !usage: generate a uniform distributed random number in [min, max] interval
    !input: - min (default 0) left limit of the interval
    !       - max (default 1) right limit of the inteval
    function uniform(min,max) result(w)
        implicit none

        real :: w                    !output 
        real, optional :: min, max   ! input (optional)
        real :: a_min,a_max ! a is short for "actual..."

        if (present(min)) then
            a_min = min
        else 
            a_min = 0
        end if 

        if (present(max)) then
            a_max = max
        else 
            a_max = 1
        end if

        call random_number(w)
        w = w*(a_max-a_min) + a_min

    end function uniform


    !##############################################
    !     Gaussian random number generator
    !##############################################
    !usage: generate a gaussian distributed random number using the BOX-MULLER method
    !input: - mu (default 0): the mu of the normal distribution
    !       - variance (default 1): the standar deviation of the normal distribution  
    double precision function gauss(sigma,mu)

        implicit none
        real*8, intent(in),optional :: sigma,mu
        real*8                  :: m_sigma,m_mu
        real*8                  :: w, v1, v2, l
        real*8                  :: s1, s2
        w = 2.d0

        if(present(sigma)) then 
            m_sigma = sigma
        else 
            m_sigma = 1. 
        end if 

        if(present(mu)) then
            m_mu = mu
        else 
            m_mu = 0.
        end if 

        do
            call random_number(s1)
            call random_number(s2)

            v1 = 2.d0*s1 - 1.d0
            v2 = 2.d0*s2 - 1.d0
            w = v1*v1 + v2*v2

            if (w.lt.1.d0) exit
        end do

        l     = v1*sqrt(-2.d0*log(w)/(w))
        l     = m_sigma*l + m_mu 
        gauss = l
    end function gauss
    !##############################################
end module