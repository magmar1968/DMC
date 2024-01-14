module array_utility
    
    private :: increase_size_di,increase_size_si
    private :: increase_size_dr,increase_size_sr

    interface increase_size
        module procedure increase_size_dr
        module procedure increase_size_sr
        module procedure increase_size_di
        module procedure increase_size_si
    end interface 

    contains 

    subroutine increase_size_dr(array,new_size)
        real*8, dimension(:),allocatable,intent(inout) :: array
        real*8, dimension(:),allocatable :: tmp
        integer,intent(in) :: new_size
        integer :: old_size

        old_size = size(array)

        allocate(tmp(new_size))
        tmp = 0.
        !copy array
        tmp(1:old_size) = array(1:old_size)
        !swap arrays 
        call move_alloc(tmp,array)
    end subroutine

    subroutine increase_size_sr(array,new_size)
        real*4, dimension(:),allocatable,intent(inout) :: array
        real*4, dimension(:),allocatable :: tmp
        integer,intent(in) :: new_size
        integer :: old_size

        old_size = size(array)

        allocate(tmp(new_size))
        tmp = 0.
        !copy array
        tmp(1:old_size) = array(1:old_size)
        !swap arrays 
        call move_alloc(tmp,array)
    end subroutine
    
    subroutine increase_size_di(array,new_size)
        integer*8, dimension(:),allocatable,intent(inout) :: array
        integer*8, dimension(:),allocatable :: tmp
        integer,intent(in) :: new_size
        integer :: old_size

        old_size = size(array)

        allocate(tmp(new_size))
        tmp = 0.
        !copy array
        tmp(1:old_size) = array(1:old_size)
        !swap arrays 
        call move_alloc(tmp,array)
    end subroutine

    subroutine increase_size_si(array,new_size)
        integer*4, dimension(:),allocatable,intent(inout) :: array
        integer*4, dimension(:),allocatable :: tmp
        integer,intent(in) :: new_size
        integer :: old_size

        old_size = size(array)

        allocate(tmp(new_size))
        tmp = 0.
        !copy array
        tmp(1:old_size) = array(1:old_size)
        !swap arrays 
        call move_alloc(tmp,array)
    end subroutine 
end module