module clock


    real,private :: T1,T2
    
    type time_res
        integer totsec
        integer seconds
        integer minutes
        integer hours
    end type 

    contains 

    subroutine start_clock()
        call cpu_time(T1)
    end subroutine

    subroutine stop_clock(delta_time) 
        type(time_res),intent(out) :: delta_time
        integer Elapsed_time
        call cpu_time(T2)
        Elapsed_time = int(T2- T1)
        delta_time%totsec  = Elapsed_time
        delta_time%hours   =  Elapsed_time/3600
        delta_time%minutes = mod(Elapsed_time,3600)/60
        delta_time%seconds = mod(Elapsed_time,60) 
    end subroutine

end module 

