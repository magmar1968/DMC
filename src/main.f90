
program main
    use DMC,Only:DMC_results,run_DMC

    implicit none
    
    type(DMC_results) :: results
    call input("./data/inputfile.dat")
    call run_DMC(results)
    print *, "######################################################"
    print *, "RESULTS:"
    print *, "  energy:    ", results%E
    print *, "  error:     ", results%error
    print *, "  cross rate:", results%cross_rate
    print *, "######################################################"


end program 