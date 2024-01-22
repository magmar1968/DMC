module DMC_print
    use DMC_parameters
    implicit none

    integer , private :: FID = 11

    contains

    !####################################################
    !#               Print TWF to File                  #
    !####################################################
    subroutine print_TWF_tofile()
        use HS_puregas,Only:get_TWF_terms
        implicit none
        integer,parameter :: NFUNC = 3 !func,funcprime,funcdoubleprime
        real*8, dimension(NFUNC,TWFNPartitions) :: TWF_corr,TWF_harm
        real*8  :: rstep,r
        integer :: i_step
        character(MAX_FILENAME_LENGHT) :: filename
        filename = trim(outfile_path)//trim(twf_filename)
        
        call get_TWF_terms(TWF_corr=TWF_corr,TWF_harm=TWF_harm,rstep=rstep,&
                           Npartitions=TWFNPartitions)

        open(FID,file=filename,status="replace")
        write(FID,"(A)") "step,r,corr,corrprime,corrdoubleprime,harm,harmprime,harmdoubleprime"
        do i_step=1,TWFNPartitions
            r = (i_step-1)*rstep
            write(FID,*) i_step,",",r,",", &
                         TWF_corr(1,i_step),",",TWF_corr(2,i_step),",",&
                         TWF_corr(3,i_step),",",TWF_harm(1,i_step),",",&
                         TWF_harm(2,i_step),",",TWF_harm(3,i_step)
        end do 
        
        close(FID)
        
    end subroutine

    !####################################################
    !#         Print Energies Evolution to file         #
    !####################################################       
    subroutine print_E_evolution_toFile(MC_step,E_avg,Energ)
        use HS_puregas,Only:Energies
        implicit none
        integer,intent(in) :: MC_step
        real*8, intent(in) :: E_avg
        type(Energies) :: Energ
        logical,save ::FirstTime = .TRUE. 
        character(MAX_FILENAME_LENGHT),save :: filename
        
        if (FirstTime) then 
            filename = trim(outfile_path)//trim(energy_evolution_filename)
            open(FID,file=filename,status="replace")
            write(FID,"(A)") "step,Eavg,EL,EL^2,Epot,Ekin,Ekinfor"
            FirstTime = .FALSE. 
            close(FID)
        end if

        open(FID,file=filename,access="append") 
        write(FID,*) MC_step,     ",", &
                     E_avg,       ",", &
                     Energ%El,    ",", &
                     Energ%El**2, ",", &
                     Energ%Epot,  ",", &
                     Energ%Ekin,  ",", & 
                     Energ%Ekinfor
        close(FID)
    end subroutine

    !##################################################
    !#           Print Density Profile to File        #
    !##################################################
    subroutine print_density_profile_toFile(density_profile)
        implicit none
        real*8,dimension(NdensProfileSteps) :: density_profile
        integer :: i_step 
        real*8 :: rmin,rmax,area,dens
        character(MAX_FILENAME_LENGHT) :: filename 
        
        filename= trim(outfile_path)//trim(dens_profile_filename)

        open(FID,file=filename,status="replace")
        write(FID,"(A)")"rmin,rmax,dens"

        !normalize for number of evaluations
        density_profile = density_profile/(NMCsteps)
        !normalize for number of particles
        density_profile = density_profile/Natoms
        !normalize for area and printing 
        do i_step = 1, NdensProfileSteps
            rmin = (i_step-1)*densProfileStep
            rmax = (i_step)  *densProfileStep
            area = PI*(rmax**2-rmin**2)
            dens = density_profile(i_step)/area

            write(FID,*) rmin,",",rmax,",",dens
        end do 
        close(FID)
    end subroutine


    !####################################################
    !#           Print Configuration to File            #
    !####################################################
    subroutine print_configuration_toFile(R,filename)
        use DMC_parameters
        implicit none
        real*8,dimension(Natoms,DIM),intent(in) :: R 
        character(*) :: filename
        integer :: i_atom

        open(FID,file=filename,status="replace")
        write(FID,"(A)") "atom,x,y"
        do i_atom=1, Natoms
            write(FID,*) i_atom,",",R(i_atom,1),",",R(i_atom,2)
        end do
        close(FID)
    end subroutine


    ! subroutine print_inital_conf_toFile(R)
    !     use DMC_parameters
    !     implicit none 
    !     real*8,intent(in),dimension(Natoms,DIM) :: R
    !     character(MAX_FILENAME_LENGHT) :: filename 
    !     filename = trim(outfile_path)//trim(init_conf_filename)
    !     call print_configuration_toFile(R,filename)
    ! end subroutine print_inital_conf_toFIle

    subroutine print_final_conf_toFile(R)
        use DMC_parameters
        implicit none
        real*8,intent(in),dimension(Natoms,DIM) :: R
        character(MAX_FILENAME_LENGHT) :: filename
        filename = trim(outfile_path)//trim(final_conf_filename)
        call print_configuration_toFile(R,filename)
    end subroutine
    !####################################################

end module 

    ! !##################################################
    ! !#           Print Atoms Path to File             #
    ! !##################################################
    ! subroutine print_atoms_path_toFile(R,step)
    !     implicit none 
    !     real*8,intent(in),dimension(Natoms,DIM) :: R
    !     integer,intent(in) :: step
    !     integer :: i_atom
    !     logical,save :: FirstTime = .TRUE.
    !     character(MAX_FILENAME_LENGHT),save :: filename 
    !     character(len=1024) :: string

    !     if (FirstTime) then 
    !         filename = trim(outfile_path)//trim(atoms_path_filename)
    !         open(FID,file=filename,status="replace")
    !         write(FID,"(A)",advance="no") "step,"
    !         do i_atom = 1, Natoms
    !             if(i_atom < 10) then 
    !                 write(string,"(A1,I1,A2,I1,A1)") "x",i_atom,",y",i_atom,","
    !             else if (i_atom >= 10 .and. i_atom < 100) then 
    !                 write(string,"(A1,I2,A2,I2,A1)") "x",i_atom,",y",i_atom,","
    !             else if (i_atom >= 100 .and. i_atom < 1000) then
    !                 write(string,"(A1,I3,A2,I3,A1)") "x",i_atom,",y",i_atom,","
    !             else 
    !                 write(string,"(A1,I3,A2,I3,A1)") "x",i_atom,",y",i_atom,","
    !             end if 

    !             write(FID,"(A)",advance="no") trim(string)
    !         end do 
    !         FirstTime = .FALSE. 
    !         write(FID,*) !escaping 
    !         close(FID)
    !     end if

    !     open(FID,file=filename,access="append")
    !     write(FID,"(I5,A1)",advance="no") step,","
    !     do i_atom = 1, Natoms
    !         write(FID,"(f15.8,A1,f15.8,A1)",advance="no") R(i_atom,1),",",R(i_atom,2),","
    !     end do 
    !     write(FID,*)
    !     close(FID)
    ! end subroutine


