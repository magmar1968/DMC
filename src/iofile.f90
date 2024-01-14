!#########################################################
!#########################################################

! This file contain the subroutines to read input files with
! the specific structure:
! "dataname"  data
! 
! To do that the string module is used 
!
! File contents:

! Subroutine                           Label               Line    
! ----------                           ---                 ----
!
! Open Inputfile                       io_open
!
! READ-DATA
!   - Read Real Data                   read_data_r
!   - Read Integer Data                read_data_i
!
!################################################################
! Last Updated 7/11/23 L.Magnoni lorenzo.magnoni@studenti.unimi.it


module iofile
    Use, intrinsic :: iso_fortran_env, Only : iostat_end
    implicit none
   

    integer,private  :: FID = 10    !file ID
    integer,private           :: IERROR
    character(len=100),private :: ioerrmsg
    
    integer,parameter,private :: max_filename_lenght = 50
    integer,parameter,private :: max_dataname_lenght = 100
    integer,parameter,private :: max_n_args = 5, max_args_len = 30
    integer,parameter,private :: max_infile_lines = 100

    character(max_dataname_lenght),private :: line
    character(len=max_args_len),private,  DIMENSION(max_n_args) :: args
    character(len=max_dataname_lenght),private, dimension(max_infile_lines) :: param_names
    character(max_filename_lenght),private  :: filename
    character,private         :: delims = ' '
    
    integer,private           :: nargs,ios
        
    
    public :: io_open,is_present
    private :: read_data_r,read_data_i
    interface read_data
        module procedure read_data_r
        module procedure read_data_i
        module procedure read_data_l
        module procedure read_data_s
    end interface
    
!###################################################################
    contains

    !####################################################
    !                    Open Inputfile
    !####################################################
    subroutine io_open(input_filename,N_lines)
        
        implicit none
        integer,optional :: N_lines
        character(*), intent(in) :: input_filename
        filename = input_filename
        open(FID,file=input_filename,status="old",action= "read")
        !count all the lines
        IERROR = 0
        N_lines = 0
        do while (IERROR == 0)
            read(FID,*,iostat=IERROR)
            Select Case(IERROR)
            Case(0)
                N_lines = N_lines + 1
            Case(iostat_end) 
                exit
            Case Default
                print *, "IERROR: ", IERROR
                print *, ioerrmsg
            End Select
        end do  
        close(FID)

        call read_parameters_names()
    end subroutine io_open

    !###################################################
    !#           Read Parameters Names                 #
    !###################################################
    subroutine read_parameters_names()
        use strings, only:parse,value,compact
        Use, intrinsic :: iso_fortran_env, Only : iostat_end
        implicit none
        integer :: param_cont

        param_cont = 0
        
        open(FID,file=filename,status='old',action='read')
        IERROR = 0
        do while(IERROR == 0)
            read(FID,"(A)",iostat= IERROR,iomsg=ioerrmsg) line

            Select Case(IERROR)
            Case(0)
                call compact(line)
                call parse(line,delims,args=args,nargs=nargs)
                param_cont = param_cont + 1
                param_names(param_cont) = args(1)
            Case(iostat_end)
                exit
            Case Default 
                print *, "IERROR: ", IERROR 
                print *, ioerrmsg
            End Select
        end do 
        close(FID)
    end subroutine

    !####################################################
    !#            dataname is Present                   #
    !####################################################
    logical function is_present(dataname)
        implicit none
        character(*),intent(in) :: dataname    
        is_present = any(param_names .eq. dataname)
        return 
    end function 
    


    !####################################################
    !              Read Real Data
    !####################################################
    subroutine read_data_r(dataname, outdata)
        use strings, only:parse,value,compact
        Use, intrinsic :: iso_fortran_env, Only : iostat_end
        implicit none
        
        character(*), intent(in) :: dataname
        real*8, intent(out)      :: outdata

        open(FID,file=filename,status="old",action= "read")
        IERROR = 0
        do while(IERROR == 0)
            read(FID,"(A)",iostat= IERROR,iomsg=ioerrmsg) line

            Select Case(IERROR)
            Case(0)
                call compact(line)
                call parse(line,delims,args=args,nargs=nargs)

                if (args(1) == dataname) then
                    call value(args(2),outdata,ios)
                    exit
                end if 
            Case(iostat_end)
                print *, dataname, " not founded"
            Case Default 
                print *, "IERROR: ", IERROR 
                print *, ioerrmsg
            End Select
        end do 
        close(FID)
    end subroutine

    !####################################################
    !                 Read Integer Data
    !####################################################
    subroutine read_data_i(dataname, outdata)
        use strings, only:parse,value,compact
        Use, intrinsic :: iso_fortran_env, Only : iostat_end
        implicit none
        
        character(*), intent(in) :: dataname
        integer, intent(out)     :: outdata  

        open(FID,file=filename,status="old",action= "read")
        IERROR = 0
        do while(IERROR == 0)
            read(FID,'(A)',iostat= IERROR,iomsg=ioerrmsg) line
            
            Select Case(IERROR)
            Case(0)
                call compact(line)
                call parse(line,delims,args=args,nargs=nargs)

                if (args(1) == dataname) then
                    call value(args(2),outdata,ios)
                    exit
                end if 
            Case(iostat_end)
                print *, dataname, " not founded"
            Case Default 
                print *, "IERROR: ", IERROR 
                print *, ioerrmsg
            End Select 
        end do 

        close(FID)
    end subroutine

    !####################################################
    !#                 Read Logical Data                #
    !####################################################
    subroutine read_data_l(dataname, outdata)
        use strings, only:parse,value,compact
        Use, intrinsic :: iso_fortran_env, Only : iostat_end
        implicit none
        
        character(*), intent(in) :: dataname
        logical,     intent(out) :: outdata
        character,parameter      :: TRUE = "T"
        character,parameter      :: FALSE = "F" 


        open(FID,file=filename,status="old",action= "read")
        IERROR = 0
        do while(IERROR == 0)
            read(FID,'(A)',iostat= IERROR,iomsg=ioerrmsg) line

            Select Case(IERROR)
            Case(0)
                call compact(line)
                call parse(line,delims,args=args,nargs=nargs)

                if (args(1) == dataname) then
                    if(args(2) == TRUE) then 
                        outdata = .TRUE.
                    else if (args(2) == FALSE) then
                        outdata = .FALSE. 
                    else 
                        print *, "IERROR: impossible to read ",args(1)
                        outdata = .FALSE.
                    end if 
                    exit 
                end if  
            Case(iostat_end)
                print *, dataname, " not founded"
                print *, IERROR
            Case Default 
                print *, "IERROR: ", IERROR 
                print *, ioerrmsg
            End Select 
        end do 

        close(FID)
    end subroutine

    !####################################################
    !#                 Read String Data                #
    !####################################################
    subroutine read_data_s(dataname, outdata)
        use strings, only:parse,value,compact
        Use, intrinsic :: iso_fortran_env, Only : iostat_end
        implicit none
        
        character(*), intent(in) :: dataname
        character(*),  intent(out) :: outdata

        open(FID,file=filename,status="old",action= "read")
        IERROR = 0
        do while(IERROR == 0)
            read(FID,'(A)',iostat= IERROR,iomsg=ioerrmsg) line

            Select Case(IERROR)
            Case(0)
                call compact(line)
                call parse(line,delims,args=args,nargs=nargs)

                if (args(1) == dataname) then
                    outdata = args(2) 
                    exit 
                end if  
            Case(iostat_end)
                print *, dataname, " not founded"
                print *, IERROR
            Case Default 
                print *, "IERROR: ", IERROR 
                print *, ioerrmsg
            End Select 
        end do 

        close(FID)
    end subroutine
end module

