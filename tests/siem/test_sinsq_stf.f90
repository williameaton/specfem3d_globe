program test_sinsq_stf

    use specfem_par
    use specfem_par_movie
    use manager_adios
    use constants, only: PI 

    implicit none
    integer :: ier,i 
    character(len=MAX_STRING_LEN) :: junk
    double precision :: time, stf, absval, realval, lochdur 

    ! initialize the MPI communicator and start the NPROCTOT MPI processes.
    call init_mpi()
    call world_rank(myrank)

    if (myrank == 0) print *,'program: test_sinsq_stf'

    ! initializes simulation parameters
    call initialize_simulation()
  
    ! sets up reference element GLL points/weights/derivatives
    call setup_GLL_points()
  
    ! starts reading the databases
    call read_mesh_databases()
  
    ! reads topography & bathymetry & ellipticity
    call read_topography_bathymetry()
  
    ! prepares sources and receivers
    call setup_sources_receivers()
  
    ! Add in test here: 

    if (myrank == 0)then 


        open(unit = IIN,file = trim(OUTPUT_FILES)//'/plot_source_time_function.txt', status = 'old', iostat=ier, form='formatted')
        if (ier /= 0 ) call exit_mpi(0,'Error opening plot_source_time_function file')

        read(IIN, '(A)')junk
        read(IIN, '(A)')junk
        read(IIN, '(A)')junk
        read(IIN, '(A)')junk

        do i = 1, NSTEP
            read(IIN,*)time, stf, absval

            ! Assuming HDUR is 20 seconds 
            lochdur = 20.0d0
            if (time < lochdur) then 
                realval =  0.5d0 + ((1/(2.0d0*PI*lochdur))*( (lochdur * sin(PI*time/lochdur)) + PI*time))
            else
                realval = 1.0d0
            endif 

            if (abs(realval - stf) > 1.e-6) then
                print *,'ERROR: stf value is incorrect'
                stop 1
            endif

        enddo
    endif 

    ! done
    if (myrank == 0) print *,'test_sinsq_stf done successfully'

    ! stop all the MPI processes, and exit
    call finalize_mpi()
end program test_sinsq_stf