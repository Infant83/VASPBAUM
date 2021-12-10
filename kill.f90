subroutine kill_job()
    use mpi_setup
    use mykind
    implicit none

#ifdef MPI
        if(COMM_KOREA%flag_split) then
            call MPI_BARRIER(COMM_KOREA%mpi_comm, mpierr)
            call mpi_finish()
        elseif(.not. COMM_KOREA%flag_split) then
            call MPI_BARRIER(mpi_comm_earth, mpierr)
            call mpi_finish()
            stop
        endif
#else

        stop

#endif        

endsubroutine

subroutine kill_job1()
    use mpi_setup
    use mykind
    implicit none

#ifdef MPI
    call MPI_ABORT(mpi_comm_earth, mpierr)
#else
    stop
#endif

endsubroutine
