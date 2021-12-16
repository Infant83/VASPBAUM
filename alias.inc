#define if_main            if(myid .eq. 0)
#define if_main_then       if(myid .eq. 0) then
#define if_main_end        endif
#define if_nsp2            if(PINPT%nspin .eq. 2)
#define if_ncol            if(PINPT%ispinor .eq. 2)
#define mpi_kill           call MPI_Abort(mpi_comm_earth, 0, mpierr)
#define if_test            if(PINPT%flag_tbfit_test)
#define write_msg          call write_log(message,3,myid)
#define write_msg_file     call write_log(message,1,myid)
#define write_msg_screen   call write_log(message,2,myid)
#define write_msg_all      call write_log(message,13,myid)
#define write_msgi         call write_log(message,print_mode,myid)
#define write_msgi_file    call write_log(message,print_mode-2,myid)
#ifdef MPI
#define kill_job           call kill_job()
#define kill_job1          call kill_job1()
#else
#define kill_job           stop
#define kill_job1          stop
#endif
#define __VBAUM_VERSION__   "0.0.1"
