module version
  use time
  use mpi_setup
  use print_io 

contains 

  function version_tag() result(ver_tag)
    implicit none
    character*132  ver_tag

    ver_tag='# VASP-BAUM version 0.0.1  (build ' // __DATE__// ' ' //__TIME__// ') '

    return
  endfunction

  subroutine version_stamp(t_start)
    implicit none
    real*8         t_start
    
    call write_log(trim( version_tag()          ),3,myid)
    call write_log(' ',3,myid)
    call write_log("#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!",3,myid)
    call write_log("#! Copyright 2018. Hyun-Jung Kim All rights reserved.!",3,myid)
    call write_log("#!           (FJZ,  h.kim@fz-juelich.de)             !",3,myid)
    call write_log("#!           (KIAS, Infant@kias.re.kr)               !",3,myid)
    call write_log("#!           (angpangmokjang@hanmail.net)            !",3,myid)
    call write_log("#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!",3,myid)
    call write_log(' ',3,myid)

    call write_log(" -------------------------------------------------------------------------",3,myid)
    call timestamp('| Program start on', t_start)
    call report_hostname()
    call write_log(" -------------------------------------------------------------------------",3,myid)

    return
  endsubroutine

endmodule
