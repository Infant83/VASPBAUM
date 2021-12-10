#include "alias.inc"

subroutine parse_very_init(PINPT)
    use parameters, only: incar
    use mpi_setup
    use mykind
    use utils, only: flag_number
    implicit none
    type(incar)         ::  PINPT
    integer(kind=sp)        narg, iarg, i, idummy
    character(len=256)      option, value, dummy
    
    PINPT%fnamelog  = 'VASPBAUM.out' ! default
    PINPT%title     = ''
    PINPT%folder_out= './'

    narg = iargc()
    do iarg = 1, narg
      call getarg(iarg, option)
      if(.not. flag_number(trim(option))) then
        if(trim(option) .eq. '-log' .or. trim(option) .eq. '-o') then
          call getarg(iarg+1, PINPT%fnamelog) ! set output file name

        elseif(trim(option) .eq. '-path') then
            call getarg(iarg+1, PINPT%folder_out)
            idummy = len_trim(PINPT%folder_out)
            if(PINPT%folder_out(idummy:idummy) .ne. '/') then
                PINPT%folder_out = trim(PINPT%folder_out)//'/'
            endif
            
        endif
      endif
    enddo

    write(dummy,'(A)')trim(PINPT%folder_out)//trim(PINPT%fnamelog)
    PINPT%fnamelog = trim(dummy)

    return
endsubroutine

subroutine parse(PINPT)
    use parameters, only: incar
    use mpi_setup
    use mykind
    use print_io
    use utils, only: flag_number
    implicit none
    type(incar  )       ::  PINPT
    integer(kind=sp)        narg, iarg, i
    integer(kind=sp)        iverbose_
    character(len=256)      option, value
    character(len=256)      dummy
    
    narg = iargc()

    write(PINPT%filenm,'(A)') trim(PINPT%folder_out)//'WAVECAR'
   !PINPT%filenm            = 'WAVECAR'    ! default
    PINPT%title             = 'VASPBAUM.OUT'
    PINPT%nelect            = 0            ! default
    PINPT%ie_init           = 1
    PINPT%ie_fina           = 999999
    PINPT%ispinor           = 2            ! default
    PINPT%icd               = 0
    PINPT%init_e            = -999.d0
    PINPT%fina_e            = -999.d0
    PINPT%sigma             = 0.01d0
    PINPT%theta             = 0.0d0
    PINPT%phi               = 0.0d0
    PINPT%nediv             = 1000
    PINPT%e_fermi           = 0.0d0
    PINPT%nkdiv             = 0 ! should be asigned via -nkdiv
    PINPT%flag_norm         = .FALSE.
    PINPT%flag_unfold       = .FALSE.
    PINPT%flag_set_unfold   = .FALSE.
    PINPT%flag_reduce       = .TRUE.


    iverbose     = 1 ! 1:full, 2:no
    print_mode   = 3 ! default verbosity
    
    do iarg = 1, narg
        call getarg(iarg, option)

        if(.not. flag_number(trim(option))) then
            if(trim(option) .eq. '-h') then
                if_main call help()

            elseif(trim(option) .eq. '-v') then
                call getarg(iarg+1, value)
                read(value, *) iverbose_
            
            elseif(trim(option) .eq. '-wf' .or. trim(option) .eq. '-f' ) then
                call getarg(iarg+1, PINPT%filenm)

            elseif(trim(option) .eq. '-s') then
                call getarg(iarg+1, value)
                read(value, *) PINPT%ispinor

            elseif(trim(option) .eq. '-soc') then
                PINPT%ispinor = 2
            
            elseif(trim(option) .eq. '-nosoc') then
                PINPT%ispinor = 1

            elseif(trim(option) .eq. '-t') then
                call getarg(iarg+1, PINPT%title)

            elseif(trim(option) .eq. '-norm') then
                call getarg(iarg+1, value)
                read(value,*) PINPT%flag_norm

            elseif(trim(option) .eq. '-nkdiv') then
                call getarg(iarg+1, value)
                read(value, *) PINPT%nkdiv

            elseif(trim(option) .eq. '-ne') then
                call getarg(iarg+1, value)
                read(value, *) PINPT%nelect

            elseif(trim(option) .eq. '-ii') then
                call getarg(iarg+1, value)
                read(value, *) PINPT%ie_init

            elseif(trim(option) .eq. '-if') then
                call getarg(iarg+1, value)
                read(value, *) PINPT%ie_fina
            
            elseif(trim(option) .eq. '-is') then
                call getarg(iarg+1, value)
                read(value, *) PINPT%ie_init
                PINPT%ie_fina = PINPT%ie_init

            elseif(trim(option) .eq. '-cd') then
                call getarg(iarg+1, value)
                read(value, *) PINPT%icd

            elseif(trim(option) .eq. '-unfold') then
                PINPT%flag_unfold = .TRUE.
            
            elseif(trim(option) .eq. '-set_unfold')then
                PINPT%flag_set_unfold = .TRUE.

            elseif(trim(option) .eq. '-no_reduce') then
                PINPT%flag_reduce = .FALSE.

            elseif(trim(option) .eq. '-ef') then
                call getarg(iarg+1, value)
                read(value, *) PINPT%e_fermi

            else if(trim(option) .eq. "-ien") then
                call getarg(iarg+1, value)
                read(value,*) PINPT%init_e
            else if(trim(option) .eq. "-fen") then
                call getarg(iarg+1, value)
                read(value,*) PINPT%fina_e
            else if(trim(option) .eq. "-theta") then
                call getarg(iarg+1, value)
                read(value,*) PINPT%theta
            else if(trim(option) .eq. "-phi") then
                call getarg(iarg+1, value)
                read(value,*) PINPT%phi
            else if(trim(option) .eq. "-nediv") then
                call getarg(iarg+1, value)
                read(value,*) PINPT%nediv
            else if(trim(option) .eq. "-sigma") then
                call getarg(iarg+1, value)
                read(value,*) PINPT%sigma

            endif
        endif

    enddo
    
!   if(PINPT%icd .ge. 2 .and. PINPT%nkdiv .eq. 0) then
!       write(6,'(A)')' !! ERROR : NKDIV = 0. Please provide NKDIV via -nkdiv option'
!       write(6,'(A)')'            Exit...'
!       kill_job
!   endif
    
    iverbose = iverbose_

    return
endsubroutine

subroutine help()
    use mpi_setup
    use mykind
    implicit none    

    write(6,'(A)')" "
    write(6,'(A)')" "
    write(6,'(A)')"          **** PROGRAM INSTRUCTION ***"
    write(6,'(A)')" "
    write(6,'(A)')" "
    write(6,'(A)')"*NOTE1 The x,y,z components of each G value are given"
    write(6,'(A)')"       in terms of the ig values and the components "
    write(6,'(A)')"       of the recip. lattice vectors according to:"
    write(6,'(A)')"        ig1*b1_x + ig2*b2_x + ig3*b3_x,"
    write(6,'(A)')"        ig1*b1_y + ig2*b2_y + ig3*b3_y, and"
    write(6,'(A)')"        ig1*b1_z + ig2*b2_z + ig3*b3_z, respectively, with"
    write(6,'(A)')"        ig1=ig(1,iplane),"
    write(6,'(A)')"        ig2=ig(2,iplane), and"
    write(6,'(A)')"        ig3=ig(3,iplane),"
    write(6,'(A)')"       where iplane=1,2,...,nplane(k,ispin) is an "
    write(6,'(A)')"       index incrementing the plane waves for specific"
    write(6,'(A)')"       k and spin values"
    write(6,'(A)')" "
    write(6,'(A)')"*NOTE2 The energy eigenvalues are complex, as provided"
    write(6,'(A)')"       in the WAVECAR file, but the imaginary part is "
    write(6,'(A)')"       zero (at least for cases investigated thus far)"
    write(6,'(A)')" "
    write(6,'(A)')"             ### POSSIBLE OPTIONS ###"
    write(6,'(A)')" -wf filename     : File name of WAVECAR to be read. Default: ./WAVECAR"
    write(6,'(A)')" -path folder     : Folder where the input/output will be read/written. Default: ./"
    write(6,'(A)')" -s   2 or 1      : for the noncollinear case, -s   2"
    write(6,'(A)')"                  : for the collinear or NM,   -s   1"
    write(6,'(A)')"                  :  Default : 2 if ISPIN 1, 1 if ISPIN 2"
    write(6,'(A)')" -soc             : equivalent to '-s 2' "
    write(6,'(A)')" -nosoc           : equivalent to '-s 1' "
    write(6,'(A)')" -t               : Title for the system   "
    write(6,'(A)')"                  : The output will be written 'title'.dat in general"
    write(6,'(A)')"                  :  Default : 'VASPBAUM'   "
    write(6,'(A)')" -norm T/F        : Whether normalize wavefunction.  "
    write(6,'(A)')"                  :  Note: It is not orthogonal due to the PAW approach"
    write(6,'(A)')" -nkdiv  nkdiv    : Number of k division between each KPATH"
    write(6,'(A)')" -ne  ne          : Specify total number of electrons. "
    write(6,'(A)')"                  : If not specified, total number of electrons will be"
    write(6,'(A)')"                  : estimated by adding up occupations at each k-point."
    write(6,'(A)')"                  : However, for semimetallic system, ne is differ for"
    write(6,'(A)')"                  : each k-points, so need to specify explicitly."
    write(6,'(A)')" -ii(if) ni(nf)   : Specify eigenvalue index ranges from ni-th to nf-th states"
    write(6,'(A)')"                  : unless '-cd' is not 1."
    write(6,'(A)')"                  :  Default -cd 0 -> -ii  1  -if VBM "
    write(6,'(A)')"                  :          -cd 1 -> -ii VBM -if CBM"
    write(6,'(A)')" -is n            : if specified, 'ni' and 'nf' will be set equal, "
    write(6,'(A)')"                  : and only this single band will be computed."
    write(6,'(A)')" -unfold          : Unfold band structure."
    write(6,'(A)')"                  : Usage: Before using this tag, you should have following files"
    write(6,'(A)')"                  :  1. Prepare POSCAR_PC  : primitive cell POSCAR to be projected"
    write(6,'(A)')"                  :  2. Prepare POSCAR_SC  : supercell      POSCAR to be unfolded "
    write(6,'(A)')"                  :  3. Prepare KPOINTS_PC : k-path of primitive BZ to be projected"
    write(6,'(A)')"                  : Note: before run with this tag, you should have WAVECAR which is"
    write(6,'(A)')"                  :       calculated with KPOINTS_SC (copied to KPOINTS for the actual run)"
    write(6,'(A)')"                  :       and generated by -set_unfold option."
    write(6,'(A)')" -ef e_fermi      : Fermi level, used in unfolding procedures. Energy shift by -e_fermi"
    write(6,'(A)')" -set_unfold      : Prepare KPOINTS for the unfolding"
    write(6,'(A)')"                  : Usage: Before using this tag, you should have following files"
    write(6,'(A)')"                  :  1. Prepare POSCAR_PC  : primitive cell POSCAR to be projected"
    write(6,'(A)')"                  :  2. Prepare POSCAR_SC  : supercell      POSCAR to be unfolded "
    write(6,'(A)')"                  :  3. Prepare KPOINTS_PC : k-path of primitive BZ to be projected"
    write(6,'(A)')"                  : Output: KPOINTS_SC file will be generated"
    write(6,'(A)')"                  :         --> copy KPOINTS_SC to KPOINTS for the calculation to get WAVECAR"
    write(6,'(A)')" -no_reduce       : do not remove duplicated K-point when generateing supercell KPOINTS with -set_unfold"
    write(6,'(A)')" -cd  1(or 0)     : Calculate spin- and k-resolved degree of circular polarization,"
    write(6,'(A)')"                  : between valence & conduction band"
    write(6,'(A)')"                  : If set to '1', Berry cuvature will not be evaluated."
    write(6,'(A)')"                  : You may set -ii and -if together, defining VBM & CBM, respectively."
    write(6,'(A)')"      2           : If -cd is set to 2, total spectrum w.r.t. the enery will be"
    write(6,'(A)')"                  : evaluated. In this case, -ien and -fen tag should be set by"
    write(6,'(A)')"                  : hand. Otherwise, default value, -ien 0 -fen 10 -nediv 1000 will"
    write(6,'(A)')"                  : be applied. The output contains degree of circular polarization"
    write(6,'(A)')"                  : with respect to photon energy for each k-point. Here, -sigma tag"
    write(6,'(A)')"                  : sets a gaussian broadning factor"
    write(6,'(A)')"      3           : If -cd is set to 3, total spectrum w.r.t. the energy will be"
    write(6,'(A)')"                  : evaluated. Same functionality with -cd 2, but instead calculate"
    write(6,'(A)')"                  : in order of kpoints list provided in lklist file with -klist option"
    write(6,'(A)')"                  : The spectral weight with unfolded band structure file (sw_file)"
    write(6,'(A)')"                  : also should be be provided with option -sw."
    write(6,'(A)')"                  :  Default : 0"
    write(6,'(A)')"  -klist kfile    : kfile lists kpoint index along k-path which will be read in sequence"
    write(6,'(A)')"  -sw sw_file     : Use sw_file generated by VaspBandUnfolding for the spectral weight"
    write(6,'(A)')"  -ien init_e     : -ien(fen) indicates initial(final) energy window to be plotted"
    write(6,'(A)')"  -fen fina_e     :              "
    write(6,'(A)')"                     NOTE: If the option '-unfold' is specified and -ien and -fen is not"
    write(6,'(A)')"                           specified, init(fina)_e will be automatically set based on "
    write(6,'(A)')"                           the energy band (min,max)"
    write(6,'(A)')"  -nediv ndiv     : How may energy grid will be divided for the spectral weight plot"
    write(6,'(A)')"  -sigma sigma    : Gaussian smearing (in eV)"
    write(6,'(A)')"  -atlist afile   : 'afile' contains information for atom_projected band structure file name."
    write(6,'(A)')"                  : Total number of atom to be highlighted, and atom indices."
    write(6,'(A)')"                  : Note that it works with -cd 2 and 3 only."
    write(6,'(A)')"                  : EX)  in afile, it read as follows"
    write(6,'(A)')"                  :   ----------------------------------------"
    write(6,'(A)')"                  :   | DOS_atom_projected.dat # LDOS file name to be read "
    write(6,'(A)')"                  :   | 6       # total atom in system   "
    write(6,'(A)')"                  :   | 3       # total atom to read     "
    write(6,'(A)')"                  :   | 1 2 5   # atom indices to read"
    write(6,'(A)')" *** Angle resolved CD calculations: *****************"
    write(6,'(A)')" * -theta theta   : angle along x-axis describing the direction of the injecting light"
    write(6,'(A)')" * -phi   phi     : angle along z-axis describing the direction of the injecting light"
    write(6,'(A)')" *                   Default: (theta,phi) = (0,0)"
    write(6,'(A)')" *                  Example:"
    write(6,'(A)')" *                   Light from z-axis(surface normal): (theta,phi) = (0.0,0.0)"
    write(6,'(A)')" *                   Light from x-axis                : (theta,phi) = (90.0,0.0)"
    write(6,'(A)')" *                   Light from y-axis                : (theta,phi) = (90.0,90.0)"
    write(6,'(A)')" *****************************************************"

    kill_job
    return
endsubroutine
