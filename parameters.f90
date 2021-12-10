#include "alias.inc"
module parameters
    use mykind
    character(len=26),public, parameter   ::alphabet='abcdefghijklmnopqrstuvwxyz'
    real(kind=dp)   , public, parameter   ::      pi=4.d0*atan(1.d0) ! Leibniz's formula for Pi
    real(kind=dp)   , public, parameter   ::     pi2=2.d0*pi
    real(kind=dp)   , public, parameter   ::    bohr=0.52917721067d0 ! meter Bohr radius
    real(kind=dp)   , public, parameter   :: hartree=27.21138602d0 ! eV   Hartree energy
    real(kind=dp)   , public, parameter   ::    hbar=4.135667662d-15/pi2  ! eV*s Plank constant
    real(kind=dp)   , public, parameter   ::       c=0.262465831d0 !! constant c = 2m/hbar**2 [1/eV Ang^2] 
    real(kind=dp)   , public, parameter   ::     k_B=8.6173303d-5  !! Boltzmann constant = 2m/hbar**2 [1/eV Ang^2] 
    real(kind=dp)   , public, parameter   :: g_elect=2.0023193043617 !! g-factor[see Mol.Phys. 98, 1597(2000) for sign]
    real(kind=dp)   , public, parameter   ::     rt2=sin( 4.d0*atan(1.d0)/4.d0 ) * 2.d0 ! sqrt(2)
    real(kind=dp)   , public, parameter   ::     rt3=sin( 4.d0*atan(1.d0)/3.d0 ) * 2.d0 ! sqrt(3)
    real(kind=dp)   , public, parameter   :: onsite_tolerance= 0.0001 !! symmetry precision
    real(kind=dp)   , public, parameter   ::     eta=epsilon(1d0) ! tiny value
    real(kind=dp)   , public              :: t1, t0 ! time
    complex(kind=dp), public, parameter   ::      zi=(0.d0,1.d0)
    complex(kind=dp), public, parameter   ::     pzi= pi*zi
    complex(kind=dp), public, parameter   ::    pzi2=2*pzi
    complex(kind=dp), public, dimension(2,2),   parameter :: pauli_0 = reshape((/  1,  0,  0,  1 /), (/2,2/))
    complex(kind=dp), public, dimension(2,2),   parameter :: pauli_x = reshape((/  0,  1,  1,  0 /), (/2,2/))
    complex(kind=dp), public, dimension(2,2),   parameter :: pauli_y = reshape((/  0,  1, -1,  0 /) * zi, (/2,2/))
    complex(kind=dp), public, dimension(2,2),   parameter :: pauli_z = reshape((/  1,  0,  0, -1 /), (/2,2/))
    integer(kind=sp), public, dimension(3,3,3), parameter :: levi_civita = reshape((/0,0,0, 0,0,-1, 0,1,0, &
                                                                       0,0,1, 0,0,0, -1,0,0, &
                                                                       0,-1,0, 1,0,0, 0,0,0/), (/3,3,3/))
    integer(kind=sp),  public, dimension(3,2),   parameter :: cyclic_axis = reshape((/2,3,1,3,1,2/), (/3,2/))

    integer(kind=sp),  public, parameter :: max_dummy       = 9999999 !! maximun number of dummy index 
    integer(kind=sp),  public, parameter :: pid_wavecar     = 10      !! wavecar id
    integer(kind=sp),  public, parameter :: pid_output      = 11      !! output id
    integer(kind=sp),  public, parameter :: pid_kpts        = 12      !! kpoint file id
    integer(kind=sp),  public, parameter :: pid_geom        = 13      !! poscar file id

    type incar !PINPT
        character(len=256)                  fnamelog            ! log file name, default = TBFIT.log
        character(len=132)                  title               ! title of the system
        character(len=256)                  filenm              ! wavecar file name
        character(len=256)                  folder_out          ! folder where the output file will be written
                                                                
        logical                             flag_cd             
        logical                             flag_bc             
        logical                             flag_unfold         ! band structure unfold
        logical                             flag_set_unfold     ! prepare unfolding calculation (KPOINTS_SC)
        logical                             flag_reduce         ! remove duplicated k-point?                                                           
        ! calculation mode                                      
        integer(kind=sp)                    icd                 ! CD calculation


        ! basic information in WAVECAR
        logical                             flag_norm           ! normalize wavefunction .TRUE. or .FALSE. Default: .FALSE.
        integer(kind=sp)                    irecl, iprec        ! record length, rtag
        real(kind=dp)                       a1(3), a2(3), a3(3) ! lattice vector
        real(kind=dp)                       b1(3), b2(3), b3(3) ! reciprocal lattice vector
        real(kind=dp)                       A(3,3)              ! a1, a2, a3
        real(kind=dp)                       B(3,3)              ! b1, b2, b3
        real(kind=dp)                       encut               ! NECUT  cutoff energy (eV)
                                                                ! ENCUT specifies the cutoff energy for the plane-wave-basis set in eV.
                                                                ! |G + k| < G_cut with ENCUT = G_cut^2 / c 
                                                                ! --> VASP convention:https://www.vasp.at/wiki/index.php/ENCUT)
        integer(kind=sp)                    ispin               ! ISPIN: 1 NM (and SOC), 2 SP(collinear)
        integer(kind=sp)                    nkpts               ! NKPTS 
        integer(kind=sp)                    nband               ! NBANDS number of eigenvalues
        integer(kind=sp)                    ispinor             ! ISPINOR: 1 NM and SP (collinear), 2 SOC
        integer(kind=sp)                    nelect              ! NELECT total number of electrons
        integer(kind=sp)                    nkdiv               ! kpoint division

        ! basic setup for calculations
        integer(kind=sp)                    ie_init, ie_fina    ! eigenvalue index for initial and final states to be calculated

        ! control parameters for CD(or unfold) calculations
        real(kind=dp)                       init_e, fina_e      ! energy window to be resolved
        integer(kind=sp)                    nediv               ! number of division of energy window
        real(kind=dp)                       theta, phi          ! incident circularly polarized light with angle (theta,phi)
                                                                ! in 'deg'
        real(kind=dp)                       sigma               ! gaussian smearing (in eV)
        real(kind=dp)                       e_fermi             ! fermi level to be shifted by -e_fermi (eV). Default: 0.0 eV
        
    endtype incar


    type eigen  !WAVEC
        real(kind=dp),      allocatable ::  kpts(:,:)           ! kpoints(3,nkpts)      ! reciprocal unit
        real(kind=dp),      allocatable ::  kpts_cart(:,:)      ! kpoints_cart(3,nkpts) ! cartesian  unit
        real(kind=dp),      allocatable ::  E(:,:,:)            ! energy(nband,ispin,nkpts)
        real(kind=dp),      allocatable ::  OCC(:,:,:)          ! occupation(nband,ispin,nkpts)
        integer(kind=sp),   allocatable ::  nplw(:)             ! nplw(nkpts) total number of plane waves for each kpoints. nplw = number of plane waves * ispinor
                                                                !  --> nplw = number of plane waves for each spinor basis * ispinor
        integer(kind=sp),   allocatable ::  nplw_screen(:)      ! nplw(nkpts) total number of plane waves screened by mapping SC BZ to PC BZ by M (should be integer)
                                                                !  --> nplw_screen = number of plane waves for each spinor
        integer(kind=sp),   allocatable ::  G(:,:,:)            ! (3,nplw_max,nkpts) G vectors that defines the plane wave for each k-point
        integer(kind=sp),   allocatable ::  Gid(:,:)            ! (  nplw_max,nkpts) G vector index (see Gid function)
                                                                ! -> each G vector will be directed by its unique id
        integer(kind=sp),   allocatable ::  GS(:,:,:)           ! (3,nplw_max,nkpts) screened G vectors of SC
        integer(kind=sp),   allocatable ::  GSid(:,:)           ! (  nplw_max,nkpts) screened G vector index of SC
        integer(kind=sp)                    nplw_max            ! maximum number of plane waves (* spinor)
        integer(kind=sp)                    plw_idx_max         ! maximum index number (Gid) for plane wave G vector indexing
        integer(kind=sp)                    nbmax(3)            ! maximum number of G vectors for each direction

        real(kind=dp),      allocatable ::  CD(:,:,:,:,:)       ! CD(3,is,nband,nband,nkpts) CD between (ie, je) at ik
                                                                ! CD(1,is,ie,je,ik) : RCD = <psi_j|A_R|phi_i>
                                                                ! CD(2,is,ie,je,ik) : LCD
                                                                ! CD(3,is,ie,je,ik) : (RCD - LCD)/(RCD+LCD)
        real(kind=dp),      allocatable ::  SW(:,:,:,:)         ! Spectral weight SW(3,is,nediv,nkpts)
                                                                ! SW(1,~~      ) : SW for RCD
                                                                ! SW(2,~~      ) : SW for LCD
                                                                ! SW(3,~~      ) : SW for  CD
                                                                ! Note: for ispinor case, one can utilize SW to 
                                                                !       compute and store magnetic CD with SW(i, ~) : sigma_i
                                                                !       where sigma_i is the Pauli matrices (i = x,y,z)
        real(kind=dp),      allocatable ::  SW_BAND(:,:,:)      ! Spectral weight for each band (nband,ispin,nkpts)
    endtype eigen       


    type poscar !PGEOM
        real(kind=dp)                       a_latt(3,3)
        real(kind=dp)                       b_latt(3,3)
        real(kind=dp)                       M(3,3)              ! transformation matrix (PGEOM_PC)

        ! if PGEOM_PC : Primitive cell info. G is not defined. kpts & kpts_cart is read from KPOINTS_BAND
        ! if PGEOM_SC : Supercell info. G is defined accordingly (see map_k_to_K_1BZ routine).
        !               kpts & kpts_cart is calculated automatically and saved.
        !               The corresponding supercell k-point index is obtained from WAVEC info.
        !#######################################################!
        integer(kind=sp),   allocatable ::  G(:,:)              ! (3,nkpts) shifting G vector w.r.t. PC
        real(kind=dp),      allocatable ::  kpts(:,:)           ! kpoints (if unfold is true) along kpath
                                                                ! k     if PGEOM_PC
                                                                ! K_1BZ if PGEOM_SC
        real(kind=dp),      allocatable ::  kpts_cart(:,:)      ! kpoints (cartesian) (if unfold is true)
        integer(kind=sp)                    nkpts               ! number of total kpoints for kpts, kpts_cart
                                                                ! nkpts = nkdiv * nline  if PGEOM_PC
                                                                !       = PINPT%nkpts    if PGEOM_SC
        integer(kind=sp),   allocatable ::  ikpt(:)             ! (PGEOM_PC%nkpts) index of K_1BZ that maps k into K_1BZ
        !#######################################################!

        integer(kind=sp)                    nkdiv               ! number of kpoint division between each high symmetry lines
        integer(kind=sp)                    nline               ! number of high symmetry lines
        real(kind=dp),      allocatable ::  kline(:,:)          ! list of high symmetry kpoints (3, nline*2)
        real(kind=dp),      allocatable ::  kline_cart(:,:)     ! list of high symmetry kpoints (cartesian)
        real(kind=dp),      allocatable ::  kdist(:)            !
        character(len=8),   allocatable ::  k_name(:)           ! name of high symmetry kpoints
        character(len=8),   allocatable ::  k_name2(:)          ! name of high symmetry kpoints (remove duplicate points)
        integer(kind=sp),   allocatable ::  k_name_index(:)     ! index of high symmetry kpoints
    endtype poscar

endmodule
