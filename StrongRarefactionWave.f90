module StrongRarefactionWaveModule
    use NumberKinds
    use Mesh
    implicit none
    
    !--------------------------------------------------
    !variables to control the simulation
    !--------------------------------------------------
    real(KREAL), parameter              :: CFL = 0.65_KREAL !CFL number
    real(KREAL)                         :: simTime = 0.0_KREAL !current simulation time
    real(KREAL), parameter              :: MAX_TIME = 10.0_KREAL !maximum simulation time
    integer(KINT), parameter            :: MAX_ITER = 1000
    integer(KINT)                       :: iter = 1_KINT !number of iteration
    real(KREAL)                         :: dt !global time step

    !GKS control parameter
    real(KREAL), parameter              :: ETA = 0.0_KREAL

    !output control
    character(len=25), parameter        :: HSTFILENAME = "StrongRarefactionWave.hst" !history file name
    character(len=25), parameter        :: RSTFILENAME = "StrongRarefactionWave.dat" !result file name
    integer, parameter                  :: HSTFILE = 20 !history file ID
    integer, parameter                  :: RSTFILE = 21 !result file ID

    !gas
    integer(KINT), parameter            :: CK = 4 !internal degree of freedom, here 4 denotes diatomic gas
    real(KREAL), parameter              :: gamma = real(k+3,KREAL)/real(k+1,KREAL) !ratio of specific heat

    !geometry
    real(KREAL), parameter              :: START_POINT = 0.0_KREAL, MID_POINT = 100.0_KREAL, END_POINT = 200.0_KREAL
    real(KREAL)                         :: xSpacing
    integer(KINT), parameter            :: POINTS_NUM = 200_KINT
    integer(KINT), parameter            :: IXMIN = 1 , IXMAX = POINTS_NUM !cell index range

    !--------------------------------------------------
    !flow field
    !--------------------------------------------------
    !index method
    !     ----------------
    !  (i)|      (i)     |(i+1)
    !     ----------------
    type(CellCenter)                    :: ctr(IXMIN-1:IXMAX+1) !cell center (with ghost cell)
    type(CellInterface)                 :: vface(IXMIN:IXMAX+1) !vertical cell interfaces

    !initial value
    real(KREAL), parameter              :: VARS_LEFT(3) = (/1.0_KREAL,-5.0_KREAL,13.5_KREAL/) !left: density, x-momentum,total energy
    real(KREAL), parameter              :: VARS_RIGHT(3) = (/1.0_KREAL,5.0_KREAL,13.5_KREAL/)!right: density, x-momentum,total energy
end module StrongRarefactionWaveModule