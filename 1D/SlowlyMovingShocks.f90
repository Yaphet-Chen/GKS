module SlowlyMovingShocksModule
    use NumberKinds
    use Mesh
    implicit none
    
    !--------------------------------------------------
    !variables to control the simulation
    !--------------------------------------------------
    real(KREAL), parameter              :: CFL = 0.65_KREAL !CFL number
    real(KREAL)                         :: simTime = 0.0_KREAL !current simulation time
    real(KREAL), parameter              :: MAX_TIME = 0.95_KREAL !maximum simulation time
    integer(KINT), parameter            :: MAX_ITER = 3000
    integer(KINT)                       :: iter = 1_KINT !number of iteration
    real(KREAL)                         :: dt !global time step

    !GKS control parameter
    real(KREAL), parameter              :: ETA = 0.5_KREAL

    !output control
    character(len=22), parameter        :: HSTFILENAME = "SlowlyMovingShocks.hst" !history file name
    character(len=22), parameter        :: RSTFILENAME = "SlowlyMovingShocks.dat" !result file name
    integer,parameter                   :: HSTFILE = 20 !history file ID
    integer,parameter                   :: RSTFILE = 21 !result file ID

    !gas
    integer(KINT), parameter            :: CK = 4 !internal degree of freedom, here 4 denotes diatomic gas
    real(KREAL), parameter              :: gamma = real(k+3,KREAL)/real(k+1,KREAL) !ratio of specific heat

    !geometry
    real(KREAL), parameter              :: START_POINT = 0.0_KREAL, MID_POINT = 0.5_KREAL, END_POINT = 1.0_KREAL
    real(KREAL)                         :: xSpacing
    integer(KINT), parameter            :: POINTS_NUM = 100_KINT
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
    real(KREAL), parameter              :: VARS_LEFT(3) = (/3.86_KREAL,-3.1266_KREAL,27.0913_KREAL/) !left: density, x-momentum,total energy
    real(KREAL), parameter              :: VARS_RIGHT(3) = (/1.0_KREAL,-3.44_KREAL,8.4168_KREAL/)!right: density, x-momentum,total energy
end module SlowlyMovingShocksModule