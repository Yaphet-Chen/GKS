module StationaryShockModule
    use NumberKinds
    use Mesh
    implicit none
    
    !--------------------------------------------------
    !variables to control the simulation
    !--------------------------------------------------
    real(KREAL), parameter              :: CFL = 0.65_KREAL !CFL number
    real(KREAL)                         :: simTime = 0.0_KREAL !current simulation time
    real(KREAL), parameter              :: MAX_TIME = 2.0_KREAL !maximum simulation time
    integer(KINT), parameter            :: MAX_ITER = 1000
    integer(KINT)                       :: iter = 1_KINT !number of iteration
    real(KREAL)                         :: dt !global time step

    !GKS control parameter
    real(KREAL), parameter              :: ETA = 0.5_KREAL

    !output control
    character(len=19), parameter        :: HSTFILENAME = "StationaryShock.hst" !history file name
    character(len=19), parameter        :: RSTFILENAME = "StationaryShock.dat" !result file name
    integer, parameter                  :: HSTFILE = 20 !history file ID
    integer, parameter                  :: RSTFILE = 21 !result file ID

    !gas
    integer(KINT), parameter            :: CK = 4 !internal degree of freedom, here 4 denotes diaatomic gas
    real(KREAL)                         :: gamma !ratio of specific heat

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
    real(KREAL), parameter              :: VARS_LEFT(3) = (/2.0_KREAL/3.0_KREAL,1.0_KREAL/sqrt(2.0_KREAL),9.0_KREAL/14.0_KREAL/) !left: density, x-momentum,total energy
    real(KREAL), parameter              :: VARS_RIGHT(3) = (/2.0_KREAL,1.0_KREAL/sqrt(2.0_KREAL),23.0_KREAL/14.0_KREAL/)!right: density, x-momentum,total energy
end module StationaryShockModule