module WoodwardColellaModule
    use NumberKinds
    implicit none
    
    !--------------------------------------------------
    !variables to control the simulation
    !--------------------------------------------------
    real(KREAL), parameter :: CFL = 0.5_KREAL !CFL number
    real(KREAL) :: simTime = 0.0_KREAL !current simulation time
    real(KREAL), parameter :: MAX_TIME = 1.0_KREAL !maximum simulation time
    integer(KINT), parameter :: MAX_ITER = 10_KINT
    integer(KINT) :: iter = 1_KINT !number of iteration
    real(KREAL) :: dt !global time step
    
    !output control
    character(len=19), parameter :: HSTFILENAME = "WoodwardColella.hst" !history file name
    character(len=19), parameter :: RSTFILENAME = "WoodwardColella.dat" !result file name
    integer,parameter :: HSTFILE = 20 !history file ID
    integer,parameter :: RSTFILE = 21 !result file ID

    !gas
    integer(KINT), parameter :: CK = 3 !internal degree of freedom, here 3 denotes monoatomic gas
    real(KREAL) :: gamma !ratio of specific heat

    !geometry
    real(KREAL), parameter :: START_POINT = 0.0_KREAL, MID_POINT_1 = 0.1_KREAL, MID_POINT_2 = 0.9_KREAL, END_POINT = 1.0_KREAL
    real(KREAL) :: xSpacing
    integer(KINT), parameter :: POINTS_NUM = 400_KINT
    integer(KINT), parameter :: IXMIN = 1 , IXMAX = POINTS_NUM !cell index range

    !initial value
    real(KREAL), parameter :: VARS_LEFT(3) = (/1.0_KREAL,0.0_KREAL,2500.0_KREAL/) !left: density, x-momentum,total energy
    real(KREAL), parameter :: VARS_MID(3) = (/1.0_KREAL,0.0_KREAL,0.025_KREAL/) !middle: density, x-momentum,total energy
    real(KREAL), parameter :: VARS_RIGHT(3) = (/1.0_KREAL,0.0_KREAL,250.0_KREAL/) !right: density, x-momentum,total energy
end module WoodwardColellaModule