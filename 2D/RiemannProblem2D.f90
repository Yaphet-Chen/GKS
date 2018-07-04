module RiemannProblem2DModule
    use NumberKinds
    use Mesh
    implicit none
    
    !--------------------------------------------------
    !variables to control the simulation
    !--------------------------------------------------
    real(KREAL), parameter              :: CFL = 0.5_KREAL !CFL number
    real(KREAL)                         :: simTime = 0.0_KREAL !current simulation time
    real(KREAL), parameter              :: MAX_TIME = 0.4_KREAL !1.6 !maximum simulation time
    integer(KINT), parameter            :: MAX_ITER = 1000 !maximum iteration numbers
    integer(KINT)                       :: iter = 1_KINT !current number of iteration
    real(KREAL)                         :: dt !global time step
    character(len=10)                   :: TAU_TYPE = "Euler"
    real(KREAL), parameter              :: C1 = 0.05_KREAL, C2 = 5.0_KREAL !parameters that control tau calculation

    !output control
    character(len=20), parameter        :: HSTFILENAME = "RiemannProblem2D.hst" !history file name
    character(len=20), parameter        :: RSTFILENAME = "RiemannProblem2D.dat" !result file name
    integer, parameter                  :: HSTFILE = 20 !history file ID
    integer, parameter                  :: RSTFILE = 21 !result file ID

    !gas properties
    real(KREAL), parameter              :: MU = 0.0_KREAL !viscosity coefficient in reference state
    integer(KINT), parameter            :: CK = 3 !internal degree of freedom, here 3 denotes diatomic gas
    real(KREAL), parameter              :: GAMMA = real(CK+4,KREAL)/real(CK+2,KREAL) !ratio of specific heat
    real(KREAL), parameter              :: PR !Prandtl number

    !geometry
    type(Grid), parameter               :: UPPER_LEFT_POINT%x = 0.0_KREAL, UPPER_LEFT_POINT%y = 2.0_KREAL
    type(Grid), parameter               :: UPPER_RIGHT_POINT%x = 2.0_KREAL, UPPER_RIGHT_POINT%y = 2.0_KREAL
    type(Grid), parameter               :: LOWER_LEFT_POINT%x = 0.0_KREAL, LOWER_LEFT_POINT%y = 0.0_KREAL
    type(Grid), parameter               :: LOWER_RIGHT_POINT%x = 2.0_KREAL, LOWER_RIGHT_POINT%y = 0.0_KREAL
    integer(KINT), parameter            :: X_POINTS_NUM = 800_KINT, Y_POINTS_NUM = 800_KINT
    integer(KINT), parameter            :: X_MID_POINT_NUM = X_POINTS_NUM/2, Y_MID_POINT_NUM = Y_POINTS_NUM/2
    !--------------------------
    ! Ghost point at boundary
    ! low order suggest 2, weno5 suggest 3
    !          UP
    !       ********
    !      |        |
    ! LEFT |        | RIGHT
    !      |        |
    !       ********
    !         DOWN
    !--------------------------
    integer(KINT), parameter            :: GHOST_NUM = 2_KINT
    integer(KINT), parameter            :: IXMIN = 1, IXMAX = X_POINTS_NUM !cell index range in x direction
    integer(KINT), parameter            :: IYMIN = 1, IYMAX = Y_POINTS_NUM !cell index range in y direction
    
    !--------------------------------------------------
    !index method
    !              hface(i,j+1)
    !            ----------------
    !            |              |
    !            |              |
    !            |              |
    ! vface(i,j) |     (i,j)    | vface(i+1,j)
    !            |              |
    !            |              |
    !            |              |
    !            ----------------
    !              hface(i,j)
    !--------------------------------------------------
    type(IntraCell)                     :: ctr(IXMIN-GHOST_NUM:IXMAX+GHOST_NUM,IYMIN-GHOST_NUM:IYMAX+GHOST_NUM) !cell center (with ghost cell)
    type(CellInterface)                 :: vface(IXMIN-GHOST_NUM+1:IXMAX+GHOST_NUM,IYMIN-GHOST_NUM:IYMAX+GHOST_NUM) !vertical interface
    type(CellInterface)                 :: hface(IXMIN-GHOST_NUM:IXMAX+GHOST_NUM,IYMIN-GHOST_NUM+1:IYMAX+GHOST_NUM) !horizontal interface

    !initial value
    real(KREAL), parameter              :: PRIMVARS_LOWER_LEFT(4) = (/1.0_KREAL,-0.75_KREAL,0.5_KREAL,1.0_KREAL/) ! x<1, y<1 zone: density, x-velocity, y-velocity, pressure
    real(KREAL), parameter              :: PRIMVARS_LOWER_RIGHT(4) = (/3.0_KREAL,-0.75_KREAL,-0.5_KREAL,1.0_KREAL/) ! x>=1, y<1 zone: density, x-velocity, y-velocity, pressure
    real(KREAL), parameter              :: PRIMVARS_UPPER_RIGHT(4) = (/1.0_KREAL,0.75_KREAL,-0.5_KREAL,1.0_KREAL/) ! x>=1, y>=1 zone: density, x-velocity, y-velocity, pressure
    real(KREAL), parameter              :: PRIMVARS_UPPER_LEFT(4) = (/2.0_KREAL,0.75_KREAL,0.5_KREAL,1.0_KREAL/) ! x<1, y>=1 zone: density, x-velocity, y-velocity, pressure

    !Initialize initial condition
    ctr(IXMIN-GHOST_NUM:IXMIN+X_MID_POINT_NUM-1,IYMIN-GHOST_NUM:IYMIN+Y_MID_POINT_NUM-1)%conVars = GetConvars(PRIMVARS_LOWER_LEFT)
    ctr(IXMIN+X_MID_POINT_NUM:IXMAX+GHOST_NUM,IYMIN-GHOST_NUM:IYMIN+Y_MID_POINT_NUM-1)%conVars = GetConvars(PRIMVARS_LOWER_RIGHT)
    ctr(IXMIN+X_MID_POINT_NUM:IXMAX+GHOST_NUM,IYMIN+Y_MID_POINT_NUM,IYMAX+GHOST_NUM)%conVars = GetConvars(PRIMVARS_UPPER_RIGHT)
    ctr(IXMIN-GHOST_NUM:IXMIN+X_MID_POINT_NUM-1,IYMIN+Y_MID_POINT_NUM,IYMAX+GHOST_NUM)%conVars = GetConvars(PRIMVARS_UPPER_LEFT)
end module RiemannProblem2DModule