module SodShockTubeModule
    use NumberKinds
    use Mesh
    implicit none
    
    !--------------------------------------------------
    !variables to control the simulation
    !--------------------------------------------------
    real(KREAL), parameter              :: CFL = 0.5_KREAL !CFL number
    real(KREAL)                         :: simTime = 0.0_KREAL !current simulation time
    real(KREAL), parameter              :: MAX_TIME = 0.15_KREAL !maximum simulation time
    integer(KINT), parameter            :: MAX_ITER = 1000
    integer(KINT)                       :: iter = 1_KINT !number of iteration
    real(KREAL)                         :: dt !global time step

    !output control
    character(len=16), parameter        :: HSTFILENAME = "SodShockTube.hst" !history file name
    character(len=16), parameter        :: RSTFILENAME = "SodShockTube.dat" !result file name
    integer,parameter                   :: HSTFILE = 20 !history file ID
    integer,parameter                   :: RSTFILE = 21 !result file ID

    !gas
    real(KREAL), parameter              :: MU = 0.0_KREAL
    integer(KINT), parameter            :: CK = 4 !internal degree of freedom, here 4 denotes diatomic gas
    real(KREAL), parameter              :: gamma = real(CK+3,KREAL)/real(CK+1,KREAL) !ratio of specific heat
    ! real(KREAL), parameter              :: gamma = 1.4_KREAL !ratio of specific heat
    ! integer(KINT), parameter            :: CK = nint((3.0_KREAL-gamma)/(gamma-1.0_KREAL)) !internal degree of freedom, here 4 denotes diatomic gas

    !geometry
    real(KREAL), parameter              :: START_POINT = 0.0_KREAL, MID_POINT = 0.5_KREAL, END_POINT = 1.0_KREAL
    integer(KINT), parameter            :: POINTS_NUM = 100_KINT
    integer(KINT), parameter            :: GHOST_NUM = 2_KINT !low order suggest 2, weno5 suggest 3
    integer(KINT), parameter            :: IXMIN = 1 , IXMAX = POINTS_NUM !cell index range

    !--------------------------------------------------
    !flow field
    !--------------------------------------------------
    !index method
    ! --------------------
    !  (i) |  (i)  | (i+1)
    ! face | cell  | face
    ! --------------------
    type(IntraCell)                     :: ctr(IXMIN-GHOST_NUM:IXMAX+GHOST_NUM) !cell center (with ghost cell)
    type(CellInterface)                 :: vface(IXMIN-GHOST_NUM+1:IXMAX+GHOST_NUM) !vertical cell interfaces

    !initial value
    real(KREAL), parameter              :: CONVARS_LEFT(3) = (/1.0_KREAL,0.0_KREAL,2.5_KREAL/) !left: density, x-momentum,total energy
    real(KREAL), parameter              :: CONVARS_RIGHT(3) = (/0.125_KREAL,0.0_KREAL,0.25_KREAL/)!right: density, x-momentum,total energy
end module SodShockTubeModule

module SjogreenSupersonicExpansionModule
    use NumberKinds
    use Mesh
    implicit none
    
    !--------------------------------------------------
    !variables to control the simulation
    !--------------------------------------------------
    real(KREAL), parameter              :: CFL = 0.5_KREAL !CFL number
    real(KREAL)                         :: simTime = 0.0_KREAL !current simulation time
    real(KREAL), parameter              :: MAX_TIME = 0.1_KREAL !maximum simulation time
    integer(KINT), parameter            :: MAX_ITER = 100
    integer(KINT)                       :: iter = 1_KINT !number of iteration
    real(KREAL)                         :: dt !global time step
    
    !output control
    character(len=31), parameter        :: HSTFILENAME = "SjogreenSupersonicExpansion.hst" !history file name
    character(len=31), parameter        :: RSTFILENAME = "SjogreenSupersonicExpansion.dat" !result file name
    integer,parameter                   :: HSTFILE = 20 !history file ID
    integer,parameter                   :: RSTFILE = 21 !result file ID

    !gas
    real(KREAL), parameter              :: MU = 0.0_KREAL
    integer(KINT), parameter            :: CK = 4 !internal degree of freedom, here 4 denotes diatomic gas
    real(KREAL), parameter              :: gamma = real(CK+3,KREAL)/real(CK+1,KREAL) !ratio of specific heat
    ! real(KREAL), parameter              :: gamma = 1.4_KREAL !ratio of specific heat
    ! integer(KINT), parameter            :: CK = nint((3.0_KREAL-gamma)/(gamma-1.0_KREAL)) !internal degree of freedom, here 4 denotes diatomic gas

    !geometry
    real(KREAL), parameter              :: START_POINT = 0.0_KREAL, MID_POINT = 0.5_KREAL, END_POINT = 1.0_KREAL
    integer(KINT), parameter            :: POINTS_NUM = 100_KINT
    integer(KINT), parameter            :: GHOST_NUM = 2_KINT !low order suggest 2, weno5 suggest 3
    integer(KINT), parameter            :: IXMIN = 1 , IXMAX = POINTS_NUM !cell index range

    type(IntraCell)                     :: ctr(IXMIN-GHOST_NUM:IXMAX+GHOST_NUM) !cell center (with ghost cell)
    type(CellInterface)                 :: vface(IXMIN-GHOST_NUM+1:IXMAX+GHOST_NUM) !vertical cell interfaces

    !initial value
    real(KREAL), parameter              :: CONVARS_LEFT(3) = (/1.0_KREAL,-2.0_KREAL,3.0_KREAL/) !left: density, x-momentum,total energy
    real(KREAL), parameter              :: CONVARS_RIGHT(3) = (/1.0_KREAL,2.0_KREAL,3.0_KREAL/)!right: density, x-momentum,total energy
end module SjogreenSupersonicExpansionModule

module WoodwardColellaModule
    use NumberKinds
    use Mesh
    implicit none
    
    !--------------------------------------------------
    !variables to control the simulation
    !--------------------------------------------------
    real(KREAL), parameter              :: CFL = 0.5_KREAL !CFL number
    real(KREAL)                         :: simTime = 0.0_KREAL !current simulation time
    real(KREAL), parameter              :: MAX_TIME = 0.038_KREAL !maximum simulation time
    integer(KINT), parameter            :: MAX_ITER = 2000_KINT
    integer(KINT)                       :: iter = 1_KINT !number of iteration
    real(KREAL)                         :: dt !global time step

    !output control
    character(len=19), parameter        :: HSTFILENAME = "WoodwardColella.hst" !history file name
    character(len=19), parameter        :: RSTFILENAME = "WoodwardColella.dat" !result file name
    integer,parameter                   :: HSTFILE = 20 !history file ID
    integer,parameter                   :: RSTFILE = 21 !result file ID

    !gas
    real(KREAL), parameter              :: MU = 0.0_KREAL
    integer(KINT), parameter            :: CK = 4 !internal degree of freedom, here 4 denotes diatomic gas
    real(KREAL), parameter              :: gamma = real(CK+3,KREAL)/real(CK+1,KREAL) !ratio of specific heat
    ! real(KREAL), parameter              :: gamma = 1.4_KREAL !ratio of specific heat
    ! integer(KINT), parameter            :: CK = nint((3.0_KREAL-gamma)/(gamma-1.0_KREAL)) !internal degree of freedom, here 4 denotes diatomic gas

    !geometry
    real(KREAL), parameter              :: START_POINT = 0.0_KREAL, MID_POINT_1 = 0.1_KREAL, MID_POINT_2 = 0.9_KREAL, END_POINT = 1.0_KREAL
    integer(KINT), parameter            :: POINTS_NUM = 400_KINT
    integer(KINT), parameter            :: GHOST_NUM = 2_KINT !low order suggest 2, weno5 suggest 3
    integer(KINT), parameter            :: IXMIN = 1 , IXMAX = POINTS_NUM !cell index range

    type(IntraCell)                     :: ctr(IXMIN-GHOST_NUM:IXMAX+GHOST_NUM) !cell center (with ghost cell)
    type(CellInterface)                 :: vface(IXMIN-GHOST_NUM+1:IXMAX+GHOST_NUM) !vertical cell interfaces

    !initial value
    real(KREAL), parameter              :: CONVARS_LEFT(3) = (/1.0_KREAL,0.0_KREAL,2500.0_KREAL/) !left: density, x-momentum,total energy
    real(KREAL), parameter              :: CONVARS_MID(3) = (/1.0_KREAL,0.0_KREAL,0.025_KREAL/) !middle: density, x-momentum,total energy
    real(KREAL), parameter              :: CONVARS_RIGHT(3) = (/1.0_KREAL,0.0_KREAL,250.0_KREAL/) !right: density, x-momentum,total energy
end module WoodwardColellaModule

module SlowlyMovingShocksModule
    use NumberKinds
    use Mesh
    implicit none
    
    !--------------------------------------------------
    !variables to control the simulation
    !--------------------------------------------------
    real(KREAL), parameter              :: CFL = 0.5_KREAL !CFL number
    real(KREAL)                         :: simTime = 0.0_KREAL !current simulation time
    real(KREAL), parameter              :: MAX_TIME = 0.95_KREAL !maximum simulation time
    integer(KINT), parameter            :: MAX_ITER = 3000
    integer(KINT)                       :: iter = 1_KINT !number of iteration
    real(KREAL)                         :: dt !global time step

    !output control
    character(len=22), parameter        :: HSTFILENAME = "SlowlyMovingShocks.hst" !history file name
    character(len=22), parameter        :: RSTFILENAME = "SlowlyMovingShocks.dat" !result file name
    integer,parameter                   :: HSTFILE = 20 !history file ID
    integer,parameter                   :: RSTFILE = 21 !result file ID

    !gas
    real(KREAL), parameter              :: MU = 0.0_KREAL
    integer(KINT), parameter            :: CK = 4 !internal degree of freedom, here 4 denotes diatomic gas
    real(KREAL), parameter              :: gamma = real(CK+3,KREAL)/real(CK+1,KREAL) !ratio of specific heat
    ! real(KREAL), parameter              :: gamma = 1.4_KREAL !ratio of specific heat
    ! integer(KINT), parameter            :: CK = nint((3.0_KREAL-gamma)/(gamma-1.0_KREAL)) !internal degree of freedom, here 4 denotes diatomic gas

    !geometry
    real(KREAL), parameter              :: START_POINT = 0.0_KREAL, MID_POINT = 0.5_KREAL, END_POINT = 1.0_KREAL
    integer(KINT), parameter            :: POINTS_NUM = 100_KINT
    integer(KINT), parameter            :: GHOST_NUM = 2_KINT !low order suggest 2, weno5 suggest 3
    integer(KINT), parameter            :: IXMIN = 1 , IXMAX = POINTS_NUM !cell index range

    type(IntraCell)                     :: ctr(IXMIN-GHOST_NUM:IXMAX+GHOST_NUM) !cell center (with ghost cell)
    type(CellInterface)                 :: vface(IXMIN-GHOST_NUM+1:IXMAX+GHOST_NUM) !vertical cell interfaces

    !initial value
    real(KREAL), parameter              :: CONVARS_LEFT(3) = (/3.86_KREAL,-3.1266_KREAL,27.0913_KREAL/) !left: density, x-momentum,total energy
    real(KREAL), parameter              :: CONVARS_RIGHT(3) = (/1.0_KREAL,-3.44_KREAL,8.4168_KREAL/)!right: density, x-momentum,total energy
end module SlowlyMovingShocksModule

module StrongRarefactionWaveModule
    use NumberKinds
    use Mesh
    implicit none
    
    !--------------------------------------------------
    !variables to control the simulation
    !--------------------------------------------------
    real(KREAL), parameter              :: CFL = 0.5_KREAL !CFL number
    real(KREAL)                         :: simTime = 0.0_KREAL !current simulation time
    real(KREAL), parameter              :: MAX_TIME = 10.0_KREAL !maximum simulation time
    integer(KINT), parameter            :: MAX_ITER = 1000
    integer(KINT)                       :: iter = 1_KINT !number of iteration
    real(KREAL)                         :: dt !global time step

    !output control
    character(len=25), parameter        :: HSTFILENAME = "StrongRarefactionWave.hst" !history file name
    character(len=25), parameter        :: RSTFILENAME = "StrongRarefactionWave.dat" !result file name
    integer, parameter                  :: HSTFILE = 20 !history file ID
    integer, parameter                  :: RSTFILE = 21 !result file ID

    !gas
    real(KREAL), parameter              :: MU = 0.0_KREAL
    integer(KINT), parameter            :: CK = 4 !internal degree of freedom, here 4 denotes diatomic gas
    real(KREAL), parameter              :: gamma = real(CK+3,KREAL)/real(CK+1,KREAL) !ratio of specific heat
    ! real(KREAL), parameter              :: gamma = 1.4_KREAL !ratio of specific heat
    ! integer(KINT), parameter            :: CK = nint((3.0_KREAL-gamma)/(gamma-1.0_KREAL)) !internal degree of freedom, here 4 denotes diatomic gas

    !geometry
    real(KREAL), parameter              :: START_POINT = 0.0_KREAL, MID_POINT = 100.0_KREAL, END_POINT = 200.0_KREAL
    integer(KINT), parameter            :: POINTS_NUM = 200_KINT
    integer(KINT), parameter            :: GHOST_NUM = 2_KINT !low order suggest 2, weno5 suggest 3
    integer(KINT), parameter            :: IXMIN = 1 , IXMAX = POINTS_NUM !cell index range
    
    type(IntraCell)                     :: ctr(IXMIN-GHOST_NUM:IXMAX+GHOST_NUM) !cell center (with ghost cell)
    type(CellInterface)                 :: vface(IXMIN-GHOST_NUM+1:IXMAX+GHOST_NUM) !vertical cell interfaces

    !initial value
    real(KREAL), parameter              :: CONVARS_LEFT(3) = (/1.0_KREAL,-5.0_KREAL,13.5_KREAL/) !left: density, x-momentum,total energy
    real(KREAL), parameter              :: CONVARS_RIGHT(3) = (/1.0_KREAL,5.0_KREAL,13.5_KREAL/)!right: density, x-momentum,total energy
end module StrongRarefactionWaveModule
