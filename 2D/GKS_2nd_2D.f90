module NumberKinds
    implicit none
    !--------------------------------------------------
    !kind selection
    !--------------------------------------------------
    integer, parameter                                  :: KREAL = kind(0.d0)
    integer, parameter                                  :: KINT = kind(1)
end module NumberKinds

!--------------------------------------------------
!>store the global constant variables
!--------------------------------------------------
module ConstantVariables
    use NumberKinds
    implicit none

    !mathematical constants
    real(KREAL), parameter                              :: PI = 4.0_KREAL*atan(1.0_KREAL) !Pi
    real(KREAL), parameter                              :: HALF = 0.5_KREAL !quantity 1/2
    real(KREAL), parameter                              :: QUATER = 0.25_KREAL !quantity 1/4

    !auxiliary constants
    real(KREAL), parameter                              :: TINYNUM = tiny(0.0_KREAL) !small value to avoid 0/0
    real(KREAL), parameter                              :: UP = 1.0_KREAL !used in sign() function

    !direction
    integer(KINT), parameter                            :: IDIRC = 1 !i direction
    integer(KINT), parameter                            :: JDIRC = 2 !j direction
    
    !rotation
    integer(KINT), parameter                            :: RN = 1 !no frame rotation
    integer(KINT), parameter                            :: RY = -1 !with frame rotation
end module ConstantVariables

!--------------------------------------------------
!basic derived type
!--------------------------------------------------
module Mesh
    use NumberKinds
    implicit none

    !cell center
    type IntraCell
        !geometry
        real(KREAL)                                     :: x, y !cell center coordinates
        real(KREAL)                                     :: area !cell area
        real(KREAL)                                     :: length(2) !length in i and j direction
        !flow field
        real(KREAL), dimension(4)                       :: conVars !density, x-momentum, y-momentum, total energy at cell center
        real(KREAL), dimension(4,2,2)                   :: faceVars
        !>faceVars represent the density, x-momentum, y-momentum, total energy at four sides of the cell
        !>face(:,m,n) 
        !>@ m: represent the selection of left/right at x-direction, up/down at y-direction.
        !>@ n: represent the selection of direction, IDIRC or JDIRC
    end type IntraCell

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
    !cell interface
    type CellInterface
        !geometry
        real(KREAL)                                     :: length !length of cell interface
        real(KREAL)                                     :: cosx,cosy !directional cosine
        !flow flux
        real(KREAL), dimension(4)                       :: flux(4) !mass flux, x momentum flux, y-momentum flux, total energy flux
    end type CellInterface

    !grid geometry(node coordinates)
    type Grid
        real(KREAL)                                     :: x, y !coordinates
    end type Grid
end module Mesh

!--------------------------------------------------
!>define some commonly used functions/subroutines
!--------------------------------------------------
module Tools
    use ConstantVariables
    implicit none
contains
    !--------------------------------------------------
    !>convert conservative variables to primary variables
    !>@param[in] w           :conservative variables
    !>@return    GetPrimary  :primitive variables
    !--------------------------------------------------
    function GetPrimary(w)
        real(KREAL), intent(in)                         :: w(4) !rho, x-momentum, y-momentum, total energy
        real(KREAL)                                     :: GetPrimary(4) !rho, x-velocity, y-velocity, pressure

        GetPrimary(1) = w(1)
        GetPrimary(2) = w(2)/w(1)
        GetPrimary(3) = w(3)/w(1)
        GetPrimary(4) = HALF*w(1)/GetLambda(w)
    end function GetPrimary

    function GetLambda(w)
        real(KREAL), intent(in)                         :: w(4) !conservative variables
        real(KREAL)                                     :: GetLambda !primitive variables

        GetLambda = (CK+2)*QUATER*w(1)/(w(4)-HALF*(w(2)**2+w(3)**2)/w(1))
    end function GetLambda

    !--------------------------------------------------
    !>convert primary variables to conservative variables
    !>@param[in] prim          :primary variables
    !>@return    GetConserved :conservative variables
    !--------------------------------------------------
    function GetConserved(prim)
        real(KREAL), intent(in)                         :: prim(4) !rho, x-velocity, y-velocity, pressure
        real(KREAL)                                     :: GetConserved(4) !rho, x-momentum, y-momentum, total energy

        GetConserved(1) = prim(1)
        GetConserved(2) = prim(1)*prim(2)
        GetConserved(3) = prim(1)*prim(3)
        GetConserved(4) = prim(4)/(GAMMA-1.0_KREAL)+HALF*prim(1)*(prim(2)**2+prim(3)**2)
    end function GetConserved

    !--------------------------------------------------
    !>obtain speed of sound
    !>@param[in] prim    :primary variables
    !>@return    GetSoundSpeed :speed of sound
    !--------------------------------------------------
    function GetSoundSpeed(prim)
        real(KREAL), intent(in)                         :: prim(4)
        real(KREAL)                                     :: GetSoundSpeed !speed of sound

        GetSoundSpeed = sqrt(GAMMA*prim(4)/prim(1))
    end function GetSoundSpeed

    !--------------------------------------------------
    !>calculate collision time
    !--------------------------------------------------
    function GetTau(preVal,midVal,nextVal) result(tau)
        real(KREAL), dimension(4), intent(in)           :: midVal,preVal,nextVal !primary variables: rho, U, V, lambda
        real(KREAL)                                     :: tau
        real(KREAL)                                     :: C1,C2
        if (TAU_TYPE .eq. "NS") then
            tau = MU/(HALF*midVal(1)/midVal(4))+C2*abs(preVal(1)/preVal(4)-nextVal(1)/nextVal(4))/(abs(preVal(1)/preVal(4)+nextVal(1)/nextVal(4))+TINYNUM)*dt
        elseif (TAU_TYPE .eq. "Euler") then
            tau = C1*dt+C2*abs(preVal(1)/preVal(4)-nextVal(1)/nextVal(4))/(abs(preVal(1)/preVal(4)+nextVal(1)/nextVal(4))+TINYNUM)*dt
        else
            stop "Error TAU_TYPE !"
        endif
    end function GetTau

    ! subroutine VanleerLimiter(midCell,leftCell,rightCell,downCell,upCell)
    !     type(IntraCell), intent(in)                     :: upCell, downCell, leftCell, rightCell
    !     type(IntraCell), intent(inout)                  :: midCell
    !     real(KREAL), dimension(4)                       :: splus, sminus

    !     !x direction
    !     splus = (rightCell%conVars-midCell%conVars)/(HALF*(rightCell%length(1)+midCell%length(1)))
    !     sminus = (midCell%conVars-leftCell%conVars)/(HALF*(midCell%length(1)+leftCell%length(1)))

    !     midCell%leftVars = midCell%conVars-HALF*midCell%length(1)*(sign(UP,splus)+sign(UP,sminus))*abs(splus)*abs(sminus)/(abs(splus)+abs(sminus)+TINYNUM)
    !     midCell%rightVars = midCell%conVars+HALF*midCell%length(1)*(sign(UP,splus)+sign(UP,sminus))*abs(splus)*abs(sminus)/(abs(splus)+abs(sminus)+TINYNUM)
        
    !     !y direction
    !     splus = (upCell%conVars-midCell%conVars)/(HALF*(upCell%length(2)+midCell%length(2)))
    !     sminus = (midCell%conVars-downCell%conVars)/(HALF*(midCell%length(2)+downCell%length(2)))

    !     midCell%downVars = midCell%conVars-HALF*midCell%length(2)*(sign(UP,splus)+sign(UP,sminus))*abs(splus)*abs(sminus)/(abs(splus)+abs(sminus)+TINYNUM)
    !     midCell%upVars = midCell%conVars+HALF*midCell%length(2)*(sign(UP,splus)+sign(UP,sminus))*abs(splus)*abs(sminus)/(abs(splus)+abs(sminus)+TINYNUM)

    !     !check whether the reconstruction failed
    !     if ( GetLambda(midCell%leftVars)<=TINYNUM .or. GetLambda(midCell%rightVars)<=TINYNUM &
    !         .or. GetLambda(midCell%downVars)<=TINYNUM .or.GetLambda(midCell%upVars)<=TINYNUM ) then
    !         midCell%leftVars = midCell%conVars
    !         midCell%rightVars = midCell%conVars
    !         midCell%downVars = midCell%conVar
    !         midCell%upVars = midCell%conVar
    !     endif
    ! end subroutine VanleerLimiter

    !--------------------------------------------------
    !>get heat flux
    !>@param[in] h,b           :distribution function
    !>@param[in] vn,vt         :normal and tangential velocity
    !>@param[in] prim          :primary variables
    !>@return    get_heat_flux :heat flux in normal and tangential direction
    !--------------------------------------------------
    ! function get_heat_flux(h,b,vn,vt,prim)
    !     real(kind=RKD),dimension(:,:),intent(in) :: h,b
    !     real(kind=RKD),dimension(:,:),intent(in) :: vn,vt
    !     real(kind=RKD),intent(in) :: prim(4)
    !     real(kind=RKD) :: get_heat_flux(2) !heat flux in normal and tangential direction

    !     get_heat_flux(1) = 0.5*(sum(weight*(vn-prim(2))*((vn-prim(2))**2+(vt-prim(3))**2)*h)+sum(weight*(vn-prim(2))*b)) 
    !     get_heat_flux(2) = 0.5*(sum(weight*(vt-prim(3))*((vn-prim(2))**2+(vt-prim(3))**2)*h)+sum(weight*(vt-prim(3))*b)) 
    ! end function get_heat_flux
end module Tools

!--------------------------------------------------
!>flux calculation
!--------------------------------------------------
module Flux
    use Tools
    implicit none
contains

end module Flux

module Solver
    use Tools
    implicit none
contains
    !--------------------------------------------------
    !>calculate time step
    !--------------------------------------------------
    subroutine TimeStep()
        real(KREAL)                                     :: prim(4) !primary variables
        real(KREAL)                                     :: sos !speed of sound
        real(KREAL)                                     :: tmax !max 1/dt allowed
        integer                                         :: i, j

        !set initial value
        tmax = 0.0_KREAL

        do j=IYMIN,IYMAX
            do i=IXMIN,IXMAX
                !convert conservative variables to primary variables
                prim = GetPrimary(ctr(i,j)%conVars)

                !sound speed
                sos = GetSoundSpeed(prim)

                !maximum velocity in x and y direction
                prim(2) = abs(prim(2))+sos
                prim(3) = abs(prim(3))+sos

                !time step
                tmax = max(tmax,(ctr(i,j)%length(2)*prim(2)+ctr(i,j)%length(1)*prim(3)+TINYNUM)/ctr(i,j)%area)
            enddo
        enddo

        !time step
        dt = CFL/tmax
    end subroutine TimeStep

    subroutine IntraCellReconstruction()
        integer(KINT)                                   :: i,j
        do i=IXMIN-1,IXMAX+1
            call VanleerLimiter(ctr(i-1),ctr(i),ctr(i+1))
            ! ctr(i)%leftVars = ctr(i)%conVars; ctr(i)%rightVars = ctr(i)%conVars !first order reconstruction, no limiter
        enddo 
    end subroutine IntraCellReconstruction

    subroutine Evolution()
        integer(KINT)                                   :: i,j

        do j=IYMIN,IYMAX
            do i=IXMIN,IXMAX+1
                call CalcFlux(ctr(i-1,j),vface(i,j),ctr(i,j),IDIRC)
            end do
        end do

        do j=IYMIN,IYMAX+1
            do i=IXMIN,IXMAX
                call CalcFlux(ctr(i,j-1),hface(i,j),ctr(i,j),JDIRC)
            end do
        end do
    end subroutine Evolution

    subroutine Update()
        integer(KINT)                                   :: i,j

        do j=IYMIN,IYMAX
            do i=IXMIN,IXMAX
                !--------------------------------------------------
                !update conservative variables
                !--------------------------------------------------
                ctr(i,j)%conVars = ctr(i,j)%conVars &
                                    +(vface(i,j)%flux*vface(i,j)%length-vface(i+1,j)%flux*vface(i+1,j)%length &
                                    +hface(i,j)%flux*hface(i,j)%length-hface(i,j+1)%flux*hface(i,j+1)%length)/ctr(i,j)%area
            end do
        end do
        end subroutine Update
end module Solver
!--------------------------------------------------
!>main program
!--------------------------------------------------
program GKS_2nd_2D
    use Solver
    use Writer
    implicit none

    !open file and write header
    open(unit=HSTFILE,file=HSTFILENAME,status="replace",action="write") !open history file
    write(HSTFILE,*) "VARIABLES = iter, simTime, dt" !write header

    !iteration
    do while(.true.)
        call TimeStep() !calculate global time step
        call Boundary() !set up Boundary condition
        call IntraCellReconstruction() !reconstruction
        call Evolution() !calculate flux across the interfaces
        call Update() !Update cell averaged value

        !check if output
        if (simTime>=MAX_TIME .or. iter>=MAX_ITER) exit

        !write iteration situation every 10 iterations
        if (mod(iter,10)==0) then
            write(*,"(A18,I15,2E15.7)") "iter,simTime,dt:",iter,simTime,dt
            write(HSTFILE,"(I15,2E15.7)") iter,simTime,dt
        end if

        iter = iter+1
        simTime = simTime+dt
    enddo

    !close history file
    close(HSTFILE)

    !output solution
    call output()

end program GKS_2nd_2D