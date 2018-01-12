module NumberKinds
    implicit none

    integer, parameter                                  :: KREAL = kind(0.d0)
    integer, parameter                                  :: KINT = kind(1)
end module NumberKinds
!--------------------------------------------------
!>store the global constant variables
!--------------------------------------------------
module ConstantVariables
    use NumberKinds
    implicit none

    real(KREAL), parameter                              :: PI = 4.0_KREAL*atan(1.0_KREAL) !Pi
    real(KREAL), parameter                              :: TINYNUM = tiny(0.0_KREAL) !small value to avoid 0/0
    real(KREAL), parameter                              :: HALF = 0.5_KREAL !quantity 1/2
    real(KREAL), parameter                              :: QUATER = 0.25_KREAL !quantity 1/4  
end module ConstantVariables

module Mesh
    use NumberKinds
    implicit none

    !--------------------------------------------------
    !basic derived type
    !--------------------------------------------------
    !cell center
    type IntraCell
        !geometry
        real(KREAL)                                     :: x !cell center coordinates
        real(KREAL)                                     :: length !length
        !flow field
        real(KREAL)                                     :: conVars(3) !density, x-momentum, total energy at cell center
        real(KREAL)                                     :: leftVars(3) !density, x-momentum, total energy at cell left end point
        real(KREAL)                                     :: rightVars(3) !density, x-momentum, total energy at cell right end point
    end type IntraCell

    !cell interface
    type CellInterface
        real(KREAL)                                     :: flux(3) !mass flux, x momentum flux, total energy flux
    end type CellInterface

    !index method
    ! --------------------
    !  (i) |  (i)  | (i+1)
    ! face | cell  | face
    ! --------------------
end module Mesh

module Initialization
    use ConstantVariables
    use SodShockTubeModule
    implicit none

contains
    !--------------------------------------------------
    !>main initialization subroutine
    !--------------------------------------------------
    subroutine Init()  
        call InitUniformMesh(START_POINT,END_POINT,POINTS_NUM) !initialize the geometry
        call InitDataOneMidPoint(CONVARS_LEFT,CONVARS_RIGHT,MID_POINT) !set the initial condition with one mid point
        ! call InitDataTwoMidPoints(CONVARS_LEFT,CONVARS_MID,CONVARS_RIGHT,MID_POINT_1,MID_POINT_2) !set the initial condition with two mid points
    end subroutine Init

    subroutine InitUniformMesh(left,right,num)
        real(KREAL), intent(in)                         :: left, right
        integer(KINT), intent(in)                       :: num
        real(KREAL)                                     :: dx
        integer(KINT)                                   :: i

        !cell length
        dx = (right-left)/(IXMAX-IXMIN+1)

        !cell center (with ghost cell)
        forall(i=IXMIN-GHOST_NUM:IXMAX+GHOST_NUM) 
            ctr(i)%x = (i-IXMIN+HALF)*dx+left
            ctr(i)%length = dx
        end forall
    end subroutine InitUniformMesh

    subroutine InitDataOneMidPoint(leftVec,rightVec,mid)
        real(KREAL), intent(in)                         :: leftVec(3), rightVec(3), mid
        integer(KINT)                                   :: index, i
        real(KREAL)                                     :: x

        index = IXMIN
        x = START_POINT
        do while(x<=mid)
            x = x+ctr(index)%length
            index = index+1
        enddo

        do i=IXMIN-GHOST_NUM,index-1
            ctr(i)%conVars = leftVec

        enddo
        do i=index,IXMAX+GHOST_NUM
            ctr(i)%conVars = rightVec
        enddo
    end subroutine InitDataOneMidPoint

    subroutine InitDataTwoMidPoints(leftVec,midVec,rightVec,mid1,mid2)
        real(KREAL), intent(in)                         :: leftVec(3), midVec(3), rightVec(3), mid1, mid2
        integer(KINT)                                   :: index1, index2, i
        real(KREAL)                                     :: x

        index1 = IXMIN
        x = START_POINT
        do while(x<=mid1)
            x = x+ctr(index1)%length
            index1 = index1+1
        enddo

        index2 = index1
        do while(x<=mid2)
            x = x+ctr(index2)%length
            index2 = index2+1
        enddo

        do i=IXMIN-GHOST_NUM,index1-1
            ctr(i)%conVars = leftVec

        enddo
        do i=index1,index2-1
            ctr(i)%conVars = midVec
        enddo
        do i=index2,IXMAX+GHOST_NUM
            ctr(i)%conVars = rightVec
        enddo
    end subroutine InitDataTwoMidPoints
end module Initialization


!--------------------------------------------------
!>KFVS solver
!--------------------------------------------------
module Solver
    use ConstantVariables
    use Initialization
    implicit none

contains
    !--------------------------------------------------
    !>convert conservative variables to primary variables
    !>@param[in] w           :conservative variables
    !>@return    GetPrimary  :primitive variables
    !--------------------------------------------------
    function GetPrimary(w)
        real(KREAL), intent(in)                         :: w(3) !conservative variables
        real(KREAL)                                     :: GetPrimary(3) !primitive variables

        GetPrimary(1) = w(1)
        GetPrimary(2) = w(2)/w(1)
        GetPrimary(3) = w(1)/(2.0_KREAL*GetLambda(w))
    end function GetPrimary

    function GetLambda(w)
        real(KREAL), intent(in)                         :: w(3)
        real(KREAL)                                     :: GetLambda

        GetLambda = (CK+1)*QUATER*w(1)/(w(3)-HALF*w(2)**2/w(1))
    end function GetLambda
    !--------------------------------------------------
    !>obtain speed of sound
    !>@param[in] prim    :primary variables
    !>@return    GetSoundSpeed :speed of sound
    !--------------------------------------------------
    function GetSoundSpeed(prim)
        real(KREAL), intent(in)                         :: prim(3)
        real(KREAL)                                     :: GetSoundSpeed !speed of sound

        GetSoundSpeed = sqrt(gamma*prim(3)/prim(1))
    end function GetSoundSpeed
    !--------------------------------------------------
    !>calculate time step
    !--------------------------------------------------
    subroutine TimeStep()
        real(KREAL)                                     :: prim(3) !primary variables
        real(KREAL)                                     :: sos !speed of sound
        real(KREAL)                                     :: vMax
        integer                                         :: i

        !set initial value
        prim = GetPrimary(ctr(IXMIN-GHOST_NUM)%conVars)
        sos = GetSoundSpeed(prim)
        vMax = max(abs(prim(2)+sos),abs(prim(2)-sos))
        dt = CFL*ctr(IXMIN-GHOST_NUM)%length/(vMax+TINYNUM)

        do i=IXMIN-GHOST_NUM+1,IXMAX+GHOST_NUM
            !convert conservative variables to primary variables
            prim = GetPrimary(ctr(i)%conVars)

            !sound speed
            sos = GetSoundSpeed(prim)

            !maximum 1/dt allowed
            vMax = max(abs(prim(2)+sos),abs(prim(2)-sos))
            !time step
            dt = min(dt,CFL*ctr(i)%length/(vMax+TINYNUM))
        end do 
    end subroutine TimeStep

    subroutine Boundary()
        integer(KINT)                                   :: i
        do i=1,GHOST_NUM
            ctr(IXMIN-i)%conVars = ctr(IXMIN-i+1)%conVars; ctr(IXMAX+i)%conVars = ctr(IXMAX+i-1)%conVars !free boundary
            ! ctr(IXMIN-i)%conVars(2) = -ctr(IXMIN-i)%conVars(2); ctr(IXMAX+i)%conVars(2) = -ctr(IXMAX+i)%conVars(2) !Reflecting boundary conditions at both sides
        enddo
    end subroutine Boundary

    subroutine IntraCellReconstruction()
        integer(KINT)                                   :: i
        do i=IXMIN-1,IXMAX+1
            call VanleerLimiter(ctr(i-1),ctr(i),ctr(i+1))
            ! ctr(i)%leftVars = ctr(i)%conVars; ctr(i)%rightVars = ctr(i)%conVars !first order reconstruction, no limiter
        enddo 
    end subroutine IntraCellReconstruction

    subroutine VanleerLimiter(leftCell,midCell,rightCell)
        type(IntraCell), intent(in)                     :: leftCell, rightCell
        type(IntraCell), intent(inout)                  :: midCell
        real(KREAL), dimension(3)                       :: splus, sminus
        real(KREAL), parameter                          :: UP = 1.0_KREAL

        splus = (rightCell%conVars-midCell%conVars)/(HALF*(rightCell%length+midCell%length))
        sminus = (midCell%conVars-leftCell%conVars)/(HALF*(midCell%length+leftCell%length))

        midCell%leftVars = midCell%conVars-HALF*midCell%length*(sign(UP,splus)+sign(UP,sminus))*abs(splus)*abs(sminus)/(abs(splus)+abs(sminus)+TINYNUM)
        midCell%rightVars = midCell%conVars+HALF*midCell%length*(sign(UP,splus)+sign(UP,sminus))*abs(splus)*abs(sminus)/(abs(splus)+abs(sminus)+TINYNUM)
        
        if ( GetLambda(midCell%leftVars)<=TINYNUM .or. GetLambda(midCell%rightVars)<=TINYNUM ) then
            midCell%leftVars = midCell%conVars
            midCell%rightVars = midCell%conVars
        endif
    end subroutine VanleerLimiter

    function GetTau(leftVal,midVal,rightVal) result(tau)
        real(KREAL), dimension(3), intent(in)           :: leftVal, midVal, rightVal !rho, U, lambda
        real(KREAL)                                     :: tau

        ! tau = MU/(HALF*midVal(1)/midVal(3))+abs(leftVal(1)/leftVal(3)-rightVal(1)/rightVal(3))/(abs(leftVal(1)/leftVal(3)+rightVal(1)/rightVal(3))+TINYNUM)*dt
        tau = 0.005_KREAL*dt+1.0_KREAL*abs(leftVal(1)/leftVal(3)-rightVal(1)/rightVal(3))/(abs(leftVal(1)/leftVal(3)+rightVal(1)/rightVal(3))+TINYNUM)*dt
    end function GetTau
    !--------------------------------------------------
    !>flux calculation
    !--------------------------------------------------
    subroutine CalcFlux(leftCell,face,rightCell)
        type(IntraCell), intent(in)                     :: leftCell, rightCell
        type(CellInterface), intent(inout)              :: face
        real(KREAL), dimension(3)                       :: al, AAL
        real(KREAL), dimension(3)                       :: ar, AAR
        real(KREAL), dimension(3)                       :: al0, ar0, AA0
        real(KREAL), dimension(3)                       :: vall, valr, val0 !rho, U, lambda
        real(KREAL), dimension(3)                       :: tempb, midConvars
        real(KREAL)                                     :: leftUtable(0:6), rightUtable(0:6)
        real(KREAL)                                     :: r(0:5), t(6), tau, eta

        vall(1) = leftCell%rightVars(1)
        vall(2) = leftCell%rightVars(2)/leftCell%rightVars(1)
        vall(3) = GetLambda(leftCell%rightVars)
        tempb = (leftCell%rightVars-leftCell%conVars)/(vall(1)*HALF*leftCell%length)
        call SolveEquMAB(al,tempb,vall)

        valr(1) = rightCell%leftVars(1)
        valr(2) = rightCell%leftVars(2)/rightCell%leftVars(1)
        valr(3) = GetLambda(rightCell%leftVars)
        tempb = (rightCell%conVars-rightCell%leftVars)/(valr(1)*HALF*rightCell%length)
        call SolveEquMAB(ar,tempb,valr)

        tempb = -matmul(MatrixWithU(vall,0),al)
        call SolveEquMAB(AAL,tempb,vall)

        tempb = -matmul(MatrixWithU(valr,0),ar)
        call SolveEquMAB(AAR,tempb,valr)

        leftUtable = ConstructUTable(vall,-1); rightUtable = ConstructUTable(valr,1)
        midConvars(1) = vall(1)*leftUtable(0)+valr(1)*rightUtable(0)
        midConvars(2) = vall(1)*leftUtable(1)+valr(1)*rightUtable(1)
        midConvars(3) = HALF*(vall(1)*(leftUtable(2)+leftUtable(0)*(HALF*CK/vall(3))) &
                            +valr(1)*(rightUtable(2)+rightUtable(0)*(HALF*CK/valr(3))))

        val0(1) = midConvars(1)
        val0(2) = midConvars(2)/midConvars(1)
        val0(3) = GetLambda(midConvars)

        tempb = (midConvars-leftCell%conVars)/(val0(1)*HALF*leftCell%length)
        call SolveEquMAB(al0,tempb,val0)
        tempb = (rightCell%conVars-midConvars)/(val0(1)*HALF*rightCell%length)
        call SolveEquMAB(ar0,tempb,val0)

        tau = GetTau(vall,val0,valr)
        eta = exp(-dt/tau)

        r(0) = dt-tau*(1.0_KREAL-eta)
        r(1) = -(1.0_KREAL-eta)/r(0)
        r(2) = (-dt+2.0_KREAL*tau*(1.0_KREAL-eta)-dt*eta)/r(0)
        r(3) = (1.0_KREAL-eta)/r(0)
        r(4) = (dt*eta-tau*(1.0_KREAL-eta))/r(0)
        r(5) = -tau*(1.0_KREAL-eta)/r(0) !for NS

        ! r(5) = 0.0_KREAL !for Euler

        tempb = r(1)*VectorWithoutU(val0,0) &
                +r(2)*(matmul(MatrixWithU(val0,-1),al0)+matmul(MatrixWithU(val0,1),ar0)) &
                +r(3)*(VectorWithoutU(vall,-1)+VectorWithoutU(valr,1)) &
                +r(4)*(matmul(MatrixWithU(vall,-1),al)+matmul(MatrixWithU(valr,1),ar)) &
                +r(5)*(matmul(MatrixWithU(vall,-1),al)+matmul(MatrixWithoutU(vall,-1),AAL) &
                        +matmul(MatrixWithU(valr,1),ar)+matmul(MatrixWithoutU(valr,1),AAR))
        call SolveEquMAB(AA0,tempb,val0)

        t(1) = dt+tau*(eta-1.0_KREAL)
        t(2) = -tau*(dt+2.0_KREAL*tau*(eta-1)+dt*eta)
        t(3) = HALF*dt**2-tau*dt-tau**2*(eta-1.0_KREAL)
        t(4) = -tau*(eta-1.0_KREAL)
        t(5) = tau*dt*eta+tau**2*(eta-1.0_KREAL)
        t(6) = tau**2*(eta-1.0_KREAL) !for NS

        ! t(6) = 0 !for Euler
        
        ! KFVS 1st
        ! we should close vanleer limiter and use 1st reconstruction, otherwise there will have oscillation.
        ! t(1) = 0.0_KREAL
        ! t(2) = 0.0_KREAL
        ! t(3) = 0.0_KREAL
        ! t(4) = dt
        ! t(5) = 0.0_KREAL
        ! t(6) = 0.0_KREAL

        ! GKS 1st
        ! we should close vanleer limiter and use 1st reconstruction, otherwise there will have oscillation.
        ! eta = 0.5_KREAL
        ! t(1) = dt*(1.0_KREAL-eta)
        ! t(2) = 0.0_KREAL
        ! t(3) = 0.0_KREAL
        ! t(4) = dt*eta
        ! t(5) = 0.0_KREAL
        ! t(6) = 0.0_KREAL

        ! KFVS 2nd
        ! t(1) = 0.0_KREAL
        ! t(2) = 0.0_KREAL
        ! t(3) = 0.0_KREAL
        ! t(4) = dt
        ! t(5) = -HALF*dt**2
        ! t(6) = 0.0_KREAL

        face%flux =  t(1)*val0(1)*VectorWithU(val0,0) &
                    +t(2)*val0(1)*(matmul(MatrixWith2U(val0,-1),al0)+matmul(MatrixWith2U(val0,1),ar0)) &
                    +t(3)*val0(1)*matmul(MatrixWithU(val0,0),AA0) &
                    +t(4)*(vall(1)*VectorWithU(vall,-1)+valr(1)*VectorWithU(valr,1)) &
                    +t(5)*(vall(1)*matmul(MatrixWith2U(vall,-1),al)+valr(1)*matmul(MatrixWith2U(valr,1),ar)) &
                    +t(6)*(vall(1)*(matmul(MatrixWith2U(vall,-1),al)+matmul(MatrixWithU(vall,-1),AAL)) &
                            +valr(1)*(matmul(MatrixWith2U(valr,1),ar)+matmul(MatrixWithU(valr,1),AAR)))

    end subroutine CalcFlux

    subroutine SolveEquMAB(a,b,val)
        real(KREAL), dimension(3), intent(in) :: b, val
        real(KREAL), dimension(3), intent(inout) :: a
        real(KREAL) :: U, lambda, AA, BB

        U = val(2); lambda = val(3)
        AA = b(2)-U*b(1); BB = 2.0_KREAL*b(3)-(U**2+HALF*(CK+1)/lambda)*b(1)

        a(3) = 4.0_KREAL*lambda**2/(CK+1)*(BB-2.0_KREAL*U*AA)
        a(2) = 2.0_KREAL*lambda*(AA-HALF*a(3)*U/lambda)
        a(1) = b(1)-a(2)*U-a(3)*(HALF*U**2+QUATER*(CK+1)/lambda)
    end subroutine SolveEquMAB

    function ConstructUTable(val,selector) result(uTable)
        real(KREAL), dimension(3), intent(in) :: val
        integer(KINT), intent(in) :: selector
        real(KREAL) :: uTable(0:6)
        real(KREAL) :: alpha, beta, U, lambda

        U = val(2); lambda = val(3)
        select case (selector)
        case (-1) !construct left table
            alpha = erfc(-sqrt(lambda)*U); beta = exp(-lambda*U**2)/sqrt(PI*lambda)
            uTable(0) = HALF*alpha
            uTable(1) = U*uTable(0)+HALF*beta
            uTable(2) = U*uTable(1)+HALF/lambda*uTable(0)
            uTable(3) = U*uTable(2)+1.0_KREAL/lambda*uTable(1)
            uTable(4) = U*uTable(3)+1.5_KREAL/lambda*uTable(2)
            uTable(5) = U*uTable(4)+2.0_KREAL/lambda*uTable(3)
            uTable(6) = U*uTable(5)+2.5_KREAL/lambda*uTable(4)
        case (0) !construct u0 table
            uTable(0) = 1
            uTable(1) = U
            uTable(2) = U*uTable(1)+HALF/lambda*uTable(0)
            uTable(3) = U*uTable(2)+1.0_KREAL/lambda*uTable(1)
            uTable(4) = U*uTable(3)+1.5_KREAL/lambda*uTable(2)
            uTable(5) = U*uTable(4)+2.0_KREAL/lambda*uTable(3)
            uTable(6) = U*uTable(5)+2.5_KREAL/lambda*uTable(4)
        case (1) !construct right table
            alpha = erfc(sqrt(lambda)*U); beta = exp(-lambda*U**2)/sqrt(PI*lambda)
            uTable(0) = HALF*alpha
            uTable(1) = U*uTable(0)-HALF*beta
            uTable(2) = U*uTable(1)+HALF/lambda*uTable(0)
            uTable(3) = U*uTable(2)+1.0_KREAL/lambda*uTable(1)
            uTable(4) = U*uTable(3)+1.5_KREAL/lambda*uTable(2)
            uTable(5) = U*uTable(4)+2.0_KREAL/lambda*uTable(3)
            uTable(6) = U*uTable(5)+2.5_KREAL/lambda*uTable(4)
            case default
            stop "something error happen in ConstructUTable"
        end select
    end function ConstructUTable

    function MatrixWithoutU(val,selector) result(uMatrix)
        real(KREAL), dimension(3), intent(in) :: val
        integer(KINT), intent(in) :: selector
        real(KINT) :: uMatrix(3,3)
        real(KREAL) :: tempUtable(0:6), lambda, xi2, xi4

        lambda = val(3); xi2 = HALF*CK/lambda; xi4 = QUATER*(CK**2+2.0_KREAL*CK)/lambda**2
        tempUtable = ConstructUTable(val,selector)
        uMatrix(1,1) = tempUtable(0); uMatrix(1,2) = tempUtable(1); uMatrix(1,3) = HALF*(tempUtable(2)+tempUtable(0)*xi2)
        uMatrix(2,1) = tempUtable(1); uMatrix(2,2) = tempUtable(2); uMatrix(2,3) = HALF*(tempUtable(3)+tempUtable(1)*xi2)
        uMatrix(3,1) = HALF*(tempUtable(2)+tempUtable(0)*xi2); uMatrix(3,2) = HALF*(tempUtable(3)+tempUtable(1)*xi2);
        uMatrix(3,3) = QUATER*(tempUtable(4)+2.0_KREAL*tempUtable(2)*xi2+tempUtable(0)*xi4)
    end function MatrixWithoutU
    function MatrixWithU(val,selector) result(uMatrix)
        real(KREAL), dimension(3), intent(in) :: val
        integer(KINT), intent(in) :: selector
        real(KINT) :: uMatrix(3,3)
        real(KREAL) :: tempUtable(0:6), lambda, xi2, xi4

        lambda = val(3); xi2 = HALF*CK/lambda; xi4 = QUATER*(CK**2+2.0_KREAL*CK)/lambda**2
        tempUtable = ConstructUTable(val,selector)
        uMatrix(1,1) = tempUtable(1); uMatrix(1,2) = tempUtable(2); uMatrix(1,3) = HALF*(tempUtable(3)+tempUtable(1)*xi2)
        uMatrix(2,1) = tempUtable(2); uMatrix(2,2) = tempUtable(3); uMatrix(2,3) = HALF*(tempUtable(4)+tempUtable(2)*xi2)
        uMatrix(3,1) = HALF*(tempUtable(3)+tempUtable(1)*xi2); uMatrix(3,2) = HALF*(tempUtable(4)+tempUtable(2)*xi2);
        uMatrix(3,3) = QUATER*(tempUtable(5)+2.0_KREAL*tempUtable(3)*xi2+tempUtable(1)*xi4)
    end function MatrixWithU
    function MatrixWith2U(val,selector) result(uMatrix)
        real(KREAL), dimension(3), intent(in) :: val
        integer(KINT), intent(in) :: selector
        real(KINT) :: uMatrix(3,3)
        real(KREAL) :: tempUtable(0:6), lambda, xi2, xi4

        lambda = val(3); xi2 = HALF*CK/lambda; xi4 = QUATER*(CK**2+2.0_KREAL*CK)/lambda**2
        tempUtable = ConstructUTable(val,selector)
        uMatrix(1,1) = tempUtable(2); uMatrix(1,2) = tempUtable(3); uMatrix(1,3) = HALF*(tempUtable(4)+tempUtable(2)*xi2)
        uMatrix(2,1) = tempUtable(3); uMatrix(2,2) = tempUtable(4); uMatrix(2,3) = HALF*(tempUtable(5)+tempUtable(3)*xi2)
        uMatrix(3,1) = HALF*(tempUtable(4)+tempUtable(2)*xi2); uMatrix(3,2) = HALF*(tempUtable(5)+tempUtable(3)*xi2);
        uMatrix(3,3) = QUATER*(tempUtable(6)+2.0_KREAL*tempUtable(4)*xi2+tempUtable(2)*xi4)
    end function MatrixWith2U

    function VectorWithoutU(val,selector) result(uVec)
        real(KREAL), dimension(3), intent(in) :: val
        integer(KINT), intent(in) :: selector
        real(KINT) :: uVec(3)
        real(KREAL) :: tempUtable(0:6), lambda, xi2

        lambda = val(3); xi2 = HALF*CK/lambda
        tempUtable = ConstructUTable(val,selector)
        uVec(1) = tempUtable(0)
        uVec(2) = tempUtable(1)
        uVec(3) = HALF*(tempUtable(2)+tempUtable(0)*xi2)
    end function VectorWithoutU
    function VectorWithU(val,selector) result(uVec)
        real(KREAL), dimension(3), intent(in) :: val
        integer(KINT), intent(in) :: selector
        real(KINT) :: uVec(3)
        real(KREAL) :: tempUtable(0:6), lambda, xi2

        lambda = val(3); xi2 = HALF*CK/lambda
        tempUtable = ConstructUTable(val,selector)
        uVec(1) = tempUtable(1)
        uVec(2) = tempUtable(2)
        uVec(3) = HALF*(tempUtable(3)+tempUtable(1)*xi2)
    end function VectorWithU

    !--------------------------------------------------
    !>calculate the flux across the interfaces
    !--------------------------------------------------
    subroutine Evolution()
        integer(KINT)                               :: i

        do i=IXMIN,IXMAX+1 !with ghost cell
            call CalcFlux(ctr(i-1),vface(i),ctr(i))
        end do
    end subroutine Evolution

    subroutine Update()
        integer(KINT)                               :: i,j

        do i=IXMIN,IXMAX
            ctr(i)%conVars = ctr(i)%conVars+1.0_KREAL/ctr(i)%length*(vface(i)%flux-vface(i+1)%flux)
            do j=1,3
                if(isnan(ctr(i)%conVars(j))) then
                    print *, ctr(i)%x, "iter:", iter
                    stop
                endif
            enddo
        end do
    end subroutine Update
end module Solver

module Writer
    use Solver
    implicit none

contains
    !--------------------------------------------------
    !>write result
    !--------------------------------------------------
    subroutine output()
        integer(KINT)                               :: i
        real(KREAL), dimension(:,:), allocatable    :: solution
        !--------------------------------------------------
        !prepare solutions
        !--------------------------------------------------
        allocate(solution(3,IXMIN-1:IXMAX+1)) !including ghost cell

        do i=IXMIN-1,IXMAX+1
            solution(:,i) = GetPrimary(ctr(i)%conVars)
        end do

        !--------------------------------------------------
        !write to file
        !--------------------------------------------------
        !open result file and write header
        open(unit=RSTFILE,file=RSTFILENAME,status="replace",action="write")
        ! write(RSTFILE,*) "VARIABLES = x, Density, U, Pressure"
        write(RSTFILE,*) "VARIABLES = x, Density, U, Pressure, Momentum"
        write(RSTFILE,*) 'ZONE  T="Time: ',simTime,'", I = ',IXMAX-IXMIN+1,', DATAPACKING=BLOCK'

        !write geometry (cell-centered)
        write(RSTFILE,"(6(ES23.16,2X))") ctr(IXMIN:IXMAX)%x

        !write solution (cell-centered)
        do i=1,3
            write(RSTFILE,"(6(ES23.16,2X))") solution(i,IXMIN:IXMAX)
        end do

        !write momentum
        write(RSTFILE,"(6(ES23.16,2X))") ctr(IXMIN:IXMAX)%conVars(2)

        !close file
        close(RSTFILE)
        deallocate(solution)
    end subroutine output
end module Writer


!--------------------------------------------------
!>main program
!--------------------------------------------------
program GKS_2nd_1D
    use Solver
    use Writer
    implicit none

    !initialization
    call Init()

    !open file and write header
    open(unit=HSTFILE,file=HSTFILENAME,status="replace",action="write") !open history file
    write(HSTFILE,*) "VARIABLES = iter, simTime, dt" !write header

    !iteration
    do while(.true.)
        call TimeStep() !calculate time step
        call Boundary() !set up Boundary condition
        call IntraCellReconstruction() !Reconstruction
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

end program GKS_2nd_1D