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
    type CellCenter
        !geometry
        real(KREAL)                                     :: x !cell center coordinates
        !flow field
        real(KREAL)                                     :: vars(3) !density, x-momentum, total energy
    end type CellCenter

    !cell interface
    type CellInterface
        real(KREAL)                                     :: flux(3) !mass flux, x momentum flux, total energy flux
    end type CellInterface
end module Mesh

module Initialization
    use ConstantVariables
    use StationaryShockModule
    implicit none
    
    contains
        !--------------------------------------------------
        !>main initialization subroutine
        !--------------------------------------------------
        subroutine Init()  
            call InitGeometry(START_POINT,END_POINT,POINTS_NUM,xSpacing) !initialize the geometry
            call initData(VARS_LEFT,VARS_RIGHT,MID_POINT) !set the initial condition with one mid point
            ! call initDataTwoPoints(VARS_LEFT,VARS_MID,VARS_RIGHT,MID_POINT_1,MID_POINT_2) !set the initial condition with two mid points
            gamma = GetGamma(CK)
        end subroutine Init

        subroutine InitGeometry(left,right,num,dx)
            real(KREAL), intent(in)                     :: left, right
            integer(KINT), intent(in)                   :: num
            real(KREAL), intent(inout)                  :: dx
            integer(KINT)                               :: i

            !cell length
            dx = (right-left)/(IXMAX-IXMIN+1)

            !cell center (with ghost cell)
            forall(i=IXMIN-1:IXMAX+1) 
                ctr(i)%x = (i-HALF)*dx+left
            end forall
        end subroutine InitGeometry

        subroutine InitData(leftVec,rightVec,mid)
            real(KREAL), intent(in)                     :: leftVec(3), rightVec(3), mid
            integer(KINT)                               :: index, i
            real(KREAL)                                 :: x

            index = IXMIN
            x = START_POINT
            do while(x<=mid)
                x = x+xSpacing
                index = index+1
            enddo

            do i=IXMIN-1,index-1
                ctr(i)%vars = leftVec
                
            enddo
            do i=index,IXMAX+1
                ctr(i)%vars = rightVec
            enddo
        end subroutine initData

        subroutine InitDataTwoPoints(leftVec,midVec,rightVec,mid1,mid2)
            real(KREAL), intent(in)                     :: leftVec(3), midVec(3), rightVec(3), mid1, mid2
            integer(KINT)                               :: index1, index2, i
            real(KREAL)                                 :: x

            index1 = IXMIN
            x = START_POINT
            do while(x<=mid1)
                x = x+xSpacing
                index1 = index1+1
            enddo

            index2 = index1
            do while(x<=mid2)
                x = x+xSpacing
                index2 = index2+1
            enddo

            do i=IXMIN-1,index1-1
                ctr(i)%vars = leftVec
                
            enddo
            do i=index1,index2-1
                ctr(i)%vars = midVec
            enddo
            do i=index2,IXMAX+1
                ctr(i)%vars = rightVec
            enddo
        end subroutine initDataTwoPoints
        !--------------------------------------------------
        !>obtain ratio of specific heat
        !>@param[in] K        :internal degree of freedom
        !>@return    GetGamma :ratio of specific heat
        !--------------------------------------------------
        function GetGamma(K)
            integer(KINT), intent(in)                   :: K
            real(KREAL) :: GetGamma

            GetGamma = real(K+3,KREAL)/real(K+1,KREAL)
        end function GetGamma
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
            real(KREAL), intent(in)                     :: w(3) !conservative variables
            real(KREAL)                                 :: GetPrimary(3) !primitive variables

            GetPrimary(1) = w(1)
            GetPrimary(2) = w(2)/w(1)
            GetPrimary(3) = (gamma-1.0)*(w(3)-HALF*w(2)**2/w(1))
        end function GetPrimary

        !--------------------------------------------------
        !>calculate time step
        !--------------------------------------------------
        subroutine TimeStep()
            real(KREAL)                                 :: prim(3) !primary variables
            real(KREAL)                                 :: sos !speed of sound
            real(KREAL)                                 :: vMax
            integer                                     :: i

            !set initial value
            vMax = 0.0_KREAL

            do i=IXMIN,IXMAX
                !convert conservative variables to primary variables
                prim = GetPrimary(ctr(i)%vars)

                !sound speed
                sos = GetSoundSpeed(prim)

                !maximum 1/dt allowed
                vMax = max(vMax,abs(prim(2))+sos)
            end do 

            !time step
            dt = CFL*xSpacing/(vMax+TINYNUM)
        end subroutine TimeStep
        !--------------------------------------------------
        !>obtain speed of sound
        !>@param[in] prim    :primary variables
        !>@return    GetSoundSpeed :speed of sound
        !--------------------------------------------------
        function GetSoundSpeed(prim)
            real(KREAL), intent(in)                     :: prim(3)
            real(KREAL)                                 :: GetSoundSpeed !speed of sound

            GetSoundSpeed = sqrt(gamma*prim(3)/prim(1))
        end function GetSoundSpeed
        !--------------------------------------------------
        !>flux calculation
        !--------------------------------------------------
        subroutine CalcFlux(cell_L,face,cell_R,locEta)
            type(CellCenter), intent(in)                :: cell_L, cell_R
            type(CellInterface), intent(inout)          :: face
            real(KREAL)                                 :: locEta
            real(KREAL)                                 :: wL(3), primL(3)
            real(KREAL)                                 :: wR(3), primR(3)
            real(KREAL)                                 :: lamL, lamR
            real(KREAL)                                 :: UL, UR, erL, erR, EL, ER
            real(KREAL)                                 :: faceVars(3), facePrim(3), faceLam, faceU, eqFlux(3), nonEqFlux(3)
            
            wL = cell_L%vars; wR = cell_R%vars
            primL = GetPrimary(wL); primR = GetPrimary(wR)
            lamL = GetLambda(wL); lamR = GetLambda(wR)
            UL = primL(2); UR = primR(2)
            erL = erfc(-sqrt(lamL)*UL); erR = erfc(sqrt(lamR)*UR)
            EL = exp(-lamL*UL**2)/sqrt(PI*lamL); ER = exp(-lamR*UR**2)/sqrt(PI*lamR)

            ! non-equilibrium flux calculation
            nonEqFlux(1) = wL(1)*(HALF*UL*erL+HALF*EL) &
                            +wR(1)*(HALF*UR*erR-HALF*ER)
            nonEqFlux(2) =  wL(1)*((HALF*UL**2+QUATER/lamL)*erL+HALF*UL*EL) &
                            +wR(1)*((HALF*UR**2+QUATER/lamR)*erR-HALF*UR*ER)
            nonEqFlux(3) =  wL(1)*(QUATER*(UL**3+HALF*(CK+3)*UL/lamL)*erL+QUATER*(UL**2+HALF*(CK+2)/lamL)*EL) &
                            +wR(1)*(QUATER*(UR**3+HALF*(CK+3)*UR/lamR)*erR-QUATER*(UR**2+HALF*(CK+2)/lamR)*ER)

            ! equilibrium satate at the cell interface
            faceVars(1) = wL(1)*HALF*erL+wR(1)*HALF*erR
            faceVars(2) = wL(1)*(HALF*UL*erL+HALF*EL)+wR(1)*(HALF*UR*erR-HALF*ER)
            faceVars(3) = wL(1)*(QUATER*(UL**2+HALF*(CK+1)/lamL)*erL+QUATER*UL*EL) &
                        +wR(1)*(QUATER*(UR**2+HALF*(CK+1)/lamR)*erR-QUATER*UR*ER)

            ! equilibrium flux at the cell interface
            facePrim = GetPrimary(faceVars); faceLam = GetLambda(faceVars)
            faceU = facePrim(2)
            
            eqFlux(1) = faceVars(1)*faceU
            eqFlux(2) = faceVars(1)*(faceU**2+HALF/faceLam)
            eqFlux(3) = faceVars(1)*HALF*faceU*(faceU**2+(CK+3)*HALF/faceLam)

            ! the final fluxes across the cell interface
            face%flux = (1.0_KREAL-locEta)*eqFlux+locEta*nonEqFlux

            contains
                function GetLambda(w)
                    real(KREAL), intent(in)             :: w(3)
                    real(KREAL)                         :: GetLambda
                    
                    GetLambda = (CK+1)*QUATER*w(1)/(w(3)-HALF*w(2)**2/w(1))
                end function GetLambda
        end subroutine CalcFlux

        !--------------------------------------------------
        !>calculate the flux across the interfaces
        !--------------------------------------------------
        subroutine Evolution()
            integer(KINT)                               :: i

            do i=IXMIN,IXMAX+1 !with ghost cell
                call CalcFlux(ctr(i-1),vface(i),ctr(i),ETA)
            end do
        end subroutine Evolution

        subroutine Update()
            integer(KINT)                               :: i
            real(KREAL)                                 :: sigma

            sigma = dt/xSpacing
            do i=IXMIN,IXMAX
                ctr(i)%vars = ctr(i)%vars+sigma*(vface(i)%flux-vface(i+1)%flux)
            end do
            ctr(IXMIN-1)%vars = ctr(IXMIN)%vars; ctr(IXMAX+1)%vars = ctr(IXMAX)%vars
            ! ctr(IXMIN-1)%vars(2) = -ctr(IXMIN-1)%vars(2); ctr(IXMAX+1)%vars(2) = -ctr(IXMAX+1)%vars(2) !Reflecting boundary conditions at both sides
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
                solution(:,i) = GetPrimary(ctr(i)%vars)
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
            write(RSTFILE,"(6(ES23.16,2X))") ctr(IXMIN:IXMAX)%vars(2)

            !close file
            close(RSTFILE)
            deallocate(solution)

        end subroutine output
end module Writer


!--------------------------------------------------
!>main program
!--------------------------------------------------
program GKS_1st_1D
    use Solver
    use Writer
    implicit none

    !initialization
    call Init()
    call output()

    !open file and write header
    open(unit=HSTFILE,file=HSTFILENAME,status="replace",action="write") !open history file
    write(HSTFILE,*) "VARIABLES = iter, simTime, dt" !write header

    !iteration
    do while(.true.)
        call TimeStep() !calculate time step
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
    end do

    !close history file
    close(HSTFILE)

    !output solution
    call output()
    
end program GKS_1st_1D