!--------------------------------------------------
!>store the global constant variables
!--------------------------------------------------
module ConstantVariables
    use NumberKinds
    implicit none
    
    real(KREAL), parameter :: PI = 4.0_KREAL*atan(1.0_KREAL) !Pi
    real(KREAL), parameter :: TINYNUM = tiny(0.0_KREAL) !small value to avoid 0/0
    real(KREAL), parameter :: HALF = 0.5_KREAL
    real(KREAL), parameter :: QUATER = 0.25_KREAL    
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
        real(KREAL) :: x !cell center xinates
        !flow field
        real(KREAL) :: vars(3) !density, x-momentum,total energy
    end type CellCenter

    !cell interface
    type CellInterface
        real(KREAL) :: flux(3) !mass flux, x momentum flux, energy flux
    end type CellInterface
end module Mesh

module Initialization
    use WoodwardColellaModule
    implicit none
    
    contains
        !--------------------------------------------------
        !>main initialization subroutine
        !--------------------------------------------------
        subroutine Init()  
            call InitGeometry(START_POINT,END_POINT,POINTS_NUM,xSpacing) !initialize the geometry
            ! call initData(VARS_LEFT,VARS_RIGHT,MID_POINT) !set the initial condition
            call initDataTwoPoints(VARS_LEFT,VARS_MID,VARS_RIGHT,MID_POINT_1,MID_POINT_2) !set the initial condition
            gamma = GetGamma(CK)
        end subroutine Init

        subroutine InitGeometry(left,right,num,dx)
            real(KREAL), intent(in) :: left, right
            integer(KINT), intent(in) :: num
            real(KREAL), intent(inout) :: dx
            integer(KINT) :: i

            !cell length
            dx = (right-left)/(IXMAX-IXMIN+1)

            !cell center (with ghost cell)
            forall(i=IXMIN-1:IXMAX+1) 
                ctr(i)%x = (i-0.5)*dx+left
            end forall
        end subroutine InitGeometry

        subroutine InitData(leftVec,rightVec,mid)
            real(KREAL), intent(in) :: leftVec(3), rightVec(3), mid
            integer(KINT) :: index, i
            real(KREAL) :: x

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
            real(KREAL), intent(in) :: leftVec(3), midVec(3), rightVec(3), mid1, mid2
            integer(KINT) :: index1, index2, i
            real(KREAL) :: x

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
        !>@param[in] ck        :internal degree of freedom
        !>@return    GetGamma :ratio of specific heat
        !--------------------------------------------------
        function GetGamma(k)
            integer(KINT), intent(in) :: k
            real(KREAL) :: GetGamma

            GetGamma = real(k+3,KREAL)/real(k+1,KREAL)
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
            real(KREAL), intent(in) :: w(3)
            real(KREAL) :: GetPrimary(3) !primary variables

            GetPrimary(1) = w(1)
            GetPrimary(2) = w(2)/w(1)
            GetPrimary(3) = (gamma-1.0)*(w(3)-HALF*w(2)**2/w(1))
        end function GetPrimary

        !--------------------------------------------------
        !>calculate time step
        !--------------------------------------------------
        subroutine TimeStep()
            real(KREAL) :: prim(3) !primary variables
            real(KREAL) :: sos !speed of sound
            real(KREAL) :: vMax
            integer :: i

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
            real(KREAL), intent(in) :: prim(3)
            real(KREAL) :: GetSoundSpeed !speed of sound

            GetSoundSpeed = sqrt(gamma*prim(3)/prim(1))
        end function GetSoundSpeed
        !--------------------------------------------------
        !>flux calculation
        !--------------------------------------------------
        subroutine CalcFlux(cell_L,face,cell_R)
            type(CellCenter), intent(in) :: cell_L, cell_R
            type(CellInterface), intent(inout) :: face
            real(KREAL) :: wL(3), primL(3), UL, erL, EL
            real(KREAL) :: wR(3), primR(3), UR, erR, ER
            real(KREAL) :: lamL, lamR

            wL = cell_L%vars; wR = cell_R%vars
            primL = GetPrimary(wL); primR = GetPrimary(wR)
            ! UL = primL(2); UR = primR(2)
            lamL = GetLambda(wL); lamR = GetLambda(wR)

            face%flux(1) = wL(1)*(HALF*primL(2)*erfc(-sqrt(lamL)*primL(2))+HALF*exp(-lamL*primL(2)**2)/sqrt(PI*lamL)) &
                        +wR(1)*(HALF*primR(2)*erfc(sqrt(lamR)*primR(2))-HALF*exp(-lamR*primR(2)**2)/sqrt(PI*lamR))
            face%flux(2) = wL(1)*(HALF*(primL(2)**2+HALF/lamL)*erfc(-sqrt(lamL)*primL(2))+HALF*primL(2)*exp(-lamL*primL(2)**2)/sqrt(PI*lamL)) &
                        +wR(1)*(HALF*(primR(2)**2+HALF/lamR)*erfc(sqrt(lamR)*primR(2))-HALF*primR(2)*exp(-lamR*primR(2)**2)/sqrt(PI*lamR))
            face%flux(3) = wL(1)*(QUATER*primL(2)*(primL(2)**2+HALF*(CK+3)/lamL)*erfc(-sqrt(lamL)*primL(2))+QUATER*(primL(2)**2+HALF*(CK+2)/lamL)*exp(-lamL*primL(2)**2)/sqrt(PI*lamL)) &
                        +wR(1)*(QUATER*primR(2)*(primR(2)**2+HALF*(CK+3)/lamR)*erfc(sqrt(lamR)*primR(2))-QUATER*(primR(2)**2+HALF*(CK+2)/lamR)*exp(-lamR*primR(2)**2)/sqrt(PI*lamR))

            ! erL = erfc(-sqrt(lamL)*UL); erR = erfc(sqrt(lamR)*UR)
            ! EL = exp(-lamL*UL**2)/sqrt(PI*lamL); ER = exp(-lamR*UR**2)/sqrt(PI*lamR)
            ! face%flux(1) = wL(1)*(UL/2.0*erL+0.5*EL) &
            !                 +wR(1)*(UR/2.0*erR-0.5*ER)
            ! face%flux(2) =  wL(1)*((UL**2/2.0+1.0/(4.0*lamL))*erL+UL/2.0*EL) &
            !                 +wR(1)*((UR**2/2.0+1/(4.0*lamR))*erR-UR/2.0*ER)
            ! face%flux(3) =  wL(1)*((UL**3/4.0+(CK+3)*UL/(8.0*lamL))*erL+(UL**2/4.0+(CK+2)/(8.0*lamL))*EL) &
            !                 +wR(1)*((UR**3/4.0+(CK+3)*UR/(8.0*lamR))*erR-(UR**2/4.0+(CK+2)/(8.0*lamR))*ER)
            contains
                function GetLambda(w)
                    real(KREAL), intent(in) :: w(3)
                    real(KREAL) :: GetLambda
                    
                    GetLambda = (CK+1)*QUATER*w(1)/(w(3)-HALF*w(2)**2/w(1))
                end function GetLambda
        end subroutine CalcFlux

        !--------------------------------------------------
        !>calculate the flux across the interfaces
        !--------------------------------------------------
        subroutine Evolution()
            integer(KINT) :: i

            do i=IXMIN,IXMAX+1 !with ghost cell
                call CalcFlux(ctr(i-1),vface(i),ctr(i))
            end do
        end subroutine Evolution

        subroutine Update()
            integer(KINT) :: i
            real(KREAL) :: sigma

            sigma = dt/xSpacing
            do i=IXMIN,IXMAX
                ctr(i)%vars = ctr(i)%vars+sigma*(vface(i)%flux-vface(i+1)%flux)
            end do
            ctr(IXMIN-1)%vars = ctr(IXMIN)%vars; ctr(IXMAX+1)%vars = ctr(IXMAX)%vars
            ctr(IXMIN-1)%vars(2) = -ctr(IXMIN-1)%vars(2); ctr(IXMAX+1)%vars(2) = -ctr(IXMAX+1)%vars(2) ! Reflecting boundary conditions at both sides
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
            integer(KINT) :: i
            real(KREAL), dimension(:,:), allocatable :: solution
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
            write(RSTFILE,*) "VARIABLES = x, Density, U, Pressure"
            write(RSTFILE,*) 'ZONE  T="Time: ',simTime,'", I = ',IXMAX-IXMIN+1,', DATAPACKING=BLOCK'

            !write geometry (cell-centered)
            write(RSTFILE,"(6(ES23.16,2X))") ctr(IXMIN:IXMAX)%x

            !write solution (cell-centered)
            do i=1,3
                write(RSTFILE,"(6(ES23.16,2X))") solution(i,IXMIN:IXMAX)
            end do
    
            !close file
            close(RSTFILE)
            deallocate(solution)

        end subroutine output
end module Writer


!--------------------------------------------------
!>main program
!--------------------------------------------------
program KFVS_1st_1D
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
    
end program KFVS_1st_1D