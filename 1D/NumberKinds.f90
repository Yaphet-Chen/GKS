! NumberKinds.f90
! Specify the kinds of real and integer.
module NumberKinds
    implicit none

    integer, parameter      :: KREAL = kind(0.d0)
    integer, parameter      :: KINT = kind(1)

end module NumberKinds