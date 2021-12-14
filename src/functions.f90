! additional functions can be added here.
! currently only has linear func


module functions
contains

real function y(x,x1,y1,m) !equation of line

    implicit none
    real, intent(in) :: x,x1,y1,m

    y = y1 + m*(x-x1)

end function y

end module functions
