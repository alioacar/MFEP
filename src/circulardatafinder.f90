module circulardatafinder
contains
subroutine circdatafinder(xp,yp,xc,yc,goodxp,goodyp)

        implicit none
        real, intent(out) :: goodxp,goodyp
        real, intent(in) :: xp,yp,xc,yc
        real,parameter :: circleradius = 0.2, tol=0.0225
        real :: badxp,badyp,d
        d=sqrt((xp-xc)**2+(yp-yc)**2)
        if (abs(d**2 - circleradius**2) .le. tol) then
                goodxp=xp
                goodyp=yp
        else
                goodxp=0
                goodyp=0
        end if

end subroutine

end module circulardatafinder
