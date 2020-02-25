module vec_solve

        implicit none

contains

subroutine solve(y,y2,y1,x2,x1,x,control,nrofpts,slope,x_in,y_out,point,newpoint)

        implicit none
        real, intent(in) :: y2,y1,x2,x1,x
        real, intent(out) :: slope,y_out,x_in
        real, intent(out),dimension(2,1) :: newpoint
        real, external :: y
        real, dimension(:,:), allocatable, intent(out) :: point
        integer, intent(out) :: nrofpts
        integer, intent(in) :: control
        real, parameter :: sampling_step = 0.1,tolerance=0.01
        integer, parameter :: nrofstates = 2
        integer :: k,i
        real :: delta_x, delta_y, x_new

        delta_x=x2-x1
        delta_y=y2-y1
        slope=delta_y/delta_x
        nrofpts=abs(delta_x/sampling_step)+1

        allocate (point(nrofstates,nrofpts))

               if (control .eq. 0) then
               y_out=y(x,x1,y1,slope)
               x_in=x
               deallocate (point)

               else if (control .eq. 1) then
               do i=1,nrofpts
               y_out=y(x_in,x1,y1,slope)
               point(1,i) = x_in
               point(2,i) = y_out
               k=k+1
               if (x1 .gt. x2) then
               x_in = x_in - sampling_step
               else if (x1 .lt. x2) then
               x_in = x_in + sampling_step
               end if

               if (abs(int(x_in*100)) -abs(int(x2*100)) .eq. 0) then !including end point
                       exit
               end if
               end do
       end if

end subroutine solve

end module vec_solve
