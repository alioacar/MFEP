module stepping
  implicit none
contains

subroutine steepestpointfinder(inputgroup,inputgroupsize,startingindex,newptindex,minpoint)

        implicit none
        integer,intent(in) :: inputgroupsize,startingindex,newptindex
        real,dimension(:,:),intent(in) :: inputgroup(3,inputgroupsize)
        integer,intent(out) :: minpoint
        integer :: i,j,k
        real :: diffx,diffy
        real,dimension(:,:) :: difftot(4,inputgroupsize)


        do j=startingindex,inputgroupsize
        diffx=inputgroup(1,j)-inputgroup(1,newptindex)
        diffy=inputgroup(2,j)-inputgroup(2,newptindex)
        difftot(1,j)=sqrt(diffx**2+diffy**2)
        difftot(2,j)=inputgroup(1,j)
        difftot(3,j)=inputgroup(2,j)
        difftot(4,j)=inputgroup(3,j)

        end do
        minpoint=minloc(difftot(1,:),dim=1,mask=(difftot(1,:) .gt. 1E-3))

end subroutine steepestpointfinder

end module stepping
