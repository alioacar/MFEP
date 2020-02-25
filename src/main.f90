program LFEP

        use vec_solve
        use circulardatafinder
        use functions
        use veclength
        use stepping


        implicit none
        real,dimension(:,:),allocatable :: fenergy
        real,dimension(:,:),allocatable :: forward,backward,forward_out,backward_out
        real,dimension(:,:),allocatable :: pointsdiff,smoothingpts,dummynewpts
        real,dimension(:,:),allocatable :: point
        real,dimension(:,:),allocatable :: goodpoints
        real,dimension(3,1) :: steepestpoint
        real,dimension(2,3) :: r,v
        real,dimension(2,1) :: firstpoint,prenewpointgrad,vecpointgrad
        real,dimension(2,1) :: newpoint
  !      real,external :: y
        real :: veclengthtofinalx,veclengthtofinaly,veclengthtofinal
        real :: veclengthtoinitialx,veclengthtoinitialy,veclengthtoinitial
        real :: pointsdiffx,pointsdiffy,minpointsdiff
        real,parameter :: leftparam=0.25,rightparam=0.75
        real,parameter :: sampling_step = 0.1, tolerance = 0.01, circle_radius = 0.2
        integer,parameter :: nrofstates = 2, dimofengfile = 3, maxiters=5000,engfileline=10
        real,dimension(3,maxiters) :: minima
        real :: delta_y,delta_x,slope,y_sub,x_in,goodxp,goodyp,slope_vec,yvec,delta_x_vec,delta_y_vec,wf
        real :: point1_y,point1_x,ynew,pre_minima,incr,dummy
        integer :: stat,k,nrofpts,i,j,ii,jj,line,dimofgoodpoints,p,weighting,iwf,percent,direction
        integer :: new_i,outputdimension,line2,outputbackdim,newdim,m,a,b
        integer :: minpoint,forwardstartingindex,backwardstartingindex,u
        integer :: forwardindex,backwardindex,sizeofnewpts,iter
        real :: diffx,diffy,diffoptfirstx,diffoptfirsty,diffoptfirst
        real :: diffoptsecx,diffoptsecy,diffoptsec,checkparambackward,checkparamforward
        real,dimension(:,:),allocatable :: difftot,forwardsmoothing,backwardsmoothing,newpts
        logical :: DEBUG ,check,testwork
        character :: bsp = char(8)
        character(len=255) :: currentpath


        DEBUG = .false. ! type false to close debugging
        call getcwd(currentpath)

!define the vector between 2 min points

testwork = .true.

do 1717 direction=1,2 !forward and backward
        if (direction .eq. 1) then
open(unit=8,file='minima',status='old')
read(8,*,iostat=stat) r
close(8)
open(33, file=TRIM(currentpath)//'/sampling/forward-withwf', status='replace')
        else
open(unit=8,file='minima',status='old')
read(8,*,iostat=stat) v
close(8)
r(2,1)=v(2,2) !backward-final-y & forward-initial-y
r(1,1)=v(1,2) !backward-final-x & forward-initial-x
r(2,2)=v(2,1) !backward-initial-y & forward-final-y
r(1,2)=v(1,1) !backward-initial-x & forward-final-x
open(34, file=TRIM(currentpath)//'/sampling/backward-withwf', status='replace')
        end if

!(r(2,1) & r(1,1) y,x final state coordinates
!(r(2,2) & r(1,2) initial state coordinates
weighting = int(r(1,3))
wf = r(2,3)
incr = wf
iwf = 1
                if (DEBUG) then
                        print *, 'inputs from minima', r(2,1),r(1,1),r(2,2),r(1,2),r(1,3),r(2,3)
                end if

                if (weighting .eq. 1) then
                2003 format (' ... WARNING!!! weighting on and weighting factor =', f5.2)
                print 2003,wf
                else if (weighting .eq. 0) then
                print *,' ... WARNING!!! weighting off'
                 else
                print *,'... ERROR!!! wrong weighting parameters'
                stop
                end if

call solve(y,r(2,1),r(2,2),r(1,1),r(1,2),r(1,2),0,nrofpts,slope,x_in,y_sub,point,newpoint) !delta_y=r(2,1)-r(2,2)
!                                             |                             !delta_x=r(1,1)-r(1,2)
!                                         control param 0                   !slope= delta_y/delta_x
!                                      (just calc. parameters)              !nrofpts=abs(delta_x/sampling_step)+1
                2004 format ('... number of sampling for initial vector =', i3)
                print 2004,nrofpts

call solve(y,r(2,1),r(2,2),r(1,1),r(1,2),x_in,1,nrofpts,slope,x_in,y_sub,point,newpoint)
!control param 1 calc vector state A to B


                print *,'... the initial vector calculated successfully'
open(25, file=TRIM(currentpath)//'/sampling/file25', status='replace')
       do i=1,nrofpts
       write (25,*) point(1,i),point(2,i)
       end do
        firstpoint(1,1)=point(1,2)
        firstpoint(2,1)=point(2,2)

       deallocate (point)
       close(25)

open(unit=9,file='freefile',status='old')
line = 0
    do
    read(9,*,end=10)
       line = line +1
    end do
10 close(9)

allocate (fenergy(dimofengfile,line))


open(unit=9,file='freefile',status='old')
read(9,*,iostat=stat) fenergy
close(9)

                print *,'... free energy file opened successfully'

17 format('Iteration: ',I5)
27 format(5A1,I5)
write(*,17,advance='no')1
do
do iter=1,maxiters


if (iter .gt. 1) then
firstpoint(1,1)=minima(1,iter-1)
firstpoint(2,1)=minima(2,iter-1)
end if
!
write(*,27,advance='no')(bsp,j=1,5),iter


dimofgoodpoints=0
do j=1,line
call circdatafinder(fenergy(1,j),fenergy(2,j),firstpoint(1,1),firstpoint(2,1),goodxp,goodyp)
        if (goodxp .ne. 0 .and. goodyp .ne. 0) then
                dimofgoodpoints=dimofgoodpoints+1
        end if
end do
allocate (goodpoints(dimofengfile,dimofgoodpoints))

p=0
do j=1,line
call circdatafinder(fenergy(1,j),fenergy(2,j),firstpoint(1,1),firstpoint(2,1),goodxp,goodyp)
        if (goodxp .ne. 0 .and. goodyp .ne. 0) then
                p=p+1
                goodpoints(1,p)=goodxp
                goodpoints(2,p)=goodyp
                goodpoints(3,p)=fenergy(3,j)
        end if
end do
open(32, file=TRIM(currentpath)//'/sampling/file32', status='replace')
do j=1,dimofgoodpoints
!write(32,'(2f9.4,203f7.1)') goodpoints(1,j),goodpoints(2,j),goodpoints(3,j)

write(32,*) goodpoints(1,j),goodpoints(2,j),goodpoints(3,j)
end do
close(32)

pre_minima=goodpoints(3,1)
do i=2,dimofgoodpoints
    if (goodpoints(3,i)<pre_minima) pre_minima=goodpoints(3,i)
end do


do j=1,dimofgoodpoints
        if(pre_minima .eq. goodpoints(3,j)) then
                minima(1,iter)=goodpoints(1,j)
                minima(2,iter)=goodpoints(2,j)
        end if
end do

        if( any(minima == pre_minima) ) then

                if (DEBUG) then
                        print *, 'WARNING !!!! back to the same point, rescaling the vector'
                end if

                !print *,'y2',minima(2,iter-1),'y1',minima(2,iter-2),'x2',minima(1,iter-1),'x1',minima(1,iter-2)
                !print *,'deltaY', minima(2,iter-1),minima(2,iter-2)
                !print *,'deltaX', minima(1,iter-1)-minima(1,iter-2)

                delta_y=minima(2,iter-1)-minima(2,iter-2)
                delta_x=minima(1,iter-1)-minima(1,iter-2)
                if (abs(delta_x) .gt. tolerance) then
                slope=delta_y/delta_x
                                else if (abs(delta_x).lt. tolerance) then                  !if slope goes infinity
                                        if (delta_y .gt. 0) then
                                                prenewpointgrad(1,1)=minima(1,iter-1)
                                                prenewpointgrad(2,1)=minima(2,iter-1)+sampling_step
                                        else if (delta_y .lt. 0) then
                                                prenewpointgrad(1,1)=minima(1,iter-1)
                                                prenewpointgrad(2,1)=minima(2,iter-1)-sampling_step
                                        end if
                end if


!            step along the gradient
                if (delta_x .gt. 0 .and. abs(delta_x) .gt. tolerance) then
                       prenewpointgrad(1,1)=minima(1,iter-1)+sampling_step
                       ynew=y(prenewpointgrad(1,1),minima(1,iter-2),minima(2,iter-2),slope)
                       prenewpointgrad(2,1)=ynew
               else if (delta_x .lt. 0 .and. abs(delta_x) .gt. tolerance) then
                       prenewpointgrad(1,1)=minima(1,iter-1)-sampling_step
                       ynew=y(prenewpointgrad(1,1),minima(1,iter-2),minima(2,iter-2),slope)
                       prenewpointgrad(2,1)=ynew
               end if



                !step along the vector
                !(r(2,1) & r(1,1) y,x final state coordinates
                if (weighting .eq. 1) then

                if (DEBUG) then
                        print *, '...applying weighting'
                end if

                       delta_x_vec=r(1,1)-minima(1,iter-1)
                       delta_y_vec=r(2,1)-minima(2,iter-1)
                       slope_vec=delta_y_vec/delta_x_vec

                if (abs(r(1,1)-prenewpointgrad(1,1))<abs(r(1,1)-minima(1,iter-1))) then
                       vecpointgrad(1,1)=minima(1,iter-1)+wf*sampling_step
                       yvec=y(vecpointgrad(1,1),minima(1,iter-1),minima(2,iter-1),slope_vec)
                       vecpointgrad(2,1)=yvec
                else
                       vecpointgrad(1,1)=minima(1,iter-1)-wf*sampling_step
                       yvec=y(vecpointgrad(1,1),minima(1,iter-1),minima(2,iter-1),slope_vec)
                       vecpointgrad(2,1)=yvec
                end if

       !         print *, 'prenew', prenewpointgrad(1,1),prenewpointgrad(2,1)
       !         print *, 'vec', vecpointgrad(1,1),vecpointgrad(2,1)

                       newpoint(1,1)=(prenewpointgrad(1,1)+vecpointgrad(1,1))/2
                       newpoint(2,1)=(prenewpointgrad(2,1)+vecpointgrad(2,1))/2
               else if (weighting .eq. 0) then

                       newpoint(1,1)=prenewpointgrad(1,1)
                       newpoint(2,1)=prenewpointgrad(2,1)
               else

                if (DEBUG) then
                        print *, 'invalid weighting, weighting should be 1 or 0'
                end if

                end if

                if (DEBUG) then
                        print *, 'WARNING !!!! checking memory'
                end if

                !print *,'ynew',newpoint(2,1),'xnew',newpoint(1,1)

                 minima(1,iter)=newpoint(1,1)
                 minima(2,iter)=newpoint(2,1)
                 !print *,'new point',newpoint(1,1),newpoint(2,1)
                 !print *, firstpoint
                 deallocate (goodpoints)
                 cycle
        end if
minima(3,iter)=pre_minima

                if (DEBUG) then
                        print *, minima(1,iter),minima(2,iter),minima(3,iter),iter
                end if

                if (direction .eq. 1) then !forward write
write(33,*) minima(1,iter),minima(2,iter),minima(3,iter),iter
                else    !backward write
write(34,*) minima(1,iter),minima(2,iter),minima(3,iter),iter
                end if

deallocate (goodpoints)

if (abs(minima(1,iter)-r(1,1)) .lt. circle_radius  .and. &
      abs(minima(2,iter)-r(2,1)) .lt. circle_radius) then
                        print *,''
                        2005 format ('Convergence on iteration =', i5)
                        print 2005,iter
                        2006 format ('Weighting factor =', f5.2)
                        print 2006,wf
                        deallocate (fenergy)
                        goto 1717
end if
                        end do

!wf=wf+0.05
!if (wf .gt. 0.5) then
!        exit
!end if
!if (weighting .eq. 0) then
!stop
!end if
!if (testwork) then
!        stop
!end if

!                        print *,''
!                        print *,'Attempt',iwf
!                        print *,'Weighting factor =',wf
!                        iwf = iwf+1

                        if (direction .eq. 1 .and. iter-1 .eq. maxiters) then
                                                print *,''
                                                print *,'!!!!!!!! No Convergence in the Forward Process'
                                                print *,'!!!!!!!! Please change the weighting factor to find'
                                                print *,'!!!!!!!! another optimization'
                                                stop
                        else if(direction .eq. 2 .and. iter-1 .eq. maxiters) then
                          print *,''
                          print *,'No Convergence in the Backward Process'
                          print *,'Please change the weighting factor to find'
                          print *,'another optimization'
                          stop
                        end if


end do
1717 continue



close(33)
close(34)

open(unit=8,file='minima',status='old')
read(8,*) v
close(8)

r(2,1)=v(2,2) !backward-final-y & forward-initial-y
r(1,1)=v(1,2) !backward-final-x & forward-initial-x
r(2,2)=v(2,1) !backward-initial-y & forward-final-y
r(1,2)=v(1,1) !backward-initial-x & forward-final-x

!(v(2,1) & v(1,1) y,x final state coordinates
!(v(2,2) & v(1,2) initial state coordinates

open(unit=9,file=TRIM(currentpath)//'/sampling/forward-withwf',status='old')
line = 0
    do
    read(9,*,end=37)
       line = line +1
    end do
37 close(9)

allocate (forward(4,line))

open(unit=10,file=TRIM(currentpath)//'/sampling/forward-withwf',status='old')
read(10,*) forward
close(10)

open(unit=11,file=TRIM(currentpath)//'/sampling/backward-withwf',status='old')
line2 = 0
    do
    read(11,*,end=13)
       line2 = line2 +1
    end do
13 close(11)

allocate (backward(4,line2))

open(unit=11,file=TRIM(currentpath)//'/sampling/backward-withwf',status='old')
read(11,*) backward
close(11)

call veclengthfinder(forward,line,v,forward_out,outputdimension)
open(unit=47,file=TRIM(currentpath)//'/sampling/file47',status='replace')
do i=1,outputdimension
write (47,*) forward_out(1,i),forward_out(2,i),forward_out(3,i)
end do

call veclengthfinder(backward,line2,r,backward_out,outputbackdim)
open(unit=46,file=TRIM(currentpath)//'/sampling/file46',status='replace')
do i=1,outputbackdim
write (46,*) backward_out(1,i),backward_out(2,i),backward_out(3,i)
end do

!eliminating bad points

newdim=outputdimension*outputbackdim
allocate (pointsdiff(3,newdim))
m=1
do i=1,outputdimension
        do j=1,outputbackdim
                pointsdiffx=forward_out(1,i)-backward_out(1,j)
                pointsdiffy=forward_out(2,i)-backward_out(2,j)
                pointsdiff(3,m)=sqrt(pointsdiffx**2+pointsdiffy**2)
                pointsdiff(1,m)=i
                pointsdiff(2,m)=j
                m=m+1
        end do
end do

!do i=1,newdim
!        print *,pointsdiff(1,i),pointsdiff(2,i),pointsdiff(3,i)
!end do

minpointsdiff=pointsdiff(3,1)
do i=2,newdim
if (pointsdiff(3,i) .lt. minpointsdiff) then
        minpointsdiff = pointsdiff(3,i)
end if
end do
!print *,minpointsdiff

do i=1,newdim
if (minpointsdiff .eq. pointsdiff(3,i)) then
        forwardindex=int(pointsdiff(1,i))
        backwardindex=int(pointsdiff(2,i))
end if
end do
open(unit=49,file=TRIM(currentpath)//'/sampling/file49',status='replace')
open(unit=50,file=TRIM(currentpath)//'/sampling/file50',status='replace')
                write (49,*) backward_out(:,backwardindex)
                write (50,*) forward_out(:,forwardindex)
                backwardstartingindex=backwardindex
                forwardstartingindex=forwardindex
                checkparambackward=backward_out(3,outputbackdim)
                checkparamforward=forward_out(3,outputdimension)

        do
        call steepestpointfinder(forward_out,outputdimension,forwardindex,forwardstartingindex,minpoint)

        diffoptfirstx=forward_out(1,forwardstartingindex)-r(1,2)
        diffoptfirsty=forward_out(2,forwardstartingindex)-r(2,2)
        diffoptfirst=sqrt(diffoptfirstx**2+diffoptfirsty**2)

        diffoptsecx=forward_out(1,minpoint)-r(1,2)
        diffoptsecy=forward_out(2,minpoint)-r(2,2)
        diffoptsec=sqrt(diffoptsecx**2+diffoptsecy**2)

        if (diffoptsec .lt. diffoptfirst) then
                write (50,*) forward_out(:,minpoint)
                !print *,backward_out(:,minpoint+backwardstartingindex)
                forwardstartingindex=minpoint
                if (checkparamforward .eq. forward_out(3,minpoint)) then
                        exit
                end if
        else
                do j=1,outputdimension
                if (j .ge. minpoint) then
                        forward_out(:,j)=forward_out(:,j+1)
                end if
                end do
                outputdimension=outputdimension-1
                if (minpoint .lt. forwardstartingindex) then
                forwardstartingindex=forwardstartingindex-1
                else
                forwardstartingindex=forwardstartingindex
                end if
        end if
        end do
        close(50)
        rewind(50)

open(unit=20,file=TRIM(currentpath)//'/sampling/file50',status='old')
line2 = 0
    do
    read(20,*,end=25)
       line2 = line2 +1
    end do
25 close(20)

allocate (forwardsmoothing(3,line2))

open(unit=20,file=TRIM(currentpath)//'/sampling/file50',status='old')
read(20,*) forwardsmoothing
close(20)
open(unit=51,file=TRIM(currentpath)//'/sampling/file51',status='replace')
do i=line2,1,-1
        write(51,*) forwardsmoothing(:,i)
end do

        do
        call steepestpointfinder(backward_out,outputbackdim,backwardindex,backwardstartingindex,minpoint)

        diffoptfirstx=backward_out(1,backwardstartingindex)-r(1,1)
        diffoptfirsty=backward_out(2,backwardstartingindex)-r(2,1)
        diffoptfirst=sqrt(diffoptfirstx**2+diffoptfirsty**2)


        diffoptsecx=backward_out(1,minpoint)-r(1,1)
        diffoptsecy=backward_out(2,minpoint)-r(2,1)
        diffoptsec=sqrt(diffoptsecx**2+diffoptsecy**2)


        if (diffoptsec .lt. diffoptfirst) then
                write (51,*) backward_out(:,minpoint)
                !print *,backward_out(:,minpoint+backwardstartingindex)
                backwardstartingindex=minpoint
                if (checkparambackward .eq. backward_out(3,minpoint)) then
                        exit
                end if
        else
                do j=1,outputbackdim
                if (j .ge. minpoint) then
                        backward_out(:,j)=backward_out(:,j+1)
                end if
                end do
                outputbackdim=outputbackdim-1
                if (minpoint .lt. backwardstartingindex) then
                backwardstartingindex=backwardstartingindex-1
                else
                backwardstartingindex=backwardstartingindex
                end if

        end if
        end do
        close(51)


                !CHAIKINS SMOOTHING ALGORITHM

open(unit=22,file=TRIM(currentpath)//'/sampling/file51',status='old')
line2 = 0
    do
    read(22,*,end=26)
       line2 = line2 +1
    end do
26 close(22)

allocate (smoothingpts(3,line2))

open(unit=22,file=TRIM(currentpath)//'/sampling/file51',status='old')
read(22,*) smoothingpts
close(22)
sizeofnewpts=2*line2

open(unit=52,file=TRIM(currentpath)//'/sampling/file52',status='replace')
!first iteration
allocate (newpts(2,sizeofnewpts))
newpts(:,1)=smoothingpts(:,1)
j=1
do i=1,sizeofnewpts-3,2
newpts(:,i+1)=smoothingpts(:,j)+(smoothingpts(:,j+1)-smoothingpts(:,j))*leftparam
newpts(:,i+2)=smoothingpts(:,j)+(smoothingpts(:,j+1)-smoothingpts(:,j))*rightparam
j=j+1
end do
newpts(:,sizeofnewpts)=smoothingpts(:,line2)

do i=1,sizeofnewpts
write (52,*) newpts(:,i)
end do
open(unit=53,file=TRIM(currentpath)//'/results/smoothpath.dat',status='replace')
!more
do iter=1,5
sizeofnewpts=2*sizeofnewpts
allocate (dummynewpts(2,sizeofnewpts))
dummynewpts(:,1)=newpts(:,1)
        j=1
        do i=1,sizeofnewpts-3,2
        dummynewpts(:,i+1)=newpts(:,j)+(newpts(:,j+1)-newpts(:,j))*leftparam
        dummynewpts(:,i+2)=newpts(:,j)+(newpts(:,j+1)-newpts(:,j))*rightparam
        j=j+1
        end do
dummynewpts(:,sizeofnewpts)=newpts(:,int(sizeofnewpts*0.5))
deallocate (newpts)
allocate (newpts(2,sizeofnewpts))
newpts(:,:)=dummynewpts(:,:)
if (iter .eq. 4) then
do i=1,sizeofnewpts
write (53,*) dummynewpts(:,i)
end do
end if
deallocate (dummynewpts)
end do

print *,"###########################################"
print *,"Path smoothed with Chaikin's corner cutting"
print *,"###########################################"
end program LFEP
