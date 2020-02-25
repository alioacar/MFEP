module veclength

        implicit none

contains

subroutine veclengthfinder(input,line,refpoints,output,outputdimension)

        implicit none
        real,intent(in),dimension(2,3) :: refpoints
        real,dimension(:,:),allocatable :: difference
        real,intent(out),dimension(:,:),allocatable :: output
        real :: veclengthtofinalx,veclengthtofinaly,veclengthtofinal
        real :: veclengthtoinitialx,veclengthtoinitialy,veclengthtoinitial
        real :: pre_minima
        integer,intent(in) :: line
        real,intent(in),dimension(4,line) :: input
        integer :: i,new_i,j,k,output_line
        integer,intent(out) :: outputdimension

        allocate (difference(1,line))

        do i=1,line
        veclengthtofinalx=input(1,i)-refpoints(1,1)
        veclengthtofinaly=input(2,i)-refpoints(2,1)
        veclengthtofinal=sqrt(veclengthtofinalx**2+veclengthtofinaly**2)

        veclengthtoinitialx=input(1,i)-refpoints(1,2)
        veclengthtoinitialy=input(2,i)-refpoints(2,2)
        veclengthtoinitial=sqrt(veclengthtoinitialx**2+veclengthtoinitialy**2)

        difference(1,i)=abs(veclengthtofinal-veclengthtoinitial)
        end do

        pre_minima=difference(1,1)
        do i=2,line
                if (difference(1,i)<pre_minima) pre_minima=difference(1,i)
        end do


do i=1,line
        if (pre_minima .eq. difference(1,i)) then
                new_i=i
                veclengthtofinalx=input(1,new_i)-refpoints(1,1)
                veclengthtofinaly=input(2,new_i)-refpoints(2,1)
                veclengthtofinal=sqrt(veclengthtofinalx**2+veclengthtofinaly**2)
                outputdimension=line-new_i+1
                k=0
                do j=new_i,line
                if (sqrt((input(1,j)-refpoints(1,1))**2+(input(2,j)-refpoints(2,1))**2) .le. veclengthtofinal) then
                        k=k+1 !counting goodpoints
                end if
                end do
                outputdimension=k
                allocate (output(3,outputdimension))

                output_line=1
                do j=new_i,line
                if (sqrt((input(1,j)-refpoints(1,1))**2+(input(2,j)-refpoints(2,1))**2) .le. veclengthtofinal) then
                        output(1,output_line)=input(1,j)
                        output(2,output_line)=input(2,j)
                        output(3,output_line)=input(3,j)
                        output_line=output_line+1
                end if
                end do

        end if
end do

end subroutine veclengthfinder
end module veclength
