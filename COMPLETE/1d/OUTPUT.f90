subroutine output(x,w,time,t,dt,Resi)
   use prmtr
!   use variable
   implicit none
   integer i,j
   integer time
   double precision t,dt,Resi
   double precision,dimension(0:CellNum):: Mach,x
   double precision,dimension(3,1:CellNum) :: w
   character*20 tmpstring
   if(time==1.or.mod(time,out_time)== 0) then
      do i=1,CellNum
         Mach(i)=sqrt((w(2,i)**2)*w(1,i)/(gamma*w(3,i)))
      enddo
      !$omp end parallel do

      write(tmpstring,'(i3.3)') int(time/out_time)
      open(10, file='data/density_'//trim(tmpstring)//'.d')
      open(11, file='data/U_Velocity_'//trim(tmpstring)//'.d')
      open(12, file='data/pressure_'//trim(tmpstring)//'.d')
      open(14, file='data/Mach_'//trim(tmpstring)//'.d')
!      !$omp parallel do private(i)
      do i=1,CellNum
         write(10, *)x(i), w(1,i)
         write(11, *)x(i), w(2,i)
         write(12, *)x(i), w(3,i)
         write(14, *)x(i), Mach(i)
      enddo
      close(10)
      close(11)
      close(12)
      close(14)
      write(*, *) tmpString, t, Resi
      open(15,file="restart.bin",form="unformatted")
         write(15) w
         write(15) time
         write(15) t
         write(15) dt
      close(15)
   end if
end subroutine output
