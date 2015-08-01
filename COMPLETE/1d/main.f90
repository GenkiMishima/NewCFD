!========================================================== 
!
!  This program is 2D implicit
!
!==========================================================
program main 
use prmtr
!use variable
   implicit none
   integer i,j,ni
   integer time_now, fl_num,time,Center
   integer*4,external::access
   double precision dt,residual,t,Resi,dx
   double precision temp0,temp1,temp2,temp3
   double precision,dimension(  0:CellNum) :: x,PreResi
   double precision,dimension(3,0:CellNum) :: w,q,Flux
   ni = CellNum

   !call set_breadth
   !!Restart{{{
   !if(access("restart.bin","r") .eq. 0) then
   !   open(15,file="restart.bin",form="unformatted")
   !      read(15) w
   !      read(15) time
   !      read(15) t
   !      read(15) dt
   !   close(15)
   !   open(15,file="check_ini.d")
   !      write(15,'(4es15.7)') temp0
   !   close(15)
   !   time_now=time+1
   !   print *,"Restart"
   !else
      t=0.d0
      dt=0.d0
      time_now=1
   !   call set_IC
   !   print *,"Initial Start"
   !end if
   !!}}}
   !call set_geojac
   !call set_conservative_variable_vector
   !call set_BC
   !call set_dt
   Center=int(dble(CellNum/2d0))

   w(1,0:Center) = 5d0
   w(2,0:Center) = 0d0
   w(3,0:Center) = 5d0
   w(1,Center:CellNum) = 1d0
   w(2,Center:CellNum) = 0d0
   w(3,Center:CellNum) = 1d0
   do i=0,ni-1
      q(1,i)=w(1,i)
      q(2,i)=w(2,i)*q(1,i)
      q(3,i)=w(3,i)/(gamma-1.d0)+0.5d0*w(1,i)*(w(2,i)**2)
      x(i) = i
   enddo
   dx = 1d-2

   open(44,file='data/Condition.txt')
   write(44,*) 'Condition'
   write(44,*) 'CFL     =',CFL
   write(44,*) 'MUSCL   =',phi
   write(44,*) 'ORDER   =',kappa
   write(44,*) 'OUTPUT  =',out_time
   close(44)
   write(*,*) 'PROCESSING' 
   write(*,*) 'CFL     =',CFL
   write(*,*) 'MUSCL   =',phi
   write(*,*) 'ORDER   =',kappa
   write(*,*) 'OUTPUT  =',out_time
   write(*,*) 'RESIDUAL=',epsilon
   print *,""
   !main loop
   open(34,file='residual.d')
   do time=time_now, time_max
   !do time=1,2000
      t=t+dt
      temp0=1.d300
      !dt
      do i=1,CellNum
         temp3=sqrt(gamma*w(3,i)/w(1,i))
         temp1=w(2,i)
         temp1=dx/(temp1+temp3)
         temp0=min(temp0,temp1)
      enddo
      dt=CFL*temp0
      PreResi(:) = q(3,:)
      !call ROE_BODYFITTED
      call SLAU_FLUX(w,Flux)
      call calc_next_step_exp(q,Flux,dt,dx)
      temp2 =0d0
      do i=1,CellNum
         temp1 =(q(3,i)-PreResi(i))**2
         temp2 = temp1+temp2
      enddo
      Resi = temp2
      !call set_w
      do i=1,ni-1
         w(1,i)=q(1,i)
         temp0=1d0/w(1,i)
         w(2,i)=q(2,i)*temp0
         w(3,i)=(gamma-1.d0)*(q(3,i)-0.5d0*w(1,i)*(w(2,i)**2))
         if(w(3,i)<0.d0)then
            write(*,*) i, t
            write(*,'(4es15.7)') w(:,i)
            call exit(1)
         endif
      enddo
      w(:,0)=w(:,1)
      w(:,CellNum)=w(:,CellNum-1)
      !call set_BC
      call output(x,w,time,t,dt,Resi)
      if(Resi<epsilon)then
         !temp_int=1
         !call output
         write(34,*) t,Resi
         write(*,*) 'convergence'
         call exit(2)
      end if
      write(34,*) dble(time),Resi
   end do
   close(34)
end program main
