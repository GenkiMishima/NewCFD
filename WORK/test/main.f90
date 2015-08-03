!========================================================== 
!
!  This program is 2D implicit
!
!==========================================================
program main 
use prmtr
use variable
   implicit none
   integer i,j
   integer time_now, fl_num
   integer*4,external::access

   call set_breadth
   !Restart{{{
   if(access("restart.bin","r") .eq. 0) then
      open(15,file="restart.bin",form="unformatted")
         read(15) w
         read(15) time
         read(15) t
         read(15) dt
      close(15)
      open(15,file="check_ini.d")
         write(15,'(4es15.7)') temp0
      close(15)
      time_now=time+1
      print *,"Restart"
   else
      t=0.d0
      time_now=1
      call set_IC
      print *,"Initial Start"
   end if
   !}}}
   call set_geojac
   call set_conservative_variable_vector
   call set_BC
   call set_dt

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
   !print *, residual,time
      t=t+dt(1,1)
      call set_dt
      !call ROE_BODYFITTED
      call SLAU_FLUX
      call set_viscous
      !call calc_next_step_exp
      !==============================
      !!TimeIntegral========================================================
      !!call calc_next_step_exp(q,Flux,Source,dt,dx,Vol)
      !call calc_next_step_inp(q,w,Flux,Source,dt,dx,Vol,sonic,A)
      call calc_next_step_inp
      !==============================
      call setRESIDUAL
      call set_w
      call set_BC
      call output
      if(residual<epsilon)then
         !temp_int=1
         !call output
         write(34,*) t,residual
         write(*,*) 'convergence'
         call exit(2)
      end if
      write(34,*) dble(time),residual
   end do
   close(34)
end program main
