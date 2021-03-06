!========================================================== 
!
!  This program is Q1D implicit
!
!==========================================================
program main 
use prmtr
!use variable
   implicit none
   integer i,j
   integer time_now, fl_num,time,Center
   integer*4,external::access
   double precision dt,residual,t,Resi,dx
   double precision,parameter :: tL   = 35d2
   double precision,parameter :: pL   = 1d7
   double precision,parameter :: uL   = 0d0
   double precision,parameter :: tR   = 35d2
   double precision,parameter :: pR   = 1d5
   double precision,parameter :: uR   = 0d0
   double precision,parameter :: xmax = 4d0
   double precision,parameter :: xmin = 0d0
   double precision rR,rL
   double precision temp0,temp1,temp2,temp3
   double precision,dimension(  0:ni) :: x,y,a,PreResi,Vol,sonic
   double precision,dimension(3,0:ni) :: w,q,Flux,Source
   rL = pL/(gas_spec*tL)
   rR = pR/(gas_spec*tR)
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
   !Quasi Term===================================
   open(40,file='Area.txt')
   read(40,'(3e16.8)') (x(i),y(i),a(i),i=1,ni)
   close(40)
   dx = (xmax-xmin)/float(ni-1)
   !open(45,file='nozzle.d')
   !do i=0,ni
   !   write(45,*) x(i),A(i)
   !end do
   !close(45)

   !call setIC(q,w,pL,uL,TL,pR,uR,tR,sonic,x,Vol,PreResi)
   do i=0,ni
      temp0 = (rR-rL)/float(ni-1)
      temp1 = (uR-uL)/float(ni-1)
      temp2 = (pR-pL)/float(ni-1)
      w(1,i)= rL+temp0*float(i-1)
      w(2,i)= uL+temp1*float(i-1)
      w(3,i)= pL+temp2*float(i-1)
      q(1,i)=w(1,i)
      q(2,i)=w(2,i)*q(1,i)
      q(3,i)=w(3,i)/(gamma-1.d0)+0.5d0*w(1,i)*(w(2,i)**2)
      temp0 = x(i+1) + 2.2d0/3d0*(x(i+1)-1.5d0)**3
      temp1 = x(i  ) + 2.2d0/3d0*(x(i  )-1.5d0)**3
      Vol(i) = temp0-temp1
      sonic(i) = sqrt(gamma*w(3,i)/w(1,i))
   end do
   PreResi(:) = q(3,:)

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
   !open(35,file='check.d')
   do time=time_now, time_max
   !do time=1,2000
      t=t+dt

      call setdt(w,dx,dt)

      !Flux================================================================
      call SLAU_FLUX(w,Flux)
      do i = 0,ni-1
         Flux(:,i) = Flux(:,i)*A(i)
         Source(2,i) = w(3,i)*(A(i+1)-A(i))/Vol(i)
         !print *,i, Flux(:,i)
      end do

      !TimeIntegral========================================================
      !call calc_next_step_exp(q,Flux,Source,dt,dx,Vol)
      call calc_next_step_inp(q,w,Flux,Source,dt,dx,Vol,sonic,A)

      !RESIDUAL===================================================
      call setRESIDUAL(q,PreResi,Resi)

      !PrimitiveValue=======================================================
      call setPrimitiveValue(w,q,sonic,t)

      !call set_BC=================================================================
      call setBC(w,pL,TL,pR,tR)

      !Vector===============================================================
      call setConservativeVector(w,q,PreResi)

      !OUTPUT==================================================================
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
   !close(35)
end program main
