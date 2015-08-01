subroutine calc_next_step_exp(q,Flux,Source,dt,dx,Vol)
   use prmtr
   implicit none
   integer i,j
   double precision,dimension(3,0:ni) :: q,Flux,Source
   double precision,dimension(  0:ni) :: Vol
   double precision dt,dx

   !$omp parallel do
   do i=1,ni-1
         q(:,i)=q(:,i)+dt/Vol(i)*(Flux(:,i  )-Flux(:,i+1))+dt*Source(:,i)
   enddo
   !$omp end parallel do

   !do j=1,ni-1
   !   temp_residual1=temp_residual1+q(3,i) 
   !enddo
   !!residual
   !residual=abs(temp_residual1-temp_residual2)

end subroutine calc_next_step_exp
