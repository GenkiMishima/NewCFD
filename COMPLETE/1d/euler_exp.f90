subroutine calc_next_step_exp(q,Flux,dt,dx)
   use prmtr
   implicit none
   integer i,j
   double precision,dimension(3,0:CellNum) :: q,Flux
   double precision dt,dx

   !$omp parallel do
   do i=1,CellNum-1
         q(:,i)=q(:,i)+dt/dx*(Flux(:,i  )-Flux(:,i+1))
   enddo
   !$omp end parallel do

   !do j=1,CellNum-1
   !   temp_residual1=temp_residual1+q(3,i) 
   !enddo
   !!residual
   !residual=abs(temp_residual1-temp_residual2)

end subroutine calc_next_step_exp
