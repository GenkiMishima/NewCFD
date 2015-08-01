subroutine calc_next_step_exp
   use prmtr
   use variable
   implicit none
   integer i,j
   double precision, dimension(4)::vis_term
   temp0=10000.d0
   temp_residual2=temp_residual1
   temp1=0.d0
   !$omp parallel do private(i) shared(j,q)
   do j=1,nj-1
      do i=1,ni-1
            q(:,i,j)=q(:,i,j)+dt(i,j)/area(i,j)*(dsj(i  ,j  )*(X_Numerical(:,i  ,j  )-vis_i(:,i  ,j  ))&
                                                -dsj(i+1,j  )*(X_Numerical(:,i+1,j  )-vis_i(:,i+1,j  ))&
                                                +dsi(i  ,j  )*(Y_Numerical(:,i  ,j  )-vis_j(:,i  ,j  ))&
                                                -dsi(i  ,j+1)*(Y_Numerical(:,i  ,j+1)-vis_j(:,i  ,j+1)))
      enddo
   enddo
   !$omp end parallel do

   do j=1,nj-1
      do i=1,ni-1
         temp_residual1=temp_residual1+q(4,i,j) 
      enddo
   enddo
   !residual
   residual=abs(temp_residual1-temp_residual2)

end subroutine calc_next_step_exp
