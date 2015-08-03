subroutine calc_next_step_rk_1
   use prmtr
   use variable
   implicit none
   integer i,j
   double precision, dimension(4)::vis_term
   !double precision temp0, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9
   temp0=10000.d0
   temp1=0.d0
   !$omp parallel do private(i) shared(j,q)
   do j=1,nj-1
      do i=1,ni-1
         !if(i<=step_backward.and.j<=step_height)then
         !   q(:,i,j)=q(:,i,j)
         !   else
            q(:,i,j)=qp(:,i,j)+0.5d0*dt(i,j)/area(i,j)&
                    *(dsj(i  ,j  )*(X_Numerical(:,i  ,j  )-vis_i(:,i  ,j  ))&
                     -dsj(i+1,j  )*(X_Numerical(:,i+1,j  )-vis_i(:,i+1,j  ))&
                     +dsi(i  ,j  )*(Y_Numerical(:,i  ,j  )-vis_j(:,i  ,j  ))&
                     -dsi(i  ,j+1)*(Y_Numerical(:,i  ,j+1)-vis_j(:,i  ,j+1)))
         !   end if
      enddo
   enddo
   !$omp end parallel do

end subroutine calc_next_step_rk_1
subroutine calc_next_step_rk_2
   use prmtr
   use variable
   implicit none
   integer i,j
   double precision, dimension(4)::vis_term
   !double precision temp0, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9
   temp0=10000.d0
   temp1=0.d0
   !$omp parallel do private(i) shared(j,q)
   do j=1,nj-1
      do i=1,ni-1
         !if(i<=step_backward.and.j<=step_height)then
         !   q(:,i,j)=q(:,i,j)
         !   else
            q(:,i,j)=qp(:,i,j)+dt(i,j)/area(i,j)*(dsj(i  ,j  )*(X_Numerical(:,i  ,j  )-vis_i(:,i  ,j  ))&
                                                 -dsj(i+1,j  )*(X_Numerical(:,i+1,j  )-vis_i(:,i+1,j  ))&
                                                 +dsi(i  ,j  )*(Y_Numerical(:,i  ,j  )-vis_j(:,i  ,j  ))&
                                                 -dsi(i  ,j+1)*(Y_Numerical(:,i  ,j+1)-vis_j(:,i  ,j+1)))
         !   end if
      enddo
   enddo
   !$omp end parallel do

   !temp_residual1(:,:)=abs(temp_residual2(:,:)-q(4,:,:))**2

   !do j=1,nj-1
   !   do i=1,ni-1
   !      temp0=temp0+temp_residual1(i,j)
   !   enddo
   !enddo
   !!residual
   !residual=sqrt(temp0)
   !temp_residual2(:,:)=q(4,:,:)

end subroutine calc_next_step_rk_2
