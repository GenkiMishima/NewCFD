subroutine set_jacobian
   use prmtr
   use variable
   implicit none
   integer i,j,m,n
   double precision rho, u, v, p, E, H, c
   double precision nuA, nuB
   double precision, parameter :: gamma_tilda=gamma-1.d0
   double precision, dimension(4,4) ::  A, B
   !$omp parallel do private(i,m,rho,u,v,p,E,H,c,&
   !$omp                     A,B)
   do j=1,nj-1
      do i=1,ni-1
         rho= w(1,i,j)
         u  = w(2,i,j)
         v  = w(3,i,j)
         p  = w(4,i,j)
         E  = p/gamma_tilda+0.5d0*rho*(u**2+v**2)
         H  = (E+p)/rho
         c  = sqrt(gamma*p/rho)
         nuA=abs(u)+c
         nuB=abs(v)+c

         A(1,1)= 0.d0 
         A(2,1)= 0.5d0*(u**2+v**2)*gamma_tilda-u**2 
         A(3,1)= -u*v
         A(4,1)= 0.5d0*u*(-2.d0*H+(u**2+v**2)*gamma_tilda)

         A(1,2)= 1.d0
         A(2,2)= u*(3.d0-gamma)
         A(3,2)= v
         A(4,2)= H-u**2*gamma_tilda

         A(1,3)= 0.d0
         A(2,3)= -v*gamma_tilda
         A(3,3)= u
         A(4,3)= -u*v*gamma_tilda

         A(1,4)= 0.d0
         A(2,4)= gamma_tilda
         A(3,4)= 0.d0
         A(4,4)= u*gamma

         B(1,1)= 0.d0 
         B(2,1)= -u*v
         B(3,1)= 0.5d0*(u**2+v**2)*gamma_tilda-v**2 
         B(4,1)= 0.5d0*v*(-2.d0*H+(u**2+v**2)*gamma_tilda)

         B(1,2)= 0.d0
         B(2,2)= v
         B(3,2)= -u*gamma_tilda
         B(4,2)= -u*v*gamma_tilda

         B(1,3)= 1.d0
         B(2,3)= u
         B(3,3)= v*(3.d0-gamma)
         B(4,3)= H-v**2*gamma_tilda

         B(1,4)= 0.d0
         B(2,4)= 0.d0
         B(3,4)= gamma_tilda
         B(4,4)= v*gamma

         Ap(:,:,i,j)=A(:,:)
         Am(:,:,i,j)=A(:,:)
         Bp(:,:,i,j)=B(:,:)
         Bm(:,:,i,j)=B(:,:)
         do m = 1,4
            Ap(m,m,i,j)=0.5d0*(Ap(m,m,i,j)+nuA)
            Am(m,m,i,j)=0.5d0*(Am(m,m,i,j)-nuA)
            Bp(m,m,i,j)=0.5d0*(Bp(m,m,i,j)+nuB)
            Bm(m,m,i,j)=0.5d0*(Bm(m,m,i,j)-nuB)
         end do
         alpha(i,j)=1.d0+dt(i,j)/dx(i,j)*nuA+dt(i,j)/dr(i,j)*nuB
      end do
   end do
   !$omp end parallel do
!   write(*,'(4es15.7)') Ap(:,:,50,41)
end subroutine set_jacobian
subroutine calc_next_step_imp
   use prmtr
   use variable
   implicit none
   integer i,j,m,n
   double precision, dimension(4,0:ni  ,0:nj  )::q_plime
   !$omp parallel do private(i,temp0)
   do j=1,nj-1
      do i=1,ni-1
         !if((i>=step_forward.and.i<=step_backward).and.j<=step_height)then
         !RHS(:,i,j)=q(:,i,j)
         !else
         temp0=1.d0/alpha(i,j)
         RHS(:,i,j)=temp0*(dt(i,j)/dx(i,j)*(X_Numerical(:,i,j)-X_Numerical(:,i+1,j))&
                          +dt(i,j)/dr(i,j)*(Y_Numerical(:,i,j)-Y_Numerical(:,i,j+1)))
         !end if
      end do
   end do
   !$omp end parallel do

   q_plime(:,0,:)=0.d0
   q_plime(:,:,0)=0.d0
   !$omp parallel do private(i,temp0)
   do j=1,nj-1
      do i=1,ni-1
         temp0=1.d0/alpha(i,j)
         q_plime(:,i,j)=RHS(:,i,j)&
                       +temp0*dt(i,j)/dx(i,j)*matmul(Ap(:,:,i-1,j),q_plime(:,i-1,j))&
                       +temp0*dt(i,j)/dr(i,j)*matmul(Bp(:,:,i,j-1),q_plime(:,i,j-1))
      end do
   end do
   !$omp end parallel do
   q_imp(:,ni,:)=0.d0
   q_imp(:,:,nj)=0.d0
   !$omp parallel do private(i,m,temp0)
   do j=1,nj-1
      n=nj-j
      do i=1,ni-1
         m=ni-i
         temp0=1.d0/alpha(m,n)
         q_imp(:,m,n)=q_plime(:,m,n)&
                     +temp0*dt(i,j)/dx(m,n)*matmul(Ap(:,:,m+1,n),q_imp(:,m+1,n))&
                     +temp0*dt(i,j)/dr(m,n)*matmul(Bp(:,:,m,n+1),q_imp(:,m,n+1))
      end do
   end do
   !$omp end parallel do

   !$omp parallel do private(i,temp0)
   do j=1,nj-1
      do i=1,ni-1
         if((i>step_forward.and.i<step_backward).and.j<step_height)then
         else
               q(:,i,j)=q(:,i,j)+q_imp(:,i,j)
         end if
      end do
   end do
   !$omp end parallel do
end subroutine calc_next_step_imp
