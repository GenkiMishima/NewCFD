subroutine SLAU_FLUX
use prmtr
use variable
implicit none
integer i,j

!   !X_MUSCL{{{
   w_left( :,:,:)=w(:,:,:)
   w_right(:,:,:)=w(:,:,:)
   !$omp parallel do private(i)
   do j=1, nj
      do i=1, ni-1
            call MUSCLlimiter(w(:,i-1,j), w(:,i,j), w(:,i+1,j), w_left(:,i,j), w_right(:,i,j))
      enddo
   enddo
   !$omp end parallel do
   !MUSCL}}}
   do j=1, nj-1
      do i=1, ni
         call SLAU(w_right(:,i-1,j),w_left(:,i,j),X_Numerical(:,i,j),nvi(:,i,j))
      end do
   end do

!   !Y_MUSCL{{{
   w_down(:,:,:)=w(:,:,:)
   w_up(  :,:,:)=w(:,:,:)
   !$omp parallel do private(i)
   do j=1, nj-1
      do i=1, ni
            call MUSCLlimiter(w(:,i,j-1), w(:,i,j), w(:,i,j+1), w_down(:,i,j), w_up(:,i,j))
      enddo
   enddo
   !$omp end parallel do
   !}}}
   do j=1, nj
      do i=1, ni-1
         call SLAU(w_up(:,i,j-1),w_down(:,i,j),Y_Numerical(:,i,j),nvj(:,i,j))
      end do
   end do
end subroutine SLAU_FLUX
