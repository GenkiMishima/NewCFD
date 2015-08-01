subroutine make_grid
   use prmtr
   use variable
   implicit none
!   integer i,j
!   double precision temp0,temp1,temp2,temp3
!   double precision, dimension(1:ni,1:nj), intent(out) :: xgrid, ygrid
   
   open(55, file='../../../grid/ideal_bl.x')
   read(55,*) temp0
   read(55,*) temp1,temp2,temp3
   if(temp0 .ne. 1)then
      print *,"temp0=",temp0
      stop "No correct"
   elseif((ni.ne.temp1).or.(nj.ne.temp2))then
      print *,"ni   =",ni
      print *,"nj   =",nj
      print *,"temp1=",temp1
      print *,"temp2=",temp2
      stop "You must adjust ni to temp1 or nj to temp2."
      call exit(1)
   endif
   read(55,*) ((xgrid(i,j),i=1,ni),j=1,nj)
   read(55,*) ((ygrid(i,j),i=1,ni),j=1,nj)
   close(55)
end subroutine make_grid

!subroutine set_breadth
!   use prmtr
!   use variable
!   implicit none
!   temp0=1.d0
!   call make_grid
!   dx(:,:)=1.d0/(ni-1)
!   dy(:,:)=1.d0/(nj-1)
!   open(44,file='grid_check.d')
!   !$omp parallel do private(i), reduction(min:temp0)
!   do j=1,nj-1
!      do i=1,ni-1
!         dx(i,j)=xgrid(i+1,j)-xgrid(i,j)
!         dy(i,j)=ygrid(i,j+1)-ygrid(i,j)
!         write(44,*) dx(i,j), dy(i,j)
!         temp0=min(dx(i,j),temp0)
!         temp0=min(dy(i,j),temp0)
!      enddo
!   enddo
!   !$omp end parallel do
!   close(44)
!end subroutine set_breadth
