
program main
use prmtr
use grid
implicit none
   integer i,j,time
   double precision, dimension(1:ni,1:nj)::xgrid,ygrid
   double precision, dimension(1:ni,1:nj)::rho_mat,u_mat,v_mat,p_mat
   character*20 tmpstring
   call make_grid(xgrid(:,:),ygrid(:,:))
   do time=0,49
   write(tmpstring,'(i2.2)') time
   open(66,file='../data/density_'//trim(tmpstring)//'.d')
   read(66,*) ((rho_mat(i,j),i=1,ni),j=1,nj)
   close(66)
   open(66,file='../data/U_Velocity_'//trim(tmpstring)//'.d')
   read(66,*) ((u_mat(i,j),i=1,ni),j=1,nj)
   close(66)
   open(66,file='../data/V_Velocity_'//trim(tmpstring)//'.d')
   read(66,*) ((v_mat(i,j),i=1,ni),j=1,nj)
   close(66)
   open(66,file='../data/pressure_'//trim(tmpstring)//'.d')
   read(66,*) ((p_mat(i,j),i=1,ni),j=1,nj)
   close(66)
   open(66,file='result/result_'//trim(tmpstring)//'.vtk')
   
   call vtk_grid(xgrid(:,:),ygrid(:,:))
   call vtk_scalar(rho_mat,"Density")
   call vtk_scalar(u_mat,"U_Velocity")
   call vtk_scalar(v_mat,"V_Velocity")
   call vtk_scalar(p_mat,"Pressure")
   close(66)
   enddo
end program main

subroutine vtk_grid(xgrid,ygrid)
use prmtr
use grid
   implicit none
   integer i,j
   double precision,dimension(ni,nj)::xgrid,ygrid
   write(66,'(a26)') '# vtk DataFile Version 2.0'
   write(66,'(a7)') '2D Data'
   write(66,'(a5)') 'ASCII'
   write(66,'(a23)') 'DATASET STRUCTURED_GRID'
   write(66,'(a11,i3,a1,i3,a1,i1)') 'DIMENSIONS ', ni, ' ', nj, ' ', 1


   write(66,'(a17)') 'FIELD FiledData 2'

   write(66,'(a15)') 'TIME 1 1 double'
   write(66,'(e14.6e3)') 0.d0
   write(66,'(a13)') 'CYCLE 1 1 int'
   write(66,'(i0)') 0


   write(66,'(a7,i7,a6)') 'POINTS ', (ni)*(nj)*1, ' float'
   do j=1,nj
      do i=1,ni
         write(66,'(e14.6e3,2(1x,e14.6e3))') xgrid(i,j), ygrid(i,j), 0.d0
      enddo
   enddo
   write(66,*)

   write(66,'(a11,i7)') 'POINT_DATA ', (ni)*(nj)*1
   write(66,'(a22)') 'VECTORS Indices float'
   do j = 1, nj
      do i = 1, ni
         write(66,'(e14.6e3,2(1x,e14.6e3))') dble(i), dble(j), 0.d0
      enddo
   enddo
   write(66,*)

end subroutine vtk_grid
subroutine vtk_scalar(arr,title)
use prmtr
   implicit none
   double precision, dimension(ni,nj), intent(in)::arr
   character(*),intent(in)::title

   integer i,j

   write(66,'(a)') 'SCALARS '//trim(title)//' float 1'
   write(66, '(a20)') 'LOOKUP_TABLE default'
   do j = 1, nj
      do i = 1, ni
         write(66,'(e14.6e3)') arr(i,j)
      enddo
   enddo
   write(66,*)
end subroutine vtk_scalar
