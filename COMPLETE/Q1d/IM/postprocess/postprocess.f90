program main
use prmtr
use variable
implicit none
   integer time_ini
   integer*4,external::access
   call make_grid
   open(66,file='result/condition.txt')
   write(66,*) 'start'
   close(66)
   if(access("restart.bin","r").eq.0)then
   open(15,file="restart.bin",form="unformatted")
      read(15) time
   close(15)
   time_ini=time+1
   else
   time_ini=1
   end if
   do time=time_ini,50000
   write(tmpstring,'(i3.3)') time
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
   open(66,file='../data/Mach_'//trim(tmpstring)//'.d')
   read(66,*) ((M_mat(i,j),i=1,ni),j=1,nj)
   close(66)
   open(66,file='../data/Vectorx_'//trim(tmpstring)//'.d')
   read(66,*) ((vecx(i,j),i=1,ni),j=1,nj)
   close(66)
   open(66,file='../data/Vectorr_'//trim(tmpstring)//'.d')
   read(66,*) ((vecr(i,j),i=1,ni),j=1,nj)
   close(66)
   open(66,file='result/result_'//trim(tmpstring)//'.vtk')
   T_mat(:,:)=p_mat(:,:)/(R_gas*rho_mat(:,:))
   call vtk_grid(xgrid(:,:),ygrid(:,:))
   call vtk_scalar(rho_mat,"Density")
   call vtk_scalar(u_mat,"U_Velocity")
   call vtk_scalar(v_mat,"V_Velocity")
   call vtk_scalar(p_mat,"Pressure")
   call vtk_scalar(M_mat,"Mach")
   call vtk_scalar(T_mat,"Temperture")
   call vtk_vector(u_mat,v_mat,"Velocity")
   call vtk_vector(vecx,vecr,"Vector")
   close(66)
   open(66,file='check.d')
   do j = 1,50
      write(66,*)(dble(j)-1)/1.2d0, u_mat(105,j)/100.d0
   end do
   close(66)
   write(*,*) time
   open(15,file="restart.bin",form="unformatted")
      write(15) time
   close(15)
   enddo
end program main

subroutine vtk_grid(xgrid,ygrid)
use prmtr
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
subroutine vtk_vector(arr1,arr2,title)
use prmtr
   implicit none
   double precision,dimension(ni,nj),intent(in)::arr1,arr2
   character(*),intent(in)::title

   integer i,j

   write(66,'(a22)') 'VECTORS '//trim(title)//' float'
   do j = 1, nj
     do i = 1, ni
       write(66,'(e14.6e3,2(1x,e14.6e3))') arr1(i,j), arr2(i,j), 0.d0
     end do
   end do
   write(66,*)
end subroutine vtk_vector
 
