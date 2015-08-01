subroutine make_grid
   use prmtr
   use variable
   implicit none
   integer i,j
   open(55, file='../../grid/ideal_bl.x')
   read(55,*) temp0
   read(55,*) temp1,temp2,temp3
   if(temp0 .ne. 1)then
      print *,"grid_number is",temp0
      stop "No correct"
   elseif((ni.ne.temp1).or.(nj.ne.temp2))then
      print *,"ni   =",ni
      print *,"nj   =",nj
      print *,"grid_x_number",temp1
      print *,"grid_r_number",temp2
      stop "You must adjust ni to x_number or nj to r_number."
      call exit(1)
   endif
   read(55,*) ((xgrid(i,j),i=1,ni),j=1,nj)
   read(55,*) ((rgrid(i,j),i=1,ni),j=1,nj)
   close(55)
   do j=1,nj
     do i=1,ni
        x(i,j)=xgrid(i,j)
        r(i,j)=rgrid(i,j)
     end do
   end do

   do i=1,ni
      x(i,   0)=2d0*x(i,   1)-x(i,   2)
      r(i,   0)=2d0*r(i,   1)-r(i,   2)
      x(i,nj+1)=2d0*x(i,nj  )-x(i,nj-1)
      r(i,nj+1)=2d0*r(i,nj  )-r(i,nj-1)
   end do

   do j=0,nj+1
      x(   0,j)=2d0*x(   1,j)-x(   2,j)
      r(   0,j)=2d0*r(   1,j)-r(   2,j)
      x(ni+1,j)=2d0*x(ni  ,j)-x(ni-1,j)
      r(ni+1,j)=2d0*r(ni  ,j)-r(ni-1,j)
   end do
end subroutine make_grid

subroutine set_breadth
   use prmtr
   use variable
   implicit none
   integer i,j
   double precision xa, xb, xc, xd
   double precision ra, rb, rc, rd
   double precision area_abd,area_bcd
   double precision temp13,temp14
   double precision ,dimension(0:ni) :: tmp,tmp23
   temp0=1.d0
   call make_grid
!   !$omp parallel do private(i), reduction(min:temp0)
   !do j=1,nj-1
   !   do i=1,ni-1
   !      dx(i,j)=xgrid(i+1,j)-xgrid(i,j)
   !      dr(i,j)=rgrid(i,j+1)-rgrid(i,j)
   !      temp0=min(dx(i,j),temp0)
   !      temp0=min(dr(i,j),temp0)
   !!      Vol=dx*dy
   !   enddo
   !enddo
!   !$omp end parallel do
   open(44,file='check_delta1.d')
   write(44,'(2es15.7)') dx
   close(44)
   !open(44,file='x.d')
   !open(45,file='r.d')
   !!$omp parallel do private(i,xa,xb,xc,xd,ra,rb,rc,rd,&
   !!$omp                     area_abd,area_bcd,area,vol,&
   !!$omp                     delta_x, delta_r,&
   !!$omp                     ab, bc, cd,da, temp0)&
   !!$omp             shared(j, nvi, nvj,dx,dr)
   do j=0,nj
      do i=0,ni
         xa=x(i  ,j  )
         xb=x(i+1,j  )
         xc=x(i+1,j+1)
         xd=x(i  ,j+1)
         ra=r(i  ,j  )
         rb=r(i+1,j  )
         rc=r(i+1,j+1)
         rd=r(i  ,j+1)
         area_abd=0.5d0*(xa*(rb-rd)+xb*(rd-ra)+xd*(ra-rb))
         area_bcd=0.5d0*(xb*(rc-rd)+xc*(rd-rb)+xd*(rb-rc))
         Area(i,j)=Area_abd+Area_bcd
         Vol (i,j)=(ra+rb+rd)/3.d0*Area_abd+(rb+rc+rd)/3.d0*Area_bcd

         delta_x(1)=x(i+1,  j)-x(i  ,j  )
         delta_x(2)=x(i+1,j+1)-x(i+1,j  )
         delta_x(3)=x(i  ,j+1)-x(i+1,j+1)
         delta_x(4)=x(i  ,j  )-x(i  ,j+1)
         delta_r(1)=r(i+1,  j)-r(i  ,j  )
         delta_r(2)=r(i+1,j+1)-r(i+1,j  )
         delta_r(3)=r(i  ,j+1)-r(i+1,j+1)
         delta_r(4)=r(i  ,j  )-r(i  ,j+1)
         AB=sqrt(delta_x(1)**2+delta_r(1)**2)
         BC=sqrt(delta_x(2)**2+delta_r(2)**2)
         CD=sqrt(delta_x(3)**2+delta_r(3)**2)
         DA=sqrt(delta_x(4)**2+delta_r(4)**2)
         
   !      write(44,*)i,j,x(i,j)!delta_
   !      write(45,*)i,j,r(i,j)!delta_
         
         dsi(i,j) =AB
         dsj(i,j) =DA
         if(i==ni)then
            dsj(ni+1,j)=BC
         else if(j==nj)then
            dsi(i,nj+1)=CD
         end if

         nvi(1,i,j) =-delta_r(4)/DA
         nvi(2,i,j) = delta_x(4)/DA
         nvj(1,i,j) =-delta_r(1)/AB
         nvj(2,i,j) = delta_x(1)/AB
         !if(i==ni-1)then
         !   nvi(1,ni,j) = delta_r(2)/BC
         !   nvi(2,ni,j) =-delta_x(2)/BC
         !end if
         !if(j==nj-1)then
         !   nvj(1,i,nj) = delta_r(3)/CD
         !   nvj(2,i,nj) =-delta_x(3)/CD
         !end if
         !if(i==ni-1.and.j==nj-1)then
         !   nvi(1,ni,nj) = delta_r(2)/BC
         !   nvi(2,ni,nj) =-delta_x(2)/BC
         !   nvj(1,ni,nj) = delta_r(3)/CD
         !   nvj(2,ni,nj) =-delta_x(3)/CD
         !end if
         if(i==ni)then
            nvi(1,ni+1,j) =-delta_r(4)/DA
            nvi(2,ni+1,j) = delta_x(4)/DA
            nvj(1,ni+1,j) = delta_r(3)/CD
            nvj(2,ni+1,j) =-delta_x(3)/CD
         end if
         if(j==nj)then
            nvi(1,i,nj+1) = delta_r(2)/BC
            nvi(2,i,nj+1) =-delta_x(2)/BC
            nvj(1,i,nj+1) =-delta_r(1)/AB
            nvj(2,i,nj+1) = delta_x(1)/AB
         end if
         !if(i==ni-1.and.j==nj-1)then
         !   nvi(1,ni+1,nj+1) = delta_r(2)/BC
         !   nvi(2,ni+1,nj+1) =-delta_x(2)/BC
         !   nvj(1,ni+1,nj+1) = delta_r(3)/CD
         !   nvj(2,ni+1,nj+1) =-delta_x(3)/CD
         !end if

         !nvi(1,i,j) =-delta_r(4)/DA
         !nvi(2,i,j) = delta_x(4)/DA
         !!nvi(2,i,j) =-delta_x(4)/DA
         !nvj(1,i,j) =-delta_r(1)/AB
         !nvj(2,i,j) = delta_x(1)/AB
         !if(i==ni-1)then
         !nvi(1,ni,j) = delta_r(2)/BC
         !nvi(2,ni,j) = delta_x(2)/BC
         !dx (  ni,j) = abs(delta_x(3))
         !else if(j==nj-1)then
         !nvj(1,i,nj) =-delta_r(3)/CD
         !tmp(i)      =-delta_x(3)/CD
         !nvj(2,i,nj) =-delta_x(3)/CD
         !dr (  i,nj) = abs(delta_r(2))
         !end if
      enddo
   enddo
   !close(44)
   !close(45)
   !open(66,file='naca.d')
   !temp14=0.d0
   !do i=1,200
   !   temp14=tmp(i)+temp14
   !enddo
   !temp13=0.d0
   !do i=1,200
   !   temp13=tmp(i)+temp13
   !   write(66,*) i, temp13/temp14
   !enddo
   !close(66)
   !stop "OK"
   !!$omp end parallel do
   !nvi(1,ni-1,nj) = 1.d0
   !nvi(2,ni-1,nj) = 0.d0
   !nvi(1,ni  ,nj) = 1.d0
   !nvi(2,ni  ,nj) = 0.d0
   !nvj(1,ni-1,nj) = 0.d0
   !nvj(2,ni-1,nj) = 1.d0
   !nvj(1,ni  ,nj) = 0.d0
   !nvj(2,ni  ,nj) = 1.d0
   !do j = 1,nj
   !   do i = 1,ni
   !      print *,i,j
   !      print '(2es15.7)',nvj(:,i,j)
   !   end do
   !end do
   !do j=1,nj-1
   !   do i=1,ni
   !      dsi(i,j)=sqrt((dx(i  ,j+1)-dx(i,j))**2+(dr(i  ,j+1)-dr(i,j))**2)
   !   end do
   !end do
   !do j=1,nj
   !   do i=1,ni-1
   !      dsj(i,j)=sqrt((dx(i+1,j  )-dx(i,j))**2+(dr(i+1,j  )-dr(i,j))**2)
   !   end do
   !end do
   !open(44,file='check_nvi1.d')
   !open(45,file='check_nvi2.d')
   !open(46,file='check_nvj1.d')
   !open(47,file='check_nvj2.d')
   !do j=1,nj
   !   do i = 1,ni
   !      write(44,*)nvi(1,i,j),i,j
   !      write(45,*)nvi(2,i,j),i,j
   !      write(46,*)nvj(1,i,j),i,j
   !      write(47,*)nvj(2,i,j),i,j
   !   end do
   !end do
   !close(44)
   !close(45)
   !close(46)
   !close(47)
end subroutine set_breadth
