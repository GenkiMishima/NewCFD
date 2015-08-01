		subroutine geomet
!
		include 'common.h'
!
! read grid data
		read(20,'(5e16.8)') ((x(i,j),r(i,j),i=1,ni), j=1,nj)
!
! calc rave_i and rave_j
!
		do j=1, nj-1
			rave_i(1:ni,j)=0.5*(r(1:ni,j)+r(1:ni,j+1))
		end do
		do i=1, ni-1
			rave_j(i,1:nj)=0.5*(r(i,1:nj)+r(i+1,1:nj))
		end do
!
! calc ds_i, ds_j
!	
		do j=1, nj-1
			ds_i(1:ni,j)=sqrt((x(1:ni,j)-x(1:ni,j+1))**2+(r(1:ni,j)-r(1:ni,j+1))**2)
		end do
		do i=1, ni-1
			ds_j(i,1:nj)=sqrt((x(i,1:nj)-x(i+1,1:nj))**2+(r(i,1:nj)-r(i+1,1:nj))**2)
		end do
!
! calc area & vol
!
		do j=1, nj-1
		do i=1, ni-1
			area_abd=0.5*(x(  i,  j)*(r(i+1,  j)-r(  i,j+1)) &
						 +x(i+1,  j)*(r(  i,j+1)-r(  i,  j)) &
						 +x(  i,j+1)*(r(  i,  j)-r(i+1,  j)) )
			area_bcd=0.5*(x(i+1,  j)*(r(i+1,j+1)-r(  i,j+1)) &
						 +x(i+1,j+1)*(r(  i,j+1)-r(i+1,  j)) &
						 +x(  i,j+1)*(r(i+1,  j)-r(i+1,j+1)) )
			area(i,j)=area_abd + area_bcd
			vol(i,j)=(r(  i,  j)+r(i+1,  j)+r(  i,j+1))/3. * area_abd &
					+(r(i+1,  j)+r(i+1,j+1)+r(  i,j+1))/3. * area_bcd
		end do
		end do
!
! calc norm_i, norm_j
!
		do i=1, ni
			do j=1, nj-1
				norm_i(1,i,j)= (r(i,j+1)-r(i,j))/ds_i(i,j)
				norm_i(2,i,j)=-(x(i,j+1)-x(i,j))/ds_i(i,j)
			end do
		end do
!
		do j=1, nj
			do i=1, ni-1
				norm_j(1,i,j)=-(r(i+1,j)-r(i,j))/ds_j(i,j)
				norm_j(2,i,j)= (x(i+1,j)-x(i,j))/ds_j(i,j)
			end do
		end do
! 
! initialize flag
		iflag(0:ni,0:nj)=0
		iflag(1:ni-1,1:nj-1)=1
!
		aburn=0.
		do j=1, nj-1
			aburn=aburn+rave_i(1,j)*ds_i(1,j)
		end do

		rmin=r(1,nj)
		do i=2, ni
			rmin=min(rmin,r(i,nj))
		end do
		athroat=0.5*rmin**2

		return
		end