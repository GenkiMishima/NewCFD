		subroutine out_field(job)
!
		include 'common.h'
		dimension dens(ni,nj)
		dimension dens_p(ni,nj), temp_p(ni,nj)
		dimension temp(ni,nj)
		include 'function.h'
!
		if(job.eq.2) then
			rewind 50
			write(50) time, ic
			write(50) q_g
			write(50) q_p
		else
!!
		do j=1, nj
		do i=1, ni
			isum=iflag(i-1,j-1)+iflag(i,j-1)+iflag(i,j)+iflag(i-1,j)
			if(isum.ne.0) then
				temp_p(i,j)=(iflag(i-1,j-1)*fnt_p(i-1,j-1)+iflag(  i,j-1)*fnt_p(  i,j-1) &
						    +iflag(  i,  j)*fnt_p(  i,  j)+iflag(i-1,  j)*fnt_p(i-1,  j)) &
						/ float(isum)
				temp(i,j)=(iflag(i-1,j-1)*fnt(i-1,j-1)+iflag(  i,j-1)*fnt(  i,j-1) &
						  +iflag(  i,  j)*fnt(  i,  j)+iflag(i-1,  j)*fnt(i-1,  j)) &
						/ float(isum)
				dens(i,j)=(iflag(i-1,j-1)*fnr(i-1,j-1)+iflag(  i,j-1)*fnr(  i,j-1) &
						  +iflag(  i,  j)*fnr(  i,  j)+iflag(i-1,  j)*fnr(i-1,  j)) &
						/ float(isum)
!				velo(1,i,j)=(iflag(i-1,j-1)*fnu(i-1,j-1)+iflag(  i,j-1)*fnu(  i,j-1) &
!						  +iflag(  i,  j)*fnu(  i,  j)+iflag(i-1,  j)*fnu(i-1,  j)) &
!						/ float(isum)
!				velo(2,i,j)=(iflag(i-1,j-1)*fnv(i-1,j-1)+iflag(  i,j-1)*fnv(  i,j-1) &
!						  +iflag(  i,  j)*fnv(  i,  j)+iflag(i-1,  j)*fnv(i-1,  j)) &
!						/ float(isum)
!				pres(i,j)=(iflag(i-1,j-1)*fnp(i-1,j-1)+iflag(  i,j-1)*fnp(  i,j-1) &
!						  +iflag(  i,  j)*fnp(  i,  j)+iflag(i-1,  j)*fnp(i-1,  j)) &
!						/ float(isum)
!				amach(i,j)=(iflag(i-1,j-1)*am(i-1,j-1)+iflag(  i,j-1)*am(  i,j-1) &
!						   +iflag(  i,  j)*am(  i,  j)+iflag(i-1,  j)*am(i-1,  j)) &
!						/ float(isum)
				dens_p(i,j)=(iflag(i-1,j-1)*fnr_p(i-1,j-1)+iflag(  i,j-1)*fnr_p(  i,j-1) &
						  +iflag(  i,  j)*fnr_p(  i,  j)+iflag(i-1,  j)*fnr_p(i-1,  j)) &
						/ float(isum)
			else
				temp_p(i,j)=0.
				dens(i,j)=0.
				temp(i,j)=0.
!				velo(1,i,j)=0.
!				velo(2,i,j)=0.
!				amach(i,j)=0.
				dens_p(i,j)=0.
			end if
		end do
		end do
!
!		write(51,'(i6)') ic
!		write(52,'(f16.8)') time
!		write(53,'(f16.8)') time
!		write(54,'(f16.8)') time
!		write(55,'(f16.8)') time*1.e3
		write(51,'(5e16.8)') ((dens(i,j),i=1,ni),j=1,nj)
		write(52,'(5e16.8)') ((temp(i,j),i=1,ni),j=1,nj)
!		write(53,'(5e16.8)') ((velo(2,i,j),i=1,ni),j=1,nj)
		write(54,'(5e16.8)') ((dens_p(i,j),i=1,ni),j=1,nj)
		write(55,'(5e16.8)') ((temp_p(i,j),i=1,ni),j=1,nj)
		end if
!
		return
		end