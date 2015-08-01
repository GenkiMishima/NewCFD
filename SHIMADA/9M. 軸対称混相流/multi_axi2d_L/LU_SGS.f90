		subroutine LU_SGS
!
		include 'common.h'
!
		dimension dedq_g(4,4), dfdq_g(4,4)
		dimension dedq_p(4,4), dfdq_p(4,4)
		dimension AP_g(4,4,ni,nj),AM_g(4,4,ni,nj),BP_g(4,4,ni,nj),BM_g(4,4,ni,nj)
		dimension AP_p(4,4,ni,nj),AM_p(4,4,ni,nj),BP_p(4,4,ni,nj),BM_p(4,4,ni,nj)
		dimension tmat_i(4,4), tmat_inv_i(4,4)
		dimension tmat_j(4,4), tmat_inv_j(4,4)
		dimension array(8,9)
!
		include 'function.h'
!
		ap_g=0.
		am_g=0.
		bp_g=0.
		bm_g=0.
!
		do j=1, nj-1
			do i=1, ni-1
!				rhs_g(1:4,i,j)=rhs_g(1:4,i,j)/alpha_g(i,j)
!				rhs_p(1:4,i,j)=rhs_p(1:4,i,j)/alpha_p(i,j)

				if(iflag(i,j).eq.1) then
				r_i = 0.5*(rave_i(i,j)+rave_i(i+1,j))
				r_j = 0.5*(rave_j(i,j)+rave_j(i,j+1))
				enx_i=0.5*(norm_i(1,i,j)+norm_i(1,i+1,j))
				enr_i=0.5*(norm_i(2,i,j)+norm_i(2,i+1,j))
				enx_j=0.5*(norm_j(1,i,j)+norm_j(1,i,j+1))
				enr_j=0.5*(norm_j(2,i,j)+norm_j(2,i,j+1))
!
				tmat_i=0.
				tmat_i(1,1)=1.
				tmat_i(4,4)=1.
				tmat_i(2,2)= enx_i
				tmat_i(2,3)=-enr_i
				tmat_i(3,2)= enr_i
				tmat_i(3,3)= enx_i
!
				tmat_j=0.
				tmat_j(1,1)=1.
				tmat_j(4,4)=1.
				tmat_j(2,2)= enx_j
				tmat_j(2,3)=-enr_j
				tmat_j(3,2)= enr_j
				tmat_j(3,3)= enx_j
!
				tmat_inv_i=0.
				tmat_inv_i(1,1)=1.
				tmat_inv_i(4,4)=1.
				tmat_inv_i(2,2)= enx_i
				tmat_inv_i(2,3)= enr_i
				tmat_inv_i(3,2)=-enr_i
				tmat_inv_i(3,3)= enx_i
!
				tmat_inv_j=0.
				tmat_inv_j(1,1)=1.
				tmat_inv_j(4,4)=1.
				tmat_inv_j(2,2)= enx_j
				tmat_inv_j(2,3)= enr_j
				tmat_inv_j(3,2)=-enr_j
				tmat_inv_j(3,3)= enx_j
!
				sos=fna(i,j)
				u=fnu(i,j)
				v=fnv(i,j)
				theta=u**2+v**2
				H=sos**2/gm+0.5*theta
                un_i = enx_i * u + enr_i * v
				un_j = enx_j * u + enr_j * v
				ut_i =-enr_i * u + enx_i * v
				ut_j =-enr_j * u + enx_j * v

				up=fnu_p(i,j)
				vp=fnv_p(i,j)
				ep=q_p(4,i,j)/q_p(1,i,j)
                un_p_i = enx_i * up + enr_i * vp
				un_p_j = enx_j * up + enr_j * vp
				ut_p_i =-enr_i * up + enx_i * vp
				ut_p_j =-enr_j * up + enx_j * vp

!
				dedq_g(1,1)=0.
				dedq_g(1,2)=1.
				dedq_g(1,3)=0.
				dedq_g(1,4)=0.
				dedq_g(2,1)=0.5*theta*gm-un_i**2
				dedq_g(2,2)=un_i*(3.-gam)
				dedq_g(2,3)=-ut_i*gm
				dedq_g(2,4)=gm
				dedq_g(3,1)=-un_i*ut_i
				dedq_g(3,2)=ut_i
				dedq_g(3,3)=un_i
				dedq_g(3,4)=0.
				dedq_g(4,1)=0.5*un_i*(-2.*H+theta*gm)
				dedq_g(4,2)=H-un_i*un_i*gm
				dedq_g(4,3)=-un_i*ut_i*gm
				dedq_g(4,4)=un_i*gam
!
				dedq_p(1,1)=0.
				dedq_p(1,2)=1.
				dedq_p(1,3)=0.
				dedq_p(1,4)=0.
				dedq_p(2,1)=-un_p_i**2
				dedq_p(2,2)=2.*un_p_i
				dedq_p(2,3)=0.
				dedq_p(2,4)=0.
				dedq_p(3,1)=-un_p_i*ut_p_i
				dedq_p(3,2)=ut_p_i
				dedq_p(3,3)=un_p_i
				dedq_p(3,4)=0.
				dedq_p(4,1)=-un_p_i*ep
				dedq_p(4,2)=ep
				dedq_p(4,3)=0.
				dedq_p(4,4)=un_p_i
!
				dfdq_g(1,1)=0.
				dfdq_g(1,2)=1.
				dfdq_g(1,3)=0.
				dfdq_g(1,4)=0.
				dfdq_g(2,1)=0.5*theta*gm-un_j**2
				dfdq_g(2,2)=un_j*(3.-gam)
				dfdq_g(2,3)=-ut_j*gm
				dfdq_g(2,4)=gm
				dfdq_g(3,1)=-un_j*ut_j
				dfdq_g(3,2)=ut_j
				dfdq_g(3,3)=un_j
				dfdq_g(3,4)=0.
				dfdq_g(4,1)=0.5*un_j*(-2.*H+theta*gm)
				dfdq_g(4,2)=H-un_j*un_j*gm
				dfdq_g(4,3)=-un_j*ut_j*gm
				dfdq_g(4,4)=un_j*gam
!
				dfdq_p(1,1)=0.
				dfdq_p(1,2)=1.
				dfdq_p(1,3)=0.
				dfdq_p(1,4)=0.
				dfdq_p(2,1)=-un_p_j**2
				dfdq_p(2,2)=2.*un_p_j
				dfdq_p(2,3)=0.
				dfdq_p(2,4)=0.
				dfdq_p(3,1)=-un_p_j*ut_p_j
				dfdq_p(3,2)=ut_p_j
				dfdq_p(3,3)=un_p_j
				dfdq_p(3,4)=0.
				dfdq_p(4,1)=-un_p_j*ep
				dfdq_p(4,2)=ep
				dfdq_p(4,3)=0.
				dfdq_p(4,4)=un_p_j
!
				dedq_g(1:4,1:4)=r_i*matmul(tmat_i(1:4,1:4),matmul(dedq_g(1:4,1:4),tmat_inv_i(1:4,1:4)))
				dfdq_g(1:4,1:4)=r_j*matmul(tmat_j(1:4,1:4),matmul(dfdq_g(1:4,1:4),tmat_inv_j(1:4,1:4)))
				dedq_p(1:4,1:4)=r_i*matmul(tmat_i(1:4,1:4),matmul(dedq_p(1:4,1:4),tmat_inv_i(1:4,1:4)))
				dfdq_p(1:4,1:4)=r_j*matmul(tmat_j(1:4,1:4),matmul(dfdq_p(1:4,1:4),tmat_inv_j(1:4,1:4)))
!
				AP_g(1:4,1:4,i,j)=0.5*dedq_g(1:4,1:4)
				AM_g(1:4,1:4,i,j)=0.5*dedq_g(1:4,1:4)
				BP_g(1:4,1:4,i,j)=0.5*dfdq_g(1:4,1:4)
				BM_g(1:4,1:4,i,j)=0.5*dfdq_g(1:4,1:4)
!
				AP_p(1:4,1:4,i,j)=0.5*dedq_p(1:4,1:4)
				AM_p(1:4,1:4,i,j)=0.5*dedq_p(1:4,1:4)
				BP_p(1:4,1:4,i,j)=0.5*dfdq_p(1:4,1:4)
				BM_p(1:4,1:4,i,j)=0.5*dfdq_p(1:4,1:4)
!
				do k=1,4
					AP_g(k,k,i,j)=AP_g(k,k,i,j)+0.5*nu_a_g(i,j)
					AM_g(k,k,i,j)=AM_g(k,k,i,j)-0.5*nu_a_g(i,j)
					BP_g(k,k,i,j)=BP_g(k,k,i,j)+0.5*nu_b_g(i,j)
					BM_g(k,k,i,j)=BM_g(k,k,i,j)-0.5*nu_b_g(i,j)

					AP_p(k,k,i,j)=AP_p(k,k,i,j)+0.5*nu_a_p(i,j)
					AM_p(k,k,i,j)=AM_p(k,k,i,j)-0.5*nu_a_p(i,j)
					BP_p(k,k,i,j)=BP_p(k,k,i,j)+0.5*nu_b_p(i,j)
					BM_p(k,k,i,j)=BM_p(k,k,i,j)-0.5*nu_b_p(i,j)
				end do
				end if
			end do
		end do
!		
! Forward Sweep
!

		do j=1, nj-1
			do i=1, ni-1
				if(j.ne.1) then
					dsj=0.5*(ds_j(i,j-1)+ds_j(i,j))
					rhs_g(1:4,i,j)=rhs_g(1:4,i,j)+dt(i,j)/vol(i,j)/alpha_g(i,j)*dsj*matmul(BP_g(1:4,1:4,i,j-1),rhs_g(1:4,i,j-1))
					rhs_p(1:4,i,j)=rhs_p(1:4,i,j)+dt(i,j)/vol(i,j)/alpha_p(i,j)*dsj*matmul(BP_p(1:4,1:4,i,j-1),rhs_p(1:4,i,j-1))
				end if
				if(i.ne.1) then
					dsi=0.5*(ds_i(i-1,j)+ds_i(i,j))
					rhs_g(1:4,i,j)=rhs_g(1:4,i,j)+dt(i,j)/vol(i,j)/alpha_g(i,j)*dsi*matmul(AP_g(1:4,1:4,i-1,j),rhs_g(1:4,i-1,j))
					rhs_p(1:4,i,j)=rhs_p(1:4,i,j)+dt(i,j)/vol(i,j)/alpha_p(i,j)*dsi*matmul(AP_p(1:4,1:4,i-1,j),rhs_p(1:4,i-1,j))
				end if

!---------------
			if(j.eq.1) then
				dsj=0.5*(ds_j(i,1)+ds_j(i,2))
				fact_g=1.+dt(i,j)/vol(i,j)/alpha_g(i,j)*dsj*nu_b_g(i,j)
				fact_p=1.+dt(i,j)/vol(i,j)/alpha_p(i,j)*dsj*nu_b_p(i,j)
				rhs_g(1:4,i,j)=rhs_g(1:4,i,j)/fact_g
				rhs_p(1:4,i,j)=rhs_p(1:4,i,j)/fact_p
			end if
			end do
!------------------
		end do
!
! backward sweep
!
		do j=nj-1,1,-1
			do i=ni-1,1,-1
				if(j.ne.nj-1) then
					dsj=0.5*(ds_j(i,j+1)+ds_j(i,j))
					rhs_g(1:4,i,j)=rhs_g(1:4,i,j)-dt(i,j)/vol(i,j)/alpha_g(i,j)*dsj*matmul(BM_g(1:4,1:4,i,j+1),rhs_g(1:4,i,j+1))
					rhs_p(1:4,i,j)=rhs_p(1:4,i,j)-dt(i,j)/vol(i,j)/alpha_p(i,j)*dsj*matmul(BM_p(1:4,1:4,i,j+1),rhs_p(1:4,i,j+1))
				end if
				if(i.ne.ni-1) then
					dsi=0.5*(ds_i(i+1,j)+ds_i(i,j))
					rhs_g(1:4,i,j)=rhs_g(1:4,i,j)-dt(i,j)/vol(i,j)/alpha_g(i,j)*dsi*matmul(AM_g(1:4,1:4,i+1,j),rhs_g(1:4,i+1,j))
					rhs_p(1:4,i,j)=rhs_p(1:4,i,j)-dt(i,j)/vol(i,j)/alpha_p(i,j)*dsi*matmul(AM_p(1:4,1:4,i+1,j),rhs_p(1:4,i+1,j))
				end if
			end do
		end do
!
		return
		end