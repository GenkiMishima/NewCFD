		subroutine J_boundary
!
		include 'common.h'
		include 'function.h'
!
		j=0
		q_g(1:4,1:ni-1,0) = q_g(1:4,1:ni-1,1)
		q_g(  3,1:ni-1,0) = - q_g(3,1:ni-1,1)
		q_p(1:4,1:ni-1,0) = q_p(1:4,1:ni-1,1)
		q_p(  3,1:ni-1,0) = - q_p(3,1:ni-1,1)
!
		j=nj
		q_g(1:4,1:ni-1,j) = q_g(1:4,1:ni-1,j-1)
		q_p(1:4,1:ni-1,j) = q_p(1:4,1:ni-1,j-1)

		do i=1, ni-1
			ux = fnu(i,j-1)
			uy = fnv(i,j-1)
			un= norm_j(1,i,j)*ux+norm_j(2,i,j)*uy
			ut=-norm_j(2,i,j)*ux+norm_j(1,i,j)*uy
			un = -un
			ux = norm_j(1,i,j)*un-norm_j(2,i,j)*ut
			uy = norm_j(2,i,j)*un+norm_j(1,i,j)*ut
			q_g(2,i,j) = q_g(1,i,j) * ux
			q_g(3,i,j) = q_g(1,i,j) * uy
!
			q_p(2,i,j) = q_p(2,i,j-1)
			q_p(3,i,j) = q_p(3,i,j-1)
		end do	
!
		return
		end