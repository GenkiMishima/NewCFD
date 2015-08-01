		fnr(i,j)=q_g(1,i,j)
		fnu(i,j)=q_g(2,i,j)/q_g(1,i,j)
		fnv(i,j)=q_g(3,i,j)/q_g(1,i,j)
		fnp(i,j)=gm*(q_g(4,i,j)-0.5*(q_g(2,i,j)**2+q_g(3,i,j)**2)/q_g(1,i,j))
		fna(i,j)=sqrt(gam*fnp(i,j)/q_g(1,i,j))
		fnt(i,j)=fnp(i,j)/(gasc*fnr(i,j))

		fnr_p(i,j)=q_p(1,i,j)
		fnu_p(i,j)=q_p(2,i,j)/q_p(1,i,j)
		fnv_p(i,j)=q_p(3,i,j)/q_p(1,i,j)
		fnt_p(i,j)=(q_p(4,i,j)/q_p(1,i,j)-0.5*(fnu_p(i,j)**2+fnv_p(i,j)**2))/cc

