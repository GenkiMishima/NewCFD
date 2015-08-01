module variable
   use prmtr
   integer i,j,time
   double precision temp0,temp1,temp2,temp3
   double precision, dimension(1:ni,1:nj):: xgrid, ygrid
   double precision, dimension(1:ni,1:nj)::rho_mat,u_mat,v_mat,p_mat,M_mat,T_mat,vecx,vecr
   character*20 tmpstring
end module variable
