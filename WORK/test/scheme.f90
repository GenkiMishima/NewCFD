
subroutine SLAU(wl,wr,flux,nv)!{{{
use prmtr
implicit none
   double precision, intent(in) , dimension(2) ::nv
   double precision, intent(in) , dimension(4) ::wl,wr
   double precision, intent(out), dimension(4) ::flux


   double precision ul,vl,rhol,gmmal,unl,phil(4),al,Ml,betal,pl
   double precision ur,vr,rhor,gmmar,unr,phir(4),ar,Mr,betar,pr

   double precision N(4)
   double precision Drho,Dp,a_bar,un_bar_abs,p_bar
   double precision g,gp,gm,chi,M_tilde,dm,un_bar_abs_p,un_bar_abs_m
   double precision p_tilde,temp0

   integer k

   !calc wlr
   rhol =wl(1)
   ul   =wl(2)
   vl   =wl(3)
   pl   =wl(4)
   phil(1)=1.d0
   phil(2)=ul
   phil(3)=vl
   temp0=pl/(gamma-1d0)+0.5d0*rhol*(ul**2+vl**2)
   phil(4)=(temp0+pl)/rhol !enthalpy
   al=sqrt(gamma*pl/rhol)
   unl=nv(1)*ul+nv(2)*vl

   rhor =wr(1)
   ur   =wr(2)
   vr   =wr(3)
   pr   =wr(4)
   phir(1)=1
   phir(2)=ur
   phir(3)=vr
   temp0=pr/(gamma-1d0)+0.5d0*rhor*(ur**2+vr**2)
   phir(4)=(temp0+pr)/rhor!enthalpy
   ar=sqrt(gamma*pr/rhor)
   unr=nv(1)*ur+nv(2)*vr

   ! '("phir")'
   ! '(es15.7)',phir
   ! '("phil")'
   ! '(es15.7)',phil

   !set average and difference
   dp  =pr-pl
   a_bar =0.5d0*(ar +al )
   p_bar =0.5d0*(pr +pl )

   ! '("dp")'
   ! '(es15.7)',dp
   ! '("a_bar")'
   ! '(es15.7)',a_bar
   ! '("p_bar")'
   ! '(es15.7)',p_bar

   ml=unl/a_bar
   mr=unr/a_bar

   !set dm
   m_tilde=min(1d0,1d0/a_bar*sqrt(0.5d0*(ul**2+vl**2+ur**2+vr**2)))
   chi=(1d0-m_tilde)**2

   gp=-max(min(ml,0d0),-1d0)
   gm= min(max(mr,0d0), 1d0)
   g=gp*gm

   un_bar_abs  =(abs(unr)*rhor+abs(unl)*rhol)/(rhor+rhol)
   un_bar_abs_p=(1d0-g)*un_bar_abs+g*abs(unl)
   un_bar_abs_m=(1d0-g)*un_bar_abs+g*abs(unr)

   dm=0.5d0*(rhol*(unl+un_bar_abs_p)+rhor*(unr-un_bar_abs_m)-chi/a_bar*dp)
   ! '("dm")'
   ! '(es15.7)',dm

  !set p_tilde
   n(:)=0d0
   n(2)=nv(1)
   n(3)=nv(2)

   if(abs(mr)<1d0) then
      betar=0.25d0*(2d0+mr)*(mr-1d0)**2
   else
      betar=0.5d0*(1d0-sign(1d0,mr))
   end if

   if(abs(ml)<1d0) then
      betal=0.25d0*(2d0-ml)*(ml+1d0)**2
   else
      betal=0.5d0*(1d0+sign(1d0,ml))
   end if

   p_tilde=p_bar+0.5d0*(betal-betar)*(pl-pr)+(1d0-chi)*(betar+betal-1d0)*p_bar

   ! '("p_tilde")'
   ! '(es15.7)',p_tilde
   !calculate tg
   flux=0.5d0*(dm+abs(dm))*phil&
       +0.5d0*(dm-abs(dm))*phir&
       +p_tilde*n
!    '("flux")'
!    '(4es15.7)',flux
       
end subroutine SLAU!}}}

subroutine ROE_BODYFITTED!{{{
use prmtr
use variable
implicit none
integer i,j
!   corner_for (:)=w(:,step_forward+1,step_height)
!   corner_back(:)=w(:,step_backward ,step_height)
!   !X_MUSCL{{{
!   !$omp parallel do
!   do j = 0, nj
!      if(j<=step_height)then
!         w(:,step_forward+1,j)= w(:,step_forward,j)
!         w(2,step_forward+1,j)=-w(2,step_forward,j)
!         w(:,step_backward,j)= w(:,step_backward+1,j)
!         w(2,step_backward,j)=-w(2,step_backward+1,j)
!      endif
!   enddo
!   !$omp end parallel do
   w_left( :,:,:)=w(:,:,:)
   w_right(:,:,:)=w(:,:,:)
   !$omp parallel do private(i)
   do j=1, nj
      do i=1, ni-1
         !if((i>step_forward.and.i<step_backward).and.j<step_height)then
         !else
            call MUSCLlimiter(w(:,i-1,j), w(:,i,j), w(:,i+1,j), w_left(:,i,j), w_right(:,i,j))
         !endif
      enddo
   enddo
   !$omp end parallel do
   !MUSCL}}}
   !X_Direction{{{
   !$omp parallel do shared(X_Numerical,j),&
   !$omp             private(i,rhol,ul,vl,pl,El,Hl,&
   !$omp                     rhor,ur,vr,pr,Er,Hr,&
   !$omp                     temp0,temp1,temp2,temp3,&
   !$omp                     temp4,uaver,vaver,Haver,&
   !$omp                     aaver,lambda,X_matrix,X_inverse,&
   !$omp                     dq,temp_arrey1,temp_arrey2,&
   !$omp                     convective_vector)
   do j=1, nj-1
      do i=1, ni
         rhol=w_right (1,i-1,j)
         ul = nvi(1,i  ,j)*w_right (2,i-1,j)+nvi(2,i  ,j)*w_right (3,i-1,j)
         vl =-nvi(2,i  ,j)*w_right (2,i-1,j)+nvi(1,i  ,j)*w_right (3,i-1,j)
         pl  =w_right (4,i-1,j)
         El  =pl/gamma_bar+0.5d0*rhol*(ul**2+vl**2)
         Hl  =(El+pl)/rhol
 
         rhor=w_left(1,i  ,j)
         ur  = nvi(1,i  ,j)*w_left(2,i  ,j)+nvi(2,i  ,j)*w_left(3,i  ,j)
         vr  =-nvi(2,i  ,j)*w_left(2,i  ,j)+nvi(1,i  ,j)*w_left(3,i  ,j)
         pr  =w_left(4,i  ,j)
         Er  =pr/gamma_bar+0.5d0*rhor*(ur**2+vr**2)
         Hr  =(Er+pr)/rhor
 
         temp0=sqrt(rhol)
         temp1=sqrt(rhor)
         temp2=1.d0/(temp0+temp1)
         temp3=temp0*temp2
         temp4=temp1*temp2
 
         uaver=temp3*ul+temp4*ur
         vaver=temp3*vl+temp4*vr
         Haver=temp3*Hl+temp4*Hr
         aaver=sqrt(gamma_bar*(Haver-0.5d0*(uaver**2+vaver**2)))
 
         lambda(:,:)=0.d0
         lambda(1,1)=abs(uaver)                 
         lambda(2,2)=abs(uaver)
         lambda(3,3)=abs(uaver-aaver) 
         lambda(4,4)=abs(uaver+aaver) 
         
         X_matrix(1,1)=0.d0
         X_matrix(2,1)=0.d0
         X_matrix(3,1)=1.d0
         X_matrix(4,1)=vaver
 
         X_matrix(1,2)=1.d0
         X_matrix(2,2)=uaver
         X_matrix(3,2)=vaver
         X_matrix(4,2)=0.5d0*(uaver**2+vaver**2)
 
         X_matrix(1,3)=1.d0
         X_matrix(2,3)=uaver-aaver
         X_matrix(3,3)=vaver
         X_matrix(4,3)=Haver-aaver*uaver
 
         X_matrix(1,4)=1.d0
         X_matrix(2,4)=uaver+aaver
         X_matrix(3,4)=vaver
         X_matrix(4,4)=Haver+aaver*uaver
 
         temp3=aaver**2
         temp0=2.d0*temp3
         temp1=gamma_bar/temp0
         temp2=1.d0/temp0
 
         X_inverse(1,1)=-vaver
         X_inverse(2,1)=temp1*(-2.d0*Haver+4.d0*temp3/gamma_bar)
         X_inverse(3,1)=temp1*(Haver+(aaver*uaver-temp3)/gamma_bar)
         X_inverse(4,1)=temp1*(Haver-(aaver*uaver+temp3)/gamma_bar)
 
         X_inverse(1,2)=0.d0
         X_inverse(2,2)=temp1*(2.d0*uaver)
         X_inverse(3,2)=temp1*(-uaver-aaver/gamma_bar)
         X_inverse(4,2)=temp1*(-uaver+aaver/gamma_bar)
 
         X_inverse(1,3)=1.d0
         X_inverse(2,3)=temp1*(2.d0*vaver)
         X_inverse(3,3)=temp1*(-vaver)
         X_inverse(4,3)=temp1*(-vaver)
 
         X_inverse(1,4)=0.d0
         X_inverse(2,4)=temp1*(-2.d0)
         X_inverse(3,4)=temp1
         X_inverse(4,4)=temp1
         
         dq(1)=   rhor-   rhol
         dq(2)=rhor*ur-rhol*ul
         dq(3)=rhor*vr-rhol*vl
         dq(4)=     Er-     El
 
         temp_arrey1(:)      =matmul(X_inverse(:,:),         dq(:))
         temp_arrey2(:)      =matmul(   lambda(:,:),temp_arrey1(:))
         convective_vector(:)=matmul( X_matrix(:,:),temp_arrey2(:))
 
         X_Numerical(1,i,j)=      rhol*ul+      rhor*ur
         X_Numerical(2,i,j)=rhol*ul**2+pl+rhor*ur**2+pr
         X_Numerical(3,i,j)=   rhol*ul*vl+   rhor*ur*vr
         X_Numerical(4,i,j)=   rhol*ul*Hl+   rhor*ur*Hr
 
         x_Numerical(:,i,j)=0.5d0*(X_Numerical(:,i,j)-convective_vector(:))
      enddo
   enddo
   !$omp end parallel do
   do j=1,nj-1
      do i=1,ni
         !T Matrix
         temp2=X_Numerical(2,i,j)
         temp3=X_Numerical(3,i,j)
         X_Numerical(2,i,j)= nvi(1,i,j)*temp2-nvi(2,i,j)*temp3
         X_Numerical(3,i,j)= nvi(2,i,j)*temp2+nvi(1,i,j)*temp3
      end do
   end do
   !X_Direction}}}
 
!   w(:,step_forward+1,step_height)=corner_for (:)
!   w(:,step_backward ,step_height)=corner_back(:)
!   !Y_MUSCL{{{
!   !$omp parallel do
!   do i = 0, ni
!      if(i>step_forward.and.i<=step_backward)then
!         w(:,i,step_height)= w(:,i,step_height+1)
!         w(3,i,step_height)=-w(3,i,step_height+1)
!      endif
!   enddo
!   !$omp end parallel do
   w_down(:,:,:)=w(:,:,:)
   w_up(  :,:,:)=w(:,:,:)
   !$omp parallel do private(i)
   do j=1, nj-1
      do i=1, ni
         !if((i>step_forward.and.i<step_backward).and.j<step_height)then
         !else
            call MUSCLlimiter(w(:,i,j-1), w(:,i,j), w(:,i,j+1), w_down(:,i,j), w_up(:,i,j))
         !endif
      enddo
   enddo
   !$omp end parallel do
   !}}}
   !Y_Direction{{{
   !$omp parallel do shared (Y_Numerical,j),&
   !$omp             private(i,rhol,ul,vl,pl,El,Hl,&
   !$omp                     rhor,ur,vr,pr,Er,Hr,&
   !$omp                     temp0,temp1,temp2,temp3,&
   !$omp                     temp4,uaver,vaver,Haver,&
   !$omp                     aaver,lambda,X_matrix,X_inverse,&
   !$omp                     dq,temp_arrey1,temp_arrey2,&
   !$omp                     convective_vector)
   do j=1, nj
      do i=1, ni-1
         rhol=w_up(1,i,j-1)
         ul  = nvj(1,i  ,j)*w_up(2,i,j-1)+nvj(2,i  ,j)*w_up(3,i,j-1)
         vl  =-nvj(2,i  ,j)*w_up(2,i,j-1)+nvj(1,i  ,j)*w_up(3,i,j-1)
         pl  =w_up(4,i,j-1)
         El  =pl/gamma_bar+0.5d0*rhol*(ul**2+vl**2)
         Hl  =(El+pl)/rhol
 
         rhor=w_down  (1,i,j  )
         ur  = nvj(1,i  ,j)*w_down  (2,i  ,j)+nvj(2,i  ,j)*w_down  (3,i  ,j)
         vr  =-nvj(2,i  ,j)*w_down  (2,i  ,j)+nvj(1,i  ,j)*w_down  (3,i  ,j)
         pr  =w_down  (4,i,j  )
         Er  =pr/gamma_bar+0.5d0*rhor*(ur**2+vr**2)            
         Hr  =(Er+pr)/rhor
 
         temp0=sqrt(rhol)
         temp1=sqrt(rhor)
         temp2=1.d0/(temp0+temp1)
         temp3=temp0*temp2
         temp4=temp1*temp2
 
         uaver=temp3*ul+temp4*ur
         vaver=temp3*vl+temp4*vr
         Haver=temp3*Hl+temp4*Hr
         aaver=sqrt(gamma_bar*(Haver-0.5d0*(uaver**2+vaver**2)))
 
         lambda(:,:)=0.d0
         lambda(1,1)=abs(uaver)                 
         lambda(2,2)=abs(uaver)
         lambda(3,3)=abs(uaver-aaver) 
         lambda(4,4)=abs(uaver+aaver) 
         
         X_matrix(1,1)=0.d0
         X_matrix(2,1)=0.d0
         X_matrix(3,1)=1.d0
         X_matrix(4,1)=vaver
 
         X_matrix(1,2)=1.d0
         X_matrix(2,2)=uaver
         X_matrix(3,2)=vaver
         X_matrix(4,2)=0.5d0*(uaver**2+vaver**2)
 
         X_matrix(1,3)=1.d0
         X_matrix(2,3)=uaver-aaver
         X_matrix(3,3)=vaver
         X_matrix(4,3)=Haver-aaver*uaver
 
         X_matrix(1,4)=1.d0
         X_matrix(2,4)=uaver+aaver
         X_matrix(3,4)=vaver
         X_matrix(4,4)=Haver+aaver*uaver
 
         temp3=aaver**2
         temp0=2.d0*temp3
         temp1=gamma_bar/temp0
         temp2=1.d0/temp0
 
         X_inverse(1,1)=-vaver
         X_inverse(2,1)=temp1*(-2.d0*Haver+4.d0*temp3/gamma_bar)
         X_inverse(3,1)=temp1*(Haver+(aaver*uaver-temp3)/gamma_bar)
         X_inverse(4,1)=temp1*(Haver-(aaver*uaver+temp3)/gamma_bar)
 
         X_inverse(1,2)=0.d0
         X_inverse(2,2)=temp1*(2.d0*uaver)
         X_inverse(3,2)=temp1*(-uaver-aaver/gamma_bar)
         X_inverse(4,2)=temp1*(-uaver+aaver/gamma_bar)
 
         X_inverse(1,3)=1.d0
         X_inverse(2,3)=temp1*(2.d0*vaver)
         X_inverse(3,3)=temp1*(-vaver)
         X_inverse(4,3)=temp1*(-vaver)
 
         X_inverse(1,4)=0.d0
         X_inverse(2,4)=temp1*(-2.d0)
         X_inverse(3,4)=temp1
         X_inverse(4,4)=temp1
         
         dq(1)=   rhor-   rhol
         dq(2)=rhor*ur-rhol*ul
         dq(3)=rhor*vr-rhol*vl
         dq(4)=     Er-     El
 
         temp_arrey1(:)      =matmul(X_inverse(:,:),         dq(:))
         temp_arrey2(:)      =matmul(   lambda(:,:),temp_arrey1(:))
         convective_vector(:)=matmul( X_matrix(:,:),temp_arrey2(:))
 
         Y_Numerical(1,i,j)=      rhol*ul+      rhor*ur
         Y_Numerical(2,i,j)=rhol*ul**2+pl+rhor*ur**2+pr
         Y_Numerical(3,i,j)=   rhol*ul*vl+   rhor*ur*vr
         Y_Numerical(4,i,j)=   rhol*ul*Hl+   rhor*ur*Hr
 
         Y_Numerical(:,i,j)=0.5d0*(Y_Numerical(:,i,j)-convective_vector(:))
      enddo
   enddo
   !$omp end parallel do
 
   do j=1,nj
      do i=1,ni-1
         !T Matrix
         temp2=Y_Numerical(2,i,j)
         temp3=Y_Numerical(3,i,j)
         Y_Numerical(2,i,j)= nvj(1,i,j)*temp2-nvj(2,i,j)*temp3
         Y_Numerical(3,i,j)= nvj(2,i,j)*temp2+nvj(1,i,j)*temp3
      enddo
   enddo
   !Y_Direction}}}
 
end subroutine ROE_BODYFITTED!}}}
subroutine SLAU_FLUX
use prmtr
use variable
implicit none
integer i,j

!   !X_MUSCL{{{
   w_left( :,:,:)=w(:,:,:)
   w_right(:,:,:)=w(:,:,:)
   !$omp parallel do private(i)
   do j=1, nj
      do i=1, ni-1
            call MUSCLlimiter(w(:,i-1,j), w(:,i,j), w(:,i+1,j), w_left(:,i,j), w_right(:,i,j))
      enddo
   enddo
   !$omp end parallel do
   !MUSCL}}}
   do j=1, nj-1
      do i=1, ni
         call SLAU(w_right(:,i-1,j),w_left(:,i,j),X_Numerical(:,i,j),nvi(:,i,j))
      end do
   end do

!   !Y_MUSCL{{{
   w_down(:,:,:)=w(:,:,:)
   w_up(  :,:,:)=w(:,:,:)
   !$omp parallel do private(i)
   do j=1, nj-1
      do i=1, ni
            call MUSCLlimiter(w(:,i,j-1), w(:,i,j), w(:,i,j+1), w_down(:,i,j), w_up(:,i,j))
      enddo
   enddo
   !$omp end parallel do
   !}}}
   do j=1, nj
      do i=1, ni-1
         call SLAU(w_up(:,i,j-1),w_down(:,i,j),Y_Numerical(:,i,j),nvj(:,i,j))
      end do
   end do
end subroutine SLAU_FLUX
