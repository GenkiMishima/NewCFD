
subroutine SLAU(rho_minus, rho_plus, u_minus, u_plus, v_minus, v_plus, p_minus, p_plus&
            &, H_minus, H_plus, flux, hou)
!{{{
use prmtr
implicit none
   integer, intent(in) :: hou
   double precision, intent(in) ::rho_minus, rho_plus, u_minus, u_plus, v_minus, v_plus, p_minus, p_plus&
            &, H_minus, H_plus
   double precision, intent(out), dimension(4) ::flux
   double precision, dimension(4) :: psi_minus, psi_plus, normal
   double precision pressure_flux, mass_flux, Mach_hat, X, Mach_minus, Mach_plus, g
   double precision sonic_plus, sonic_minus, sonic_ave, Vel_plus, Vel_minus
   double precision Veln_ave, Vel_ave_plus, Vel_ave_minus, beta_plus, beta_minus
   double precision temp0, temp1

   psi_minus(1)=1.d0
   psi_minus(2)=u_minus
   psi_minus(3)=v_minus
   psi_minus(4)=H_minus
   psi_plus (1)=1.d0
   psi_plus (2)=u_plus
   psi_plus (3)=v_plus
   psi_plus (4)=H_plus
   if(hou==1)then
      normal(1)=0.d0
      normal(2)=1.d0
      normal(3)=0.d0
      normal(4)=0.d0
   elseif(hou==2) then
      normal(1)=0.d0
      normal(2)=0.d0
      normal(3)=1.d0
      normal(4)=0.d0
   endif

   sonic_plus =sqrt(gamma*p_plus /rho_plus )
   sonic_minus=sqrt(gamma*p_minus/rho_minus)
   sonic_ave  =0.5d0*(sonic_plus+sonic_minus)
   
   Vel_plus =u_plus *normal(2)+v_plus *normal(3)
   Vel_minus=u_minus*normal(2)+v_minus*normal(3)
   Mach_plus =Vel_plus /sonic_ave
   Mach_minus=Vel_minus/sonic_ave
   
   if(abs(Mach_plus)<1)then
      beta_plus=0.25d0*(2.d0-Mach_plus)*(Mach_plus+1)**2
   else
      beta_plus=0.5d0*(1.d0+Mach_plus)
   endif
   if(abs(Mach_minus)<1)then
      beta_minus=0.25d0*(2.d0+Mach_minus)*(Mach_minus+1)**2
   else
      beta_minus=0.5d0*(1.d0-Mach_minus)
   endif
   
   g=-max(min(Mach_plus,0.d0),-1.d0)*min(max(Mach_minus,0.d0),1.d0)
   Veln_ave=(rho_plus*abs(Vel_plus)+rho_minus*abs(Vel_minus))/(rho_plus+rho_minus)
   
   Vel_ave_plus =(1.d0-g)*Veln_ave+g*abs(Vel_plus )
   Vel_ave_minus=(1.d0-g)*Veln_ave+g*abs(Vel_minus)
   Mach_hat=min(1.d0, sqrt(0.5d0*(u_plus**2+v_plus**2+u_minus**2+v_minus**2)))
   X=(1.d0-Mach_hat)**2
   
   pressure_flux=0.5d0*(p_plus+p_minus)+0.5d0*(beta_plus-beta_minus)*(p_plus-p_minus)&
   &+0.5d0*(1.d0-X)*(beta_plus+beta_minus-1.d0)*(p_plus+p_minus)
   mass_flux=0.5d0*(rho_plus*(Vel_plus+Vel_ave_plus)+rho_minus*(Vel_minus-Vel_ave_minus)&
   &-X/sonic_ave*abs(p_plus-p_minus))
   
   temp0=(Mass_flux+abs(Mass_flux))
   temp1=(Mass_flux-abs(Mass_flux))
   flux(:)=0.5d0*temp0*psi_plus(:)+0.5d0*temp1*psi_minus(:)+pressure_flux*normal(:)
   !}}}
end subroutine SLAU
subroutine ROE_BODYFITTED
use prmtr
use variable
implicit none
integer i,j
      corner_for (:)=w(:,step_forward+1,step_height)
      corner_back(:)=w(:,step_backward ,step_height)
      !X_MUSCL{{{
      !!$omp parallel do
      !do j = 0, nj
      !   if(j<=step_height)then
      !      w(:,step_forward+1,j)= w(:,step_forward,j)
      !      w(2,step_forward+1,j)=-w(2,step_forward,j)
      !      w(:,step_backward,j)= w(:,step_backward+1,j)
      !      w(2,step_backward,j)=-w(2,step_backward+1,j)
      !   endif
      !enddo
      !!$omp end parallel do
      w_left( :,:,:)=w(:,:,:)
      w_right(:,:,:)=w(:,:,:)
      !$omp parallel do private(i)
      do j=1, nj-1
         do i=1, ni
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
      !$omp             private(i,rhol,vnl,vtl,pl,El,Hl,&
      !$omp                     rhor,vnr,vtr,pr,Er,Hr,&
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
            !ul  =w_right (2,i-1,j)
            !vl  =w_right (3,i-1,j)
            pl  =w_right (4,i-1,j)
            El  =pl/gamma_bar+0.5d0*rhol*(ul**2+vl**2)
            Hl  =(El+pl)/rhol

            rhor=w_left(1,i  ,j)
            vnr  = nvi(1,i  ,j)*w_left(2,i  ,j)+nvi(2,i  ,j)*w_left(3,i  ,j)
            vtr  =-nvi(2,i  ,j)*w_left(2,i  ,j)+nvi(1,i  ,j)*w_left(3,i  ,j)
            ur  =vnr
            vr  =vtr
            !ur  =w_left(2,i  ,j)
            !vr  =w_left(3,i  ,j)
            pr  =w_left(4,i  ,j)
            Er  =pr/gamma_bar+0.5d0*rhor*(ur**2+vr**2)
            Hr  =(Er+pr)/rhor
            !rhol=w_left (1,i-1,j)
            !ul = nvi(1,i  ,j)*w_left (2,i-1,j)+nvi(2,i  ,j)*w_left (3,i-1,j)
            !vl =-nvi(2,i  ,j)*w_left (2,i-1,j)+nvi(1,i  ,j)*w_left (3,i-1,j)
            !!ul  =w_left (2,i-1,j)
            !!vl  =w_left (3,i-1,j)
            !pl  =w_left (4,i-1,j)
            !El  =pl/gamma_bar+0.5d0*rhol*(ul**2+vl**2)
            !Hl  =(El+pl)/rhol

            !rhor=w_right(1,i  ,j)
            !vnr  = nvi(1,i  ,j)*w_right(2,i  ,j)+nvi(2,i  ,j)*w_right(3,i  ,j)
            !vtr  =-nvi(2,i  ,j)*w_right(2,i  ,j)+nvi(1,i  ,j)*w_right(3,i  ,j)
            !ur  =vnr
            !vr  =vtr
            !!ur  =w_right(2,i  ,j)
            !!vr  =w_right(3,i  ,j)
            !pr  =w_right(4,i  ,j)
            !Er  =pr/gamma_bar+0.5d0*rhor*(ur**2+vr**2)
            !Hr  =(Er+pr)/rhor

            !temp0=sqrt(rhol)
            !temp1=sqrt(rhor)
            !temp2=1.d0/(temp0+temp1)
            !temp3=temp0*temp2
            !temp4=temp1*temp2

            !uaver=temp3*ul+temp4*ur
            !vaver=temp3*vl+temp4*vr
            !Haver=temp3*Hl+temp4*Hr
            !aaver=sqrt(gamma_bar*(Haver-0.5d0*(uaver**2+vaver**2)))


            !lambda(:,:)=0.d0
            !lambda(1,1)=abs(uaver)                 
            !lambda(2,2)=abs(uaver)
            !lambda(3,3)=abs(uaver-aaver) 
            !lambda(4,4)=abs(uaver+aaver) 
            !
            !X_matrix(1,1)=0.d0
            !X_matrix(2,1)=0.d0
            !X_matrix(3,1)=1.d0
            !X_matrix(4,1)=vaver

            !X_matrix(1,2)=1.d0
            !X_matrix(2,2)=uaver
            !X_matrix(3,2)=vaver
            !X_matrix(4,2)=0.5d0*(uaver**2+vaver**2)

            !X_matrix(1,3)=1.d0
            !X_matrix(2,3)=uaver-aaver
            !X_matrix(3,3)=vaver
            !X_matrix(4,3)=Haver-aaver*uaver

            !X_matrix(1,4)=1.d0
            !X_matrix(2,4)=uaver+aaver
            !X_matrix(3,4)=vaver
            !X_matrix(4,4)=Haver+aaver*uaver

            !temp3=aaver**2
            !temp0=2.d0*temp3
            !temp1=gamma_bar/temp0
            !temp2=1.d0/temp0

            !X_inverse(1,1)=-vaver
            !X_inverse(2,1)=temp1*(-2.d0*Haver+4.d0*temp3/gamma_bar)
            !X_inverse(3,1)=temp1*(Haver+(aaver*uaver-temp3)/gamma_bar)
            !X_inverse(4,1)=temp1*(Haver-(aaver*uaver+temp3)/gamma_bar)

            !X_inverse(1,2)=0.d0
            !X_inverse(2,2)=temp1*(2.d0*uaver)
            !X_inverse(3,2)=temp1*(-uaver-aaver/gamma_bar)
            !X_inverse(4,2)=temp1*(-uaver+aaver/gamma_bar)

            !X_inverse(1,3)=1.d0
            !X_inverse(2,3)=temp1*(2.d0*vaver)
            !X_inverse(3,3)=temp1*(-vaver)
            !X_inverse(4,3)=temp1*(-vaver)

            !X_inverse(1,4)=0.d0
            !X_inverse(2,4)=temp1*(-2.d0)
            !X_inverse(3,4)=temp1
            !X_inverse(4,4)=temp1
            !
            !dq(1)=   rhor-   rhol
            !dq(2)=rhor*ur-rhol*ul
            !dq(3)=rhor*vr-rhol*vl
            !dq(4)=     Er-     El

            !temp_arrey1(:)      =matmul(X_inverse(:,:),         dq(:))
            !temp_arrey2(:)      =matmul(   lambda(:,:),temp_arrey1(:))
            !convective_vector(:)=matmul( X_matrix(:,:),temp_arrey2(:))

            !X_Numerical(1,i,j)=      rhol*ul+      rhor*ur
            !X_Numerical(2,i,j)=rhol*ul**2+pl+rhor*ur**2+pr
            !X_Numerical(3,i,j)=   rhol*ul*vl+   rhor*ur*vr
            !X_Numerical(4,i,j)=   rhol*ul*Hl+   rhor*ur*Hr

            !x_Numerical(:,i,j)=0.5d0*(X_Numerical(:,i,j)-convective_vector(:))

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

      w(:,step_forward+1,step_height)=corner_for (:)
      w(:,step_backward ,step_height)=corner_back(:)
      !Y_MUSCL{{{
      !!$omp parallel do
      !do i = 0, ni
      !   if(i>step_forward.and.i<=step_backward)then
      !      w(:,i,step_height)= w(:,i,step_height+1)
      !      w(3,i,step_height)=-w(3,i,step_height+1)
      !   endif
      !enddo
      !!$omp end parallel do
      w_down(:,:,:)=w(:,:,:)
      w_up(  :,:,:)=w(:,:,:)
      !$omp parallel do private(i)
      do j=1, nj
         do i=1, ni-1
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
      !$omp             private(i,rhol,vnl,vtl,pl,El,Hl,&
      !$omp                     rhor,vnr,vtr,pr,Er,Hr,&
      !$omp                     temp0,temp1,temp2,temp3,&
      !$omp                     temp4,uaver,vaver,Haver,&
      !$omp                     aaver,lambda,X_matrix,X_inverse,&
      !$omp                     dq,temp_arrey1,temp_arrey2,&
      !$omp                     convective_vector)
      do j=1, nj
         do i=1, ni-1
            rhol=w_up(1,i,j-1)
            vnl  = nvj(1,i  ,j)*w_up(2,i,j-1)+nvj(2,i  ,j)*w_up(3,i,j-1)
            vtl  =-nvj(2,i  ,j)*w_up(2,i,j-1)+nvj(1,i  ,j)*w_up(3,i,j-1)
            ul  =vnl
            vl  =vtl
            !ul  =w_up(2,i,j-1)
            !vl  =w_up(3,i,j-1)
            pl  =w_up(4,i,j-1)
            El  =pl/gamma_bar+0.5d0*rhol*(ul**2+vl**2)
            Hl  =(El+pl)/rhol

            rhor=w_down  (1,i,j  )
            vnr  = nvj(1,i  ,j)*w_down  (2,i  ,j)+nvj(2,i  ,j)*w_down  (3,i  ,j)
            vtr  =-nvj(2,i  ,j)*w_down  (2,i  ,j)+nvj(1,i  ,j)*w_down  (3,i  ,j)
            ur  =vnr
            vr  =vtr
            !ur  =w_down  (2,i,j  )
            !vr  =w_down  (3,i,j  )
            pr  =w_down  (4,i,j  )
            Er  =pr/gamma_bar+0.5d0*rhor*(ur**2+vr**2)            
            Hr  =(Er+pr)/rhor
            !rhol=w_down(1,i,j-1)
            !vnl  = nvj(1,i  ,j)*w_down(2,i,j-1)+nvj(2,i  ,j)*w_down(3,i,j-1)
            !vtl  =-nvj(2,i  ,j)*w_down(2,i,j-1)+nvj(1,i  ,j)*w_down(3,i,j-1)
            !ul  =vnl
            !vl  =vtl
            !!ul  =w_down(2,i,j-1)
            !!vl  =w_down(3,i,j-1)
            !pl  =w_down(4,i,j-1)
            !El  =pl/gamma_bar+0.5d0*rhol*(ul**2+vl**2)
            !Hl  =(El+pl)/rhol

            !rhor=w_up  (1,i,j  )
            !vnr  = nvj(1,i  ,j)*w_up  (2,i  ,j)+nvj(2,i  ,j)*w_up  (3,i  ,j)
            !vtr  =-nvj(2,i  ,j)*w_up  (2,i  ,j)+nvj(1,i  ,j)*w_up  (3,i  ,j)
            !ur  =vnr
            !vr  =vtr
            !!ur  =w_up  (2,i,j  )
            !!vr  =w_up  (3,i,j  )
            !pr  =w_up  (4,i,j  )
            !Er  =pr/gamma_bar+0.5d0*rhor*(ur**2+vr**2)            
            !Hr  =(Er+pr)/rhor

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

            !temp0=sqrt(rhol)
            !temp1=sqrt(rhor)
            !temp2=1.d0/(temp0+temp1)
            !temp0=temp0*temp2
            !temp1=temp1*temp2

            !uaver=temp0*ul+temp1*ur
            !vaver=temp0*vl+temp1*vr
            !Haver=temp0*Hl+temp1*Hr
            !aaver=sqrt(gamma_bar*(Haver-0.5d0*(uaver**2+vaver**2)))

            !lambda(:,:)=0.d0
            !lambda(1,1)=abs(vaver)                 
            !lambda(2,2)=abs(vaver)
            !lambda(3,3)=abs(vaver-aaver) 
            !lambda(4,4)=abs(vaver+aaver) 
            !
            !Y_matrix(1,1)=0.d0
            !Y_matrix(2,1)=1.d0
            !Y_matrix(3,1)=0.d0
            !Y_matrix(4,1)=uaver

            !Y_matrix(1,2)=1.d0
            !Y_matrix(2,2)=uaver
            !Y_matrix(3,2)=vaver
            !Y_matrix(4,2)=0.5d0*(uaver**2+vaver**2)

            !Y_matrix(1,3)=1.d0
            !Y_matrix(2,3)=uaver
            !Y_matrix(3,3)=vaver-aaver
            !Y_matrix(4,3)=Haver-aaver*vaver

            !Y_matrix(1,4)=1.d0
            !Y_matrix(2,4)=uaver
            !Y_matrix(3,4)=vaver+aaver
            !Y_matrix(4,4)=Haver+aaver*vaver

            !temp3=aaver**2
            !temp0=2.d0*temp3
            !temp1=gamma_bar/temp0
            !temp2=1.d0/temp0
            !Y_inverse(1,1)=-uaver
            !Y_inverse(2,1)=temp1*(-2.d0*Haver+(4.d0*temp3)/gamma_bar)
            !Y_inverse(3,1)=temp1*(Haver+(aaver*(vaver-aaver)/gamma_bar))
            !Y_inverse(4,1)=temp1*(Haver-(aaver*(vaver+aaver)/gamma_bar))

            !Y_inverse(1,2)=1.d0
            !Y_inverse(2,2)=temp1*(2.d0*uaver)
            !Y_inverse(3,2)=temp1*(-uaver)
            !Y_inverse(4,2)=temp1*(-uaver)

            !Y_inverse(1,3)=0.d0
            !Y_inverse(2,3)=temp1*(2.d0*vaver)
            !Y_inverse(3,3)=temp1*(-vaver-aaver/gamma_bar)
            !Y_inverse(4,3)=temp1*(-vaver+aaver/gamma_bar)

            !Y_inverse(1,4)=0.d0
            !Y_inverse(2,4)=temp1*(-2.d0)
            !Y_inverse(3,4)=temp1
            !Y_inverse(4,4)=temp1

            !dq(1)=   rhor-   rhol
            !dq(2)=rhor*ur-rhol*ul
            !dq(3)=rhor*vr-rhol*vl
            !dq(4)=     Er-     El

            !temp_arrey1(:)      =matmul(Y_inverse(:,:),         dq(:))
            !temp_arrey2(:)      =matmul(   lambda(:,:),temp_arrey1(:))
            !convective_vector(:)=matmul( Y_matrix(:,:),temp_arrey2(:))

            !Y_Numerical(1,i,j)=      rhol*vl+      rhor*vr
            !Y_Numerical(2,i,j)=   rhol*ul*vl+   rhor*ur*vr
            !Y_Numerical(3,i,j)=rhol*vl**2+pl+rhor*vr**2+pr
            !Y_Numerical(4,i,j)=   rhol*vl*Hl+   rhor*vr*Hr

            !Y_Numerical(:,i,j)=0.5d0*(Y_Numerical(:,i,j)-convective_vector(:))

            !temp0=sqrt(rhol)
            !temp1=sqrt(rhor)
            !temp2=1.d0/(temp0+temp1)
            !temp3=temp0*temp2
            !temp4=temp1*temp2

            !uaver=temp3*ul+temp4*ur
            !vaver=temp3*vl+temp4*vr
            !Haver=temp3*Hl+temp4*Hr
            !aaver=sqrt(gamma_bar*(Haver-0.5d0*(uaver**2+vaver**2)))

            !lambda(:,:)=0.d0
            !lambda(1,1)=abs(vaver)                 
            !lambda(2,2)=abs(vaver)
            !lambda(3,3)=abs(vaver-aaver) 
            !lambda(4,4)=abs(vaver+aaver) 
            !
            !X_matrix(1,1)=0.d0
            !X_matrix(2,1)=1.d0
            !X_matrix(3,1)=0.d0
            !X_matrix(4,1)=uaver

            !X_matrix(1,2)=1.d0
            !X_matrix(2,2)=uaver
            !X_matrix(3,2)=vaver
            !X_matrix(4,2)=0.5d0*(uaver**2+vaver**2)

            !X_matrix(1,3)=1.d0
            !X_matrix(2,3)=uaver
            !X_matrix(3,3)=vaver-aaver
            !X_matrix(4,3)=Haver-aaver*vaver

            !X_matrix(1,4)=1.d0
            !X_matrix(2,4)=uaver
            !X_matrix(3,4)=vaver+aaver
            !X_matrix(4,4)=Haver+aaver*vaver

            !temp3=aaver**2
            !temp0=2.d0*temp3
            !temp1=gamma_bar/temp0
            !temp2=1.d0/temp0
            !X_inverse(1,1)=-uaver
            !X_inverse(2,1)=temp1*(-2.d0*Haver+(4.d0*temp3)/gamma_bar)
            !X_inverse(3,1)=temp1*(Haver+(aaver*(vaver-aaver)/gamma_bar))
            !X_inverse(4,1)=temp1*(Haver-(aaver*(vaver+aaver)/gamma_bar))

            !X_inverse(1,2)=1.d0
            !X_inverse(2,2)=temp1*(2.d0*uaver)
            !X_inverse(3,2)=temp1*(-uaver)
            !X_inverse(4,2)=temp1*(-uaver)

            !X_inverse(1,3)=0.d0
            !X_inverse(2,3)=temp1*(2.d0*vaver)
            !X_inverse(3,3)=temp1*(-vaver-aaver/gamma_bar)
            !X_inverse(4,3)=temp1*(-vaver+aaver/gamma_bar)

            !X_inverse(1,4)=0.d0
            !X_inverse(2,4)=temp1*(-2.d0)
            !X_inverse(3,4)=temp1
            !X_inverse(4,4)=temp1

            !dq(1)=   rhor-   rhol
            !dq(2)=rhor*ur-rhol*ul
            !dq(3)=rhor*vr-rhol*vl
            !dq(4)=     Er-     El

            !temp_arrey1(:)      =matmul(X_inverse(:,:),         dq(:))
            !temp_arrey2(:)      =matmul(   lambda(:,:),temp_arrey1(:))
            !convective_vector(:)=matmul( X_matrix(:,:),temp_arrey2(:))

            !Y_Numerical(1,i,j)=      rhol*vl+      rhor*vr
            !Y_Numerical(2,i,j)=   rhol*ul*vl+   rhor*ur*vr
            !Y_Numerical(3,i,j)=rhol*vl**2+pl+rhor*vr**2+pr
            !Y_Numerical(4,i,j)=   rhol*vl*Hl+   rhor*vr*Hr

            !Y_Numerical(:,i,j)=0.5d0*(Y_Numerical(:,i,j)-convective_vector(:))

            !Y_Numerical(1,i,j)=      rhol*ul+      rhor*ur
            !Y_Numerical(2,i,j)=   rhol*ul*vl+   rhor*ur*vr
            !Y_Numerical(3,i,j)=rhol*ul**2+pl+rhor*ur**2+pr
            !Y_Numerical(4,i,j)=   rhol*ul*Hl+   rhor*ur*Hr


            !Y_Numerical(:,i,j)=0.5d0*(Y_Numerical(:,i,j)-convective_vector(:))

            !temp2=Y_Numerical(2,i,j)
            !temp3=Y_Numerical(3,i,j)
            !
            !Y_Numerical(2,i,j)= nvj(1,i,j)*temp2+nvj(2,i,j)*temp3
            !Y_Numerical(3,i,j)=-nvj(2,i,j)*temp2+nvj(1,i,j)*temp3
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

end subroutine ROE_BODYFITTED
