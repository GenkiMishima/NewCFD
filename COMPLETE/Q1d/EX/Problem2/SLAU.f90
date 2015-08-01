subroutine SLAU(rho_minus, rho_plus, u_minus, u_plus, v_minus, v_plus, p_minus, p_plus&
            &, H_minus, H_plus, flux, hou)
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
end subroutine SLAU
