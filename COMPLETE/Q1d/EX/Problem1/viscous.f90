subroutine set_geojac !{{{
use prmtr
use variable
implicit none
double precision dxdxi,drdxi,dxdeta,drdeta
double precision det
integer i,j

do j=1,nj-1
   do i=1,ni
      dxdxi=x(i,j)-x(i-1,j)
      drdxi=r(i,j)-r(i-1,j)
      
      dxdeta=(x(i,j+1)+x(i-1,j+1)&
             -x(i,j-1)-x(i-1,j-1))*0.25d0
      drdeta=(r(i,j+1)+r(i-1,j+1)&
             -r(i,j-1)-r(i-1,j-1))*0.25d0
      
      det=dxdxi*drdeta-drdxi*dxdeta
      
      !geojaci(1,1,i,j)=dxdxi
      !geojaci(1,2,i,j)=dxdeta
      !geojaci(2,1,i,j)=drdxi
      !geojaci(2,2,i,j)=drdeta

      geojaci(1,1,i,j)= drdeta/det
      geojaci(1,2,i,j)=-dxdeta/det
      geojaci(2,1,i,j)=-drdxi /det
      geojaci(2,2,i,j)= dxdxi /det
   end do
end do

do j=1,nj
   do i=1,ni-1
      dxdeta=x(i,j)-x(i,j-1)
      drdeta=r(i,j)-r(i,j-1)

      dxdxi=(x(i+1,j)+x(i+1,j-1)&
            -x(i-1,j)-x(i-1,j-1))*0.25d0
      drdxi=(r(i+1,j)+r(i+1,j-1)&
            -r(i-1,j)-r(i-1,j-1))*0.25d0
      
      det=dxdxi*drdeta-drdxi*dxdeta
      
      !geojacj(1,1,i,j)=dxdxi
      !geojacj(1,2,i,j)=dxdeta
      !geojacj(2,1,i,j)=drdxi
      !geojacj(2,2,i,j)=drdeta

      geojacj(1,1,i,j)= drdeta/det
      geojacj(1,2,i,j)=-dxdeta/det
      geojacj(2,1,i,j)=-drdxi /det
      geojacj(2,2,i,j)= dxdxi /det
   end do
end do
end subroutine set_geojac  !}}}
subroutine set_viscous
use prmtr
use variable
implicit none
double precision, dimension(4) :: Ev, Fv
double precision, dimension(4) :: dwdxi,dwdata,dwdx,dwdr
double precision, dimension(0:ni+1,0:nj+1) :: Tdeg 
double precision k_heat
double precision dudxi,dudeta
double precision dudx,dudr
double precision dvdxi,dvdeta
double precision dvdx,dvdr
double precision dTdxi,dTdeta
double precision dTdx,dTdr
double precision rho,u,v,Temper
double precision Rgas
double precision tau_xx,tau_xr,tau_rr
double precision qx,qr
double precision divu

double precision D
double precision, parameter::Sc=1d0
double precision, parameter::Prentl=7d-1
double precision temp10, temp11, temp12, temp13

integer i,j,k

!$omp parallel do shared(Tdeg,j),&
!$omp             private(i)
do j=0,nj
   do i=0,ni
      Tdeg(i,j)=w(4,i,j)/w(1,i,j)/gas_specific
   end do
end do
!$omp end parallel do

!set TGvi!{{{
!$omp parallel do shared(vis_i,j),&
!$omp             private(rho,u,v,Temper,&
!$omp             k_heat,dudxi,dvdxi,dTdxi,&
!$omp             dudeta,dvdeta,dTdeta,&
!$omp             dudx,dvdx,dTdx,dudr,dvdr,dTdr,&
!$omp             divu,tau_xx,tau_xr,tau_rr,qx,qr,&
!$omp             Ev,Fv,i)
do j=1,nj-1
   do i=1,ni
      rho    = (w(1 ,i,j)+w(1 ,i-1,j))*0.5d0
      u      = (w(2 ,i,j)+w(2 ,i-1,j))*0.5d0
      v      = (w(3 ,i,j)+w(3 ,i-1,j))*0.5d0
      Temper = (Tdeg(i,j)+Tdeg(i-1,j))*0.5d0 !total energy

      k_heat=mu*gamma/(gamma-1d0)*gas_specific/Prentl               !k from mu by Prandtl number

      dudxi  =  w(2,i,j) -w(2,i-1,j)
      dvdxi  =  w(3,i,j) -w(3,i-1,j)
      dTdxi  = Tdeg(i,j)-Tdeg(i-1,j)

      dudeta =  (w(2,i-1,j+1) +w(2,i,j+1)&
                -w(2,i-1,j-1) -w(2,i,j-1))*0.25d0
      dvdeta =  (w(3,i-1,j+1) +w(3,i,j+1)&
                -w(3,i-1,j-1) -w(3,i,j-1))*0.25d0
      dTdeta = (Tdeg(i-1,j+1)+Tdeg(i,j+1)&
               -Tdeg(i-1,j-1)-Tdeg(i,j-1))*0.25d0

      dudx=dudxi*geojaci(1,1,i,j)+dudeta*geojaci(2,1,i,j)
      dvdx=dvdxi*geojaci(1,1,i,j)+dvdeta*geojaci(2,1,i,j)
      dTdx=dTdxi*geojaci(1,1,i,j)+dTdeta*geojaci(2,1,i,j)
      dudr=dudxi*geojaci(1,2,i,j)+dudeta*geojaci(2,2,i,j)
      dvdr=dvdxi*geojaci(1,2,i,j)+dvdeta*geojaci(2,2,i,j)
      dTdr=dTdxi*geojaci(1,2,i,j)+dTdeta*geojaci(2,2,i,j)

      divu=dudx+dvdr

      tau_xx=-2d0/3d0*mu*divu+2d0*mu*dudx
      tau_rr=-2d0/3d0*mu*divu+2d0*mu*dvdr
      tau_xr= mu*(dvdx+dudr)
      qx  = -k_heat*dTdx
      qr  = -k_heat*dTdr

      Ev(1)=0d0
      Ev(2)=tau_xx
      Ev(3)=tau_xr
      Ev(4)=u*tau_xx+v*tau_xr-qx

      Fv(1)=0d0
      Fv(2)=tau_xr
      Fv(3)=tau_rr
      Fv(4)=u*tau_xr+v*tau_rr-qr

      vis_i(:,i,j)=Ev*nvi(1,i,j)+Fv*nvi(2,i,j)
   end do
end do
!$omp end parallel do
!}}}

!set TGvj!{{{
!$omp parallel do shared(vis_j,j),&
!$omp             private(rho,u,v,Temper,&
!$omp             k_heat,dudxi,dvdxi,dTdxi,&
!$omp             dudeta,dvdeta,dTdeta,&
!$omp             dudx,dvdx,dTdx,dudr,dvdr,dTdr,&
!$omp             divu,tau_xx,tau_xr,tau_rr,qx,qr,&
!$omp             Ev,Fv,i)
do j=1,nj
   do i=1,ni-1
      rho    = (w(1    ,i,j-1)+w(1    ,i,j))*0.5d0
      u      = (w(2    ,i,j-1)+w(2    ,i,j))*0.5d0
      v      = (w(3    ,i,j-1)+w(3    ,i,j))*0.5d0
      Temper =    (Tdeg(i,j-1)+   Tdeg(i,j))*0.5d0

      k_heat=mu*gamma/(gamma-1d0)*gas_specific/Prentl               !k from mu by Prandtl number

      dudxi  =  (w(2,i+1,j) +w(2,i+1,j-1)&
                -w(2,i-1,j) -w(2,i-1,j-1))*0.25d0
      dvdxi  =  (w(3,i+1,j) +w(3,i+1,j-1)&
                -w(3,i-1,j) -w(3,i-1,j-1))*0.25d0
      dTdxi  = (Tdeg(i+1,j)+Tdeg(i+1,j-1)&
               -Tdeg(i-1,j)-Tdeg(i-1,j-1))*0.25d0

      dudeta =  w(2,i,j) -w(2,i,j-1)
      dvdeta =  w(3,i,j) -w(3,i,j-1)
      dTdeta = Tdeg(i,j)-Tdeg(i,j-1)

      dudx=dudxi*geojacj(1,1,i,j)+dudeta*geojacj(2,1,i,j)
      dvdx=dvdxi*geojacj(1,1,i,j)+dvdeta*geojacj(2,1,i,j)
      dTdx=dTdxi*geojacj(1,1,i,j)+dTdeta*geojacj(2,1,i,j)
      dudr=dudxi*geojacj(1,2,i,j)+dudeta*geojacj(2,2,i,j)
      dvdr=dvdxi*geojacj(1,2,i,j)+dvdeta*geojacj(2,2,i,j)
      dTdr=dTdxi*geojacj(1,2,i,j)+dTdeta*geojacj(2,2,i,j)

      !!2 dimentional-plane
      divu=dudx+dvdr

      tau_xx=-2d0/3d0*mu*divu+2d0*mu*dudx
      tau_rr=-2d0/3d0*mu*divu+2d0*mu*dvdr
      tau_xr= mu*(dvdx+dudr)
      qx  = -k_heat*dTdx
      qr  = -k_heat*dTdr

      Ev(1)=0d0
      Ev(2)=tau_xx
      Ev(3)=tau_xr
      Ev(4)=u*tau_xx+v*tau_xr-qx

      Fv(1)=0d0
      Fv(2)=tau_xr
      Fv(3)=tau_rr
      Fv(4)=u*tau_xr+v*tau_rr-qr

      vis_j(:,i,j)=Ev*nvj(1,i,j)+Fv*nvj(2,i,j)
   end do
end do
!$omp end parallel do
!}}}

end subroutine set_viscous
