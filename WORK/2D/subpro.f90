subroutine setdt(w,dx,dt)
use prmtr
implicit none
integer i
double precision,dimension(1:3,0:ni)::w
double precision dx,dt
double precision temp(0:3)
temp(0)=1.d300
do i=1,ni-1
   temp(3)=sqrt(gamma*w(3,i)/w(1,i))
   temp(1)=w(2,i)
   temp(1)=dx/(temp(1)+temp(3))
   temp(0)=min(temp(0),temp(1))
enddo
dt=CFL*temp(0)
end subroutine setdt

subroutine setIC(q,w,pL,uL,TL,pR,uR,tR,sonic,x,Vol,PreResi)
use prmtr
implicit none
integer i
double precision,dimension(1:3,0:ni)::w,q
double precision,dimension(    0:ni)::sonic,PreResi,x,Vol
double precision pL,uL,rL,TL,pR,uR,rR,tR
double precision temp(0:3)
do i=0,ni
   temp(0) = (rR-rL)/float(ni-1)
   temp(1) = (uR-uL)/float(ni-1)
   temp(2) = (pR-pL)/float(ni-1)
   w(1,i)= rL+temp(0)*float(i-1)
   w(2,i)= uL+temp(1)*float(i-1)
   w(3,i)= pL+temp(2)*float(i-1)
   q(1,i)=w(1,i)
   q(2,i)=w(2,i)*q(1,i)
   q(3,i)=w(3,i)/(gamma-1.d0)+0.5d0*w(1,i)*(w(2,i)**2)
   temp(0) = x(i+1) + 2.2d0/3d0*(x(i+1)-1.5d0)**3
   temp(1) = x(i  ) + 2.2d0/3d0*(x(i  )-1.5d0)**3
   Vol(i) = temp(0)-temp(1)
   sonic(i) = sqrt(gamma*w(3,i)/w(1,i))
end do
PreResi(:) = q(3,:)
end subroutine setIC
subroutine setBC(w,pL,TL,pR,tR)
use prmtr
implicit none
integer i
double precision,dimension(1:3,0:ni)::w
double precision pL,TL,pR,tR
double precision temp(0:3)
w(1,0)=pL/(TL*gas_spec)
w(2,0)=1.5d0*w(2,1)-0.5d0*w(2,2)
w(3,0)=pL

!w(1,ni)=1d0
!w(2,ni)=w(2,1)
!w(3,ni)=pR
w(:,ni)=1.5d0*w(:,ni-1)-0.5d0*w(:,ni-2)

end subroutine setBC

subroutine setRESIDUAL(q,PreResi,Resi)
use prmtr
implicit none
integer i
double precision,dimension(1:3,0:ni)::q
double precision,dimension(    0:ni)::PreResi
double precision temp(0:3)
double precision Resi
temp(2) =0d0
do i=1,ni-1
   temp(1) =(q(3,i)-PreResi(i))**2
   temp(2) = temp(1)+temp(2)
   !write(35,'(i4,f15.4,f15.4,f15.4)') i,q(:,i)
enddo
Resi = temp(2)
end subroutine setRESIDUAL

subroutine setPrimitiveValue(w,q,sonic,t)
use prmtr
implicit none
integer i
double precision,dimension(1:3,0:ni)::q,w
double precision,dimension(    0:ni)::sonic
double precision temp(0:3)
double precision t
do i=1,ni-1
   w(1,i)=q(1,i)
   temp(0)=1d0/w(1,i)
   w(2,i)=q(2,i)*temp(0)
   w(3,i)=(gamma-1.d0)*(q(3,i)-0.5d0*w(1,i)*(w(2,i)**2))
   sonic(i) = sqrt(gamma*w(3,i)/w(1,i))
   if(w(3,i)<0.d0)then
      write(*,*) i, t
      write(*,'(4es15.7)') w(:,i)
      call exit(1)
   endif
enddo
end subroutine setPrimitiveValue
subroutine setConservativeVector(w,q,PreResi)
use prmtr
implicit none
integer i
double precision,dimension(1:3,0:ni)::q,w
double precision,dimension(    0:ni)::PreResi
double precision temp(0:3)
do i=0,ni
   q(1,i)=w(1,i)
   q(2,i)=w(2,i)*q(1,i)
   q(3,i)=w(3,i)/(gamma-1.d0)+0.5d0*w(1,i)*(w(2,i)**2)
enddo
PreResi(:) = q(3,:)
end subroutine setConservativeVector
