!==========================================================
!
!  This program is 2Dimentional CFD
!
!==========================================================
module global
   double precision, parameter :: pi=2.d0*acos(0.d0)
   double precision, parameter :: gamma=1.4d0
!MUSCL INCREMENT
   double precision, parameter :: phi=1.d0
   double precision, parameter :: kappa=1.d0
   double precision, parameter :: dx=0.01d0
   double precision, parameter :: dy=0.01d0
   double precision, parameter :: dt=0.001d0
!X direction initial number
   double precision, parameter :: rho_left_initial=5.d0
   double precision, parameter :: u_left_initial=0.d0
   double precision, parameter :: v_left_initial=0.d0
   double precision, parameter :: p_left_initial=5.d0
   double precision, parameter :: rho_right_initial=1.d0
   double precision, parameter :: u_right_initial=0.d0
   double precision, parameter :: v_right_initial=0.d0
   double precision, parameter :: p_right_initial=1.d0
!Y direction initial number
   double precision, parameter :: rho_down_initial=5.d0
   double precision, parameter :: u_down_initial=0.d0
   double precision, parameter :: p_down_initial=5.d0
   double precision, parameter :: rho_up_initial=1.d0
   double precision, parameter :: u_up_initial=0.d0
   double precision, parameter :: p_up_initial=1.d0
end module global

program main 
use global
   implicit none
   integer i,j,k,m,n,o,time
   integer X_max, Y_max, T_max
   double precision rho(-2:1003,-2:1003), rhorho(-2:1003,-2:1003),&
   & rho_right(-1:1002,-1:1002), rho_left(-1:1002,-1:1002), rho_up(-1:1002,-1:1002), rho_down(-1:1002,-1:1002)
   double precision u(-2:1003,-2:1003), uu(-2:1003,-2:1003),&
   &u_right(-1:1002,-1:1002), u_left(-1:1002,-1:1002), u_up(-1:1002,-1:1002) ,u_down(-1:1002,-1:1002)
   double precision v(-2:1003,-2:1003), vv(-2:1003,-2:1003),&
   &v_right(-1:1002,-1:1002), v_left(-1:1002,-1:1002), v_up(-1:1002,-1:1002) ,v_down(-1:1002,-1:1002)
   double precision p(-2:1003,-2:1003), pp(-2:1003,-2:1003),&
   &p_right(-1:1002,-1:1002), p_left(-1:1002,-1:1002), p_up(-1:1002,-1:1002) ,p_down(-1:1002,-1:1002)
   double precision E(-2:1003,-2:1003), EE(-2:1003,-2:1003),&
   &E_right(-1:1002,-1:1002), E_left(-1:1002,-1:1002), E_up(-1:1002,-1:1002) ,E_down(-1:1002,-1:1002)
   double precision H(-2:1003,-2:1003), HH(-2:1003,-2:1003),&
   &H_right(-1:1002,-1:1002), H_left(-1:1002,-1:1002), H_up(-1:1002,-1:1002) ,H_down(-1:1002,-1:1002)
   double precision x(-2:1003,-2:1003), y(-2:1003,-2:1003)
   double precision X_matrix(4,4), X_inverse(4,4), lamda(4,4), Y_matrix(4,4), Y_inverse(4,4)
   double precision u_ave, v_ave, H_ave, a_ave, t
   double precision Xlamda_matrix(4,4), XlamdaX_matrix(4,4),Ylamda_matrix(4,4), YlamdaY_matrix(4,4)
   double precision q_minus(4), q_plus(4), convective_vector(4), X_Numerical_minus(4), X_Numerical_plus(4)
   double precision X_Flux(4), Y_Flux(4), Y_Numerical_minus(4), Y_Numerical_plus(4)
   double precision, Dimension(4,4) :: A_matrix, B_matrix
   double precision q(4), qq(4)
   double precision Mach_number(-1:1002,-1:1002)
   double precision temp0, temp1, temp2, temp3, temp4, temp5, place
   character*20 tmpstring
   X_max=int(dx**(-1))
   Y_max=int(dy**(-1))
   T_max=int(dt**(-1))


   open(10,file='output1.d')
         do i = -2, X_max+3
            place = dx*0.5d0+dble(i)*dx
            if(place < 0.5d0) then
               rho(i,:)=rho_left_initial
               u(i,:)=u_left_initial
               v(i,:)=v_left_initial
               p(i,:)=p_left_initial
            else
               rho(i,:)=rho_right_initial
               u(i,:)=u_right_initial
               v(i,:)=v_right_initial
               p(i,:)=p_right_initial
            endif
            E(i,:)=p(i,1)/(gamma-1)+rho(i,1)/2.d0*(u(i,1)**2+v(i,1)**2)
            H(i,:)=(E(i,1)+p(i,1))/rho(i,1)
            write(10, *) place, rho(i,1)
         enddo
   close(10)

   do time=1, 5*T_max
      t=dble(time) * dt
      !MUSCL Parallel
      
      do j=-1, Y_max+2
         do i=-1, X_max+2

            call MUSCLlimiter( rho(i-1,j), rho(i,j), rho(i+1,j),  rho_left(i,j),  rho_right(i,j))
            call MUSCLlimiter(   u(i-1,j),   u(i,j),   u(i+1,j),    u_left(i,j),    u_right(i,j))
            call MUSCLlimiter(   v(i-1,j),   v(i,j),   v(i+1,j),    v_left(i,j),    v_right(i,j))
            call MUSCLlimiter(   p(i-1,j),   p(i,j),   p(i+1,j),    p_left(i,j),    p_right(i,j))
            call MUSCLlimiter(   E(i-1,j),   E(i,j),   E(i+1,j),    E_left(i,j),    E_right(i,j))
            call MUSCLlimiter(   H(i-1,j),   H(i,j),   H(i+1,j),    H_left(i,j),    H_right(i,j))
 

   !        write (*, *) rho_left(i), rho_right(i), u_left(i), u_right(i)
         enddo
      enddo
      !MUSCL Vertical
      do j=-1, Y_max+2
         do i=-1, X_max+2

            call MUSCLlimiter( rho(i,j-1), rho(i,j), rho(i,j+1),  rho_down(i,j),  rho_up(i,j))
            call MUSCLlimiter(   u(i,j-1),   u(i,j),   u(i,j+1),    u_down(i,j),    u_up(i,j))
            call MUSCLlimiter(   v(i,j-1),   v(i,j),   v(i,j+1),    v_down(i,j),    v_up(i,j))
            call MUSCLlimiter(   p(i,j-1),   p(i,j),   p(i,j+1),    p_down(i,j),    p_up(i,j))
            call MUSCLlimiter(   E(i,j-1),   E(i,j),   E(i,j+1),    E_down(i,j),    E_up(i,j))
            call MUSCLlimiter(   H(i,j-1),   H(i,j),   H(i,j+1),    H_down(i,j),    H_up(i,j))

   !        write (*, *) rho_down(i), rho_up(i), u_down(i), u_up(i)
         enddo
      enddo

      do j=1, Y_max
         do i=1, X_max
!===============================================================================================
!start parallel=================================================================================
            !Between i and i-1
            u_ave=u_average(rho_left(i-1,j),rho_right(i,j),u_left(i-1,j),u_right(i,j))
            v_ave=v_average(rho_left(i-1,j),rho_right(i,j),v_left(i-1,j),v_right(i,j))
            H_ave=H_average(rho_left(i-1,j),rho_right(i,j),H_left(i-1,j),H_right(i,j))
            a_ave=a_average(u_ave,v_ave,H_ave)
            lamda(:,:)=0.d0
            
            lamda(1,1)=abs(u_ave)
            lamda(2,2)=abs(u_ave)
            lamda(3,3)=abs(u_ave-a_ave)
            lamda(4,4)=abs(u_ave+a_ave)
            X_matrix=X_coefficience(u_ave,v_ave,H_ave,a_ave)
            X_inverse=X_inverse_co(u_ave,v_ave,H_ave,a_ave)

            A_matrix=A_coefficience(u_ave,v_ave,H_ave,a_ave)

            do m=1,4
               do n=1,4
                  temp1=0.d0
                  do o=1,4
                     temp1=X_matrix(m,o)*lamda(o,n)+temp1
                  enddo
                  Xlamda_matrix(m,n)=temp1
                  !write(*,*) Ans(i,j)
               enddo
            enddo

            do m=1,4
               do n=1,4
                  temp1=0.d0
                  do o=1,4
                     temp1=Xlamda_matrix(m,o)*X_inverse(o,n)+temp1
                  enddo
                  XlamdaX_matrix(m,n)=temp1
                  !write(*,*) Ans(i,j)
               enddo
            enddo
            
            if(time==1 .and. i==1 .and. j==1) then
            open(11, file='test.d')
               do m=1,4
                  do n=1,4
                    write(11, *) XlamdaX_matrix(m,n), A_matrix(m,n) 
                  enddo
               enddo
            close(11)
            endif

            q_minus(1)=rho_right(i,j)-rho_left(i-1,j)
            q_minus(2)=u_right(i,j)-u_left(i-1,j)
            q_minus(3)=v_right(i,j)-v_left(i-1,j)
            q_minus(4)=p_right(i,j)-p_left(i-1,j)

            do m=1,4
               temp1=0.d0
               do o=1,4
                  temp1=XlamdaX_matrix(m,o)*q_minus(o)+temp1
               enddo
               convective_vector(m)=temp1
               !write(*,*) Anp(i)
            enddo
            
            temp1=rho_left(i-1,j)*u_left(i-1,j)+rho_right(i,j)*u_right(i,j)
            temp2=rho_left(i-1,j)*u_left(i-1,j)**2+p_left(i-1,j)+rho_right(i,j)*u_right(i,j)**2+p_right(i,j)
            temp3=rho_left(i-1,j)*u_left(i-1,j)*v_left(i-1,j)+rho_right(i,j)*u_right(i,j)*v_right(i,j)
            temp4=rho_left(i-1,j)*u_left(i-1,j)*H_left(i-1,j)+rho_right(i,j)*u_right(i,j)*H_right(i,j)

            X_Numerical_minus(1)=0.5d0*(temp1-convective_vector(1))
            X_Numerical_minus(2)=0.5d0*(temp2-convective_vector(2))
            X_Numerical_minus(3)=0.5d0*(temp3-convective_vector(3))
            X_Numerical_minus(4)=0.5d0*(temp4-convective_vector(4))

            !Between i and i+1

            u_ave=u_average(rho_left(i,j),rho_right(i+1,j),u_left(i,j),u_right(i+1,j))
            v_ave=v_average(rho_left(i,j),rho_right(i+1,j),v_left(i,j),v_right(i+1,j))
            H_ave=H_average(rho_left(i,j),rho_right(i+1,j),H_left(i,j),H_right(i+1,j))
            a_ave=a_average(u_ave,v_ave,H_ave)
            lamda(:,:)=0.d0
            lamda(1,1)=abs(u_ave)
            lamda(2,2)=abs(u_ave)
            lamda(3,3)=abs(u_ave-a_ave)
            lamda(4,4)=abs(u_ave+a_ave)
            X_matrix=X_coefficience(u_ave,v_ave,H_ave,a_ave)
            X_inverse=X_inverse_co(u_ave,v_ave,H_ave,a_ave)


            A_matrix=A_coefficience(u_ave,v_ave,H_ave,a_ave)

            do m=1,4
               do n=1,4
                  temp1=0.d0
                  do o=1,4
                     temp1=X_matrix(m,o)*lamda(o,n)+temp1
                  enddo
                  Xlamda_matrix(m,n)=temp1
                  !write(*,*) Ans(i,j)
               enddo
            enddo


            do m=1,4
               do n=1,4
                  temp1=0.d0
                  do o=1,4
                     temp1=Xlamda_matrix(m,o)*X_inverse(o,n)+temp1
                  enddo
                  XlamdaX_matrix(m,n)=temp1
                  !write(*,*) Ans(i,j)
               enddo
            enddo
            

            q_plus(1)=rho_right(i+1,j)-rho_left(i,j)
            q_plus(2)=u_right(i+1,j)-u_left(i,j)
            q_plus(3)=v_right(i+1,j)-v_left(i,j)
            q_plus(4)=p_right(i+1,j)-p_left(i,j)

            do m=1,4
               temp1=0.d0
               do o=1,4
                  temp1=XlamdaX_matrix(m,o)*q_plus(o)+temp1
               enddo
               convective_vector(m)=temp1
               !write(*,*) Anp(i)
            enddo
            
            temp1=rho_left(i,j)*u_left(i,j)+rho_right(i+1,j)*u_right(i+1,j)
            temp2=rho_left(i,j)*u_left(i,j)**2+p_left(i,j)+rho_right(i+1,j)*u_right(i+1,j)**2+p_right(i+1,j)
            temp3=rho_left(i,j)*u_left(i,j)*v_left(i,j)+rho_right(i+1,j)*u_right(i+1,j)*v_right(i+1,j)
            temp4=rho_left(i,j)*u_left(i,j)*H_left(i,j)+rho_right(i+1,j)*u_right(i+1,j)*H_right(i+1,j)

            X_Numerical_plus(1)=0.5d0*(temp1-convective_vector(1))
            X_Numerical_plus(2)=0.5d0*(temp2-convective_vector(2))
            X_Numerical_plus(3)=0.5d0*(temp3-convective_vector(3))
            X_Numerical_plus(4)=0.5d0*(temp4-convective_vector(4))

            X_Flux(1)=X_Numerical_minus(1)-X_Numerical_plus(1)
            X_Flux(2)=X_Numerical_minus(2)-X_Numerical_plus(2)
            X_Flux(3)=X_Numerical_minus(3)-X_Numerical_plus(3)
            X_Flux(4)=X_Numerical_minus(4)-X_Numerical_plus(4)
!end parallel===============================================================================================
!============================================================================================================
!==========================================================================================================
!start vertical============================================================================================

            u_ave=u_average(rho_down(i,j-1),rho_up(i,j),u_down(i,j-1),u_up(i,j))
            v_ave=v_average(rho_down(i,j-1),rho_up(i,j),v_down(i,j-1),v_up(i,j))
            H_ave=H_average(rho_down(i,j-1),rho_up(i,j),H_down(i,j-1),H_up(i,j))
            a_ave=a_average(u_ave,v_ave,H_ave)
            lamda(:,:)=0.d0
            lamda(1,1)=abs(v_ave)
            lamda(2,2)=abs(v_ave)
            lamda(3,3)=abs(v_ave-a_ave)
            lamda(4,4)=abs(v_ave+a_ave)
            Y_matrix=Y_coefficience(u_ave,v_ave,H_ave,a_ave)
            Y_inverse=Y_inverse_co(u_ave,v_ave,H_ave,a_ave)

            B_matrix=B_coefficience(u_ave,v_ave,H_ave,a_ave)
            
            do m=1,4
               do n=1,4
                  temp1=0.d0
                  do o=1,4
                     temp1=Y_matrix(m,o)*lamda(o,n)+temp1
                  enddo
                  Ylamda_matrix(m,n)=temp1
                  !write(*,*) Ans(i,j)
               enddo
            enddo

            do m=1,4
               do n=1,4
                  temp1=0.d0
                  do o=1,4
                     temp1=Ylamda_matrix(m,o)*Y_inverse(o,n)+temp1
                  enddo
                  YlamdaY_matrix(m,n)=temp1
                  !write(*,*) Ans(i,j)
               enddo
            enddo

            q_minus(1)=rho_up(i,j)-rho_down(i,j-1)
            q_minus(2)=u_up(i,j)-u_down(i,j-1)
            q_minus(3)=v_up(i,j)-v_down(i,j-1)
            q_minus(4)=p_up(i,j)-p_down(i,j-1)

            do m=1,4
               temp1=0.d0
               do o=1,4
                  temp1=YlamdaY_matrix(m,o)*q_minus(o)+temp1
               enddo
               convective_vector(m)=temp1
               !write(*,*) Anp(i)
            enddo
            
            temp1=rho_down(i,j-1)*u_down(i,j-1)+rho_up(i,j)*u_up(i,j)
            temp2=rho_down(i,j-1)*u_down(i,j-1)**2+p_down(i,j-1)+rho_up(i,j)*u_up(i,j)**2+p_up(i,j)
            temp3=rho_down(i,j-1)*u_down(i,j-1)*v_down(i,j-1)+rho_up(i,j)*u_up(i,j)*v_up(i,j)
            temp4=rho_down(i,j-1)*u_down(i,j-1)*H_down(i,j-1)+rho_up(i,j)*u_up(i,j)*H_up(i,j)

            Y_Numerical_minus(1)=0.5d0*(temp1-convective_vector(1))
            Y_Numerical_minus(2)=0.5d0*(temp2-convective_vector(2))
            Y_Numerical_minus(3)=0.5d0*(temp3-convective_vector(3))
            Y_Numerical_minus(4)=0.5d0*(temp4-convective_vector(4))

            !Between j and j+1

            u_ave=u_average(rho_down(i,j),rho_up(i,j+1),u_down(i,j),u_up(i,j+1))
            v_ave=v_average(rho_down(i,j),rho_up(i,j+1),v_down(i,j),v_up(i,j+1))
            H_ave=H_average(rho_down(i,j),rho_up(i,j+1),H_down(i,j),H_up(i,j+1))
            a_ave=a_average(u_ave,v_ave,H_ave)
            lamda(:,:)=0.d0
            lamda(1,1)=abs(v_ave)
            lamda(2,2)=abs(v_ave)
            lamda(3,3)=abs(v_ave-a_ave)
            lamda(4,4)=abs(v_ave+a_ave)
            Y_matrix=Y_coefficience(u_ave,v_ave,H_ave,a_ave)
            Y_inverse=Y_inverse_co(u_ave,v_ave,H_ave,a_ave)


            B_matrix=B_coefficience(u_ave,v_ave,H_ave,a_ave)

            do m=1,4
               do n=1,4
                  temp1=0.d0
                  do o=1,4
                     temp1=Y_matrix(m,o)*lamda(o,n)+temp1
                  enddo
                  Ylamda_matrix(m,n)=temp1
                  !write(*,*) Ans(i,j)
               enddo
            enddo

            do m=1,4
               do n=1,4
                  temp1=0.d0
                  do o=1,4
                     temp1=Ylamda_matrix(m,o)*Y_inverse(o,n)+temp1
                  enddo
                  YlamdaY_matrix(m,n)=temp1
                  !write(*,*) Ans(i,j)
               enddo
            enddo

            q_plus(1)=rho_up(i,j+1)-rho_down(i,j)
            q_plus(2)=u_up(i,j+1)-u_down(i,j)
            q_plus(3)=v_up(i,j+1)-v_down(i,j)
            q_plus(4)=p_up(i,j+1)-p_down(i,j)

            do m=1,4
               temp1=0.d0
               do o=1,4
                  temp1=YlamdaY_matrix(m,o)*q_plus(o)+temp1
               enddo
               convective_vector(m)=temp1
               !write(*,*) Anp(i)
            enddo
            
            temp1=rho_down(i,j)*u_down(i,j)+rho_up(i,j+1)*u_up(i,j+1)
            temp2=rho_down(i,j)*u_down(i,j)**2+p_down(i,j)+rho_up(i,j+1)*u_up(i,j+1)**2+p_up(i,j+1)
            temp3=rho_down(i,j)*u_down(i,j)*v_down(i,j)+rho_up(i,j+1)*u_up(i,j+1)*v_up(i,j+1)
            temp4=rho_down(i,j)*u_down(i,j)*H_down(i,j)+rho_up(i,j+1)*u_up(i,j+1)*H_up(i,j+1)

            Y_Numerical_plus(1)=0.5d0*(temp1-convective_vector(1))
            Y_Numerical_plus(2)=0.5d0*(temp2-convective_vector(2))
            Y_Numerical_plus(3)=0.5d0*(temp3-convective_vector(3))
            Y_Numerical_plus(4)=0.5d0*(temp4-convective_vector(4))

            Y_Flux(1)=Y_Numerical_minus(1)-Y_Numerical_plus(1)
            Y_Flux(2)=Y_Numerical_minus(2)-Y_Numerical_plus(2)
            Y_Flux(3)=Y_Numerical_minus(3)-Y_Numerical_plus(3)
            Y_Flux(4)=Y_Numerical_minus(4)-Y_Numerical_plus(4)
!end vertical==============================================================================================
!==========================================================================================================
            
            qq(1)=rho(i,j)+dt/dx*X_Flux(1)+dt/dy*Y_Flux(1)
            qq(2)=rho(i,j)*u(i,j)+dt/dx*X_Flux(2)+dt/dy*Y_Flux(2)
            qq(3)=rho(i,j)*v(i,j)+dt/dx*X_Flux(3)+dt/dy*Y_Flux(3)
            qq(4)=E(i,j)+dt/dx*X_Flux(4)+dt/dy*Y_Flux(4)

            rhorho(i,j)=qq(1)
            uu(i,j)=qq(2)/qq(1)
            vv(i,j)=qq(3)/qq(1)
            pp(i,j)=(gamma-1.d0)*(qq(4)-qq(1)/2.d0*(uu(i,j)**2+vv(i,j)**2))
            EE(i,j)=qq(4)
            HH(i,j)=(EE(i,j)+pp(i,j))/rhorho(i,j)
         enddo
      enddo
      do j=1,Y_max
         do i=1,X_max
            rho(i,j)=rhorho(i,j)
            u(i,j)=uu(i,j)
            v(i,j)=vv(i,j)
            p(i,j)=pp(i,j)
            E(i,j)=EE(i,j)
            H(i,j)=HH(i,j)
            Mach_number(i,j)=u(i,j)/(gamma*p(i,j)/rho(i,j))
         enddo
      enddo
      if(mod(t,0.01d0)<= dt .and. t< 0.5d0) then
      write(tmpstring,'(i2.2)') int(t/0.01d0+0.5d0)
      open(10, file='data/density_'//trim(tmpstring)//'.d')
      open(11, file='data/U_velocity_'//trim(tmpstring)//'.d')
      open(12, file='data/V_velocity_'//trim(tmpstring)//'.d')
      open(13, file='data/pressure_'//trim(tmpstring)//'.d')
      open(14, file='data/Mach_'//trim(tmpstring)//'.d')
      do i=1,X_max
      place=dx*0.5+dble(i-1)*dx
      write(10, *) place, rho(i,50)
      write(11, *) place, u(i,50)
      write(12, *) place, v(i,50)
      write(13, *) place, p(i,50)
      write(14, *) place, Mach_number(i,50)
      enddo
      close(10)
      close(11)
      close(12)
      close(13)
      close(14)
      write(*, *) tmpString
      end if
 

   enddo
contains
function u_average(rho_left, rho_right, u_left, u_right)
   double precision rho_left, rho_right, u_left, u_right
   double precision u_average
   u_average=(rho_left**(0.5)*u_left+rho_right**(0.5)*u_right)/(rho_left**(0.5)+rho_right**(0.5))
end function u_average

function v_average(rho_left, rho_right, v_left, v_right)
   double precision rho_left, rho_right, v_left, v_right
   double precision v_average
   v_average=(rho_left**(0.5)*v_left+rho_right**(0.5)*v_right)/(rho_left**(0.5)+rho_right**(0.5))
end function v_average

function H_average(rho_left, rho_right, H_left, H_right)
   double precision rho_left, rho_right, H_left, H_right
   double precision H_average
   H_average=(rho_left**(0.5)*H_left+rho_right**(0.5)*H_right)/(rho_left**(0.5)+rho_right**(0.5))
end function H_average

function a_average(u_average, v_average, H_average)
   double precision u_average, v_average, H_average
   double precision a_average
   double precision temp0
   temp0=gamma-1
   a_average=(temp0*(H_average-0.5d0*(u_average**2+v_average**2)))**(0.5)
end function a_average

function X_coefficience(u,v,H,a)
   double precision X_coefficience(4,4)
   double precision u,v,H,a
   X_coefficience(1,1)=0.d0
   X_coefficience(1,2)=1.d0
   X_coefficience(1,3)=1.d0
   X_coefficience(1,4)=1.d0
   X_coefficience(2,1)=0.d0
   X_coefficience(2,2)=u
   X_coefficience(2,3)=u-a
   X_coefficience(2,4)=u+a
   X_coefficience(3,1)=1.d0
   X_coefficience(3,2)=v
   X_coefficience(3,3)=v
   X_coefficience(3,4)=v
   X_coefficience(4,1)=v
   X_coefficience(4,2)=0.5d0*(u**2+v**2)
   X_coefficience(4,3)=H-a*u
   X_coefficience(4,4)=H+a*u
end function X_coefficience

function Y_coefficience(u,v,H,a)
   double precision Y_coefficience(4,4)
   double precision u,v,H,a
   Y_coefficience(1,1)=0.d0
   Y_coefficience(1,2)=1.d0
   Y_coefficience(1,3)=1.d0
   Y_coefficience(1,4)=1.d0
   Y_coefficience(2,1)=1.d0
   Y_coefficience(2,2)=u
   Y_coefficience(2,3)=u
   Y_coefficience(2,4)=u
   Y_coefficience(3,1)=0.d0
   Y_coefficience(3,2)=v
   Y_coefficience(3,3)=v-a
   Y_coefficience(3,4)=v+a
   Y_coefficience(4,1)=u
   Y_coefficience(4,2)=0.5d0*(u**2+v**2)
   Y_coefficience(4,3)=H-a*v
   Y_coefficience(4,4)=H+a*v
end function Y_coefficience

function X_inverse_co(u,v,H,a)
   double precision X_inverse_co(4,4)
   double precision u,v,H,a
   double precision temp0,temp1,temp2,temp3
   temp0=gamma-1.d0
   temp1=temp0/(2.d0*a**2)
   X_inverse_co(1,1)=temp1*(-(2.d0*a**2*v)/temp0)
   X_inverse_co(1,2)=temp1*(0.d0)
   X_inverse_co(1,3)=temp1*((2.d0*a**2)/temp0)
   X_inverse_co(1,4)=temp1*(0.d0)
   X_inverse_co(2,1)=temp1*(-2.d0*H+(4.d0*a**2)/temp0)
   X_inverse_co(2,2)=temp1*(2.d0*u)
   X_inverse_co(2,3)=temp1*(2.d0*v)
   X_inverse_co(2,4)=temp1*(-2.d0)
   X_inverse_co(3,1)=temp1*(H+(a*(u-a))/temp0)
   X_inverse_co(3,2)=temp1*(-u-a/temp0)
   X_inverse_co(3,3)=temp1*(-v)
   X_inverse_co(3,4)=temp1*(1.d0)
   X_inverse_co(4,1)=temp1*(H-(a*(u+a))/temp0)
   X_inverse_co(4,2)=temp1*(-u+a/temp0)
   X_inverse_co(4,3)=temp1*(-v)
   X_inverse_co(4,4)=temp1*(1.d0)
end function X_inverse_co

function Y_inverse_co(u,v,H,a)
   double precision Y_inverse_co(4,4)
   double precision u,v,H,a
   double precision temp0,temp1,temp2,temp3
   temp0=gamma-1.d0
   temp1=temp0/(2.d0*a**2)
   Y_inverse_co(1,1)=temp1*(-(2.d0*a**2*u)/temp0)
   Y_inverse_co(1,2)=temp1*((2.d0*a**2)/temp0)
   Y_inverse_co(1,3)=temp1*(0.d0)
   Y_inverse_co(1,4)=temp1*(0.d0)
   Y_inverse_co(2,1)=temp1*(-2.d0*H+(4.d0*a**2)/temp0)
   Y_inverse_co(2,2)=temp1*(2.d0*u)
   Y_inverse_co(2,3)=temp1*(2.d0*v)
   Y_inverse_co(2,4)=temp1*(-2.d0)
   Y_inverse_co(3,1)=temp1*(H+(a*(v-a))/temp0)
   Y_inverse_co(3,2)=temp1*(-u)
   Y_inverse_co(3,3)=temp1*(-v-a/temp0)
   Y_inverse_co(3,4)=temp1*(1.d0)
   Y_inverse_co(4,1)=temp1*(H-(a*(v+a))/temp0)
   Y_inverse_co(4,2)=temp1*(-u)
   Y_inverse_co(4,3)=temp1*(-v+a/temp0)
   Y_inverse_co(4,4)=temp1*(1.d0)
end function Y_inverse_co

function A_coefficience(u,v,H,a)
   double precision A_coefficience(4,4)
   double precision u,v,H,a
   double precision temp0,temp1,temp2
   temp0=gamma-1
   A_coefficience(1,1)=0.d0
   A_coefficience(1,2)=1.d0
   A_coefficience(1,3)=0.d0
   A_coefficience(1,4)=0.d0
   A_coefficience(2,1)=0.5d0*(u**2+v**2)*temp0-u**2
   A_coefficience(2,2)=u*(3.d0-gamma)
   A_coefficience(2,3)=-v*temp0
   A_coefficience(2,4)=temp0
   A_coefficience(3,1)=-u*v
   A_coefficience(3,2)=v
   A_coefficience(3,3)=u
   A_coefficience(3,4)=0.d0
   A_coefficience(4,1)=0.5d0*u*(-2.d0*H+(u**2+v**2)*temp0)
   A_coefficience(4,2)=H-u**2*temp0
   A_coefficience(4,3)=-u*v*temp0
   A_coefficience(4,4)=u*gamma
end function A_coefficience

function B_coefficience(u,v,H,a)
   double precision B_coefficience(4,4)
   double precision u,v,H,a
   double precision temp0,temp1,temp2
   temp0=gamma-1
   B_coefficience(1,1)=0.d0
   B_coefficience(1,2)=0.d0
   B_coefficience(1,3)=1.d0
   B_coefficience(1,4)=0.d0
   B_coefficience(2,1)=-u*v
   B_coefficience(2,2)=v
   B_coefficience(2,3)=u
   B_coefficience(2,4)=0.d0
   B_coefficience(3,1)=0.5d0*(u**2+v**2)*temp0-v**2
   B_coefficience(3,2)=-u*temp0
   B_coefficience(3,3)=v*(3.d0-gamma)
   B_coefficience(3,4)=temp0
   B_coefficience(4,1)=0.5d0*v*(-2*H+(u**2+v**2)*temp0)
   B_coefficience(4,2)=-u*v*temp0
   B_coefficience(4,3)=H-v**2*temp0
   B_coefficience(4,4)=v*gamma
end function B_coefficience


end program main

subroutine MUSCLlimiter(a,b,c,dl,dr)
  use global
  implicit none
  double precision a,b,c,dl,dr
  double precision dm, dp, dlim, dright,dleft, theta, Rr, Rl, dlimr, dliml

  double precision s

  dm=b-a
  dp=c-b

  !dright=0.5d0*((1d0+kappa)*dm+(1d0-kappa)*dp)
  !dleft =0.5d0*((1d0-kappa)*dm+(1d0+kappa)*dp)
  !theta=dp/dm

  !Rr=min(4d0      /((1d0+kappa)+(1d0-kappa)*theta), 1d0)
  !Rl=min(4d0*theta/((1d0-kappa)+(1d0+kappa)*theta), 1d0)

  !dlimr=Rr*dright
  !dliml=Rl*dleft
  !dlim=0.5d0*(sign(1d0,dlimr)+sign(1d0,dliml))*min(abs(dlimr),abs(dliml))

  !dl=b+phi*0.5d0*dlim
  !dr=b-phi*0.5d0*dlim

  s=(2d0*dp*dm)/(dp**2+dm**2+1d-300)

  dl=b+0.25d0*s*((1d0+kappa*s)*dp+(1d0-kappa*s)*dm)
  dr=b-0.25d0*s*((1d0+kappa*s)*dm+(1d0-kappa*s)*dp)
end subroutine MUSCLlimiter

subroutine SLAU
   use global
   implicit none


end subroutine SLAU
