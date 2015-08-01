module prmtr
   implicit none
   double precision, parameter :: pi=2.d0*acos(0.d0)
   double precision, parameter :: gamma=1.4d0
   double precision, parameter :: R_gas=287.1d0
   !Grid number
   integer         , parameter :: ni=224
   integer         , parameter :: nj=100
   !MUSCL INCREMENT
   double precision, parameter :: phi=1.d0
   double precision, parameter :: kappa=1.d0/3.d0
   !Division
   double precision, parameter :: dx=0.01d0
   double precision, parameter :: dy=0.01d0
   double precision, parameter :: dt=0.001d0
   integer         , parameter :: X_max=1d0/dx
   integer         , parameter :: Y_max=1d0/dy
   integer         , parameter :: T_max=1d0/dt
end module prmtr

