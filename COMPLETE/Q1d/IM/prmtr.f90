module prmtr
   implicit none
   double precision, parameter :: pi=2.d0*acos(0.d0)
   double precision, parameter :: mu=5d-2 !viscous term
   double precision, parameter :: gamma=1.4d0
   double precision, parameter :: gas_constant=8.3144621d0
   double precision, parameter :: molar_mass=28.96d0
   double precision, parameter :: gas_spec=gas_constant/molar_mass*1.d3
   integer,parameter :: CellNum = 151
   !Grid number
   !integer         , parameter :: ni=224
   !integer         , parameter :: nj=100
   !MUSCL INCREMENT
   double precision, parameter :: phi=1.d0
   double precision, parameter :: kappa=1.d0/3.d0
   !Division
   double precision, parameter :: CFL= 1d1
   !STEP PARAMETER
   !integer, parameter :: step_forward =101
   !integer, parameter :: step_backward=103
   !integer, parameter :: step_height  =181
   !OUTPOT TIME
   integer, parameter :: out_time = 300
   integer, parameter :: time_max = 500000000
   !RESIDUAL
   double precision, parameter :: epsilon= 1d-10
end module prmtr
