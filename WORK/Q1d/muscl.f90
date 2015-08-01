subroutine MUSCLlimiter(a,b,c,dl,dr)
   use prmtr
   implicit none

   double precision,dimension(4)::a,b,c,dl,dr
   double precision,dimension(4)::dm, dp, dlim, dright,dleft, theta, Rr, Rl, dlimr, dliml

   double precision,dimension(4)::s
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

   !dr=b+phi*0.5d0*dlim
   !dl=b-phi*0.5d0*dlim

   s=(2d0*dp*dm+1d-30)/(dp**2+dm**2+1d-30)

   dr=b+0.25d0*s*((1d0+kappa*s)*dp+(1d0-kappa*s)*dm)
   dl=b-0.25d0*s*((1d0+kappa*s)*dm+(1d0-kappa*s)*dp)
end subroutine MUSCLlimiter
