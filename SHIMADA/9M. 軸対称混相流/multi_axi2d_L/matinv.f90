!----------
      subroutine matinv(array)
	  implicit real*8 (a-h,o-z)
!
!     version 2.0  Jan.22, 1992  (Vectorized)
!
      parameter (nv=8)
      dimension array(nv,nv+1)
      integer pivr
!
      m=nv
      n=nv+1
!
! initialization
!
      do k=1, m
!
!        search for pivot element
!
         piv=0.
         pivr = 0
         do i= k, m
            w=abs(array(i,k))
            if(piv.lt.w) then
               piv=w
               pivr=i
            end if
         end do
         if(pivr.eq.0) return
!
!        interchange rows to put pivot element on diagonal
!
         do j=1, n
            if(pivr.ne.k) then
               w=array(pivr,j)
               array(pivr,j)=array(k,j)
               array(k,j)=w
            end if
         end do
!
         piv=array(k,k)
         array(k,k)=1.
!
!        divide pivot row by pivot element
!
         do j=1,n
            array(k,j)=array(k,j)/piv
         end do
!
!        reduce non-pivot rows
!
         do i=1, m
            if(i.ne.k) then
               w=array(i,k)
               array(i,k)=0.
               do j=1, n
                  array(i,j)=array(i,j)-array(k,j)*w
               end do
            end if
         end do
      end do
!
      return
      end
