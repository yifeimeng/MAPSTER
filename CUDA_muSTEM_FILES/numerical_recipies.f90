!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    SUBROUTINE spline(x,y,n,yp1,ypn,y2)
    use m_precision
    integer(4) n
    real(fp_kind) yp1,ypn,x(n),y(n),y2(n)
    PARAMETER (NMAX=641)
    integer(4) i,k
    real(fp_kind) p,qn,sig,un,u(n)
    
	if (yp1.gt..99e30_fp_kind) then
       y2(1)=0.0_fp_kind
       u(1)=0.0_fp_kind
    else
        y2(1)=-0.5_fp_kind
        u(1)=(3.0_fp_kind/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
    endif
    do i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.0_fp_kind
        y2(i)=(sig-1.0_fp_kind)/p
        u(i)=(6.0_fp_kind*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
    enddo
	continue
      if (ypn.gt..99e30_fp_kind) then
        qn=0.0_fp_kind
        un=0.0_fp_kind
      else
        qn=0.5_fp_kind
        un=(3.0_fp_kind/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.0_fp_kind)
      do k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
	  enddo
    continue
    return
    END
!  (C) Copr. 1986-92 Numerical Recipes Software ]2+18Z9.

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    SUBROUTINE splint(xa,ya,y2a,n,x,y)
    use m_precision
    integer(4) n
    real(fp_kind) x,y,xa(n),y2a(n),ya(n)
    integer(4) k,khi,klo
    real(fp_kind) a,b,h
    klo=1
    khi=n
 1  if (khi-klo.gt.1) then
      k=(khi+klo)/2
      if(xa(k).gt.x)then
        khi=k
      else
        klo=k
      endif
    goto 1
    endif
    h=xa(khi)-xa(klo)
    if (h.eq.0.0_fp_kind) pause 'bad xa input in splint'
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.0_fp_kind
    return
    END
!  (C) Copr. 1986-92 Numerical Recipes Software ]2+18Z9.

!--------------------------------------------------------------------------------------
      function gasdev(idum)
      !	Numerical Recipies function,
      !	Returns a normally distributed deviate with zero mean
      !	and unit variance. Uses ran1 as a source of uniform deviates 
      use m_precision

	implicit none

	integer(4) iset, idum
	real(fp_kind) v1, v2, ran1, r, fac, gset
      real(fp_kind) gasdev
	save iset, gset
	data iset/0/

1   if( iset .eq. 0 ) then
       v1 = 2.0_fp_kind*ran1(idum)-1.0_fp_kind
	   v2=2.0_fp_kind*ran1(idum)-1.0_fp_kind
	   r=v1*v1 + v2*v2
	   if(r .ge. 1.0_fp_kind) goto 1
	   fac=sqrt(-2.0_fp_kind*log(r)/r)
	   gset=v1*fac
	   gasdev=v2*fac
	   iset=1
	else
	   gasdev=gset
	   iset=0
	endif

	return
	end

!--------------------------------------------------------------------------------------
      FUNCTION ran1(idum)
      !	random number generator from numerical recipies
      !	Returns UNIFORM random number
      !	Pg 271
      use m_precision
      implicit none

      INTEGER(4) idum,IA,IM,IQ,IR,NTAB,NDIV

      real(fp_kind) AM,EPS,RNMX,ran1
      PARAMETER (IA=16807,IM=2147483647,IQ=127773,IR=2836,NTAB=32)
      INTEGER*4 j,k,iv(NTAB),iy
      SAVE iv,iy
      DATA iv /NTAB*0/, iy /0/
      
      NDIV=1+(IM-1)/NTAB
      AM=1.0_fp_kind/IM
      EPS=3.0e-16_fp_kind
      RNMX=1.e0_fp_kind-EPS

      if (idum.le.0.or.iy.eq.0) then
        idum=max(-idum,1)
        do j=NTAB+8,1,-1
          k=idum/IQ
          idum=IA*(idum-k*IQ)-IR*k
          if (idum.lt.0) idum=idum+IM
          if (j.le.NTAB) iv(j)=idum
        enddo
        iy=iv(1)
      endif
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      ran1=min(AM*iy,RNMX)
      return
      END
      
      
!c----------------------------------------------------------------------
!c     Subroutine trapzd taken from numerical recipes in F77 2nd ED
!c     p. 131
!c     calculates the integral of a function using the simple trapezoidal method
!c     func is the function to be integrated
!C     a and b and the limits of integration
!c     s is the result of the integration
!c     When n=1 gives the crudest estimate of the integral, subsequent
!c     calls with n=2,3,... in that order will improve the accuracy
!c     of s by adding 2^(n-2) additional interior points.
!c     s should not be modified between sequential calls
!c----------------------------------------------------------------------

	recursive subroutine trapzd(func,a,b,s,n)
	implicit none

	integer(4) n
	real(8) a,b,s,func
	
	external func

	integer*4 it, j
	real(8) del, sum, tnm, x
	
	if (n.eq.1) then
	   s = 0.5d0*(b-a)*(func(a)+func(b))
	else
	   it = 2**(n-2)
	   tnm=it
	   del=(b-a)/tnm
	   x=a+0.5d0*del
	   sum=0.0d0
	   do j=1, it
	      sum=sum+func(x)
	      x=x+del
	   enddo
	   s=0.5d0*(s+(b-a)*sum/tnm)
	endif
	return
	end



!c----------------------------------------------------------------------
!c     Subroutine qromb taken from numerical recipes in F77 2nd ED
!c     p. 134
!c     calculates the integral of the function func using romberg integration
!c     a and b are the limits of the integral and s is the result
!c 
!c     Originally the test for convergence was
!c     if (abs(dss) .le. rel_tol*abs(ss)) 
!c     but this fails if the integral goes to zero
!c     so the test was changed to have an absolute and a relative tolerance
!c----------------------------------------------------------------------

	recursive subroutine qromb(func,a,b,ssum)
	implicit none

	integer(4) JMAX, JMAXP, K, KM
	real(8)  func, ssum, abs_tol, rel_tol
    real(8) :: a, b
	external func
	!parameter(abs_tol=1.0e-10_fp_kind,rel_tol=1.0e-6_fp_kind,JMAX=40,JMAXP=JMAX+1,K=5,KM=K-1)
    !parameter(abs_tol=1.0e-10_fp_kind,rel_tol=1.0e-5_fp_kind,JMAX=40,JMAXP=JMAX+1,K=5,KM=K-1)

	integer(4) :: j
	
	real(8) :: dss
	real(8),allocatable :: h(:),s(:)
	
	!if(fp_kind.eq.Double) then
	!abs_tol=1.0e-10
	!rel_tol=1.0e-6
	!else
	    abs_tol=1.0e-10
	    rel_tol=5.0e-6
	!endif
	
	JMAX=40
	JMAXP=JMAX+1
	K=5
	KM=K-1
	
	if(allocated(h)) deallocate(h)
	allocate(h(JMAXP))
	if(allocated(s)) deallocate(s)
	allocate(s(JMAXP))
	!real(fp_kind) :: dss, h(JMAXP), s(JMAXP)

	h(1) = 1.0d0
	do j=1,JMAX
	   call trapzd(func,a,b,s(j),j)
	   if (j .ge. k) then
	      call polint(h(j-KM), s(j-KM), K, 0.0d0, ssum, dss)
	      if (abs(dss) .le. max(abs_tol, rel_tol*abs(ssum))) then
	           return
	      endif
	   endif
	   s(j+1)=s(j)
	   h(j+1)=0.25d0*h(j)
	enddo
!c	pause 'too many steps in qromb'
	END

!c----------------------------------------------------------------------
!c     Subroutine polint taken from numerical recipes in F77 2nd ED
!c     p. 104
!c     calculates polynomial interpolation
!c----------------------------------------------------------------------

	recursive subroutine polint(xa,ya,n,x,y,dy)
    implicit none

	integer(4) n, NMAX
	real(8) dy,x,y,xa(n),ya(n)
	parameter (NMAX=10)  ! Largest anticipated value of n

!c	Given arrays xa nd ya, each of length n, and given a value of x,
!c     this routine returns a value y, and an error estimate dy. If P(x)
!c     is the polynomial of degree N-1 such that P(xa_i)=ya_i, i=1,...,n,
!c     then the returned value y=P(x)
	integer(4) i,m,ns
	real(8) den, dif, dift, ho,hp,w,c(NMAX),d(NMAX)
	ns=1
	dif=abs(x-xa(1))

	do i=1, n  ! Here we find the index ns of the closest table entry
	   dift=abs(x-xa(i))
	   if (dift.lt.dif) then
	      ns=i
	      dif=dift
	   endif
	   c(i)=ya(i)   !and initialize the tableau of c's and d's
	   d(i)=ya(i)
	enddo

	y=ya(ns)
	ns=ns-1        !this is the initial approximation to y.
	do m=1, n-1    !for each column of the tableau, we loop over the currect
	               !c's and d's and update them
	   do i=1, n-m
	      ho=xa(i)-x
	      hp=xa(i+m)-x
	      w=c(i+1)-d(i)
	      den=ho-hp
	      if(den .eq. 0.0d0) pause 'failure in polint'
!              This error can occur only if two input xa's are (to within roudoff) identical.
	      den=w/den
	      d(i)=hp*den   !here the c's and d's are updated.
	      c(i)=ho*den
	   enddo
	   if(2*ns .lt. n-m) then
            dy=c(ns+1)
	   else
	      dy=d(ns)
	      ns=ns-1
	   endif
	   y=y+dy
	enddo
	return
	end

      SUBROUTINE splin2(x1a,x2a,ya,y2a,m,n,x1,x2,y)
      use m_precision
      implicit none
      INTEGER(4) m,n,NN
      REAL(fp_kind) x1,x2,y,x1a(m),x2a(n),y2a(m,n),ya(m,n)
      PARAMETER (NN=1000)
!CU    USES spline,splint
      INTEGER(4) j,k
      REAL(fp_kind) y2tmp(n),ytmp(n),yytmp(n)
      
      
      !$OMP PARALLEL DO PRIVATE(ytmp,y2tmp) SHARED(x2a,yytmp)
      do j=1,m
        do k=1,n
          ytmp(k)=ya(j,k)
          y2tmp(k)=y2a(j,k)
        enddo
        call splint(x2a,ytmp,y2tmp,n,x2,yytmp(j))
      enddo
      !$OMP END PARALLEL DO
      
      call spline(x1a,yytmp,m,1.e30_fp_kind,1.e30_fp_kind,y2tmp)
      call splint(x1a,yytmp,y2tmp,m,x1,y)
      return
      END

      SUBROUTINE splie2(x1a,x2a,ya,m,n,y2a)
      use m_precision
      implicit none
      INTEGER(4) m,n,NN
      REAL(fp_kind) x1a(m),x2a(n),y2a(m,n),ya(m,n)
      PARAMETER (NN=10000)
!CU    USES spline
      INTEGER(4) j,k
      REAL(fp_kind) y2tmp(n),ytmp(n)
      do j=1,m
        do k=1,n
          ytmp(k)=ya(j,k)
11      enddo
        call spline(x2a,ytmp,n,1.0e31_fp_kind,1.e31_fp_kind,y2tmp)
        
        do k=1,n
          y2a(j,k)=y2tmp(k)
12      enddo
13    enddo
      return
      END
      

      SUBROUTINE polin2(x1a,x2a,ya,m,n,x1,x2,y,dy)
      implicit none
      INTEGER(4) m,n,NMAX,MMAX
      REAL(8) dy,x1,x2,y,x1a(m),x2a(n),ya(m,n)
      PARAMETER (NMAX=20,MMAX=20)
      !USES polint
      INTEGER(4) j,k
      REAL(8) ymtmp(MMAX),yntmp(NMAX)
      do j=1,m
        do k=1,n
          yntmp(k)=ya(j,k)
        enddo
        call polint(x2a,yntmp,n,x2,ymtmp(j),dy)
      enddo
      call polint(x1a,ymtmp,m,x1,y,dy)
      return
      END