module m_elsa 
    
    implicit none
      
    interface elsa_ext
        module procedure single_elsa_ext
        module procedure double_elsa_ext
    end interface
    
    
    contains
    
    
      
!--------------------------------------------------------------------------------------
!     This function returns the electron scattering factors in
!     ANGSTROM units. The atom type is given by the index K, and
!     the scattering vector squared (in A ^ (-2)) by S2. This
!     has no relativistic mass correction at this stage. This
!     correction is carried out elsewhere in the program (in CCD
!     or RELM).
!
!     Uses 5 gaussian summation coefficients from Waasmaier and Kirfel
!     modified 4/6/96 for ions
!	elsa_ext() determines the electron scattering factor using the
!	shielded coulomb extrapolation for scattering vectors beyond
!	the limits for which the parameterisation is fitted.
!	
!	see Elsa for further descriptions, as they are essentially the
!	same function (AJD)

	
    function double_elsa_ext(nt,k,atomf,s2)
    
        implicit none
      
        integer(4) nt,m,n,k
        real(8) atomf(13,nt),s2,deltak,alpha,total,xbs
	    real(8) c1,fxray,s2l,a2
        real(8) double_elsa_ext
        data c1/2.393367d-02/

        ! alpha = inverse yukawa range in A (only used for small s2 in ions)      
        alpha=0.02d0

        if(s2.lt.1.0d-03) then
                ! Small scattering vector limit
                total=atomf(11,k)
                do n=1,5
                      total=total+atomf(n,k)
                enddo
                ! if deltak = 0 then we have a neutral atom else deltak = Z - sum a_i - c      
                deltak=float(nint(atomf(12,k)-total))
                total = 0.0d0
                do n= 1, 5
	                m = n + 5
	                total = total + atomf(n,k) * atomf(m,k)* (1.0 - atomf(m,k)/2.0d0*s2)
                enddo
                total=total+deltak/(s2+alpha**2d0)
                double_elsa_ext = c1 * total
        else
                if (s2.gt.36.0d0) then
		          ! Large scattering vector limit
	                s2l = s2 + atomf(13,k)
                      double_elsa_ext = c1 * atomf(12,k) / s2l
		          else
       	          ! scattering vector in 
		          ! parameterisation range
                      fxray = atomf(11,k)
                      do n = 1, 5
	                      m = n + 5
	                      xbs = - atomf(m,k) * s2
	                      xbs = exp(xbs)
	                      fxray = fxray + atomf(n,k) * xbs
                      enddo
                      double_elsa_ext = c1 * (atomf(12,k) - fxray) / s2
                endif
        endif

    end function double_elsa_ext
    
    
    
    function single_elsa_ext(nt,k,atomf,s2)
    
        implicit none
      
        integer(4) nt,m,n,k
        real(4) atomf(13,nt),s2,deltak,alpha,total,xbs
	    real(4) c1,fxray,s2l,a2
        real(4) single_elsa_ext
        data c1/2.393367e-02/

        ! alpha = inverse yukawa range in A (only used for small s2 in ions)      
        alpha=0.02

        if(s2.lt.1.0e-03) then
                ! Small scattering vector limit
                total=atomf(11,k)
                do n=1,5
                      total=total+atomf(n,k)
                enddo
                ! if deltak = 0 then we have a neutral atom else deltak = Z - sum a_i - c      
                deltak=float(nint(atomf(12,k)-total))
                total = 0.0
                do n= 1, 5
	                m = n + 5
	                total = total + atomf(n,k) * atomf(m,k)* (1.0 - atomf(m,k)/2.0*s2)
                enddo
                total=total+deltak/(s2+alpha**2)
                single_elsa_ext = c1 * total
        else
                if (s2.gt.36.0) then
		          ! Large scattering vector limit
	                s2l = s2 + atomf(13,k)
                      single_elsa_ext = c1 * atomf(12,k) / s2l
		          else
       	          ! scattering vector in 
		          ! parameterisation range
                      fxray = atomf(11,k)
                      do n = 1, 5
	                      m = n + 5
	                      xbs = - atomf(m,k) * s2
	                      xbs = exp(xbs)
	                      fxray = fxray + atomf(n,k) * xbs
                      enddo
                      single_elsa_ext = c1 * (atomf(12,k) - fxray) / s2
                endif
        endif

    end function single_elsa_ext
    
    
    
end module m_elsa
