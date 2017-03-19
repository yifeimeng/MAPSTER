module m_electron
    
    implicit none
    
    contains
    
    
    
      function wavev(e)
      !  this function returns the wavevector in one over lambda, in a-1,
      !  for an input electron energy e in ev.
      
      use m_precision
      
      implicit none
      
      real(fp_kind) c1,c2,e
      real(fp_kind) wavev
      data c1, c2 / 9.7846113e-07_fp_kind, 12.263868_fp_kind /
      
      wavev = sqrt( e + c1 *e ** 2.0_fp_kind ) / c2
      
      end
   
    
    
    
   end module
   