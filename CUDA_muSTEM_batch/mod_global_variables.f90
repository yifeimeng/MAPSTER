    module global_variables
    
    use m_precision, only: fp_kind
    
    implicit none
    
    integer(4) :: nt,nm,i_xtl                   !number of atom types, number of atoms, xtl file read flag
    integer(4) :: nopiy,nopix,npixels           !supercell
    integer(4) :: nopix_ucell,nopiy_ucell       !unit cell
    integer(4) :: ifactorx,ifactory             !unit cell tilings
    real(fp_kind) :: deltay,deltax              !real space spacing between pixels in x and y
    real(fp_kind) :: normalisation
    
    real(fp_kind), allocatable :: bwl_mat(:,:)            !bandwidth limiting matrix
                                          
    integer(4), allocatable :: nat(:)                     !number of each atom type in the unit cell
    real(fp_kind), allocatable :: tau(:,:,:)              !position of the atoms in the unit cell
    real(fp_kind), allocatable :: atf(:,:)                !atomic number, occupancy and DWF (urms)
    real(fp_kind), allocatable :: atomf(:,:),fx(:)        !electron scattering factor parameterisation from elsa
    real(fp_kind)  :: a0(3),deg(3),ekv,ss(7)              !a b c unit cell lengths, angle between a b c, accelerating voltage, tricyclinc info
    real(fp_kind)  :: thetad,surfn(3),orthog(3,3)
    real(fp_kind)  :: volts,ak                            !mean inner potential, wavevector (corrected for refraction)
    real(fp_kind)  :: ak1,relm                            !wavevector in freespace, relativistically corrected mass
    real(fp_kind)  :: claue(3),Kz                         !Specimen tilt vector and z component of incident wave vector
    real(fp_kind)  :: bvec(3)                             !Beam tilt vector
    
    !sample thickness and slicing variables
    real(fp_kind) :: thickness                                   !sample thickness
    integer(4) :: n_cells                                 !number of unit cells

    complex(fp_kind), allocatable :: prop(:,:,:)          !propagator
    complex(fp_kind), allocatable :: fz(:,:,:)            !the scattering factors, in reciprocal space, calculated on the grid (supercell)
    complex(fp_kind), allocatable :: fz_DWF(:,:,:)        !the DWF smear_array, in reciprocal space, calculated on the grid (supercell)
    complex(fp_kind), allocatable :: sinc(:,:)            !sinc function to correct for pixelation in the potential construction
    complex(fp_kind), allocatable :: inverse_sinc(:,:)    !1/sinc function to correct for pixelation in the potential construction

    real(fp_kind)  :: uvw1(3),uvw2(3)                     !real space scan vectors that are parallel 
    integer(4) :: ig1(3),ig2(3),izone(3)                  !to the reciprocal space vectors ig1,ig2.
        
    character*120 :: substance                            !label for the crystal substance
    character*10, allocatable :: substance_atom_types(:)

    !output variables
    integer(4) :: ndet                                    !number of integrating detectors
    real(fp_kind), allocatable :: outer(:),inner(:)       !detector ranges (in inverse angstrom)
    
    !interpolation variables
    integer(4) :: output_nopiy,output_nopix               !output number of pixels in the interpolated output image
    integer(4) :: tilex,tiley                             !Interpolated image tiling
    real(fp_kind)  ::  bwl_rad                            !band width limiting radius (in q space)
                                                
    !Constants data
    real(fp_kind),parameter :: pi = atan(1.0_fp_kind)*4.0_fp_kind 
    real(fp_kind),parameter :: tp = atan(1.0_fp_kind)*8.0_fp_kind 
    real(fp_kind),parameter :: hsq_on_twom = 150.4132_fp_kind  !h^2/2m
    complex(fp_kind) :: ci = cmplx(0.0_fp_kind,1.0_fp_kind)
    
    real(fp_kind),parameter :: const1 = 9.7846113e-07_fp_kind       ! const1 = 1/(2mc^2) in eV^(-1)
    real(fp_kind),parameter :: const2 = 12.263868_fp_kind           ! const2 = h  in eV.A
    real(fp_kind),parameter :: bohr = 0.529177_fp_kind              ! bohr radius
    real(fp_kind),parameter :: ryd = 13.60535_fp_kind               ! rydberg constant
    real(fp_kind),parameter :: fsc = 7.29735e-3_fp_kind             ! fsc = fine structure constant (dimensionless)
    real(fp_kind),parameter :: hbarc = 1973.26_fp_kind              ! hbarc = hbar * c in eV A units
    
    !logical types to pick inelastic calculations
    logical :: ionization = .false.
    logical :: adf 
    logical :: EELS = .false.
    
    logical :: qep
    
    logical :: on_the_fly = .false.
    logical :: high_accuracy
    
	contains
     
    
    
      subroutine constants()
      
          use m_precision
          use m_electron
      
          implicit none
      
      
          real(fp_kind) con1,c2,delk
          data con1, c2 / 511.00657_fp_kind, 1.9569222e-03_fp_kind /

          relm = ( con1 + ekv + 0.001_fp_kind * volts ) * c2
          
    ! ak is the incident wave vector in the solid corrected for refraction
          ak   = wavev( ekv * 1000_fp_kind + volts)
          
    ! ak1 is the incident wave vector without corrections for refraction
          ak1   = wavev( ekv * 1000_fp_kind )
    !Initialise tilt to zero (on-axis)
          Kz = ak1
          claue = 0_fp_kind
          
          delk = ak - ak1

	    write(6,*) 'Pertinent quantities for the crystal:'
          write(6,111) ekv,ak,char(143),ak1,char(143),delk,char(143),relm
      111 format('          E = ',F12.3,' keV',/,&
         &' Crystal Ko = ',G20.10,1x,a1,'-1 (1/lambda)',/,&
         &'  Vacuum Ko = ',G20.10,1x,a1,'-1 (1/lambda)',/,&
         &'   Delta Ko = ',G20.10,1x,a1,'-1 (1/lambda)',/,&
         &'     m* / m = ',g12.5,&
         &' (relativistic mass increment)',/,/)

      end
      
    
      end module global_variables
