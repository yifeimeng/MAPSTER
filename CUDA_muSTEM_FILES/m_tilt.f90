module m_tilt

    use m_precision, only: fp_kind
    use global_variables, only: ifactory,ifactorx
    implicit none

    logical :: tilt_illumination = .false.
    real(fp_kind) :: alpha, beta
    
    contains
    
    
    
    subroutine prompt_tilt
        
        use m_user_input, only: get_input
        
        implicit none
        
        integer :: i_tilt, i_cshift
        
10      write(*,*) 'The illumination can be tilted off the beam axis.'
        write(*,*) '<1> Beam tilt'
        write(*,*) '<2> No beam tilt'
        call get_input('<1> Beam tilt <2> No beam tilt', i_tilt)
        write(*,*)
        
        if (i_tilt.eq.1) then
            tilt_illumination = .true.
            
            call setup_tilt
            
        elseif (i_tilt.eq.2) then
            tilt_illumination = .false.
        
        else
            goto 10
            
        endif

20      write(*,*) 'The specimen can be tilted off the specimen axis.'
        write(*,*) '<1> Specimen tilt'
        write(*,*) '<2> No specimen tilt'
        call get_input('<1> Specimen tilt <2> No specimen tilt', i_tilt)
        write(*,*)
        
        if (i_tilt.eq.1) call setup_specimen_tilt
        
    end subroutine
    
    
    
    subroutine setup_tilt
    
        use m_user_input, only: get_input
        use global_variables, only: bvec,ak1,ig1,ig2,ss
        use m_crystallography, only:trimi
        
        implicit none
        
        real(fp_kind)::tilt_theta,tilt_phi,tilt_phi_,tilt_theta_
        
        write(*,*) 'Please enter the beam tilt in mrad:'
        write(*,*)
        
        call get_input('Beam tilt in mrad', tilt_theta)
        
        write(*,*)
        write(*,*) 'Please enter the azimuth of the beam tilt, measured from'
        write(*,*) '[100] ("East") clockwise to [010] ("South") in mrad:'
        write(*,*) 
        call get_input('Beam tilt azimuth in mrad', tilt_phi)
        
        write(*,*)

        bvec  = ak1 * [ sin(tilt_theta*1e-3_fp_kind)*cos(tilt_phi*1e-3_fp_kind)/trimi(ig1,ss), sin(tilt_theta*1e-3_fp_kind)*sin(tilt_phi*1e-3_fp_kind)/trimi(ig2,ss),0_fp_kind]
        bvec = float(nint(bvec*[ifactory,ifactorx,0_fp_kind]))/[ifactory,ifactorx,1_fp_kind]
        tilt_theta_ = asin(sqrt(sum((bvec*[trimi(ig1,ss),trimi(ig1,ss),0_fp_kind])**2))/ak1)*1e3
        tilt_phi_ = atan2(bvec(2),bvec(1))*1e3
        write(6,30) tilt_theta,tilt_theta_,tilt_phi,tilt_phi_
        write(*,*)
                
30     format(&
        &1x,'To ensure that the illumination wave function is continuous',/,&
        &1x,'at the boundary of the simulation grid, the beam tilt has',/,&
        &1x,'been rounded from ',f5.1,' to ',f5.1,' mrad and the azimuth has ',/,&
        &1x,'been rounded from ',f7.1,' to ',f7.1,' mrad. If you would prefer ',/,&
        &1x,'that the tilt used in simulation was closer to the value ',/,&
        &1x,'inputted, please consider increasing the dimensions of the',/,&
        &1x,'supercell used in the simulation by increasing the unit cell',/,&
        &1x,'tiling.')
        
        
    end subroutine
     
    subroutine setup_specimen_tilt
        
        use global_variables, only: claue,ak1,Kz,ss,ig1,ig2
        use m_crystallography
        use m_user_input, only: get_input
        
        implicit none
        
        real(fp_kind)::tilt_theta,tilt_phi
        
        
    
        write(*,*) 'Please enter the specimen tilt in mrad:'
        write(*,*)
        
        call get_input('Specimen tilt in mrad', tilt_theta)
        
        write(*,*)
        write(*,*) 'Please enter the azimuth of the specimen tilt, measured from'
        write(*,*) '[100] ("East") clockwise to [010] ("South") in mrad:'
        write(*,*) 
        call get_input('Specimen tilt azimuth in mrad', tilt_phi)
        
        write(*,*)

        Kz   = ak1 * cos( tilt_theta*1e-3_fp_kind )
        claue  = ak1 * [ sin(tilt_theta*1e-3_fp_kind)*cos(tilt_phi*1e-3_fp_kind)/trimi(ig1,ss), sin(tilt_theta*1e-3_fp_kind)*sin(tilt_phi*1e-3_fp_kind)/trimi(ig2,ss),0_fp_kind]
    end subroutine   
    
    subroutine tilt_wave_function(psi)
        implicit none
        
        complex(fp_kind) :: psi(:,:)
        
        if (tilt_illumination) then
            call tilt_wave_function_phase_ramp(psi)
        endif
        
    end subroutine
    
    
    
    subroutine tilt_wave_function_phase_ramp(psi)
        
        use global_variables, only: ifactory, ifactorx, pi,bvec
        use output
        implicit none
        
        complex(fp_kind) :: psi(:,:)
        
        complex(fp_kind),allocatable:: psi2(:,:)
        integer :: nopiy, nopix
        real(fp_kind) :: shift(3), shift_frac(3),bvec_(3)
        integer :: ny, nx
        
        nopiy = size(psi, 1)
        nopix = size(psi, 2)
        allocate(psi2(nopiy,nopix))
        shift = [alpha*ifactorx, beta*ifactory, 0.0_fp_kind]
        shift_frac = shift / [nopiy, nopix, 1]
        bvec_ = bvec/[nopiy,nopix,1.0_fp_kind]*[ifactory,ifactorx,0.0_fp_kind]
        bvec_ = [bvec_(2),bvec_(1),0]
        !write(*,*) bvec_
        !$OMP PARALLEL DO PRIVATE(nx, ny)
        do nx = 1, nopix
        do ny = 1, nopiy
                
            psi(ny,nx) = psi(ny,nx) * exp(cmplx(0.0_fp_kind, 2*pi*dot_product([ny-1, nx-1, 0], bvec_), fp_kind))
        enddo
        enddo
        !$OMP END PARALLEL DO
        !call binary_out(nopiy,nopix,atan2(imag(psi2),real(psi2)),'phase_ramp')                
    end subroutine
    
    subroutine tilt_wave_function_shift(psi)
    
        use global_variables, only: ifactory, ifactorx, pi,ig1,ig2,bvec
        use output
        use CUFFT_wrapper
        
        implicit none
        
        complex(fp_kind) :: psi(:,:)
        
        integer :: nopiy, nopix,shiftx,shifty
        real(fp_kind) :: shift(3), shift_frac(3),ky(3),kx(3)!,bvec(3)
        integer :: ny, nx,m1,m2
        real(fp_kind),allocatable::shift_array(:,:)
        complex(fp_kind),allocatable :: shift_grate(:,:)
        integer,allocatable,dimension(:)::fftfreqy,fftfreqx
        
        nopiy = size(psi, 1)
        nopix = size(psi, 2)

        allocate(shift_array(nopiy,nopix),fftfreqy(nopiy),fftfreqx(nopix),shift_grate(nopiy,nopix))
        
        shifty = (nopiy-1)/2-1
        shiftx = (nopix-1)/2-1
        shift_array = 0
        fftfreqy = fftfreq(nopiy)
        fftfreqx = fftfreq(nopix)
        
        do ny = 1, nopiy
            write(*,*) abs(fftfreqy(ny)-bvec(1)*ifactory)
            if(abs(fftfreqy(ny)-bvec(1)*ifactory)<1) then
                do nx =1,nopix
                    write(*,*) abs(fftfreqx(nx)-bvec(2)*ifactorx)
                    if(abs(fftfreqx(nx)-bvec(2)*ifactorx)<1) then
                        shift_array(ny,nx)=1                       
                        exit
                    endif
                enddo
            exit
            endif
        enddo
        call binary_out_unwrap(nopiy,nopix,shift_array,'unwrapped_shift_array')
        call binary_out(nopiy,nopix,quad_shift(shift_array),'shift_array')
        call fft2(nopiy,nopix,cmplx(shift_array,kind=fp_kind),nopiy,shift_grate,nopiy)
        call binary_out(nopiy,nopix,atan2(imag(shift_grate),real(shift_grate)),'phase_grate')
        
    end subroutine
    
    function fftfreq(n)
    integer*4,intent(in)::n
    integer*4::fftfreq(n),i
    
    if (mod(n,2)==0) then
        !Even case
        do i=1,n
            fftfreq(i) = -n/2+i
        enddo
    else
        !Odd case
        do i=1,n
            fftfreq(i) = -n/2+i-1
        enddo
    endif
    end function
    
    function quad_shift(array)
        real(fp_kind)::array(:,:)
        integer*4::nopiy,nopix
        real(fp_kind)::quad_shift(size(array, 1),size(array, 2))
    
        nopiy = size(array, 1)
        nopix = size(array, 2)
        
        quad_shift = cshift(cshift(array, nopix/2,dim = 2),nopiy/2,dim = 1)
    
    end function
    !
end module
