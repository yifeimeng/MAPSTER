function qep_stem_GPU_memory() result(required_memory)
    
    use m_precision, only: fp_kind
    use global_variables, only: nopiy, nopix, ifactory, ifactorx, on_the_fly, ndet, ionization, eels
    use m_lens, only: imaging
    use m_slicing, only: n_slices
    use m_qep, only: n_qep_grates, quick_shift, phase_ramp_shift
    
    implicit none
    
    real(fp_kind) :: required_memory
    
    real(fp_kind) :: array_size
    integer :: array_count
    
    array_size = 4.0_fp_kind * nopiy * nopix
    
    array_count = 2*(4 + n_slices + n_slices*n_qep_grates) + 2 + ndet
        
    if (on_the_fly.or.quick_shift.or.phase_ramp_shift) array_count = array_count + 2
    if (phase_ramp_shift) array_count = array_count + 2
    if (ionization) array_count = array_count + 1 + n_slices
    if (eels) array_count = array_count + 1
    
    ! Temporary array used in cuda_stem_detector()
    array_count = array_count + 1
    
    required_memory = array_count * array_size
    
    if (phase_ramp_shift) required_memory = required_memory + 8.0_fp_kind*(nopiy*ifactory + nopix*ifactorx)
    
end function

    

subroutine qep_stem
    
    use global_variables
    use m_lens
    use m_user_input
    use m_qep
    use m_precision
    use output
    use cuda_array_library
    use cudafor
    use cuda_ms
    use CUFFT
	use cufft_wrapper
    use local_ionization
    use cuda_potential
    use m_slicing
    use cuda_setup
    use m_probe_scan, only: nysample, nxsample, probe_positions, scan_quarter
    use m_tilt, only: tilt_wave_function
    use m_string, only: to_string
    use m_multislice, only: make_qep_grates, output_probe_intensity, output_cell_list, cell_map, output_thickness_list, setup_propagators
    use m_potential, only: precalculate_scattering_factors
    
    implicit none
    
    !dummy variables
    integer(4) ::  i,j,l,m,i_qep_pass
    integer(4) ::  count
    integer(4) ::  shifty,shiftx
    integer(4) ::  ny,nx,i_df,idet
    integer(4) ::  total_slices
    
    !random variables
    integer(4) :: idum
    real(fp_kind) :: ran1
    complex(fp_kind) :: temp_transf(nopiy,nopix)
       
    !probe variables
    complex(fp_kind) :: psi(nopiy,nopix)
    complex(fp_kind) :: trans(nopiy,nopix), psi_initial(nopiy,nopix)
     
    !output
    real(fp_kind) :: cbed(nopiy,nopix)
    real(fp_kind) :: image(nopiy,nopix)
    
    !STEM image variables
    real(fp_kind) :: masks(nopiy,nopix,ndet)                       !detector masks
    real(fp_kind) :: stem_image(nysample,nxsample,n_df,ndet)
    real(fp_kind) :: stem_elastic_image(nysample,nxsample,n_df,ndet)
    real(fp_kind) :: stem_inelastic_image(nysample,nxsample,n_df)
    real(fp_kind) :: stem_ion_image(nysample,nxsample,n_df)
    real(fp_kind) :: eels_correction_image(nysample,nxsample,n_df)
    real(fp_kind) :: temp_image(nopiy,nopix)
    
    !diagnostic variables
    real(fp_kind) :: intensity
    real(fp_kind) :: t1, delta
    
    !output variables
    character(120) :: fnam, fnam_temp, fnam_det

    !device variables
	integer :: plan
    complex(fp_kind),device,dimension(nopiy,nopix) :: psi_d,psi_elastic_d,psi_out_d,psi_initial_d
    real(fp_kind),device,dimension(nopiy,nopix) :: cbed_d,temp_d
	complex(fp_kind),device,allocatable :: prop_d(:,:,:),transf_d(:,:,:,:)
	complex(fp_kind),device,allocatable :: trans_d(:,:), shift_arrayx_d(:,:),shift_arrayy_d(:,:),shift_array_d(:,:)
    real(fp_kind),device,allocatable :: masks_d(:,:,:),ion_potential_d(:,:,:)
    real(fp_kind),device,allocatable :: ion_image_d(:,:),eels_correction_detector_d(:,:)
    
    !device variables for on the fly potentials
    complex(fp_kind),device,allocatable,dimension(:,:) :: bwl_mat_d,inverse_sinc_d,fz_mu_d,fz_dwf_d
    complex(fp_kind),device,allocatable,dimension(:,:,:) :: fz_d
    real(fp_kind),device,allocatable :: inelastic_potential_d(:,:)

    real(fp_kind),allocatable :: probe_intensity(:,:,:)
    character(1024) :: filename
    
    real(fp_kind)::qep_stem_GPU_memory
    
    
    
    call GPU_memory_message(qep_stem_GPU_memory(), on_the_fly)
    
    
    
    write(*,*) '|----------------------------------|'
	write(*,*) '|      Pre-calculation setup       |'
	write(*,*) '|----------------------------------|'
    write(*,*)
            
    ! Precalculate the scattering factors on a grid
    call precalculate_scattering_factors()
    
    if (on_the_fly) then
        call cuda_setup_many_phasegrate()               !setup the atom co-ordinate for on the fly potentials (the arguments are simply so that
        if(ionization) call make_fz_EELS_EDX()          !setup the scattering factor grid for ionization
        
    else
        idum = seed_rng()
        call make_qep_grates(idum)
        call make_local_inelastic_potentials()          !setup the REAL space inelastic potential (ionization and adf) for QUEP ADF is disabled
    
    endif
    
 	call setup_propagators
    
    do i = 1, ndet
        call make_detector_mask(inner(i),outer(i),masks(:,:,i))
    enddo

    if (on_the_fly.or.quick_shift.or.phase_ramp_shift) then
        allocate(trans_d(nopiy,nopix))
    endif
    
	t1 = secnds(0.0)
    
    
    
    write(*,*) '|--------------------------------|'
	write(*,*) '|      Calculation running       |'
	write(*,*) '|--------------------------------|'
    write(*,*)
    
	! Plan the fourier transforms
    if(fp_kind.eq.8)then
	    call cufftPlan(plan,nopix,nopiy,CUFFT_Z2Z)
	else
        call cufftPlan(plan,nopix,nopiy,CUFFT_C2C)
    endif
    
    ! Copy host arrays to the device
    if (ndet.gt.0) then
        allocate(masks_d(nopiy,nopix,ndet))
        masks_d = masks
    endif
    
    if (ionization) allocate(ion_image_d(nopiy,nopix))
    
    if(EELS) then
        allocate(eels_correction_detector_d(nopiy,nopix))
        eels_correction_detector_d=eels_correction_detector
    endif
    
    allocate(prop_d(nopiy,nopix,n_slices))
    prop_d = prop
    
    if (on_the_fly) then
        allocate(bwl_mat_d(nopiy,nopix))
        allocate(inverse_sinc_d(nopiy,nopix))
        allocate(fz_d(nopiy,nopix,nt))
        fz_d = fz 
        if(ionization) then
            allocate(inelastic_potential_d(nopiy,nopix))
            allocate(fz_mu_d(nopiy,nopix))
            allocate(fz_dwf_d(nopiy,nopix))
            fz_mu_d = fz_mu
            fz_dwf_d = fz_dwf(:,:,kval)
        endif
        inverse_sinc_d = inverse_sinc
        bwl_mat_d = bwl_mat
    else
        allocate(transf_d(nopiy,nopix,n_qep_grates,n_slices))
  	    transf_d = qep_grates
        if (ionization) then
            allocate(ion_potential_d(nopiy,nopix,n_slices))
            ion_potential_d = ionization_potential
        endif
    
        if (phase_ramp_shift) then    
            allocate(shift_array_d(nopiy,nopix))
            allocate(shift_arrayx_d(nopix,ifactorx))
            allocate(shift_arrayy_d(nopiy,ifactory))
            shift_arrayx_d = shift_arrayx
            shift_arrayy_d = shift_arrayy
        endif
    endif
    
    if (output_probe_intensity) then
        allocate(probe_intensity(nopiy,nopix,size(output_thickness_list)))
    endif
    
    stem_image = 0.0_fp_kind
    stem_elastic_image = 0.0_fp_kind
    if (ionization) stem_ion_image = 0.0_fp_kind
    if (eels) eels_correction_image = 0.0_fp_kind
    
    intensity = 1.0d0
    
	do i_df = 1, n_df
	do ny = 1, nysample
    do nx = 1, nxsample
        write(6,903,ADVANCE='NO') achar(13), i_df, n_df, ny, nysample, nx, nxsample, intensity
903     format(a1,' df:',i3,'/',i3,' y:',i3,'/',i3,' x:',i3,'/',i3,'  Intensity:', f6.3, ' (to monitor BWL)')	
    
        if (scan_quarter) then
            if (probe_positions(1,ny,nx).gt.0.501_fp_kind .or. probe_positions(2,ny,nx).gt.0.501_fp_kind) cycle
        endif
    

        cbed_d=0.0_fp_kind
        if (ionization) ion_image_d=0.0_fp_kind
        psi_elastic_d=0.0_fp_kind
        call make_stem_wfn(psi_initial,defoci(i_df),probe_positions(:,ny,nx))
        
        call tilt_wave_function(psi_initial)
        
        psi_initial_d = psi_initial
        
        if (output_probe_intensity) then
            probe_intensity = 0.0_fp_kind
        endif
        
        do i_qep_pass = 1, n_qep_passes 
        
            ! Reset wavefunction
            psi_d = psi_initial_d
            
            do i = 1, n_cells
	            do j = 1, n_slices
                    ! Accumulate ionization cross section
                    if(ionization) then
                        call cuda_mod<<<blocks,threads>>>(psi_d,temp_d,1.0_fp_kind,nopiy,nopix)
                        if(on_the_fly) then
                            call cuda_make_ion_potential(inelastic_potential_d,tau_slice(:,:,:,j),nat_slice(:,j),plan,fz_mu_d,fz_dwf_d,inverse_sinc_d,Volume_array(j),kval)
                            call cuda_multiplication<<<blocks,threads>>>(temp_d,inelastic_potential_d, temp_d,prop_distance(j),nopiy,nopix)
                        else
                            call cuda_multiplication<<<blocks,threads>>>(temp_d,ion_potential_d(:,:,j), temp_d,prop_distance(j),nopiy,nopix)
                        endif
                        call cuda_addition<<<blocks,threads>>>(temp_d,ion_image_d,ion_image_d,1.0_fp_kind,nopiy,nopix)
                    endif
                    
                    ! Phase grate
				    nran = floor(n_qep_grates*ran1(idum)) + 1
                    if(on_the_fly) then
                        call cuda_fph_make_potential(trans_d,ccd_slice_array(j),tau_slice,nat_slice(:,j),j,prop_distance(j),idum,plan,fz_d,inverse_sinc_d,bwl_mat_d)
                        call cuda_multiplication<<<blocks,threads>>>(psi_d,trans_d, psi_out_d,1.0_fp_kind,nopiy,nopix)
                    elseif(quick_shift) then
                        shiftx = floor(ifactorx*ran1(idum)) * nopix_ucell
                        shifty = floor(ifactory*ran1(idum)) * nopiy_ucell
                        call cuda_cshift<<<blocks,threads>>>(transf_d(:,:,nran,j),trans_d,nopiy,nopix,shifty,shiftx)
                        call cuda_multiplication<<<blocks,threads>>>(psi_d,trans_d, psi_out_d,1.0_fp_kind,nopiy,nopix)
                    elseif(phase_ramp_shift) then                       !randomly shift phase grate
                        shiftx = floor(ifactorx*ran1(idum)) + 1
                        shifty = floor(ifactory*ran1(idum)) + 1
                        call cuda_make_shift_array<<<blocks,threads>>>(shift_array_d,shift_arrayy_d(:,shifty),shift_arrayx_d(:,shiftx),nopiy,nopix)     !make the qspace shift array
                        call cuda_multiplication<<<blocks,threads>>>(transf_d(:,:,nran,j),shift_array_d, trans_d,1.0_fp_kind,nopiy,nopix) !multiply by the qspace shift array
                        call cufftExec(plan,trans_d,trans_d,CUFFT_INVERSE)                                                                    !inverse fourier transform
                        call cuda_multiplication<<<blocks,threads>>>(psi_d,trans_d, psi_out_d,sqrt(normalisation),nopiy,nopix)              !do the phase grate multiplication
                    else
                        call cuda_multiplication<<<blocks,threads>>>(psi_d,transf_d(:,:,nran,j), psi_out_d,1.0_fp_kind,nopiy,nopix)
                    endif
                    
                    ! Propagate
			    	call cufftExec(plan,psi_out_d,psi_d,CUFFT_FORWARD)
                    call cuda_multiplication<<<blocks,threads>>>(psi_d,prop_d(:,:,j), psi_out_d,normalisation,nopiy,nopix)
                    call cufftExec(plan,psi_out_d,psi_d,CUFFT_INVERSE)
                    
		        enddo ! End loop over slices
                
                if (output_probe_intensity) then
                    if (output_cell_list(i)) then
                        psi = psi_d
                        probe_intensity(:,:,cell_map(i)) = probe_intensity(:,:,cell_map(i)) + abs(psi)**2
                    endif
                endif
                
            enddo ! End loop over cells
            
            ! Accumulate elastic wave function
	        call cuda_addition<<<blocks,threads>>>(psi_elastic_d,psi_d,psi_elastic_d,1.0_fp_kind,nopiy,nopix)
            
            ! Accumulate diffaction pattern
            call cufftExec(plan,psi_d,psi_out_d,CUFFT_FORWARD)
            call cuda_mod<<<blocks,threads>>>(psi_out_d,temp_d,normalisation,nopiy,nopix)
            call cuda_addition<<<blocks,threads>>>(cbed_d,temp_d,cbed_d,1.0_fp_kind,nopiy,nopix)
	    enddo ! End loop over QEP passes
        
        intensity = get_sum(psi_d)
        
        ! Integrate the diffraction pattern
	    do idet = 1, ndet
            stem_image(ny,nx,i_df,idet) = cuda_stem_detector(cbed_d,masks_d(:,:,idet))
	    enddo
        
        if (ionization) stem_ion_image(ny,nx,i_df) = get_sum(ion_image_d)
        
        if (eels) eels_correction_image(ny,nx,i_df) = cuda_stem_detector(cbed_d,eels_correction_detector_d)
        
        ! Integrate the elastic diffraction pattern
        call cufftExec(plan,psi_elastic_d,psi_out_d,CUFFT_FORWARD)
        call cuda_mod<<<blocks,threads>>>(psi_out_d,temp_d,normalisation,nopiy,nopix) 
        do idet = 1, ndet
            stem_elastic_image(ny,nx,i_df,idet)=cuda_stem_detector(temp_d,masks_d(:,:,idet))
	    enddo
        
        
        if (output_probe_intensity) then
            probe_intensity = probe_intensity / n_qep_passes
            
            filename = trim(adjustl(output_prefix)) // '_ProbeIntensity'
            if (n_df.gt.1) filename = trim(filename) // '_df' // to_string(i_df)
            if (nysample.gt.1) filename = trim(filename) // '_ny' // to_string(ny)
            if (nxsample.gt.1) filename = trim(filename) // '_nx' // to_string(nx)
            filename = trim(filename) // '.bin'
            
            open(4985, file=filename, form='binary', convert='big_endian')
            write(4985) probe_intensity
            close(4985)
            
        endif

	enddo ! End loop over x probe positions
	enddo ! End loop over y probe positions
	enddo ! End loop over defocus series
	
    ! QEP normalisation
    
    if(ionization) stem_ion_image = stem_ion_image/float(n_qep_passes)
    
    if(EELS) eels_correction_image = eels_correction_image/float(n_qep_passes)

    if (ndet.gt.0) then
	    stem_image = stem_image/float(n_qep_passes)
        stem_elastic_image = stem_elastic_image/float(n_qep_passes*n_qep_passes)
    endif

    delta = secnds(t1)
    
    write(*,*) 
    write(*,*) 
    
    write(*,*) 'Calculation is finished.'
    write(*,*) 
    write(*,*) 'Time elapsed ', delta, ' seconds.'
    write(*,*)    
    
    open(unit=9834, file=trim(adjustl(output_prefix))//'_timing.txt', access='append')
    write(9834, '(a, g, a, /)') 'The multislice calculation took ', delta, 'seconds.'
    close(9834)
    
    if (fp_kind.eq.8) then
        write(*,*) 'The following files were outputted (as 64-bit big-endian floating point):'
	else
        write(*,*) 'The following files were outputted (as 32-bit big-endian floating point):'
	endif
    
    write(*,*)
    
    do idet = 1, ndet
        fnam_temp = trim(adjustl(output_prefix)) // '_DiffPlaneTotal_Detector'
        call add_zero_padded_int(fnam_temp, fnam_det, idet, 2)
        call output_stem_image(stem_image(:,:,:,idet),fnam_det)
        
        fnam_temp = trim(adjustl(output_prefix)) // '_DiffPlaneElastic_Detector'
        call add_zero_padded_int(fnam_temp, fnam_det, idet, 2)
        call output_stem_image(stem_elastic_image(:,:,:,idet),fnam_det)
        
        stem_inelastic_image = stem_image(:,:,:,idet) - stem_elastic_image(:,:,:,idet)
        fnam_temp = trim(adjustl(output_prefix)) // '_DiffPlaneTDS_Detector'
        call add_zero_padded_int(fnam_temp, fnam_det, idet, 2)
        call output_stem_image(stem_inelastic_image,fnam_det)
        

    enddo

    !ionization
    if (ionization) then
        if (EELS) then
            fnam = trim(adjustl(output_prefix)) // '_EELS'
            call output_stem_image(stem_ion_image,fnam)
            
            fnam = trim(adjustl(output_prefix)) // '_EELS_CorrectionMap'
            call output_stem_image(eels_correction_image,fnam)
            
            stem_ion_image = stem_ion_image*eels_correction_image
            fnam = trim(adjustl(output_prefix)) // '_EELS_Corrected'
            call output_stem_image(stem_ion_image,fnam)
            
        else 
            fnam = trim(adjustl(output_prefix)) // '_EDX'
            call output_stem_image(stem_ion_image,fnam)
        
        endif
    endif

end subroutine qep_stem
