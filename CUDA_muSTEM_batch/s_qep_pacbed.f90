function qep_pacbed_GPU_memory() result(required_memory)
    
    use m_precision, only: fp_kind
    use global_variables, only: nopiy, nopix, ifactory, ifactorx, on_the_fly
    use m_lens, only: imaging
    use m_slicing, only: n_slices
    use m_qep, only: n_qep_grates, quick_shift, phase_ramp_shift
    
    implicit none
    
    real(fp_kind) :: required_memory
    
    real(fp_kind) :: array_size
    integer :: array_count
    
    array_size = 4.0_fp_kind * nopiy * nopix
    
    array_count = 2*(3 + n_slices + n_slices*n_qep_grates) + 2
        
    if (on_the_fly.or.quick_shift.or.phase_ramp_shift) array_count = array_count + 2
    if (phase_ramp_shift) array_count = array_count + 2
    
    required_memory = array_count * array_size
    
    if (phase_ramp_shift) required_memory = required_memory + 8.0_fp_kind*(nopiy*ifactory + nopix*ifactorx)
    
end function

    

subroutine qep_pacbed
    
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
	use m_string
    use m_probe_scan, only: nysample, nxsample, probe_positions
    use m_tilt, only: tilt_wave_function
    use m_multislice, only: make_qep_grates, setup_propagators
    use m_potential, only: precalculate_scattering_factors
    
    implicit none
    
    !dummy variables
    integer(4) ::  i, j, i_qep_pass
    integer(4) ::  shifty,shiftx
    integer(4) ::  ny, nx
    
    !random variables
    integer(4) :: idum
    real(fp_kind) :: ran1
       
    !probe variables
    complex(fp_kind) :: psi(nopiy,nopix)
    complex(fp_kind) :: psi_initial(nopiy,nopix)
         
    real(fp_kind),dimension(nopiy, nopix) :: pacbed_pattern,fourDSTEM_pattern
    
    !diagnostic variables
    real(fp_kind) :: intensity,t1, delta
    
    !output variables
    character(120) :: filename

    !device variables
	integer :: plan
    complex(fp_kind),device,dimension(nopiy,nopix) :: psi_initial_d, psi_d, psi_out_d
    real(fp_kind),device,dimension(nopiy,nopix) :: temp_cbed_d, pacbed_pattern_d,fourDSTEM_pattern_d
	complex(fp_kind),device,allocatable :: prop_d(:,:,:), transf_d(:,:,:,:)
	complex(fp_kind),device,allocatable :: trans_d(:,:), shift_arrayx_d(:,:), shift_arrayy_d(:,:), shift_array_d(:,:)
    
    !device variables for on the fly potentials
    complex(fp_kind),device,allocatable,dimension(:,:) :: bwl_mat_d, inverse_sinc_d
    complex(fp_kind),device,allocatable,dimension(:,:,:) :: fz_d

	logical:: fourDSTEM

    real(fp_kind) :: qep_pacbed_GPU_memory
    
    
    
    call GPU_memory_message(qep_pacbed_GPU_memory(), on_the_fly)
    
    
    
    write(*,*) '|----------------------------------|'
	write(*,*) '|      Pre-calculation setup       |'
	write(*,*) '|----------------------------------|'
    write(*,*)
    
	write(*,*) 'Output diffraction patterns for each scan position?'
	call get_input('<1> Diffraction pattern for each probe position',idum)
	fourDSTEM = idum == 1

    call precalculate_scattering_factors
    
    if (on_the_fly) then
        call cuda_setup_many_phasegrate
        
    else
        idum = seed_rng()
        call make_qep_grates(idum)
    
    endif
    
 	call setup_propagators
    
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
    allocate(prop_d(nopiy,nopix,n_slices))
    prop_d = prop
    
    if(on_the_fly) then
        allocate(bwl_mat_d(nopiy,nopix))
        allocate(inverse_sinc_d(nopiy,nopix))
        allocate(fz_d(nopiy,nopix,nt))
        fz_d = fz 
        inverse_sinc_d = inverse_sinc
        bwl_mat_d = bwl_mat
    else
        allocate(transf_d(nopiy,nopix,n_qep_grates,n_slices))
  	    transf_d = qep_grates
        
        if(phase_ramp_shift) then 
            allocate(shift_array_d(nopiy,nopix))
            allocate(shift_arrayx_d(nopix,ifactorx))
            allocate(shift_arrayy_d(nopiy,ifactory))
            shift_arrayx_d = shift_arrayx
            shift_arrayy_d = shift_arrayy
        endif
    endif
	
    intensity = 1.0_fp_kind
    pacbed_pattern_d = 0.0_fp_kind
    
	do ny = 1, nysample
    do nx = 1, nxsample
        write(6,903,ADVANCE='NO') achar(13), ny, nysample, nx, nxsample, intensity
903     format(a1, ' y:', i3, '/', i3, ' x:', i3, '/', i3, '  Intensity:', f6.3, ' (to monitor BWL)')	
    
        call make_stem_wfn(psi_initial, probe_df, probe_positions(:,ny,nx))
        
        call tilt_wave_function(psi_initial)
        
        psi_initial_d = psi_initial
        
		!Reset 4DSTEM pattern to zero
		fourDSTEM_pattern_d = 0_fp_kind

        do i_qep_pass = 1, n_qep_passes 
        
            ! Reset wavefunction
            psi_d = psi_initial_d
            
			

            do i = 1, n_cells
	            do j = 1, n_slices
                    
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
            enddo ! End loop over cells
            
            ! Accumulate PACBED pattern (and 4DSTEM diffraction pattern)
            call cufftExec(plan, psi_d, psi_out_d, CUFFT_FORWARD)
            call cuda_mod<<<blocks,threads>>>(psi_out_d, temp_cbed_d, normalisation, nopiy, nopix)
            if (fourDSTEM) call cuda_addition<<<blocks,threads>>>(fourDSTEM_pattern_d, temp_cbed_d, fourDSTEM_pattern_d, 1.0_fp_kind, nopiy, nopix)
			call cuda_addition<<<blocks,threads>>>(pacbed_pattern_d, temp_cbed_d, pacbed_pattern_d, 1.0_fp_kind, nopiy, nopix)

			
            
	    enddo ! End loop over QEP passes
        
        intensity = get_sum(psi_d)

        !Output 4D STEM diffraction pattern
		if(fourDSTEM) then
				filename = trim(adjustl(output_prefix)) //'_pp_'//to_string(nx)//'_'//to_string(ny)//'_Diffraction_pattern'
				fourDSTEM_pattern = fourDSTEM_pattern_d
				call binary_out_unwrap(nopiy, nopix, fourDSTEM_pattern/n_qep_passes, filename,write_to_screen=.false.)
		endif
		!end

	enddo ! End loop over x probe positions
	enddo ! End loop over y probe positions
	
    ! Copy PACBED pattern to CPU
    pacbed_pattern = pacbed_pattern_d
    
    ! QEP normalisation
    pacbed_pattern = pacbed_pattern / n_qep_passes
    
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

    filename = trim(adjustl(output_prefix)) // '_PACBED_pattern'
    call binary_out_unwrap(nopiy, nopix, pacbed_pattern, filename)
    
end subroutine qep_pacbed
