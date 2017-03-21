module m_multislice
    
    use m_precision, only: fp_kind
    
    implicit none

    logical :: save_grates = .false.
    logical :: load_grates = .false.
    character(1024) :: grates_filename
    
    logical :: output_probe_intensity = .false.
    logical,allocatable :: output_cell_list(:)
    real(fp_kind),allocatable :: output_thickness_list(:)
    integer,allocatable :: cell_map(:)
    
    contains
    
    
    
    subroutine prompt_output_probe_intensity
    
        use m_user_input, only: get_input
        use global_variables, only: thickness, nopiy, nopix
        use output, only: output_prefix
        use m_string, only: to_string
        
        implicit none
        
        integer :: i_output
        real(fp_kind) :: thickness_interval
        
        write(*,*) '|-----------------------------------|'
        write(*,*) '|      Output probe intensity       |'
        write(*,*) '|-----------------------------------|'
        write(*,*)
        write(*,*) 'The probe intensity as a function of thickness'
        write(*,*) 'can be outputted to file at each probe position.'
        write(*,*) 'The user is advised that the outputted dataset'
        write(*,*) 'may be very large.'
        write(*,*) 
10      write(*,*) '<0> Proceed'
        write(*,*) '<1> Output probe intensity'
        call get_input('<0> Proceed <1> Output probe intensity', i_output)
        write(*,*)
  
        select case (i_output)
            case (0)
                return
                
            case (1)
                output_probe_intensity = .true.
                
20              write(*,*) 'At what thickness interval (in Angstroms)'
                write(*,*) 'should intensities be outputted?'
                call get_input('At what thickness interval should intensities be outputted?', thickness_interval)
                write(*,*)
                
                if (thickness_interval.le.0.0_fp_kind .or. thickness_interval.gt.thickness) then
                    write(*,*) 'ERROR: invalid thickness.'
                    goto 20
                endif
                
                call generate_cell_list(thickness_interval)
                
                call write_thicknesss_to_file
                
                write(*,*) 'The probe intensities will be written to the files'
                write(*,*)
                write(*,*) '    ' // trim(adjustl(output_prefix)) // '_ProbeIntensity*.bin'
                write(*,*)
                if (fp_kind.eq.4) then
                    write(*,*) 'as 32-bit big-endian floats.'
                elseif (fp_kind.eq.8) then
                    write(*,*) 'as 64-bit big-endian floats.'
                endif
                write(*,*) 'Each file contains a sequence of ' // to_string(nopiy) // 'x' // to_string(nopix) // ' arrays.'
                write(*,*)
                
            case default
                goto 10
                
        end select
        
    end subroutine
    
    
    
    subroutine generate_cell_list(thickness_interval)
    
        use global_variables, only: thickness, n_cells, a0
        
        implicit none
        
        real(fp_kind) :: thickness_interval
        
        integer :: i_cell, count
        real(fp_kind) :: t
        
        allocate(output_cell_list(n_cells))
        output_cell_list = .false.
        
        allocate(cell_map(n_cells))
        cell_map = 0
        
        t = 0.0_fp_kind
        count = 0
        
        do while (t.le.thickness)
            i_cell = nint(t/a0(3))
            
            if (i_cell.ge.1) then
                output_cell_list(i_cell) = .true.
                count = count + 1
            endif
            
            t = t + thickness_interval
        enddo
        
        allocate(output_thickness_list(count))
        
        count = 0
        
        do i_cell = 1, n_cells
            if (output_cell_list(i_cell)) then
                count = count + 1
                output_thickness_list(count) = i_cell * a0(3)
                cell_map(i_cell) = count
            endif            
        enddo
        
    end subroutine
    
        
    
    subroutine write_thicknesss_to_file
    
        use output, only: output_prefix
        
        implicit none
        
        integer :: i
        character(1024) :: filename
        
        filename = trim(adjustl(output_prefix))//'_probe_intensity_thicknesss.txt'
        
        write(*,*) 'The thicknesses at which the probe intensity'
        write(*,*) 'is being outputted have been written to'
        write(*,*)
        write(*,*) '    ' // trim(filename)
        write(*,*)
        
        open(unit=8734, file=filename)
        
        do i = 1, size(output_thickness_list)
            write(8734, *) output_thickness_list(i)
        enddo
        
        close(8734)        
        
    end subroutine
    
    
    
    subroutine prompt_save_load_grates
    
        use m_user_input, only: get_input
        use output, only: output_prefix
        
        implicit none
        
        integer :: i_save_load, i_retry
        logical :: exists
        
        write(*,*) '|---------------------------------------------|'
        write(*,*) '|      Save/load transmission functions       |'
        write(*,*) '|---------------------------------------------|'
        write(*,*)
        write(*,*) 'Warning: the files outputted when saving may be very large.'
        write(*,*)
1       write(*,*) '<0> Proceed without saving or loading'
        write(*,*) '<1> Save transmission functions'
        write(*,*) '<2> Load transmission functions'
        call get_input('<0> continue <1> save <2> load', i_save_load)
        write(*,*)
        
        select case (i_save_load)
            case (0)
                return
                
            case (1)
                save_grates = .true.
                grates_filename = trim(adjustl(output_prefix)) // '_transmission_functions.bin'
                
                write(*,*) 'The transmission functions will be saved to the file'
                write(*,*)
                write(*,*) '    ' // trim(grates_filename)
                write(*,*)
                write(*,*) 'They can be loaded for later calculations provided'
                write(*,*) 'the following parameters are identical:'
                write(*,*) '  - The xtl file'
                write(*,*) '  - The slicing of the unit cell'
                write(*,*) '  - The choice of thermal scattering model (QEP vs. absorptive)'
                write(*,*) '  - The tiling of the unit cell'
                write(*,*) '  - The number of pixels'
                write(*,*) '  - (For absorptive model: whether absorption is included)'
                write(*,*) '  - (For QEP model: the number of distinct transmission functions)'
                write(*,*) '  - (For QEP model: phase ramp shift choice)'
                write(*,*)
                
            case (2)
                write(*,*) 'It is up to the user to ensure that the parameters used'
                write(*,*) 'to create the loaded transmission functions are consistent'
                write(*,*) 'with those of the current calculation:'
                write(*,*) '  - The xtl file'
                write(*,*) '  - The slicing of the unit cell'
                write(*,*) '  - The choice of thermal scattering model (QEP vs. absorptive)'
                write(*,*) '  - The tiling of the unit cell'
                write(*,*) '  - The number of pixels'
                write(*,*) '  - (For absorptive model: whether absorption is included)'
                write(*,*) '  - (For QEP model: the number of distinct transmission functions)'
                write(*,*) '  - (For QEP model: phase ramp shift choice)'
                write(*,*)
2               write(*,*) 'Enter filename of transmission functions:'
                call get_input('filename of transmission functions', grates_filename)
                write(*,*)
                
                inquire(file=grates_filename, exist=exists)
                
                if (.not.exists) then
                    write(*,*) 'ERROR: cannot find this file.'
3                   write(*,*) '<1> Enter again'
                    write(*,*) '<2> Proceed without loading'
                    call get_input('<1> Enter again <2> Proceed without loading', i_retry)
                    write(*,*)
                    
                    select case (i_retry)
                        case (1)
                            goto 2
                        case (2)
                            return
                        case default
                            goto 3
                    end select
                    
                endif
                
                load_grates = .true.
            case default
                goto 1
                
        end select
        
    end subroutine
    
    
    
    subroutine make_qep_grates(idum)
    
        use m_precision, only: fp_kind
	    use cufft_wrapper, only: fft2, ifft2
        use global_variables, only: nopiy, nopix, nt, relm, tp, ak, ak1, atf, high_accuracy, ci, pi, bwl_mat,Kz
        use m_qep, only: displace, qep_grates,phase_ramp_shift, n_qep_grates
        use m_slicing, only: n_slices, nat_slice, a0_slice, tau_slice, maxnat_slice, ss_slice
        use m_string, only: to_string
        use output, only: output_prefix
        use m_potential, only: make_qep_potential, make_site_factor_generic, make_site_factor_cuda, make_site_factor_hybrid
        
        implicit none
    
        integer(4) :: idum
    
        integer(4) :: i, j, m, n
        complex(fp_kind) :: projected_potential(nopiy,nopix)
        complex(fp_kind) :: temp(nopiy,nopix)
        real(fp_kind) :: ccd_slice
        complex(fp_kind) :: scattering_pot(nopiy,nopix,nt)
        real(fp_kind),allocatable :: mod_tau(:,:,:,:,:)
    
        real(fp_kind) :: t1, delta
    
        procedure(make_site_factor_generic),pointer :: make_site_factor
        
	    if(allocated(mod_tau)) deallocate(mod_tau)
	    if(allocated(qep_grates)) deallocate(qep_grates)
	
        allocate(qep_grates(nopiy,nopix,n_qep_grates,n_slices))
        allocate(mod_tau(3,nt,maxnat_slice,n_slices,n_qep_grates))
 	
        if (load_grates) then
            write(*,*) 'Loading transmission functions from file...'
            write(*,*)
            open(unit=3984, file=grates_filename, form='binary', convert='big_endian')
            read(3984) qep_grates
            close(3984)
            
            ! Make sure the random number sequence is as if grates were calculated
            do j = 1, n_slices
	        do i = 1, n_qep_grates
 	        do m = 1, nt
	        do n = 1, nat_slice(m,j)
			    call displace(tau_slice(1:3,m,n,j),mod_tau(1:3,m,n,j,i),sqrt(atf(3,m)),a0_slice,idum)
            enddo
            enddo
            enddo
            enddo
            
            return
        endif
        
        t1 = secnds(0.0_fp_kind)
        
        if (high_accuracy) then
            make_site_factor => make_site_factor_cuda
        else
            make_site_factor => make_site_factor_hybrid
        endif

        do j = 1, n_slices
	        write(*,'(1x, a, a, a, a, a)') 'Calculating transmission functions for slice ', to_string(j), '/', to_string(n_slices), '...'
        
    198	    write(6,199) to_string(sum(nat_slice(:,j)))
    199     format(1x, 'Number of atoms in this slice: ', a) 

		    ccd_slice = relm / (tp * ak * ss_slice(7,j))

	        do i = 1, n_qep_grates
    200	        format(a1, 1x, i3, '/', i3)
	            write(6,200, advance='no') achar(13), i, n_qep_grates
          
                ! Randomly displace the atoms
 	            do m = 1, nt
	                do n = 1, nat_slice(m,j)
			            call displace(tau_slice(1:3,m,n,j),mod_tau(1:3,m,n,j,i),sqrt(atf(3,m)),a0_slice,idum)
	                enddo
                enddo
	      
                call make_qep_potential(projected_potential, mod_tau(:,:,:,j,i), nat_slice(:,j), ccd_slice, make_site_factor)

                projected_potential = real(projected_potential)
                
                qep_grates(:,:,i,j) = exp(ci*pi*a0_slice(3,j)/Kz*projected_potential)
            
                ! Bandwith limit the phase grate
                call fft2(nopiy, nopix, qep_grates(:,:,i,j), nopiy, temp, nopiy)
	            qep_grates(:,:,i,j) = temp * bwl_mat
            
                if (.not.phase_ramp_shift) then
                    ! If not using phase_ramp_shifts, get the phase grates in real space
                    call ifft2(nopiy, nopix, qep_grates(:,:,i,j), nopiy, qep_grates(:,:,i,j), nopiy)
                endif
	        enddo ! End loop over grates
        
            write(*,*)
            write(*,*)
        
	    enddo ! End loop over slices
	
	    delta = secnds(t1)
        
        write(*,*) 'The calculation of transmission functions for the QEP model took ', delta, 'seconds.'
        write(*,*)
        
        open(unit=9834, file=trim(adjustl(output_prefix))//'_timing.txt', access='append')
        write(9834, '(a, g, a, /)') 'The calculation of transmission functions for the QEP model took ', delta, 'seconds.'
        close(9834)
        
        if (save_grates) then
            write(*,*) 'Saving transmission functions to file...'
            write(*,*)
            open(unit=3984, file=grates_filename, form='binary', convert='big_endian')
            write(3984) qep_grates
            close(3984)
        endif        
    
    end subroutine make_qep_grates


    
    subroutine make_absorptive_grates
    
        use m_precision, only: fp_kind
	    use cufft_wrapper, only: fft2, ifft2
        use global_variables, only: nopiy, nopix, nt, relm, tp, ak, Kz, atf, high_accuracy, ci, pi, bwl_mat
        use m_absorption, only: transf_absorptive
        use m_slicing, only: n_slices, nat_slice, a0_slice, tau_slice, maxnat_slice, ss_slice
        use m_string, only: to_string
        use output, only: output_prefix
        use m_potential, only: make_absorptive_potential, make_site_factor_generic, make_site_factor_cuda, make_site_factor_hybrid
        
        implicit none
    
        integer(4) :: j, m, n
        complex(fp_kind) :: projected_potential(nopiy,nopix)
        complex(fp_kind) :: temp(nopiy,nopix)
        real(fp_kind) :: ccd_slice
        complex(fp_kind) :: scattering_pot(nopiy,nopix,nt)
    
        real(fp_kind) :: t1, delta
    
        procedure(make_site_factor_generic),pointer :: make_site_factor
        	
        if(allocated(transf_absorptive)) deallocate(transf_absorptive)   
        allocate(transf_absorptive(nopiy,nopix,n_slices))
        
        if (load_grates) then
            write(*,*) 'Loading transmission functions from file...'
            write(*,*)
            open(unit=3984, file=grates_filename, form='binary', convert='big_endian')
            read(3984) transf_absorptive
            close(3984)
            return
        endif
        
        t1 = secnds(0.0_fp_kind)
        
        if (high_accuracy) then
            make_site_factor => make_site_factor_cuda
        else
            make_site_factor => make_site_factor_hybrid
        endif

        do j = 1, n_slices
	        write(*,'(1x, a, a, a, a, a)') 'Calculating transmission functions for slice ', to_string(j), '/', to_string(n_slices), '...'
        
    198	    write(6,199) to_string(sum(nat_slice(:,j)))
    199     format(1x, 'Number of atoms in this slice: ', a, /) 

		    ccd_slice = relm / (tp * ak * ss_slice(7,j))

            call make_absorptive_potential(projected_potential, tau_slice(:,:,:,j), nat_slice(:,j), ccd_slice, ss_slice(7,j), make_site_factor)

            transf_absorptive(:,:,j) = exp(ci*pi*a0_slice(3,j)/Kz*projected_potential)
            
            ! Bandwith limit the phase grate
            call fft2(nopiy, nopix, transf_absorptive(:,:,j), nopiy, temp, nopiy)
            temp = temp * bwl_mat
            call ifft2(nopiy, nopix, temp, nopiy, transf_absorptive(:,:,j), nopiy)
                
	    enddo ! End loop over slices
	
	    delta = secnds(t1)
        
        write(*,*) 'The calculation of transmission functions for the absorptive model took ', delta, 'seconds.'
        write(*,*)
        
        open(unit=9834, file=trim(adjustl(output_prefix))//'_timing.txt', access='append')
        write(9834, '(a, g, a, /)') 'The calculation of transmission functions for the absorptive model took ', delta, 'seconds.'
        close(9834)
    
        if (save_grates) then
            write(*,*) 'Saving transmission functions to file...'
            write(*,*)
            open(unit=3984, file=grates_filename, form='binary', convert='big_endian')
            write(3984) transf_absorptive
            close(3984)
        endif        
        
    end subroutine make_absorptive_grates



    subroutine setup_propagators()
      
        use global_variables, only: nopiy, nopix, bwl_mat, prop, ak1, ss, ig1, ig2, ifactory, ifactorx,Kz,claue
        use m_precision, only: fp_kind
        use m_slicing, only: n_slices, prop_distance
        use output
      
        implicit none
      
        integer(4) :: i
        
        if(allocated(prop)) deallocate(prop)
        allocate(prop(nopiy,nopix,n_slices))
        
        do i = 1, n_slices
	        call make_propagator(nopiy,nopix,prop(:,:,i),prop_distance(i),Kz,ss,ig1,ig2,claue,ifactorx,ifactory)
	        prop(:,:,i) = prop(:,:,i) * bwl_mat
            !call binary_out_unwrap(nopiy,nopix, atan2(real(prop(:,:,i)),imag(prop(:,:,i)))* bwl_mat,'Propagator')
        enddo

    end subroutine

      
      
    subroutine make_propagator(nopiy,nopix,prop,dz,ak1,ss,ig1,ig2,claue,ifactorx,ifactory)

        use m_precision, only: fp_kind
        use m_crystallography, only: trimr
            
        implicit none

        integer(4) :: nopiy,nopix
        complex(fp_kind) :: prop(nopiy,nopix)        
        real(fp_kind) :: ak1, ss(7), claue(3), dz
        integer(4) :: ifactorx, ifactory, ig1(3), ig2(3)
        
        integer(4) :: m1, m2, kx(3), ky(3), shifty, shiftx
        real(fp_kind) :: kr(3), q_sq

        real(fp_kind),parameter :: pi = atan(1.0d0)*4.0d0
        integer(4) :: ny, nx
        
        shifty = (nopiy-1)/2-1
        shiftx = (nopix-1)/2-1

        do ny = 1, nopiy

	        m2 = mod( ny+shifty, nopiy) - shifty -1
		
	        ky = m2 * ig2	
	  
	        do nx = 1, nopix

	            m1 = mod( nx+shiftx, nopix) - shiftx -1
	 
	            kx = m1 * ig1
                
	            kr = kx/float(ifactorx) + ky/float(ifactory) - claue
	            q_sq = trimr(kr,ss)**2

	            prop(ny,nx) = exp(cmplx(0.0d0, -pi*dz*q_sq/ak1, fp_kind ))

	        enddo
        enddo

	end subroutine
    
    
	    
end module

