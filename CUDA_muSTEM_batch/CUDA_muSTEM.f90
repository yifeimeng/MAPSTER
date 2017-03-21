!--------------------------------------------------------------------------------
!   Program: MU_STEM (GPU VERSION)
!    
!   Description:    Calculate (S)TEM images and diffraction patterns using the
!                   multislice approach.
!                   Includes the contributions of TDS using either the Hall & Hirsch
!                   absorptive model or the Quantum Excitation of Phonons model.
!                   Both plane wave and convergent beam illumination may be used.
!                   STEM EDX/EELS images can be calculated within the local approximation.
!
!                   Outputs are big-endian floating point numbers in either
!                   32-bit or 64-bit precision, depending on which precision
!                   is chosen in mod_precision.f90.
!    
!   Maintainer:     Les Allen
!   Email:          lja@unimelb.edu.au
!   Date:           May 2016
!   Requirements:   PGI Fortran
!
!   version:        4.5  
!                       
!--------------------------------------------------------------------------------

    program MU_STEM
    
        use m_user_input, only: init_input, get_input
        use global_variables, only: high_accuracy, nt, atf, nat, atomf, volts, ss, qep, adf, constants, nopiy, nopix
        use m_lens, only: pw_illum, cb_illum, setup_imaging_lens_parameters, setup_probe_lens_parameters, setup_defocus_series
        use local_ionization, only: setup_inelastic_ionization_types
        use cuda_setup, only: setup_GPU
        use cuda_array_library, only: set_blocks
        use m_slicing, only: setup_slicing_depths, calculate_slices
        use output, only: setup_output_prefix
        use m_probe_scan, only: setup_probe_scan, setup_probe_scan_pacbed, place_probe, probe_initial_position
        use m_tilt, only: prompt_tilt
        use m_qep, only: setup_qep_parameters
        use m_absorption, only: complex_absorption, prompt_include_absorption, setup_local_diffraction_plane_geometry
        use m_potential, only: prompt_high_accuracy
        use m_multislice, only: prompt_save_load_grates, prompt_output_probe_intensity
        
        implicit none
        
        integer :: i_illum, i_tds_model, i_cb_menu, i_cb_calc_type
        
        logical :: nopause = .false.
        character(512)::command_argument      
     
        integer :: i_arg
        
108     write(6,109)
109     format(&
       &1x,'|----------------------------------------------------------------------------|',/,&
       &1x,'|              Melbourne University (scanning) transmission electron         |',/,&
       &1x,'|                            microscopy computing suite                      |',/,&
       &1x,'|      __       __  __    __   ______  ________  ________  __       __       |',/,&
       &1x,'|     |  \     /  \|  \  |  \ /      \|        \|        \|  \     /  \      |',/,&
       &1x,'|     | $$\   /  $$| $$  | $$|  $$$$$$\\$$$$$$$$| $$$$$$$$| $$\   /  $$      |',/,&
       &1x,'|     | $$$\ /  $$$| $$  | $$| $$___\$$  | $$   | $$__    | $$$\ /  $$$      |',/,&
       &1x,'|     | $$$$\  $$$$| $$  | $$ \$$    \   | $$   | $$  \   | $$$$\  $$$$      |',/,&
       &1x,'|     | $$\$$ $$ $$| $$  | $$ _\$$$$$$\  | $$   | $$$$$   | $$\$$ $$ $$      |',/,&
       &1x,'|     | $$ \$$$| $$| $$__/ $$|  \__| $$  | $$   | $$_____ | $$ \$$$| $$      |',/,&
       &1x,'|     | $$  \$ | $$ \$$    $$ \$$    $$  | $$   | $$     \| $$  \$ | $$      |',/,&
       &1x,'|      \$$      \$$  \$$$$$$   \$$$$$$    \$$    \$$$$$$$$ \$$      \$$      |',/,&
       &1x,'|                                                                            |',/,&
       &1x,'|    This software may only be used for non-commercial academic purposes.    |',/,&
       &1x,'|   Copyright 2016 A.J. Dalfonso, S.D. Findlay, B. D. Forbes and L.J. Allen. |',/,&
       &1x,'|      The University of Melbourne, Parkville 3010, Victoria, Australia.     |',/,&
       &1x,'|           No warranties of any kind are given or implied.                  |',/,&
       &1x,'|                       email: a.j@dalfonso.com.au                           |',/,&
       &1x,'|                                                                            |',/,&
       &1x,'|                                                                            |',/,&
       &1x,'|    Version 4.5                                                             |',/,&
       &1x,'|                                                                            |',/,&
       &1x,'|    Note: pass the argument "nopause" (without quotation marks)             |',/,&
       &1x,'|          e.g. muSTEM.exe nopause                                           |',/,&
       &1x,'|          to avoid pauses.                                                  |',/,&
       &1x,'|                                                                            |',/,&
       &1x,'|----------------------------------------------------------------------------|',/)
       
        ! Process command line arguments
        do i_arg = 1, command_argument_count()
            call get_command_argument(i_arg, command_argument)
            select case (trim(adjustl(command_argument)))
            case ('nopause')
                nopause = .true.
            end select
        enddo
       
        if (.not. nopause) then
            write(*,*) ' Press enter to continue.'
            read(*,*)
        endif
       
        
        
        ! Set up user input routines	
        call init_input
       
        ! Set up CPU multithreading
        call setup_threading
        
        ! Set up GPU
        call setup_GPU

        write(*,*) '|------------------------------|'
        write(*,*) '|        Dataset output        |'
        write(*,*) '|------------------------------|'
        write(*,*)
        
        ! Ask for the prefix for dataset output
        call setup_output_prefix
        
        write(*,*) '|--------------------------------|'
        write(*,*) '|        Specimen details        |'
        write(*,*) '|--------------------------------|'
        write(*,*)
        
        ! Read in the xtl file
        call set_xtl_global_params 
        call validate_xtl
        
      
        ! Calculate the mean inner potential
        call set_volts(nt, atf, nat, atomf, volts, ss)       
        
        ! Set electron quantities
        call constants           
        
        ! Ask for slicing scheme
        call setup_slicing_depths
        
        ! Ask for thickness
        call setup_specimen_thickness
       
        ! Set up the potential calculation method
        call prompt_high_accuracy    
        
        
        
        write(*,*) '|----------------------------|'
        write(*,*) '|        Illumination        |'
        write(*,*) '|----------------------------|'
        write(*,*)
110     write(*,*) '<1> Plane wave      (including HRTEM images and diffraction patterns)'
        write(*,*) '<2> Convergent-beam (including STEM images and CBED patterns)'
                                                                                               
        call get_input('Calculation type', i_illum) 
        write(*,*)
       
        ! Validate illumination choice
        if (i_illum.ne.1 .and. i_illum.ne.2 .and. i_illum.ne.99) goto 110
     
        
        ! Set illumination flags
        pw_illum = (i_illum == 1)
        cb_illum = .not. pw_illum
        

        
        ! For convergent-beam, set up the probe forming lens
        if (cb_illum) call setup_probe_lens_parameters

        
        
        write(*,*) '|----------------------------------------|'
        write(*,*) '|        Thermal scattering model        |'
        write(*,*) '|----------------------------------------|'
        write(*,*)
112     write(*,*) '<1> Quantum Excitation of Phonons (QEP) model'
        write(*,*) '     (includes thermal intensity component)'
        write(*,*) '<2> Absorptive model'
        write(*,*) '     (only includes absorption from elastic wave function)'

        call get_input('<1> QEP <2> ABS', i_tds_model)
        write(*,*)
        
        ! Validate choice
        if (i_tds_model.ne.1 .and. i_tds_model.ne.2) goto 112
        
        ! Set TDS model flags
        complex_absorption = (i_tds_model == 2)
        qep = .not. complex_absorption
        
        ! Prompt for including absorptive potential
        if (complex_absorption) call prompt_include_absorption
                
        ! Set the unit cell tiling and grid size
        call set_tiling_grid

        ! Prompt user for a tilt for either kind of illumination
        call prompt_tilt        
        
        ! Setup the CUDA thread hierachy for nopiy, nopix arrays
        call set_blocks(nopiy, nopix)
        
        ! Calculate the slices of the supercell
        call calculate_slices
        
        ! Calculate the bandwidth limiting matrix
        call make_bwl_mat                           
        
        ! Ask for QEP parameters
        if (qep) call setup_qep_parameters
        
        ! Save/load transmission functions
        call prompt_save_load_grates
        
        ! Set up the imaging lens
        if (pw_illum) call setup_imaging_lens_parameters
        
        ! Choose the convergent-beam calculation type
        if (cb_illum) then
            write(*,*) '|--------------------------------|'
            write(*,*) '|        Calculation type        |'
            write(*,*) '|--------------------------------|'
            write(*,*)
            write(*,*) 'Choose a calculation type:'
115         write(*,*) '<1> CBED pattern'
            write(*,*) '<2> PACBED pattern (can also output diffraction patterns for each probe position)'
            write(*,*) '<3> STEM image (BF/ABF/ADF/EELS/EDX)'
            
            call get_input('<1> CBED <2> PACBED <3> STEM', i_cb_calc_type)
            write(*,*)
            
            select case (i_cb_calc_type)
                case (1)
                    ! CBED pattern
                    call place_probe(probe_initial_position)
                    
                case (2)
                    ! PACBED pattern
                    call setup_probe_scan_pacbed
                    
                case (3)
                    ! STEM images
                    call setup_defocus_series
                    call setup_probe_scan
                    
                    if (qep) then
                        adf = .false. 
                        call setup_integration_measurements
                        
                    else
                        adf = .true.
                        call setup_local_diffraction_plane_geometry
                    
                    endif

                    call setup_inelastic_ionization_types
                    
                case default
                    goto 115
                    
            end select
            
            
        endif
        
        
        
        !-----------------------------
        ! Start simulation
        !-----------------------------
            
        if (pw_illum .and. qep) then
            ! Plane wave QEP
            call qep_tem

        elseif (pw_illum .and. complex_absorption) then
            ! Plane wave absorptive
            call absorptive_tem
           
        elseif (cb_illum .and. qep .and. i_cb_calc_type.eq.1) then
            ! Convergent-beam QEP CBED
            call qep_tem

        elseif (cb_illum .and. qep .and. i_cb_calc_type.eq.2) then
            ! Convergent-beam QEP PACBED
            call qep_pacbed

        elseif (cb_illum .and. qep .and. i_cb_calc_type.eq.3) then
            ! Convergent-beam QEP STEM
            call qep_stem
            
        elseif (cb_illum .and. complex_absorption .and. i_cb_calc_type.eq.1) then
            ! Convergent-beam absorptive CBED
            call absorptive_tem

        elseif (cb_illum .and. complex_absorption .and. i_cb_calc_type.eq.2) then
            ! Convergent-beam absorptive PACBED
            call absorptive_pacbed

        elseif (cb_illum .and. complex_absorption .and. i_cb_calc_type.eq.3) then
            ! Convergent-beam absorptive STEM
            call absorptive_stem
            
        endif
         
        write(*,*)
  
        if (.not. nopause) then
            write(*,*) ' Press enter to exit.'
            read(*,*) 
        endif
    end program

   
   
    subroutine setup_threading()
    
        use m_string, only: to_string
    
        implicit none
      
        integer*4 :: num_threads
        integer*4 :: omp_get_max_threads, omp_get_num_procs
    
        num_threads = omp_get_num_procs() 
        
        call omp_set_num_threads(num_threads)
        
        write(*,*) '|----------------------------------|'
        write(*,*) '|        CPU multithreading        |'
        write(*,*) '|----------------------------------|'
        write(*,*)
        write(*,*) 'The number of threads being used on the CPU is: ' // to_string(num_threads)
        write(*,*)
    
    end subroutine         
