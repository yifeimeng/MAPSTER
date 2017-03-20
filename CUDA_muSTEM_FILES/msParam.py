class Param:
	device = 0
	output_name = 'N/A'
	crystal_file_name = 'N/A'
	num_slice = 1
	thickness = 1
	potential_method = 1 # 1 is reciprocal space method (accuracy), 2 is hybrid method (speed)
	beam_type = 1
	cs = 0
	c5 = 0
	convergence_angle = 10
	defocus = 10
	tds_method = 2 # 1 is qep, 2 is absorption model
	absorp_choice = 1
	tile_x = 1
	tile_x = 1
	pixel_x = 256
	pixel_y = 256
	beam_tilt = 0
	specimen_tile = 0
	output_type = 2 # convergent beam only, 1 is CBED, 2 is PACBED, 3 is STEM
    otf_potential = 0 # on-the-fly calculation of the potential, 1 is otf, 0 is precalculated (faster)

	def writeParam(self):
        fo = open('input_param.txt', 'w+')
		fo.write('Device used for calculation\n')
		fo.write(str(device) + '\n')
		fo.write('Output filename\n')
		fo.write(output_name + '\n')
		fo.write('Input crystal file name\n')
		fo.write(crystal_file_name + '\n')
		fo.write('Slice unit cell <1> yes <2> no\n')
		fo.write('1\n')
		fo.write('Number of distinct potentials\n')
		fo.write(str(num_slice) + '\n')
		fo.write('<1> manual or <2> auto slicing\n')
		fo.write('2\n')
		fo.write('Thickness\n')
		fo.write(str(thickness) + '\n')
		fo.write('Scattering factor accuracy\n')
		fo.write(str(potential_method) + '\n')
		fo.write('Calculation type\n')
		fo.write(str(beam_type) + '\n')
	
		if param['beam_type'] == 2:
			fo.write('Cs coefficient\n')
			fo.write(str(cs) + '\n')
			fo.write('C5 coefficient\n')
			fo.write(str(c5) + '\n')
			fo.write('aperture cutoff\n')
			fo.write(str(convergence_angle) + '\n')
			fo.write('Defocus\n')
			fo.write(str(defocus) + '\n')
			fo.write('Change probe lens parameters\n')
			fo.write('0\n')

		fo.write('<1> QEP <2> ABS\n')
		fo.write(str(tds_method) + '\n')

		if param['tds_method'] == 2:
			fo.write('<1> Absorption <2> No absorption\n')
			fo.write(str(absorp_choice) + '\n')

		fo.write('Tile supercell x\n')
		fo.write(str(tile_x) + '\n')
		fo.write('Tile supercell y\n')
		fo.write(str(tile_y) + '\n')
		fo.write('Number of pixels in x\n')
		fo.write(str(pixel_x) + '\n')
		fo.write('Number of pixels in y\n')
		fo.write(str(pixel_y) + '\n')
		fo.write('<1> Continue <2> Change\n')
		fo.write('1\n')
		fo.write('<1> Beam tilt <2> No beam tilt\n')
		fo.write('2\n')
		fo.write('<1> Specimen tilt <2> No specimen tilt\n')
		fo.write('2\n')
	
		if param['tds_method'] == 1:
			#add qep parameters
			pass

		fo.write('<0> continue <1> save <2> load\n')
		fo.write('0\n')
	
        # parallel beam
		if param['beam_type'] == 1:
			fo.write('Cs coefficient\n')
			fo.write(str(cs) + '\n')
			fo.write('C5 coefficient\n')
			fo.write(str(c5) + '\n')
			fo.write('aperture cutoff\n')
			fo.write(str(covergence_angle) + '\n')
			fo.write('Defocus\n')
			fo.write(str(defocus) + '\n')
			fo.write('Change probe lens parameters\n')
			fo.write('0\n')
			fo.write('<0> Precalculated potentials <1> On-the-fly calculation\n')
			fo.write(str(otf_potential) + '\n')
        
        # convergent beam
		if param['beam_type'] == 2:
			fo.write('<1> CBED <2> PACBED <3> STEM\n')
			fo.write(str(output_type) + '\n')

			if param['output_type'] == 2:
				fo.write('<1> Diffraction pattern for each probe position\n')
				fo.write('0\n')

		#add statements for single cbed and stem 

        fo.close()
		return;
