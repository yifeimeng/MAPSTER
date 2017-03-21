from sqlalchemy_declarative import Base, Sample, Parameter

class Param:
	device = 0
	output_name = 'N/A'
	crystal_file_name = 'N/A'
	num_slice = 2
	thickness = 200
	potential_method = 1 # 1 is reciprocal space method (accuracy), 2 is hybrid method (speed)
	beam_type = 2
	cs = 0
	c5 = 0
	convergence_angle = 10
	defocus = 0
	tds_method = 2 # 1 is qep, 2 is absorption model
	absorp_choice = 1
	tile_x = 10
	tile_y = 10
	pixel_x = 256
	pixel_y = 256
	beam_tilt = 0
	specimen_tilt = 0
	output_type = 2 # convergent beam only, 1 is CBED, 2 is PACBED, 3 is STEM
	otf_potential = 0 # on-the-fly calculation of the potential, 1 is otf, 0 is precalculated (faster)

	def writeParam(self):
        	fo = open('input_param.txt', 'w+')
		fo.write('Device used for calculation\n')
		fo.write(str(self.device) + '\n')
		fo.write('Output filename\n')
		fo.write(self.output_name + '\n')
		fo.write('Input crystal file name\n')
		fo.write(self.crystal_file_name + '\n')
		fo.write('Slice unit cell <1> yes <2> no\n')
		fo.write('1\n')
		fo.write('Number of distinct potentials\n')
		fo.write(str(self.num_slice) + '\n')
		fo.write('<1> manual or <2> auto slicing\n')
		fo.write('2\n')
		fo.write('Thickness\n')
		fo.write(str(self.thickness) + '\n')
		fo.write('Scattering factor accuracy\n')
		fo.write(str(self.potential_method) + '\n')
		fo.write('Calculation type\n')
		fo.write(str(self.beam_type) + '\n')
	
		if self.beam_type == 2:
			fo.write('Cs coefficient\n')
			fo.write(str(self.cs) + '\n')
			fo.write('C5 coefficient\n')
			fo.write(str(self.c5) + '\n')
			fo.write('aperture cutoff\n')
			fo.write(str(self.convergence_angle) + '\n')
			fo.write('Defocus\n')
			fo.write(str(self.defocus) + '\n')
			fo.write('Change probe lens parameters\n')
			fo.write('0\n')

		fo.write('<1> QEP <2> ABS\n')
		fo.write(str(self.tds_method) + '\n')

		if self.tds_method == 2:
			fo.write('<1> Absorption <2> No absorption\n')
			fo.write(str(self.absorp_choice) + '\n')

		fo.write('Tile supercell x\n')
		fo.write(str(self.tile_x) + '\n')
		fo.write('Tile supercell y\n')
		fo.write(str(self.tile_y) + '\n')
		fo.write('Number of pixels in x\n')
		fo.write(str(self.pixel_x) + '\n')
		fo.write('Number of pixels in y\n')
		fo.write(str(self.pixel_y) + '\n')
		fo.write('<1> Continue <2> Change\n')
		fo.write('1\n')
		fo.write('<1> Beam tilt <2> No beam tilt\n')
		fo.write('2\n')
		fo.write('<1> Specimen tilt <2> No specimen tilt\n')
		fo.write('2\n')
	
		if self.tds_method == 1:
			#add qep parameters
			pass

		fo.write('<0> continue <1> save <2> load\n')
		fo.write('0\n')
	
       		# parallel beam
		if self.beam_type == 1:
			fo.write('Cs coefficient\n')
			fo.write(str(self.cs) + '\n')
			fo.write('C5 coefficient\n')
			fo.write(str(self.c5) + '\n')
			fo.write('aperture cutoff\n')
			fo.write(str(self.convergence_angle) + '\n')
			fo.write('Defocus\n')
			fo.write(str(self.defocus) + '\n')
			fo.write('Change image lens parameters\n')
			fo.write('0\n')
			fo.write('<0> Precalculated potentials <1> On-the-fly calculation\n')
			fo.write(str(self.otf_potential) + '\n')
        
       		# convergent beam
		if self.beam_type == 2:
			fo.write('<1> CBED <2> PACBED <3> STEM\n')
			fo.write(str(self.output_type) + '\n')

			if self.output_type == 2:
				fo.write('<1> Diffraction pattern for each probe position\n')
				fo.write('0\n')

		#add statements for single cbed and stem 

        	fo.close()
		return;
	
	def copy2Parameter(self):
		new_parameter = Parameter()
		new_parameter.num_slice = self.num_slice
		new_parameter.thickness = self.thickness
		new_parameter.potential_method = self.potential_method
		new_parameter.beam_type = self.beam_type
		new_parameter.tds_method = self.tds_method
		new_parameter.tile_x = self.tile_x
		new_parameter.tile_y = self.tile_y
		new_parameter.pixel_x = self.pixel_x
		new_parameter.pixel_y = self.pixel_y
		new_parameter.beam_tilt = self.beam_tilt
		new_parameter.specimen_tilt = self.specimen_tilt
		new_parameter.cs = self.cs
		new_parameter.c5 = self.c5
		new_parameter.convergence_angle = self.convergence_angle
		new_parameter.defocus = self.defocus
		return new_parameter
		
