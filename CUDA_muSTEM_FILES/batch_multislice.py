import subprocess


param = {'Device': 0}
param['output_name'] = 'test_output'
param['crystal_file_name'] = "SrTiO3_001_300keV.xtl"
param['num_slices'] = 2
param['thickness'] = 200
param['potential_method'] = 1 # 1 is reciprocal space method (accuracy), 2 is hybrid method (speed)
param['beam_type'] = 2 # 1 is parallel beam, 2 is convergent beam
param['Cs'] = 0
param['C5'] = 0
param['convergence_angle'] = 9.5
param['defocus'] = 0
param['tds_method'] = 2 # 1 is qep, 2 is absorption model
param['absorp_choice'] = 1 # 1 is absorption, 2 is no absorption
param['tile_x'] = 8
param['tile_y'] = 8
param['pixel_x'] = 256
param['pixel_y'] = 256
param['output_type'] = 2 # convergent beam only, 1 is CBED, 2 is PACBED, 3 is STEM


print('Start batch job...')
#create the parameter file, define the file name here
file_param = open('input_parameter.txt', 'w+')

writeParameterFile(file_param, param)

file_param.close()


subprocess.call("./mu_STEM nopause")



