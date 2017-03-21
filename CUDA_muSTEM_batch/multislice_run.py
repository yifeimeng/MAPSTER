import subprocess
from multislice_param import Param

# initilize an instance of Param
currParam = Param()

print('Start batch job...')
#create the parameter file, define the file name here

convergence_angle_list = [3.5, 5, 10, 17]
currParam.convergence_angle = convergence_angle_list[0]
currParam.tile_x = 20
currParam.tile_y = 20
currParam.crystal_file_name = 'BTO_GSO_100' + '_' + '0' + '.xtl'
currParam.output_name = 'BTO_GSO_test'
currParam.writeParam()
subprocess.call(['./mu_STEM', 'nopause'])
