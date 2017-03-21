import subprocess
import shutil
from multislice_crystal import writeXtl
from multislice_param import Param
from genfind import gen_find
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from sqlalchemy_declarative import Base, Sample, Parameter

engine = create_engine('sqlite:///sqlalchemy_test.db')
Base.metadata.bind = engine

DBSession = sessionmaker(bind = engine)
session = DBSession()

# initilize an instance of Param
currParam = Param()

print('Start batch job...')
#create the parameter file, define the file name here

sampleName = 'BTO_GSO_100'
convergence_angle_list = [3.5, 5, 10, 17]
num_xtl_file = 21

#load the sample information into the Sample table
for i in range(0, num_xtl_file):
	crystalFileName = sampleName + '_var_' + str(i) + '.xtl'
	new_sample = Sample(name = sampleName, crystal_file = crystalFileName, var_1 = i)
	session.add(new_sample)
	# create the relavant crystal file
	writeXtl(crystalFileName, 0)

session.commit()

# update any simulation parameters here
currParam.thickness = 500

# calculate the number of thickness layers
# this design is implemented in the fortran codes
# output one thickness result for every 8 slices
n_thickness = int(int(currParam.thickness/4.0760 + 0.5)*currParam.num_slice/8)

for i in range(len(convergence_angle_list)):
    # update the parameter
    currParam.convergence_angle = convergence_angle_list[i]
    for j in range(0, num_xtl_file):
        # update the xtl file name
	crystalFileName = sampleName + '_var_' + str(i) + '.xtl'
        currParam.crystal_file_name = crystalFileName
        
        # define the output prefix
        currParam.output_name = 'BTO_GSO' + '_p' + str(i) + '_s' + str(j)

        # perform the calculation
        currParam.writeParam()
        subprocess.call(['./mu_STEM', 'nopause'])
        
        # load the parameter information into the Parameter table, considering the automatic thickness calculation
        for i_thickness in range(0, n_thickness):
			new_parameter = currParam.copy2Parameter()
			new_parameter.sample_id = j
			new_parameter.img_name = currParam.output_name + '_z' + str(j).zfill(3) + '_' + str(currParam.pixel_x) + 'x' + str(currParam.pixel_y) + '.bin'
			session.add(new_parameter)
			session.commit()


src = '/home/yifei/Documents/MAPSTER/CUDA_muSTEM_batch' # input
dst = '/home/yifei/Documents/muSTEM_library' # desired     location

filesToMove = gen_find("*.bin",src)
for name in filesToMove:
	shutil.move(name, dst)

