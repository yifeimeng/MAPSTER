import subprocess
import msParam.py
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

convergence_angle_list = [3.5, 5, 10, 17]
num_xtl_file = 20

#load the sample information into the Sample table
for i in range(0, num_xtl_file)
     new_sample = Sample(id = i, crystal_file = 'BTO_GSO_100' + '_' + str(i) + '.xtl', var_1 = j)
     session.add(new_sample)
     session.commit()

for i in range(len(convergence_angle_list))
    # update the parameter
    currParam.convergence_angle = convergence_angle_list[i]
    for j in range(0, num_xtl_file)
        # update the sample
        currParam.crystal_file_name = 'BTO_GSO_100' + '_' + str(i) + '.xtl'
        
        # define the output prefix
        currParam.output_name = 'BTO_GSO' + '_p' + str(i) + '_s' + str(j)
        # perform the calculation
        currParam.writeParam()
        subprocess.call("./mu_STEM nopause")
        
        # load the parameter information into the Parameter table, considering the automatic thickness calculation
        for i_thickness in range(0, n_thickness)
            currThickness =
            new_parameter = currParam.createParameter()
            session.add(new_parameter)
            session.commit()





