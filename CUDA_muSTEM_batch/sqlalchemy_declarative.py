import os
import sys
from sqlalchemy import Column, ForeignKey, Integer, Float, String
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship
from sqlalchemy import create_engine

Base = declarative_base()

class Sample(Base):
	__tablename__ = 'sample'
	id = Column(Integer, primary_key = True)
	name = Column(String(20))
	var_1 = Column(Integer)
	var_2 = Column(Integer)
	var_3 = Column(Integer)
	crystal_file = Column(String(50))

class Parameter(Base):
	__tablename__ = 'parameter'	
	id = Column(Integer, primary_key = True)
	# Simulation model parameter
	num_slice = Column(Integer)
	thickness = Column(Integer)
	potential_method = Column(Integer)
	beam_type = Column(Integer)
	tds_method = Column(Integer)
	tile_x = Column(Integer)
	tile_y = Column(Integer)
	pixel_x = Column(Integer)
	pixel_y = Column(Integer)
	beam_tilt = Column(Integer)
	specimen_tilt = Column(Integer)
	# Instrument parameters
	cs = Column(Float)
	c5 = Column(Float)
	convergence_angle = Column(Float)
	defocus = Column(Float)
	# Sample parameters
	sample_id = Column(Integer, ForeignKey('sample.id'))
	sample = relationship(Sample)
    	# Image file name
    	img_name = Column(String(20))

engine = create_engine('sqlite:///sqlalchemy_test.db', echo = True)

Base.metadata.create_all(engine)


	
	
