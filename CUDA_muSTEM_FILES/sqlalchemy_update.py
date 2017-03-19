from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from sqlalchemy_declarative import Base, Sample, Parameter

engine = create_engine('sqlite:///sqlalchemy_test.db')
Base.metadata.bind = engine

DBSession = sessionmaker(bind = engine)
session = DBSession()

new_sample = Sample(name = 'BTO')
session.add(new_sample)
session.commit()

new_parameter = Parameter(crystal_file = 'BTO.xtl')
session.add(new_parameter)
session.commit()
