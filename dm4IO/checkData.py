import h5py
from PIL import Image
import numpy as np

fh5 = h5py.File('/home/yifei/Documents/test_data/Capture61_0000.h5', 'r')
dset = fh5['DPStack']
array_1d = dset[0,:]
array_2d = np.reshape(array_1d, (-1, 896))
print(array_2d)
img = Image.fromarray(array_2d)
img.show()
