# genfind.py
#
# A function that generates files that match a given filename pattern

import os
import shutil
import fnmatch

def gen_find(filepat,top):
    for path, dirlist, filelist in os.walk(top):
        for name in fnmatch.filter(filelist,filepat):
            yield os.path.join(path,name)



 
