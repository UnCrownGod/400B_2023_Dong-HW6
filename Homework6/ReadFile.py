## ASTR400B Homework 2
## This script is designed to open and read a data file and return specific data as result

# import all packages we need in the following code
import numpy as np
import astropy.units as u

# define read function
def Read(file_name):
    file = open(file_name, 'r') # open the file
    
    line1 = file.readline()# read line by line
    label1, value1 = line1.split()# split label and value
    line2 = file.readline()# read line by line
    label2, value2 = line2.split()#split label and value
    
    time = float(value1) * u.Myr# apply unit
    total_particles = float(value2)# total num of particles
    file.close()# close file
    
    data = np.genfromtxt(file_name, dtype=None, names=True, skip_header=3)#skipping header
    
    return time, total_particles, data # return data we need in particle_properties