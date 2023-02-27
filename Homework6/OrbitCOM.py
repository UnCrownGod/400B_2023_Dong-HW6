# Homework 6
# Ezekiel Dong




# import modules
import numpy as np
import astropy.units as u
from astropy.constants import G

# import plotting modules
import matplotlib.pyplot as plt
import matplotlib

# my modules
from ReadFile import Read
# Step 1: modify CenterOfMass so that COM_P now takes a parameter specifying 
# by how much to decrease RMAX instead of a factor of 2
from CenterOfMass2 import CenterOfMass




def OrbitCOM(galaxy, start, end, snap, n)
    """
    Goal:
        function that loops over all the desired snapshots to compute the COM pos and vel as a function of time.
        
    Inputs:
        galaxy(str)--name of galaxy
        start(int)--number of first snapshot
        end(int)--number of last snapshot
        n(int)--the intervial 
        
    Outputs: 
        file(file)--file of center of mass and velocity at each snapshot by step
    """
    
    # compose the filename for output
    fileout = 'orbit_' + galaxy + '.txt'
    
    #  set tolerance and VolDec for calculating COM_P in CenterOfMass
    # for M33 that is stripped more, use different values for VolDec
    delta = 0.1
    volDec = 2
    if galaxy == 'M33':
        volDec = 4
    
    # generate the snapshot id sequence 
    # it is always a good idea to also check if the input is eligible (not required)
    snap_ids = np.arange(start, end)
    if len(snap_ids) == 0:# if empty
        os.exit()# code will stop while empty
    
    # initialize the array for orbital info: t, x, y, z, vx, vy, vz of COM
    orbit = np.zeros([len(snap_ids), 7])# empty array
    
    # a for loop 
    for i, snap_id in enumerate(snap_ids):# loop over files
        
        # compose the data filename (be careful about the folder)
        
        # add a string of the filenumber to the value “000” , code from HW5
        ilbl = '000' + str(snap)
        # remove all but the last 3 digits
        ilbl = ilbl[-3:]
        filename = '%s_' % (galaxy) + ilbl + '.txt'
                
        # read data in the given file using Read
        time, total_particles, data = Read(filename)
        
        # Initialize an instance of CenterOfMass class, using disk particles
        
        COM = CenterOfMass(filename, 2)
        
        # Store the COM pos and vel. Remember that now COM_P required VolDec
       
        COM_Position = COM.COM_P(0.1, volDec)
        COM_Velocity = COM.COM_V(COM_p[0], COM_p[1], COM_p[2])
        
        # store the time, pos, vel in ith element of the orbit array,  without units (.value)
        # note that you can store 
        # a[i] = var1, *tuple(array1)
        orbit[i] = COM.time.value / 1000, * tuple(COM_p.value), * tuple(COM_v.value)

        
        # print snap_id to see the progress
        print(snap_id)
        
    # write the data to a file
    # we do this because we don't want to have to repeat this process 
    # this code should only have to be called once per galaxy.
    np.savetxt(fileout, orbit, fmt = "%11.3f"*7, comments='#',
               header="{:>10s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}"\
                      .format('t', 'x', 'y', 'z', 'vx', 'vy', 'vz'))




# Recover the orbits and generate the COM files for each galaxy
# read in 800 snapshots in intervals of n=5
# Note: This might take a little while - test your code with a smaller number of snapshots first! 
OrbitCOM('MW', 0, 800, 5)
OrbitCOM('M31', 0, 800, 5)
OrbitCOM('M33', 0, 800, 5)



# Read in the data files for the orbits of each galaxy that you just created
# headers:  t, x, y, z, vx, vy, vz
# using np.genfromtxt
MW_Orbit = np.genfromtxt('Orbit_MW.txt')
M31_Orbit = np.genfromtxt('Orbit_M31.txt')
M33_Orbit = np.genfromtxt('Orbit_M33.txt')



# function to compute the magnitude of the difference between two vectors 
# You can use this function to return both the relative position and relative velocity for two 
# galaxies over the entire orbit  
def magnitude_diff(vec1, vec2)
    '''
    Goal:
        find the diff of two vec
    
    Input:
        vec1(array)--first vec
        vec2(array)--second vec
    
    Ourput:
        magnitude_pos(array)--position magnitude
        magnitude-vel(array)--velocity magnitude
    '''
    magnitude_pos = np.sqrt((vec1['x'] - vec2['x']) ** 2 + (vec1['y'] - vec2['y']) ** 2 + (vec1['z'] - vec2['z']) ** 2)
    magnitude_vel= np.sqrt((vec1['vx'] - vec2['vx']) ** 2 + (vec1['vy'] - vec2['vy']) ** 2 + (vec1['vz'] - vec2['vz']) ** 2)
    return magnitude_pos, magnitude_vel



# Determine the magnitude of the relative position and velocities 
# of MW and M31
MW_M31_mag_pos, MW_M31_mag_vel = magnitude_diff(MW_Orbit, M31_Orbit)

# Determine the magnitude of the relative position and velocities 
# of M33 and M31
M31_M33_mag_pos, M31_M33_mag_vel = magnitude_diff(M31_Orbit, M33_Orbit)

# Plot the magnitude of the separation
plt.title('Separation')
plt.plot(MW_Orbit['t'], MW_M31_mag_pos, 'b', label = 'M31 - MW')
plt.plot(M33_Orbit['t'], M31_M33_mag_pos, 'c', label = 'M33 - M31')
plt.xlabel('time(Gyrs)')
plt.ylabel('position sepration(kpc)')

# Plot the magnitude of the relative velocity
plt.title('Relative velocity')
plt.plot(MW_Orbit['t'], MW_M31_mag_vel, 'b', label = 'M31 - MW')
plt.plot(M33_Orbit['t'], M31_M33_mag_vel, 'c', label = 'M33 - M31')
plt.xlabel('time(Gyrs)')
plt.ylabel('relative velocity(km/s)')