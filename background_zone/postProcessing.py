
import numpy as np
import matplotlib.pyplot as plt
import os
import subprocess

# Run the sequence of commands to get motion data
subprocess.run("cat log.simulation | grep -w 'Time =' | cut -d' ' -f3 | tr -d '(' > time", shell=True)
subprocess.run("cat log.simulation | grep 'of mass' | cut -d' ' -f8 | tr -d '(' > xcenter", shell=True)
subprocess.run("cat log.simulation | grep 'of mass' | cut -d' ' -f9 | tr -d '(' > ycenter", shell=True)
subprocess.run("cat log.simulation | grep 'of mass' | cut -d' ' -f10 | tr -d ')' > zcenter", shell=True)
subprocess.run("cat log.simulation | grep 'Linear' | cut -d' ' -f7 | tr -d '(' > vx", shell=True)
subprocess.run("cat log.simulation | grep 'Linear' | cut -d' ' -f8 | tr -d ' ' > vy", shell=True)
subprocess.run("cat log.simulation | grep 'Linear' | cut -d' ' -f9 | tr -d ')' > vz", shell=True)
subprocess.run("cat log.simulation | grep 'Angular' | cut -d' ' -f7 | tr -d '(' > wx", shell=True)
subprocess.run("cat log.simulation | grep 'Angular' | cut -d' ' -f8 | tr -d '(' > wy", shell=True)
subprocess.run("cat log.simulation | grep 'Angular' | cut -d' ' -f9 | tr -d ')' > wz", shell=True)
subprocess.run("cat log.simulation | grep 'Courant Number' | awk '{print $4}' > Co_mean", shell=True)
subprocess.run("cat log.simulation | grep 'Courant Number' | awk '{print $6}' > Co_max", shell=True)

subprocess.run("cat log.simulation | grep -w 'Time =' | cut -d' ' -f3 | tr -d '(' > time", shell=True)
subprocess.run("cat log.simulation | grep 'Orientation' | awk '{print $2}' | tr -d '()' | awk -F' ' '{print $1}' > r11", shell=True)
subprocess.run("cat log.simulation | grep 'Orientation' | awk '{print $3}' | tr -d '()' | awk -F' ' '{print $1}' > r12", shell=True)
subprocess.run("cat log.simulation | grep 'Orientation' | awk '{print $4}' | tr -d '()' | awk -F' ' '{print $1}' > r13", shell=True)
subprocess.run("cat log.simulation | grep 'Orientation' | awk '{print $5}' | tr -d '()' | awk -F' ' '{print $1}' > r21", shell=True)
subprocess.run("cat log.simulation | grep 'Orientation' | awk '{print $6}' | tr -d '()' | awk -F' ' '{print $1}' > r22", shell=True)
subprocess.run("cat log.simulation | grep 'Orientation' | awk '{print $7}' | tr -d '()' | awk -F' ' '{print $1}' > r23", shell=True)
subprocess.run("cat log.simulation | grep 'Orientation' | awk '{print $8}' | tr -d '()' | awk -F' ' '{print $1}' > r31", shell=True)
subprocess.run("cat log.simulation | grep 'Orientation' | awk '{print $9}' | tr -d '()' | awk -F' ' '{print $1}' > r32", shell=True)
subprocess.run("cat log.simulation | grep 'Orientation' | awk '{print $10}' | tr -d '()' | awk -F' ' '{print $1}' > r33", shell=True)

# Create plotfiles with data to plot
subprocess.run("paste time r11 r12 r13 r21 r22 r23 r31 r32 r33 > plotfile2", shell=True)
subprocess.run("paste time xcenter ycenter zcenter vx vy vz wx wy wz Co_mean Co_max > plotfile", shell=True)

# Function to load the generated data from the first plotfile
def load_data_():
    plotfile_path = "/home/dell/Desktop/fsi_6dof_validation/background_zone/plotfile"
    
    # Attempt to load the data with error handling for inconsistent rows
    try:
        # Use np.genfromtxt with 'invalid_raise=False' to skip rows with inconsistent column count
        data = np.genfromtxt(plotfile_path, delimiter=None, invalid_raise=False)  # Skip problematic rows
        
        # Check if data has the correct number of columns (10 columns expected)
        if data.ndim == 1:
            print("Data appears to be one-dimensional. Checking file format...")
            with open(plotfile_path, 'r') as f:
                for i, line in enumerate(f.readlines()[:5]):
                    print(f"Line {i}: {line.strip()}")
            raise ValueError("Data is not in the expected multi-column format.")
        
        # Verify that each row has exactly 10 columns, otherwise skip
        expected_columns = 12
        data = [row for row in data if len(row) == expected_columns]
        
        if len(data) == 0:
            raise ValueError("No valid data found with the correct number of columns.")
        
        data = np.array(data)  # Convert to numpy array for further processing
        
        # Extract columns for time, center of mass coordinates, linear velocities, and angular velocities
        time = data[:, 0]
        xcenter = data[:, 1]
        ycenter = data[:, 2]
        zcenter = data[:, 3]
        vx = data[:, 4]
        vy = data[:, 5]
        vz = data[:, 6]
        wx = data[:, 7]
        wy = data[:, 8]
        wz = data[:, 9]
        Co_mean = data[:, 10]
        #Co_max = data[:, 11]              
        return time, xcenter, ycenter, zcenter, vx, vy, vz, wx, wy, wz, Co_mean
    
    except Exception as e:
        print(f"Error loading data from {plotfile_path}: {e}")
        return None

# Function to load the generated data from the plotfile2
def load_data():
    plotfile_path = "/home/dell/Desktop/fsi_6dof_validation/background_zone/plotfile2"
    
    # Attempt to load the data with error handling for inconsistent rows
    try:
        # Use np.genfromtxt with 'invalid_raise=False' to skip rows with inconsistent column count
        data = np.genfromtxt(plotfile_path, delimiter=None, invalid_raise=False)  # Skip problematic rows
        
        # Check if data has the correct number of columns (10 columns expected)
        if data.ndim == 1:
            print("Data appears to be one-dimensional. Checking file format...")
            with open(plotfile_path, 'r') as f:
                for i, line in enumerate(f.readlines()[:5]):
                    print(f"Line {i}: {line.strip()}")
            raise ValueError("Data is not in the expected multi-column format.")
        
        # Verify that each row has exactly 10 columns, otherwise skip
        expected_columns = 10
        data = [row for row in data if len(row) == expected_columns]
        
        if len(data) == 0:
            raise ValueError("No valid data found with the correct number of columns.")
        
        data = np.array(data)  # Convert to numpy array for further processing
        print(data)
        # Extract columns for time, center of mass coordinates, linear velocities, and angular velocities
        time = data[:, 0]
        r11 = data[:, 1]
        r12 = data[:, 2]
        r13 = data[:, 3]
        r21 = data[:, 4]
        r22 = data[:, 5]
        r23 = data[:, 6]
        r31 = data[:, 7]
        r32 = data[:, 8]
        #r33 = data[:, 9]
 
        return time, r11, r12, r13, r21, r22, r23, r31, r32, r33 
    
    except Exception as e:
        print(f"Error loading data from {plotfile_path}: {e}")
        return None


# Load data from the combined plotfile2
data = np.loadtxt("plotfile2", unpack=True)
time = data[0]
r11, r12, r13 = data[1], data[2], data[3]
r21, r22, r23 = data[4], data[5], data[6]
r31, r32, r33 = data[7], data[8], data[9]

# Compute Euler angles (roll, pitch, yaw)
roll = np.arctan2(r32, r33)  # coning angle 
pitch = np.arcsin(-r31)      # Pitch angle
yaw = np.arctan2(r21, r11)   # Yaw 

# Convert angles to degrees for readability
roll_deg = np.degrees(roll)
pitch_deg = np.degrees(pitch)
yaw_deg = np.degrees(yaw)


# Paths to the force and moment data files
force_file_path = "/home/dell/Desktop/fsi_6dof_validation/background_zone/postProcessing/forces1/0/force.dat"
moment_file_path = "/home/dell/Desktop/fsi_6dof_validation/background_zone/postProcessing/forces1/0/moment.dat"

# Load the force data
force_data = np.genfromtxt(force_file_path, skip_header=4, invalid_raise=False)

# Extract force components from the data
time_force = force_data[:, 0]
total_fx, total_fy, total_fz = force_data[:, 1], force_data[:, 2], force_data[:, 3]

# Load the moment data
moment_data = np.genfromtxt(moment_file_path, skip_header=4, invalid_raise=False)

# Extract moment components from the moment data
time_moment = moment_data[:, 0]
total_mx, total_my, total_mz = moment_data[:, 1], moment_data[:, 2], moment_data[:, 3]

# Load the motion data from plotfile
motion_data = load_data_()

if motion_data is None:
    print("Error: Motion data could not be loaded.")
else:
    # Unpack motion data
    time, xcenter, ycenter, zcenter, vx, vy, vz, wx, wy, wz, Co_mean= motion_data

    # Create subplots for motion data and force & moment data
    fig, axs = plt.subplots(3, 1, figsize=(10, 16))

    # Plot Center of Mass position
    axs[0].plot(time, xcenter, label='x', color='r')
    axs[0].plot(time, ycenter, label='y', color='g')
    axs[0].plot(time, zcenter, label='z', color='b')
    axs[0].set_title("Centre of Mass [m]", fontsize=12)
    axs[0].set_xlabel('Time', fontsize=12)
    axs[0].legend(loc="upper right", fontsize=10)

    # Plot velocity of Centre of Mass
    axs[1].plot(time, vx, label='vx', color='r')
    axs[1].plot(time, vy, label='vy', color='g')
    axs[1].plot(time, vz, label='vz', color='b')
    axs[1].set_title("Linear Velocity [m/s]", fontsize=12)
    axs[1].set_xlabel('Time', fontsize=12)
    axs[1].legend(loc="upper right", fontsize=10)

    # Plot Angular Velocity of Centre of Mass
    axs[2].plot(time, wx, label='wx', color='r')
    axs[2].plot(time, wy, label='wy', color='g')
    axs[2].plot(time, wz, label='wz', color='b')
    axs[2].set_title("Angular Velocity [rad/s]", fontsize=12)
    axs[2].set_xlabel('Time', fontsize=12)
    axs[2].legend(loc="upper right", fontsize=10)

    # Total force plot
    fig_force, axs_force = plt.subplots(1, 1, figsize=(10, 5))
    axs_force.plot(time_force, total_fx, label="Total Fx", color="r")
    axs_force.plot(time_force, total_fy, label="Total Fy", color="g")
    axs_force.plot(time_force, total_fz, label="Total Fz", color="b")
    axs_force.set_xlabel("Time", fontsize=12)
    axs_force.set_ylabel("Total Force [N]", fontsize=12)
    axs_force.legend(loc="upper right", fontsize=10)

    # Total moment plot
    fig_moment, axs_moment = plt.subplots(1, 1, figsize=(10, 5))
    axs_moment.plot(time_moment, total_mx, label="Total Mx", color="r")
    axs_moment.plot(time_moment, total_my, label="Total My", color="g")
    axs_moment.plot(time_moment, total_mz, label="Total Mz", color="b")
    #axs_moment.set_ylim([-6e-8, 1e-5])
    axs_moment.set_xlabel("Time", fontsize=12)
    axs_moment.set_title("Total Moment [Nm]", fontsize=12)
    axs_moment.legend(loc="upper right", fontsize=12)

# Plot Euler angles
fig, axs = plt.subplots(3, 1, figsize=(10, 16))

# Roll
axs[0].plot(time, pitch_deg, label="coning angle", color='r')
axs[0].set_title("Coning angle time evolution", fontsize=14)
axs[0].set_xlabel("Time [s]", fontsize=12)
axs[0].set_title("Angle [degrees]", fontsize=12)
axs[0].legend(fontsize=10)
axs[0].grid(False)

# Pitch
axs[1].plot(time, roll_deg, label="Pitch", color='g')
axs[1].set_title("Pitch angle time evolution", fontsize=14)
axs[1].set_xlabel("Time [s]", fontsize=12)
axs[1].set_ylabel("Angle [degrees]", fontsize=12)
axs[1].legend(fontsize=10)
axs[1].grid(False)

# Yaw
axs[2].plot(time, yaw_deg, label="Yaw", color='b')
axs[2].set_title("Yaw angle time evolution", fontsize=14)
axs[2].set_xlabel("Time [s]", fontsize=12)
axs[2].set_ylabel("Angle [degrees]", fontsize=12)
axs[2].legend(fontsize=10)
axs[2].grid(False)

plt.tight_layout()
plt.show()

