import numpy as np
import matplotlib.pyplot as plt
import os

# Function to run the makeFiles script
def run_make_files():
    try:
        # Check if the script exists and is executable
        make_files_path = '/home/lollilol/Desktop/fsi_rigid/background_zone/makeFiles'
        if os.path.isfile(make_files_path) and os.access(make_files_path, os.X_OK):
            # Run the script
            subprocess.run(['./makeFiles'], check=True)
            print("makeFiles script executed successfully.")
        else:
            print("Error: makeFiles script not found or not executable.")
    except Exception as e:
        print(f"Error executing makeFiles script: {e}")


# Function to load the generated data from the plotfile
def load_data():
    plotfile_path = "/home/lollilol/Desktop/fsi_rigid/background_zone/plotfile"
    
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
        
        return time, xcenter, ycenter, zcenter, vx, vy, vz, wx, wy, wz
    
    except Exception as e:
        print(f"Error loading data from {plotfile_path}: {e}")
        return None

# Paths to the force and moment data files
force_file_path = "/home/lollilol/Desktop/fsi_rigid/background_zone/postProcessing/forces1/0/force.dat"
moment_file_path = "/home/lollilol/Desktop/fsi_rigid/background_zone/postProcessing/forces1/0/moment.dat"

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

# Load the motion data (center of mass and velocities)
motion_data = load_data()

if motion_data is None:
    print("Error: Motion data could not be loaded.")
else:
    # Unpack motion data
    time, xcenter, ycenter, zcenter, vx, vy, vz, wx, wy, wz = motion_data

    # Create subplots for motion data and force & moment data
    fig, axs = plt.subplots(4, 1, figsize=(10, 16))

    # Plot Center of Mass position
    axs[0].plot(time, xcenter, label='x', color='r')
    axs[0].plot(time, ycenter, label='y', color='g')
    axs[0].plot(time, zcenter, label='z', color='b')
    axs[0].set_title("Centre of Mass")
    axs[0].set_xlabel('Time')
    axs[0].set_ylabel('Centre of Mass [m]')
    axs[0].legend(loc="upper right")
    axs[0].grid()

    # Plot velocity of Centre of Mass
    axs[1].plot(time, vx, label='vx', color='r')
    axs[1].plot(time, vy, label='vy', color='g')
    axs[1].plot(time, vz, label='vz', color='b')
    axs[1].set_title("Velocity of Centre of Mass")
    axs[1].set_xlabel('Time')
    axs[1].set_ylabel('Velocity [m/s]')
    axs[1].legend(loc="upper right")
    axs[1].grid()

    # Plot Angular Velocity of Centre of Mass
    axs[2].plot(time, wx, label='wx', color='r')
    axs[2].plot(time, wy, label='wy', color='g')
    axs[2].plot(time, wz, label='wz', color='b')
    axs[2].set_title("Angular Velocity of Centre of Mass")
    axs[2].set_xlabel('Time')
    axs[2].set_ylabel('Angular Velocity [rad/s]')
    axs[2].legend(loc="upper right")
    axs[2].grid()

    # Total force plot
    fig_force, axs_force = plt.subplots(1, 1, figsize=(10, 5))
    axs_force.plot(time_force, total_fx, label="Total Fx", color="r")
    axs_force.plot(time_force, total_fy, label="Total Fy", color="g")
    axs_force.plot(time_force, total_fz, label="Total Fz", color="b")
    axs_force.set_ylim([0, 6e-4])
    axs_force.set_xlabel("Time")
    axs_force.set_ylabel("Total Force [N]")
    axs_force.legend(loc="upper right")
    axs_force.grid()

    # Total moment plot
    fig_moment, axs_moment = plt.subplots(1, 1, figsize=(10, 5))
    axs_moment.plot(time_moment, total_mx, label="Total Mx", color="r")
    axs_moment.plot(time_moment, total_my, label="Total My", color="g")
    axs_moment.plot(time_moment, total_mz, label="Total Mz", color="b")
    axs_moment.set_ylim([-6e-8, 1e-5])
    axs_moment.set_xlabel("Time")
    axs_moment.set_ylabel("Total Moment [Nm]")
    axs_moment.legend(loc="upper right")
    axs_moment.grid()

    # Adjust layout to prevent overlap
    plt.tight_layout()

    # Show all plots
    plt.show()
