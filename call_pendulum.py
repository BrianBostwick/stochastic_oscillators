import ctypes
import numpy as np
import matplotlib.pyplot as plt
import os
import datetime

# Load the shared library (adjust if needed for Windows)
lib = ctypes.CDLL("./libpendulum.so")  # Linux/macOS

# Define function prototype: void simulate_pendulum(double, double, int, double*, double*, double*, int)
lib.simulate_pendulum.argtypes = [
    ctypes.c_double,  # gamma
    ctypes.c_double,  # amplitude
    ctypes.c_int,     # num_kicks
    ctypes.POINTER(ctypes.c_double),  # t_out
    ctypes.POINTER(ctypes.c_double),  # theta_out
    ctypes.POINTER(ctypes.c_double),  # omega_out
    ctypes.c_int      # num_steps
]


# Create a timestamped directory under 'data/YYYY-MM-DD/'
def create_plots_directory():
    date_folder = datetime.datetime.now().strftime("%Y-%m-%d")  # e.g., "2025-03-14"
    base_dir = "plots"
    plots_dir = os.path.join(base_dir, date_folder)
    os.makedirs(plots_dir, exist_ok=True)  # Create the directory if it doesn't exist
    return plots_dir

# Create a timestamped filename for plots
def generate_plot_filename(plots_dir, scan_number, prefix, extension = "png"):
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S") # e.g., 2024-03-14_15-30-00
    filename = f"{plots_dir}/scan{scan_number}_{prefix}_{timestamp}.{extension}"
    return filename

# Create a timestamped directory under 'data/YYYY-MM-DD/'
def create_data_directory():
    date_folder = datetime.datetime.now().strftime("%Y-%m-%d")  # e.g., "2025-03-14"
    base_dir = "data"
    data_dir = os.path.join(base_dir, date_folder)
    os.makedirs(data_dir, exist_ok=True)  # Create the directory if it doesn't exist
    return data_dir

# Automatically determine the next available scan number
def get_next_scan_number(data_dir):
    existing_scans = [
        int(f.split("_")[0][4:])  # Extract scan number from filenames like "scan1_..."
        for f in os.listdir(data_dir) if f.startswith("scan") and f.endswith(".txt")
    ]
    return max(existing_scans, default=0) + 1  # Increment the last scan number

# Create a timestamped filename
def generate_filename(data_dir, scan_number, prefix, extension):
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")  # e.g., 2024-03-14_15-30-00
    filename = f"{data_dir}/scan{scan_number}_{prefix}_{timestamp}.{extension}"
    return filename

# Ensure the file does not exist and create a new file
def get_unique_filename(data_dir, scan_number, prefix, extension):
    filename = generate_filename(data_dir, scan_number, prefix, extension)
    while os.path.exists(filename):  # Check if the file already exists
        filename = generate_filename(prefix, extension)  # Generate a new timestamped name
    return filename

# Create the data directory for saving the simulations.
data_dir = create_data_directory()

# Simulation parameters
gamma_values = [1.0, 10.0, 100.0]  # Damping values
amplitude_values = [10.0, 50.0, 100.0]  # Kick amplitudes
num_kicks_values = [100020, 100000, 10000000]  # Number of kicks
num_steps = int(30.0 / 0.01)  # 20s simulation with dt=0.01

# Allocate memory for output
t_out = (ctypes.c_double * num_steps)()
theta_out = (ctypes.c_double * num_steps)()
omega_out = (ctypes.c_double * num_steps)()

# Create a figure with multiple subplots
fig, axes = plt.subplots(len(gamma_values), len(amplitude_values), figsize=(15, 10), sharex=True, sharey=True)

run_index = 0

# Determine the next available scan number
scan_number = get_next_scan_number(data_dir)

for i, gamma in enumerate(gamma_values):
    for j, amplitude in enumerate(amplitude_values):
        num_kicks = num_kicks_values[j]  # Adjust kicks based on amplitude

        print(f"Running: γ={gamma}, A={amplitude}, Kicks={num_kicks}")
        
        # Call the C function
        lib.simulate_pendulum(gamma, amplitude, num_kicks, t_out, theta_out, omega_out, num_steps)

        # Convert to NumPy arrays
        t = np.array(t_out)
        theta = np.array(theta_out)
        omega = np.array(omega_out)
        
        # Save data to a uniquely named file
        data_filename = get_unique_filename(data_dir, scan_number, f"run{run_index}_pendulum_data", "txt")
        
        np.savetxt(
            data_filename, 
            np.column_stack((t, theta, omega)), 
            header="Time (s)   Theta (radians)   Omega (rad/s)", 
            fmt="%.6f"
            )
        
        # Save metadata file with simulation parameters
        metadata_filename = data_filename.replace(".txt", "_metadata.txt")
        
        with open(metadata_filename, "w") as meta_file:
            meta_file.write(f"Simulation Timestamp: {datetime.datetime.now()}")
            meta_file.write(f"Data File: {data_filename}")
            meta_file.write(f"Gamma (Damping): {i}")
            meta_file.write(f"Amplitude of Kicks: {j}")
            meta_file.write(f"Number of Kicks: {num_kicks}")
            meta_file.write(f"Number of Steps: {num_steps}")
        
        print(f"Simulation data saved to {data_filename}")
        print(f"Metadata saved to {metadata_filename}")
        run_index = run_index + 1
        
        # Load data from file (for verification)
        loaded_data = np.loadtxt(data_filename)
        t_loaded, theta_loaded, omega_loaded = loaded_data[:, 0], loaded_data[:, 1], loaded_data[:, 2]

        # Wrap the data around 2Pi for plotting. 
        theta_wrapped = (theta_loaded + np.pi) % (2 * np.pi) - np.pi

        # Select subplot position
        ax = axes[i, j]

        # Plot angular displacement vs. time
        ax.plot(
            t_loaded, 
            theta_wrapped, 
            '.',
            label=f"γ={gamma}, A={amplitude}, Kicks={num_kicks}", 
            alpha=0.8
            )
        ax.set_xlabel("Time (s)")
        ax.set_ylabel("θ (radians)")
        ax.axhline(np.pi, color="gray", linestyle="--", alpha=1.0)
        ax.axhline(-np.pi, color="gray", linestyle="--", alpha=1.0)
        ax.legend()
        ax.set_ylim(-np.pi - 0.5, np.pi + 0.5)
        ax.grid()

# Adjust layout for better appearance
data_filename = get_unique_filename(data_dir, scan_number, f"run{run_index}_pendulum_data", "txt")
plt.suptitle("Nonlinear Pendulum with Stochastic Kicks", fontsize=16)
plt.tight_layout(rect=[0, 0, 1, 0.97])

plots_dir = create_plots_directory()
prefix_plots  = "kicked_pendulum"
plt.savefig(generate_plot_filename(plots_dir, scan_number, prefix_plots))
plt.show()
