#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 16 14:48:17 2025

@author: brianbostwick
"""

import ctypes
import numpy as np
import matplotlib.pyplot as plt
import os
import datetime

# Create a timestamped directory under 'data/YYYY-MM-DD/'
def create_data_directory():
    date_folder = datetime.datetime.now().strftime("%Y-%m-%d")  # e.g., "2025-03-14"
    base_dir = "data"
    data_dir = os.path.join(base_dir, date_folder)
    os.makedirs(data_dir, exist_ok=True)  # Create the directory if it doesn't exist
    return data_dir

# Load the shared library
lib = ctypes.CDLL("./liblatticependulum.so")  # Adjust for Windows if needed

# Define function prototype
lib.simulate_pendulums.argtypes = [
    ctypes.c_double,  # gamma
    ctypes.c_double,  # amplitude
    ctypes.c_int,     # num_kicks
    ctypes.POINTER(ctypes.c_double),  # t_out
    ctypes.POINTER(ctypes.c_double),  # theta_out
    ctypes.POINTER(ctypes.c_double),  # omega_out
    ctypes.c_int      # num_steps
]

# Simulation parameters
gamma = 0.5
amplitude = 10.0
num_kicks = 10
num_steps = int(120.0 / 0.01)  # 30s simulation with dt=0.01
N = 100  # Number of pendulums in the lattice, this must mach the c code.

# Allocate memory
t_out = (ctypes.c_double * num_steps)()
theta_out = (ctypes.c_double * (num_steps * N))()
omega_out = (ctypes.c_double * (num_steps * N))()

# Call C function
lib.simulate_pendulums(gamma, amplitude, num_kicks, t_out, theta_out, omega_out, num_steps)

# Convert to NumPy arrays
t = np.array(t_out)
theta = np.array(theta_out).reshape((num_steps, N))
omega = np.array(omega_out).reshape((num_steps, N))

theta_wrapped = (theta + np.pi) % (2 * np.pi) - np.pi


# # Save data
# data_dir = create_data_directory()
# os.makedirs("data", exist_ok=True)
# date_str = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
# data_filename = f"{data_dir}/lattice_pendulum_{date_str}.txt"
# np.savetxt(data_filename, np.column_stack((t, theta, omega)), 
#             header="Time (s) " + " ".join([f"Theta_{i}" for i in range(N)]) + " " + " ".join([f"Omega_{i}" for i in range(N)]),
#             fmt="%.6f")
# print(f"Simulation data saved to {data_filename}")


# Plot results
fig, ax = plt.subplots(2, 1, figsize=(10, 8))

# Plot angular displacement
for i in range(N):
    ax[0].plot(t, theta[:, i], label=f"Pendulum {i+1}")
ax[0].set_xlabel("Time (s)")
ax[0].set_ylabel("Theta (radians)")
# ax[0].legend()
ax[0].grid()
ax[0].set_title("Angular Displacement of Coupled Pendulums")

# Plot angular velocity
for i in range(N):
    ax[1].plot(t, omega[:, i], label=f"Pendulum {i+1}")
ax[1].set_xlabel("Time (s)")
ax[1].set_ylabel("Omega (rad/s)")
# ax[1].legend()
ax[1].grid()
ax[1].set_title("Angular Velocity of Coupled Pendulums")

plt.tight_layout()
plt.show()


plt.hist(theta[1,:])
plt.show()