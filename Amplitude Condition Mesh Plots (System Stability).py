import numpy as np
import matplotlib.pyplot as plt
import datetime
import time
import os

# Parameters
g = 9.81  # Acceleration due to gravity (m/s^2)
l = 0.25  # Length of the pendulum (m)
a0 = 1  # Coefficient a0
b = 0     # Coefficient b
m = 1.0   # Mass of the pendulum bob (kg)

# Time span
t_span = np.linspace(0, 3, 1500)

# Function to apply the constraint
def apply_constraint(theta, theta_dot):
    if theta >= np.pi:
        theta_dot = -theta_dot  
        # Reverse theta_dot when theta reaches pi
    elif theta <= 0:
        theta_dot = -theta_dot 
        # Reverse theta_dot when theta reaches -pi
    return theta_dot

# Function to solve the differential equation for a given omega
def solve_pendulum(omega, theta0, theta_dot0):
    # Initialize arrays to store results
    theta_values = np.zeros_like(t_span)
    theta_dot_values = np.zeros_like(t_span)
    
    # Iterate through time steps 
    theta = theta0
    theta_dot = theta_dot0
    
    for i, t in enumerate(t_span):
        # Apply the constraint
        theta_dot = apply_constraint(theta, theta_dot)
        
        # Update theta_dot and theta 
        theta_double_dot = (-(g / l) + a0 * omega**2 
        * np.cos(omega * t)) * np.sin(theta)
        + (b / m) * theta_dot
        theta_dot += theta_double_dot * (t_span[1] - t_span[0])
        theta += theta_dot * (t_span[1] - t_span[0])
        
        # Store the results
        theta_values[i] = theta
        theta_dot_values[i] = theta_dot
    
    return np.max(theta_values)

# Range of angular frequencies (omega)
omega_range = np.linspace(1,10, 12)  # Use three different 
values of omega

# Range of initial conditions
theta0_range = np.linspace(0, np.pi, 100)
theta_dot0_range = np.linspace(-2, 2, 100)

# Create mesh grids for initial conditions
Theta0, ThetaDot0 = np.meshgrid(theta0_range, theta_dot0_range)

# Create arrays to store amplitude categories
amplitude_categories = np.zeros((len(theta0_range), 
len(theta_dot0_range), len(omega_range)))

# Start a timer before the loop
start_time = time.time()

# Calculate amplitude categories 
for i, omega in enumerate(omega_range):
    for j in range(len(theta0_range)):
        for k in range(len(theta_dot0_range)):
            theta0 = Theta0[j, k]
            theta_dot0 = ThetaDot0[j, k]
            max_amplitude = solve_pendulum(omega=omega, 
            theta0=theta0, theta_dot0=theta_dot0)
            
            if max_amplitude < np.pi:
                amplitude_categories[j, k, i] = 1
            else:
                amplitude_categories[j, k, i] = 2

    # Calculate the remaining processing time
    elapsed_time = time.time() - start_time
    average_time_per_plot = elapsed_time / (i + 1) 
    if i > 0 else elapsed_time
    remaining_plots = len(omega_range) - (i + 1)
    remaining_time = average_time_per_plot * remaining_plots

    # Print remaining processing time
    print(f"Remaining processing time: {remaining_time:.2f} seconds")

    # Create 2D mesh plot for the current omega
    plt.figure(figsize=(8, 6))
    plt.pcolormesh(theta0_range, theta_dot0_range, 
    amplitude_categories[:, :, i], cmap='coolwarm', 
    shading='auto', vmin=1, vmax=2)
    plt.colorbar(label='Amplitude Category', 
    ticks=[1.5, 2.5], boundaries=[1, 2, 3], format='%1i')
    plt.xlabel('Initial Angle (theta0)')
    plt.ylabel('Initial Angular Velocity (theta_dot0)')
    plt.title(f'Omega = {omega}')
    plt.grid(True)
