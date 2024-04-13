import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Parameters
g = 9.81  # Acceleration due to gravity (m/s^2)
l = 0.25   # Length of the pendulum (m)
a0 = 0.5  # Coefficient a0
omega = 10  # Angular frequency (rad/s)
b = 0   # Coefficient b
m = 1.0   # Mass of the pendulum bob (kg) (does not affect equation of motion, set 1 for convenience)

# Initial conditions
theta0 = 1.5       # Initial angle (rad)
theta_dot0 = 0.0  # Initial angular velocity (rad/s)

# Time span
t_span = np.linspace(0, 10, 300000)

# Initialize arrays to store results
theta_values = np.zeros_like(t_span)
theta_dot_values = np.zeros_like(t_span)
energy_values = np.zeros_like(t_span)

# Define the differential equation
def diff_eq(t, y):
    theta, theta_dot = y
    theta_double_dot = (g / l - a0 * omega**2 * np.cos(omega * t)
    * np.sin(theta)
    + b / m * theta_dot) / (-l)
    return [theta_dot, theta_double_dot]

# Apply boundary constraint
def apply_constraint(theta, theta_dot):
    if theta <= 0:
        theta_dot = -theta_dot
        theta = 0
    elif theta >= np.pi:
        theta_dot = -theta_dot
        theta = np.pi
    return theta, theta_dot

# Define the energy function
def energy_function(theta, theta_dot, t):
    kinetic_energy = 0.5 * m * (l**2 * theta_dot**2 +a0**2 
    * omega**2 * (np.sin(omega * t))**2  - 2* theta_dot * l * a0 * omega 
    * np.sin(omega * t)) 
    potential_energy =  m * g * l * (1 - np.cos(theta))
    return kinetic_energy + potential_energy

# Solve the differential equation and track energy
y0 = [theta0, theta_dot0]

for i, t in enumerate(t_span):
    sol = solve_ivp(diff_eq, [t, t + (t_span[1] - t_span[0])], 
    y0, t_eval=[t + (t_span[1] - t_span[0])])
    y0 = sol.y[:, -1]
    
    # Apply boundary constraint
    y0[0], y0[1] = apply_constraint(y0[0], y0[1])
    
    # Store the results
    theta_values[i] = y0[0]
    theta_dot_values[i] = y0[1]
    
    # Calculate and store the energy
    energy = energy_function(y0[0], y0[1], t)
    energy_values[i] = energy

# Plot energy as a function of time
plt.figure(figsize=(10, 4))
plt.plot(t_span, energy_values)
plt.xlabel('Time (s)')
plt.ylabel('Energy (E)')
plt.title('Energy E as a Function of Time')
plt.grid()
plt.show()

#Plot energy as a function of theta
plt.figure(figsize=(10, 4))
plt.plot(theta_values, energy_values)
plt.xlabel('theta (rad)')
plt.ylabel('Energy (E)')
plt.title('Energy E as a Function of theta')
plt.grid()
plt.show()
