from turbie_mod import (read_turbie_parameters, build_matrices, read_ct_data,
interpolate_ct, read_wind_data, aerodynamic_force, build_A_matrix, dydt)

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

def get_wind_speed(file):
    """
    Description to be determined soon...
    """
    parts = file.name.split("_")
    filename_int_wind_speed = int(parts[1])

    return filename_int_wind_speed

# Input file paths
turbie_filepath = "inputs/turbie_inputs/turbie_parameters.txt"
ct_filepath = "inputs/turbie_inputs/CT.txt"
wind_filepath = "inputs/wind_files/wind_TI_0.1/wind_10_ms_TI_0.1.txt"
wind_basepath = Path("inputs/wind_files")

wind_folders = sorted(wind_basepath.iterdir())
for folder in wind_folders:

    wind_files = sorted(folder.glob("*.txt"), key=get_wind_speed)

    for wind_file in wind_files:
        print(wind_file)


# Read Turbie parameters and build model matrices
Turbie_params = read_turbie_parameters(turbie_filepath)
M, C, K = build_matrices(Turbie_params)

# Build first order system matrix
A = build_A_matrix(M, C, K)

# Read CT vs wind speeds lookup table
wind_speeds, ct_values = read_ct_data(ct_filepath)

# Test CT interpolation function
test_wind_speed = 12.4
test_ct = interpolate_ct(test_wind_speed, wind_speeds, ct_values)

# Read wind time series and calculate mean wind speed and TI
time, wind_speed = read_wind_data(wind_filepath)
mean_wind_speed = np.mean(wind_speed)
TI = np.std(wind_speed) / mean_wind_speed

# Calculating CT for mean wind speed of input file
ct = interpolate_ct(mean_wind_speed, wind_speeds, ct_values)

# Test aerodynamic force at initial conditions
test_blade_speed = 0
test_Faero = aerodynamic_force(wind_speed[0],test_blade_speed,ct,Turbie_params)

# Initial state vector
y0 = np.array([0, 0, 0, 0])

# Test dydt at initial conditions
test_dydt = dydt(0, y0, A, M, ct, Turbie_params, time, wind_speed)

# Use of IVP for solving ydot
t_span = (time[0], time[-1])
t_eval = time
solution = solve_ivp(dydt, t_span, y0, t_eval=t_eval,
                     args=(A, M, ct, Turbie_params, time, wind_speed),
                     method="RK45")
# Extract simulation results
x1 = solution.y[0]
x2 = solution.y[1]
xdot1 = solution.y[2]
xdot2 = solution.y[3]
blade_displacement = x1 - x2
tower_displacement = x2
# Print results for testing

print("A =")
print(A)

print("Interpolated CT =", test_ct)

print(wind_speeds)
print(ct_values)

print("Mean wind speed =", mean_wind_speed)
print("TI =", TI)

print("CT =", ct)

print("Initial wind speed =", wind_speed[0], "m/s")
print("Initial blade speed =", test_blade_speed, "m/s")
print("Initial aerodynamic force =", test_Faero, "N")

print("dy/dt at t = 0:")
print(test_dydt)

print("Simulation successful =", solution.success)
print("Solver message =", solution.message)
print("Number of output times =", len(solution.t))
print("Solution shape =", solution.y.shape)

print("Blade displacement range:",
      np.min(blade_displacement),
      np.max(blade_displacement))

print("Tower displacement range:",
      np.min(tower_displacement),
      np.max(tower_displacement))

# Plot CT lookup table for verification
plt.plot(wind_speeds, ct_values, "o-")
plt.xlabel("Wind speed [m/s]")
plt.ylabel("CT [-]")
plt.title("CT lookup table")
plt.grid()
plt.show()

# Plot wind speed time series for verification
plt.plot(time, wind_speed)
plt.axhline(mean_wind_speed, color="red", linestyle="--",
            label="Mean wind speed")
plt.xlabel("Time [s]")
plt.ylabel("Wind speed [m/s]")
plt.title("Wind speed time series")
plt.grid()
plt.legend()
plt.show()

plt.plot(solution.t, blade_displacement, label="Blade system")
plt.plot(solution.t, tower_displacement, label="Tower system")
plt.xlabel("Time [s]")
plt.ylabel("Displacement [m]")
plt.title("Turbie displacement response")
plt.grid()
plt.legend()
plt.show()
