from turbie_mod import (read_turbie_parameters, build_matrices, read_ct_data,
interpolate_ct, read_wind_data, aerodynamic_force)

import numpy as np
import matplotlib.pyplot as plt

# Input file paths
turbie_filepath = "inputs/turbie_inputs/turbie_parameters.txt"
ct_filepath = "inputs/turbie_inputs/CT.txt"
wind_filepath = "inputs/wind_files/wind_TI_0.1/wind_10_ms_TI_0.1.txt"

# Read Turbie parameters and build model matrices
Turbie_params = read_turbie_parameters(turbie_filepath)
M, C, K = build_matrices(Turbie_params)

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

# Print results for testing
print("Interpolated CT =", test_ct)

print(wind_speeds)
print(ct_values)

print("Mean wind speed =", mean_wind_speed)
print("TI =", TI)

print("CT =", ct)

print("Initial wind speed =", wind_speed[0], "m/s")
print("Initial blade speed =", test_blade_speed, "m/s")
print("Initial aerodynamic force =", test_Faero, "N")

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
