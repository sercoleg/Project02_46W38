from turbie_mod import (read_turbie_parameters, build_matrices, read_ct_data,
interpolate_ct, read_wind_data, aerodynamic_force, build_A_matrix, dydt)

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

def get_wind_speed(file):
    """
    Extracts the nominal wind speed from wind input filename.

    Parameter
    ---------
    file : Path
            Path of the wind input file.

    Returns
    -------
    filename_int_wind_speed : int
            Nominal wind speed from filename [m/s].
    """
    parts = file.name.split("_")
    filename_int_wind_speed = int(parts[1])

    return filename_int_wind_speed

# Input file paths
turbie_filepath = "inputs/turbie_inputs/turbie_parameters.txt"
ct_filepath = "inputs/turbie_inputs/CT.txt"
wind_basepath = Path("inputs/wind_files")

# Read Turbie parameters and build model matrices
Turbie_params = read_turbie_parameters(turbie_filepath)
M, C, K = build_matrices(Turbie_params)
# Build first order system matrix
A = build_A_matrix(M, C, K)

# Read CT vs wind speeds lookup table
wind_speeds, ct_values = read_ct_data(ct_filepath)

# Initial state vector for IVP solver
y0 = np.array([0, 0, 0, 0])

# Loading wind speed data and solving dynamic system for each wind speed file
wind_folders = sorted(wind_basepath.iterdir())
for folder in wind_folders:

    wind_files = sorted(folder.glob("*.txt"), key=get_wind_speed)
    TI_name = folder.name.replace("wind_", "")
    output_folder = Path("outputs") / TI_name
    output_folder.mkdir(parents=True, exist_ok=True)
    # Empty list for storing statistics of current TI category
    statistics = []

    for wind_file in wind_files:
            # Read wind time series
            time, wind_speed = read_wind_data(wind_file)

            # Calculate mean wind speed
            mean_wind_speed = np.mean(wind_speed)

            # Calculate CT for current simulation
            ct = interpolate_ct(mean_wind_speed, wind_speeds, ct_values)

            # Define simulation time
            t_span = (time[0],time[-1])

            # Solve dynamic system
            solution = solve_ivp(dydt, t_span, y0, t_eval=time,
                                 args=(A, M, ct, Turbie_params, time,
                                       wind_speed), method="RK45")

            # Extract simulation results
            x1 = solution.y[0]
            x2 = solution.y[1]
            xdot1 = solution.y[2]
            xdot2 = solution.y[3]

            # Calculate relative displacements
            blade_displacement = x1 - x2
            tower_displacement = x2

            # Temporary check
            print(folder.name, wind_file.name, mean_wind_speed, solution.success)
            results = np.column_stack((solution.t, wind_speed, x1, x2,
                                       xdot1, xdot2, blade_displacement,
                                       tower_displacement))

            nominal_wind_speed = get_wind_speed(wind_file)

            output_filename = f"results_{nominal_wind_speed}_ms.txt"
            output_filepath = output_folder / output_filename

            np.savetxt(output_filepath, results, fmt="%.6f", header=(
                 "Time[s] Wind[m/s] x1[m] x2[m] "
                 "xdot1[m/s] xdot2[m/s] "
                 "Blade_rel_disp[m] Tower_rel_disp[m]"))


            blade_mean = np.mean(blade_displacement)
            blade_std = np.std(blade_displacement)

            tower_mean = np.mean(tower_displacement)
            tower_std = np.std(tower_displacement)

            # Store statistics for current wind speed case
            statistics.append([nominal_wind_speed, mean_wind_speed, blade_mean,
                               blade_std, tower_mean, tower_std])

    statistics = np.array(statistics)

    statistics_filename = f"statistics_{TI_name}.txt"
    statistics_filepath = output_folder / statistics_filename

    np.savetxt(statistics_filepath, statistics, fmt="%.6f", header=(
         "Nominal_Wind[m/s] Mean_Wind[m/s] "
         "Blade_Mean[m] Blade_Std[m] "
         "Tower_Mean[m] Tower_Std[m]"))

#-----------------------------------------------------------------------------
# Create folder for output figures and saving figures before plotting
figures_folder = Path("outputs") / "figures"
figures_folder.mkdir(parents=True, exist_ok=True)
#-----------------------------------------------------------------------------
# Plot time series for selected wind case : TI = 10%; 10 m/s.
selected_results = np.loadtxt("outputs/TI_0.1/results_10_ms.txt", skiprows=1)
plot_time = selected_results[:, 0]
plot_wind = selected_results[:, 1]
plot_blade = selected_results[:, 6]
plot_tower = selected_results[:, 7]

fig, ax = plt.subplots(2, 1, sharex=True)

ax[0].plot(plot_time, plot_wind)
ax[0].set_ylabel("Wind speed [m/s]")
ax[0].set_title("10 m/s wind case - TI 0.10")
ax[0].grid()

ax[1].plot(plot_time, plot_blade, label="Blade system")
ax[1].plot(plot_time, plot_tower, label="Tower system")
ax[1].set_xlabel("Time [s]")
ax[1].set_ylabel("Displacement [m]")
ax[1].grid()
ax[1].legend()

plt.tight_layout()
fig.savefig(figures_folder / "time_series_10_ms_TI_0.1.png", dpi=300)
plt.show()
#-----------------------------------------------------------------------------
# Load displacement statistics for each TI category
stats_005 = np.loadtxt("outputs/TI_0.05/statistics_TI_0.05.txt", skiprows=1)
stats_01 = np.loadtxt("outputs/TI_0.1/statistics_TI_0.1.txt", skiprows=1)
stats_015 = np.loadtxt("outputs/TI_0.15/statistics_TI_0.15.txt", skiprows=1)

# Plot displacement statistics for all TI categories
fig, ax = plt.subplots(1, 2, figsize=(12, 5))


# Mean displacement - TI 0.05
ax[0].plot(stats_005[:, 0], stats_005[:, 2], "o-",
           color="blue", label="Blade - TI 0.05")
ax[0].plot(stats_005[:, 0], stats_005[:, 4], "o--",
           color="blue", label="Tower - TI 0.05")

# Mean displacement - TI 0.10
ax[0].plot(stats_01[:, 0], stats_01[:, 2], "o-",
           color="green", label="Blade - TI 0.10")
ax[0].plot(stats_01[:, 0], stats_01[:, 4], "o--",
           color="green", label="Tower - TI 0.10")

# Mean displacement - TI 0.15
ax[0].plot(stats_015[:, 0], stats_015[:, 2], "o-",
           color="red", label="Blade - TI 0.15")
ax[0].plot(stats_015[:, 0], stats_015[:, 4], "o--",
           color="red", label="Tower - TI 0.15")

ax[0].set_title("Mean displacement")
ax[0].set_xlabel("Wind speed [m/s]")
ax[0].set_ylabel("Mean displacement [m]")
ax[0].grid()
ax[0].legend()


# Standard deviation - TI 0.05
ax[1].plot(stats_005[:, 0], stats_005[:, 3], "o-",
           color="blue", label="Blade - TI 0.05")
ax[1].plot(stats_005[:, 0], stats_005[:, 5], "o--",
           color="blue", label="Tower - TI 0.05")

# Standard deviation - TI 0.10
ax[1].plot(stats_01[:, 0], stats_01[:, 3], "o-",
           color="green", label="Blade - TI 0.10")
ax[1].plot(stats_01[:, 0], stats_01[:, 5], "o--",
           color="green", label="Tower - TI 0.10")

# Standard deviation - TI 0.15
ax[1].plot(stats_015[:, 0], stats_015[:, 3], "o-",
           color="red", label="Blade - TI 0.15")
ax[1].plot(stats_015[:, 0], stats_015[:, 5], "o--",
           color="red", label="Tower - TI 0.15")

ax[1].set_title("Displacement standard deviation")
ax[1].set_xlabel("Wind speed [m/s]")
ax[1].set_ylabel("Standard deviation [m]")
ax[1].grid()
ax[1].legend()


fig.suptitle("Turbie displacement statistics")
plt.tight_layout()
fig.savefig(figures_folder / "displacement_statistics.png",  dpi=300)
plt.show()
#----------------------------------------------------------------------

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
