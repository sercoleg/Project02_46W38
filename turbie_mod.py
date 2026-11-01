import numpy as np

def read_turbie_parameters(filepath):
    """
    Reads the Turbie model parameters from a text file.

    Parameter
    ---------
    filepath : str
        Path to the Turbie parameter file.

    Returns
    -------
    Turbie_params : dict
        Turbie parameters with parameter names and units as keys.
    """
    Turbie_params = {}
    with open(filepath, "r") as file:
        lines = file.readlines()

    for line in lines:
        if not line.startswith("#"):
            parts = line.split("#")
            value = float(parts[0])
            name = parts[1].strip()
            Turbie_params[name] = value

    return Turbie_params

def build_matrices(Turbie_params):
    """
    Builds the mass M, damping C and stiffness K matrices of the Turbie model.

    Parameter
    ---------
    Turbie_params : dict
        Turbie model parameters.

    Returns
    -------
    M : numpy.ndarray
        Mass matrix.

    C : numpy.ndarray
        Damping matrix.

    K : numpy.ndarray
        Stiffness matrix.
    """
    # Calculate equivalent masses and build mass matrix
    m1 = 3 * Turbie_params["mb [kg]"]
    m2 = Turbie_params["mn [kg]"] + Turbie_params["mh [kg]"] + Turbie_params["mt [kg]"]
    M = np.array([[m1,0],[0,m2]])
    # Build damping matrix
    c1 = Turbie_params["c1 [N/(m/s)]"]
    c2 = Turbie_params["c2 [N/(m/s)]"]
    C = np.array([[c1, -c1], [-c1, c1 + c2]])
    # Build stiffness matrix
    k1 = Turbie_params["k1 [N/(m/s)^2]"]
    k2 = Turbie_params["k2 [N/(m/s)^2]"]
    K = np.array([[k1, -k1], [-k1, k1 + k2]])

    return M, C ,K

def read_ct_data(filepath):
    """
    Reads wind speed and thrust coefficient data from a text file.

    Parameter
    ---------
    filepath : str
        Path to the CT input file.

    Returns
    -------
    wind_speeds : list
        Wind speeds [m/s].

    ct_values : list
        Thrust coefficient values [-].
    """
    wind_speeds = []
    ct_values = []

    # Read all lines from the CT input file
    with open(filepath,"r") as file:
        lines = file.readlines()

    # Extract wind speed and CT values, ignoring comment lines
    for line in lines:
        if not line.startswith("#"):
            parts = line.split()
            wind_speeds.append(float(parts[0]))
            ct_values.append(float(parts[1]))

    return wind_speeds, ct_values

def interpolate_ct(mean_wind_speed, wind_speeds, ct_values):
    """
    Interpolates the thrust coefficient CT for a given mean wind speed.

    Parameters
    ----------
    mean_wind_speed : float
        Mean wind speed of the simulation [m/s].

    wind_speeds : list
        Wind speeds from the CT lookup table [m/s].

    ct_values : list
        Thrust coefficient values from the CT lookup table [-].

    Returns
    -------
    ct : float
        Interpolated thrust coefficient [-].
    """

    ct = np.interp(mean_wind_speed, wind_speeds, ct_values)

    return ct

def read_wind_data(filepath):
    """
    Reads time and wind speed data from a wind file.

    Parameter
    ---------
    filepath : str
        Path to the wind input file.

    Returns
    -------
    time : numpy.ndarray
        Time [s].

    wind_speed : numpy.ndarray
        Wind speed [m/s].
    """

    data = np.loadtxt(filepath, skiprows=1)

    time = data[:, 0]
    wind_speed = data[:, 1]

    return time, wind_speed

def aerodynamic_force(wind_speed, blade_speed, ct, Turbie_params):
    """
    Calculates the aerodynamic force acting on the blades.

    Parameters
    ----------
    wind_speed : float
        Instant wind speed acting on the rotor blades [m/s].

    blade_speed : float
        Instant blades speed in the wind direction [m/s].

    ct : float
        Thrust coefficient constant for the entire wind speed data set [-].

    Turbie_params : dict
        Set of parameters of Turbie model

    Returns
    -------
    Faero : float
        Aerodynamic force on the blade [N]
    """
    # Extracting constant parameters from Turbie model.
    rho = Turbie_params["rho [kg/m3]"]
    Dr = Turbie_params["Dr [m]"]

    # Calculating swept area of the wind turbine rotor
    A = np.pi*Dr**2/4

    # Auxiliar variable u(t) - xdot1(t) for compacting Faero equation
    speed_diff = wind_speed - blade_speed

    # Calculating aerodynamic force on the rotor blades.
    Faero = 0.5*rho*ct*A*speed_diff*abs(speed_diff)

    return Faero
