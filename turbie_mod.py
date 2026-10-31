import numpy as np

def read_turbie_parameters(filepath):
    params = {}
    with open(filepath, "r") as file:
        lines = file.readlines()

    for line in lines:
        if not line.startswith("#"):
            parts = line.split("#")
            value = float(parts[0])
            name = parts[1].strip()
            params[name] = value

    return params

def build_matrices(params):
    m1 = 3 * params["mb [kg]"]
    m2 = params["mn [kg]"] + params["mh [kg]"] + params["mt [kg]"]
    M = np.array([[m1,0],[0,m2]])
    c1 = params["c1 [N/(m/s)]"]
    c2 = params["c2 [N/(m/s)]"]
    C = np.array([[c1, -c1], [-c1, c1 + c2]])
    k1 = params["k1 [N/(m/s)^2]"]
    k2 = params["k2 [N/(m/s)^2]"]
    K = np.array([[k1, -k1], [-k1, k1 + k2]])

    return M, C ,K
