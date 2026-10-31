from turbie_mod import read_turbie_parameters, build_matrices, read_ct_data, interpolate_ct

turbie_filepath = "inputs/turbie_inputs/turbie_parameters.txt"
ct_filepath = "inputs/turbie_inputs/CT.txt"

Turbie_params = read_turbie_parameters(turbie_filepath)
M, C, K = build_matrices(Turbie_params)

wind_speeds, ct_values = read_ct_data(ct_filepath)

test_wind_speed = 12.4
test_ct = interpolate_ct(test_wind_speed, wind_speeds, ct_values)

print("Interpolated CT =", test_ct)

print(wind_speeds)
print(ct_values)
