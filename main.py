from turbie_mod import read_turbie_parameters, build_matrices

filepath = "inputs/turbie_inputs/turbie_parameters.txt"

params = read_turbie_parameters(filepath)
M, C, K = build_matrices(params)
print(params)
print("M =")
print(M)

print("C =")
print(C)

print("K =")
print(K)
