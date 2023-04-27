UDI: User Defined Input
Version = 1.1
Date = 2023/04/27
Function:
UDI interpolates the scatter data from other sources on the inlet boundary, using the scipy.interpolate.CubicSpline
Now supports char-riem-inv BC and dirichlet BC for linear solver
Usage:
The user-defined data is called 'q_{variable name}_{real/imag}' in the boundary condition of .ini file for complex input.

'q_{variable name}' for real input  
variable name = rho, u ,v, p, w
real/imag = r, i
It is linked with coordinates data on the inlet boundary, sorted as the same order as the inlet coordinates.

Usage example:
[soln-bcs-inlet]
type = dirichlet
rho = A_rho*(q_rho_r*cos(w*t) - q_rho_i*sin(w*t))
u = A_u*(q_u_r*cos(w*t) - q_u_i*sin(w*t))
v = A_v*(q_v_r*cos(w*t) - q_v_i*sin(w*t))
p = A_p*(q_p_r*cos(w*t) - q_p_i*sin(w*t))

[soln-bcs-inlet]
type = char-riem-inv
rho = q_rho
u = q_u
v = q_v
p = P #(use a constant)

Requirement:
Input Data:
reads .mat data from MATLAB, .mat version <= v7.
only supports input of 4(2D) and 5(3D) variables, no single variable can be missing, otherwise may crash.
currently only supports complex input.
the input data of this version is ordered from large value to small value, so the data is flipped.
the variable name saved in the .mat file should be same as the file name

Input Format (.ini file):
[soln-UDI]
UDI-datatype = real / complex
UDI-dirname = 
UDI-filenames = rho_q, u_q, v_q, p_q # the variables must have the order rho, u, v, p, w
UDI-format = .mat
UDI-coords-filename = # coordinates for the input raw data