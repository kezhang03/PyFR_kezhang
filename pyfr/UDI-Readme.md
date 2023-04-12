UDI: User Defined Input
Version = 1.0 
Date = 2023/04/12
Function:
UDI interpolates the scatter data from other sources on the inlet boundary for the Dirichlet boundary condition of linearised solver, using the scipy.interpolate.CubicSpline

Usage:
The user-defined data is called 'q_{variable name}_{real/imag}' in the boundary condition of .ini file. 
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

Requirement:
Input Data:
reads .mat data from MATLAB, .mat version <= v7.
only supports input of 4(2D) and 5(3D) variables, no single variable can be missing, otherwise may crash.
currently only supports complex input.
the input data of this version is ordered from large value to small value, so the data is flipped.

Input Format (.ini file):
[soln-UDI]
UDI-dirname = 
UDI-filenames = rho_q, u_q, v_q, p_q # the variables must have the order rho, u, v, p, w
UDI-format = .mat
UDI-coords-filename = # coordinates for the input raw data