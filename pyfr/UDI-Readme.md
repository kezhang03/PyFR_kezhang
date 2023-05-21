UDI: User Defined Input
Version = 1.2
Date = 2023/05/21
Function:
UDI interpolates the scatter data from other sources on the inlet boundary, using the scipy.interpolate.interp1d and constant extrapolation based on the last value.
Supported BC:
Non-liner solver: char-riem-inv
Linear Solver: dirichlet(LNS)
Usage:
The user-defined data is called 'q_{variable name}_{real/imag}' in the boundary condition of .ini file for complex input.
The 2nd input should be q2_rho or q2_rho_i.

Update:
V 1.2: now supports 2 input simultaneously, one input is same as before so old .ini file can be used without any modification

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

[soln-bcs-inlet]
type = char-riem-inv
rho = q2_rho + A*(q_rho_r*cos(w*(t-82800)) + q_rho_i*sin(w*(t-82800)))
u = q2_u + A*(q_u_r*cos(w*(t-82800)) + q_u_i*sin(w*(t-82800)))
v = q2_v + A*(q_v_r*cos(w*(t-82800)) + q_v_i*sin(w*(t-82800)))
p = P + A*(q_p_r*cos(w*(t-82800)) + q_p_i*sin(w*(t-82800)))

Requirement:
Input Data:
reads .mat data from MATLAB, .mat version <= v7.
only supports input of 4(2D) and 5(3D) variables, no single variable can be missing, otherwise may crash.
currently only supports complex input.
input data must be ordered from small value to large 
the variable name saved in the .mat file should be same as the file name

Input Format (.ini file):
[soln-UDI]
UDI-datatype = real / complex
UDI-dirname = 
UDI-filenames = rho_q, u_q, v_q, p_q # the variables must have the order rho, u, v, p, w
UDI-format = .mat
UDI-coords-filename = # coordinates for the input raw data

UDI-dirname2 = .
UDI-filenames2 = rho_inlet, u_inlet, v_inlet, p_inlet
UDI-coords-filename2 = y_inlet
UDI-format2 = .mat
UDI-datatype2 = real