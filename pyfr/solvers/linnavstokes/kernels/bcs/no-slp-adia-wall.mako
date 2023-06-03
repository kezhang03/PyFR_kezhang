# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.linnavstokes.kernels.bcs.common'/>

<%pyfr:macro name='bc_rsolve_state' params='ul, nl, ur' externs='ploc, t'>
    ur[0] = ul[0];
% for i in range(ndims):
    // ur[${i + 1}] = -ul[${i + 1}];
    ur[${i + 1}] = 0.0;
% endfor
    ur[${bnvars - 1}] = ul[${bnvars - 1}];
    ${pyfr.expand('bc_common_base_var_copy', 'ul', 'ur')};
</%pyfr:macro>

<%pyfr:macro name='bc_ldg_state' params='ul, nl, ur' externs='ploc, t'>
    ur[0] = ul[0];
% for i in range(ndims):
    ur[${i + 1}] = 0.0;
% endfor
    ur[${bnvars - 1}] = ul[${bnvars - 1}];
    ${pyfr.expand('bc_common_base_var_copy', 'ul', 'ur')};
</%pyfr:macro>

<%pyfr:macro name='bc_ldg_grad_state' params='ul, nl, grad_ul, grad_ur'>
    // 1/rho_base
    fpdtype_t invrhob = 1.0/ul[${bnvars}];
    // 1/(gamma-1)
    fpdtype_t invgmo = ${1/(c['gamma']-1)};
% if ndims == 2:
    //baseflow variables
    fpdtype_t rhob = ul[4];
    fpdtype_t pb = ul[7];

    //baseflow derivatives
    fpdtype_t rhob_x = grad_ul[0][4];
    fpdtype_t rhob_y = grad_ul[1][4];
    fpdtype_t pb_x = grad_ul[0][7];
    fpdtype_t pb_y = grad_ul[1][7];

    //perturbation variables
    fpdtype_t rho = ul[0];
    fpdtype_t p = ul[3];

    //perturbation derivatives
    fpdtype_t rho_x = grad_ul[0][0];
    fpdtype_t rho_y = grad_ul[1][0];

    fpdtype_t p_x = grad_ul[0][3];
    fpdtype_t p_y = grad_ul[1][3];

    // Compute temperature derivatives (c_v*dT/d[x,y,z])
    fpdtype_t Tb = invgmo*pb/rhob;
    fpdtype_t Tb_x = invgmo*(invrhob*pb_x-invrhob*invrhob*pb*rhob_x);
    fpdtype_t Tb_y = invgmo*(invrhob*pb_y-invrhob*invrhob*pb*rhob_y);

    fpdtype_t T_x = Tb_x * (p/pb-rho/rhob) + Tb * (p_x/pb - p/pb/pb*pb_x - rho_x*invrhob + rho*invrhob*invrhob*rhob_x);
    fpdtype_t T_y = Tb_y * (p/pb-rho/rhob) + Tb * (p_y/pb - p/pb/pb*pb_y - rho_y*invrhob + rho*invrhob*invrhob*rhob_y);

    // Copy all fluid-side gradients across to wall-side gradients
    // We have already copied the baseflow gradient before
    ${pyfr.expand('bc_common_grad_copy', 'ul', 'nl', 'grad_ul', 'grad_ur')};

    // Correct copied across in-fluid temp gradients to in-wall gradients
    grad_ur[0][3] -= rhob/Tb*(nl[0]*nl[0]*(T_x-Tb_x*(p/pb-rho/rhob)) + nl[0]*nl[1]*(T_y-Tb_y*(p/pb-rho/rhob)));
    grad_ur[1][3] -= rhob/Tb*(nl[1]*nl[0]*(T_x-Tb_x*(p/pb-rho/rhob)) + nl[1]*nl[1]*(T_y-Tb_y*(p/pb-rho/rhob)));

% elif ndims == 3:
    //baseflow variables
    fpdtype_t rhob = u_in[5];
    fpdtype_t pb = u_in[9];

    //baseflow derivatives
    fpdtype_t rhob_x = grad_uin[0][5];
    fpdtype_t rhob_y = grad_uin[1][5];
    fpdtype_t rhob_z = grad_uin[2][5];

    fpdtype_t pb_x = grad_uin[0][9];
    fpdtype_t pb_y = grad_uin[1][9];
    fpdtype_t pb_z = grad_uin[2][9];

    //perturbation variables
    fpdtype_t rho = u_in[0];
    fpdtype_t p = u_in[4];

    //perturbation derivatives
    fpdtype_t rho_x = grad_ul[0][0];
    fpdtype_t rho_y = grad_ul[1][0];
    fpdtype_t rho_z = grad_ul[2][0];

    fpdtype_t p_x = grad_ul[0][4];
    fpdtype_t p_y = grad_ul[1][4];
    fpdtype_t p_z = grad_ul[2][4];

    // Compute temperature derivatives (c_v*dT/d[x,y,z])
    fpdtype_t Tb = invgmo*pb/rhob;
    fpdtype_t Tb_x = invgmo*(invrhob*pb_x-pb*invrhob*invrhob*rhob_x);
    fpdtype_t Tb_y = invgmo*(invrhob*pb_y-pb*invrhob*invrhob*rhob_y);
    fpdtype_t Tb_z = invgmo*(invrhob*pb_z-pb*invrhob*invrhob*rhob_z);

    fpdtype_t T_x = Tb_x * (p/pb-rho/rhob) + Tb * (p_x/pb - p/pb/pb*pb_x - rho_x*invrhob + rho*invrhob*invrhob*rhob_x);
    fpdtype_t T_y = Tb_y * (p/pb-rho/rhob) + Tb * (p_y/pb - p/pb/pb*pb_y - rho_y*invrhob + rho*invrhob*invrhob*rhob_y);
    fpdtype_t T_z = Tb_z * (p/pb-rho/rhob) + Tb * (p_z/pb - p/pb/pb*pb_z - rho_z*invrhob + rho*invrhob*invrhob*rhob_z);

    // Copy all fluid-side gradients across to wall-side gradients
    ${pyfr.expand('bc_common_grad_copy', 'ul', 'nl', 'grad_ul', 'grad_ur')};

    // Correct copied across in-fluid temp gradients to in-wall gradients
    grad_ur[0][4] -= rhob/Tb*(nl[0]*nl[0]*(T_x-Tb_x*(p/pb-rho/rhob)) + nl[0]*nl[1]*(T_y-Tb_y*(p/pb-rho/rhob)));
    grad_ur[1][4] -= rhob/Tb*(nl[1]*nl[0]*(T_x-Tb_x*(p/pb-rho/rhob)) + nl[1]*nl[1]*(T_y-Tb_y*(p/pb-rho/rhob)));
    grad_ur[2][4] -= rhob/Tb*(nl[2]*nl[0]*(T_x-Tb_x*(p/pb-rho/rhob)) + nl[2]*nl[1]*(T_y-Tb_y*(p/pb-rho/rhob)));
    grad_ur[0][4] -= rhob/Tb*(nl[0]*nl[2]*(T_z-Tb_z*(p/pb-rho/rhob)));
    grad_ur[1][4] -= rhob/Tb*(nl[1]*nl[2]*(T_z-Tb_z*(p/pb-rho/rhob)));
    grad_ur[2][4] -= rhob/Tb*(nl[2]*nl[2]*(T_z-Tb_z*(p/pb-rho/rhob)));

% endif
</%pyfr:macro>

<%pyfr:alias name='bc_ldg_base_grad_state' func='bc_common_base_grad_copy'/>