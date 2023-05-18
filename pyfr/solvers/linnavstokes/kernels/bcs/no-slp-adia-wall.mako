# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.linnavstokes.kernels.bcs.common'/>

<%pyfr:macro name='bc_rsolve_state' params='ul, nl, ur' externs='ploc, t'>
    // fpdtype_t invgmo = ${1/(c['gamma']-1)};
    // fpdtype_t rhob = ul[${bnvars}];
    // fpdtype_t pb = ul[${nvars-1}];
    // fpdtype_t Tb = invgmo*pb/rhob;
    // fpdtype_t p = ul[${bnvars-1}];
    // fpdtype_t rho = ul[0];
    // fpdtype_t T = (p/pb-rho/rhob)*Tb;

    // ur[0] = ul[0]+rhob*T/Tb;
    ur[0] = 0;
% for i in range(ndims):
    ur[${i + 1}] = -ul[${i + 1}];
% endfor
    // ur[${bnvars - 1}] = ul[${bnvars - 1}]-pb*T/Tb;
    ur[${bnvars - 1}] = 0;
    ${pyfr.expand('bc_common_base_var_copy', 'ul', 'ur')};
</%pyfr:macro>

<%pyfr:macro name='bc_ldg_state' params='ul, nl, ur' externs='ploc, t'>
    // fpdtype_t invgmo = ${1/(c['gamma']-1)};
    // fpdtype_t rhob = ul[${bnvars}];
    // fpdtype_t pb = ul[${nvars-1}];
    // fpdtype_t Tb = invgmo*pb/rhob;
    // fpdtype_t p = ul[${bnvars-1}];
    // fpdtype_t rho = ul[0];
    // fpdtype_t T = (p/pb-rho/rhob)*Tb;

    // ur[0] = ul[0]+rhob*T/Tb;
    ur[0] = 0;
% for i in range(ndims):
    ur[${i + 1}] = 0.0;
% endfor
    // ur[${bnvars - 1}] = ul[${bnvars - 1}]-pb*T/Tb;
    ur[${bnvars - 1}] = 0;
    ${pyfr.expand('bc_common_base_var_copy', 'ul', 'ur')};
</%pyfr:macro>

<%pyfr:macro name='bc_ldg_grad_state' params='ul, nl, grad_ul, grad_ur'>
    ${pyfr.expand('bc_common_grad_copy', 'ul', 'nl', 'grad_ul', 'grad_ur')};
</%pyfr:macro>

<%pyfr:alias name='bc_ldg_base_grad_state' func='bc_common_base_grad_copy'/>
