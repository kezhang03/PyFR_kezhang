# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.linnavstokes.kernels.bcs.common'/>

<%pyfr:macro name='bc_rsolve_state' params='ul, nl, ur' externs='ploc, t'>
    fpdtype_t pb = ul[${nvars-1}];
    fpdtype_t rhob = ul[${bnvars}];
    ur[0] = rhob*ul[${bnvars - 1}]/pb;
% for i in range(ndims):
    ur[${i + 1}] = -ul[${i + 1}];
% endfor
    ur[${bnvars - 1}] = pb*ul[0]/rhob;
    ${pyfr.expand('bc_common_base_var_copy', 'ul', 'ur')};
</%pyfr:macro>

<%pyfr:macro name='bc_ldg_state' params='ul, nl, ur' externs='ploc, t'>
    fpdtype_t pb = ul[${nvars-1}];
    fpdtype_t rhob = ul[${bnvars}];
    ur[0] = rhob*ul[${bnvars - 1}]/pb;
% for i in range(ndims):
    ur[${i + 1}] = 0.0;
% endfor
    ur[${bnvars - 1}] = pb*ul[0]/rhob;
    ${pyfr.expand('bc_common_base_var_copy', 'ul', 'ur')};
</%pyfr:macro>

<%pyfr:alias name='bc_ldg_base_grad_state' func='bc_common_base_grad_copy'/>
<%pyfr:alias name='bc_ldg_grad_state' func='bc_grad_copy'/>