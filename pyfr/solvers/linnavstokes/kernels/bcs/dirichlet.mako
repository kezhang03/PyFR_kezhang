# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.linnavstokes.kernels.bcs.common'/>

<%pyfr:macro name='bc_rsolve_state' params='ul, nl, ur' externs='ploc, t, q'>
    ur[0] = ${c['rho']};
% for i, v in enumerate('uvw'[:ndims]):
    ur[${i+1}] = (${c['rho']})*(${c[v]});
% endfor
    ur[${bnvars - 1}] = ${c['p']};

// copy the baseflow
    ${pyfr.expand('bc_common_base_var_copy', 'ul', 'ur')};
</%pyfr:macro>

<%pyfr:alias name='bc_ldg_state' func='bc_rsolve_state'/>

<%pyfr:alias name='bc_ldg_grad_state' func='bc_common_grad_zero'/>
