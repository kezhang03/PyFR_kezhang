<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.lineuler.kernels.rsolvers.${rsolver}'/>
<%include file='pyfr.solvers.linnavstokes.kernels.flux'/>

<% tau = c['ldg-tau'] %>

<%pyfr:macro name='bc_common_flux_state' params='ul, gradul, nl, magnl'>
    // Viscous states
    fpdtype_t ur[${nvars}], gradur[${ndims}][${nvars}];
    // MODIFICATION FOR LINEAR SOLVER
    ${pyfr.expand('bc_ldg_base_grad_state', 'ul', 'nl', 'gradul', 'gradur')};
    // MODIFICATION FOR LINEAR SOLVER
    ${pyfr.expand('bc_ldg_state', 'ul', 'nl', 'ur')};
    ${pyfr.expand('bc_ldg_grad_state', 'ul', 'nl', 'gradul', 'gradur')};

    fpdtype_t fvr[${ndims}][${bnvars}] = {{0}};
    ${pyfr.expand('viscous_flux_add', 'ur', 'gradur', 'fvr')};

    // Inviscid (Riemann solve) state
    ${pyfr.expand('bc_rsolve_state', 'ul', 'nl', 'ur')};

    // Perform the Riemann solve
    fpdtype_t ficomm[${bnvars}], fvcomm;
    ${pyfr.expand('rsolve', 'ul', 'ur', 'nl', 'ficomm')};

% for i in range(bnvars):
    fvcomm = ${' + '.join(f'nl[{j}]*fvr[{j}][{i}]' for j in range(ndims))};
% if tau != 0.0:
    fvcomm += ${tau}*(ul[${i}] - ur[${i}]);
% endif

    ul[${i}] = magnl*(ficomm[${i}] + fvcomm);
% endfor
</%pyfr:macro>
