# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.lineuler.kernels.flux'/>

<%pyfr:macro name='rsolve' params='ul, ur, n, nf'>
    // Compute the left and right fluxes + velocities and pressures
    fpdtype_t fl[${ndims}][${bnvars}], fr[${ndims}][${bnvars}];
    fpdtype_t vl[${ndims}], vr[${ndims}];
    fpdtype_t pl, pr;
    fpdtype_t pbl = ul[${nvars-1}];
    fpdtype_t pbr = ur[${nvars-1}];

    ${pyfr.expand('inviscid_flux', 'ul', 'fl', 'pl', 'vl')};
    ${pyfr.expand('inviscid_flux', 'ur', 'fr', 'pr', 'vr')};

    // Sum the left and right velocities and take the normal
    fpdtype_t nv = ${pyfr.dot('n[{i}]', 'vl[{i}] + vr[{i}]', i=ndims)};

    // Estimate the maximum wave speed / 2
    // 2023/05/08: remove the perturbation variables in estimating the wave speed
    fpdtype_t a = sqrt(${0.25*c['gamma']}*(pbl + pbr)/(ul[${bnvars}] + ur[${bnvars}]))
                + 0.25*fabs(nv);

    // Output
% for i in range(bnvars):
    nf[${i}] = 0.5*(${' + '.join('n[{j}]*(fl[{j}][{i}] + fr[{j}][{i}])'
                                 .format(i=i, j=j) for j in range(ndims))})
             + a*(ul[${i}] - ur[${i}]);

% endfor
</%pyfr:macro>
