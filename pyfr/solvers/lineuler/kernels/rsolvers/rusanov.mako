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
    //fpdtype_t vbl[${ndims}], vbr[${ndims}];

    fpdtype_t c0l = sqrt(${c['gamma']}*pbl/ul[${bnvars}]);
    fpdtype_t c0r = sqrt(${c['gamma']}*pbr/ur[${bnvars}]);
    fpdtype_t invrhobl = ul[${bnvars}];
    fpdtype_t invrhobr = ur[${bnvars}];

    // Estimate the maximum wave speed / 2
    // Use the Nektar++ way, explained in SPECFEM2D-DG-LNS paper.
    // fpdtype_t a = max(abs(ul[${bnvars+1}]+ul[1]*invrhobl-c0l),abs(ul[${bnvars+1}]+ul[1]*invrhobl+c0l));
    // a = max(abs(ur[${bnvars+1}]+ur[1]*invrhobr-c0r), a);
    // a = max(abs(ur[${bnvars+1}]+ur[1]*invrhobr+c0r), a);

    fpdtype_t a = max(abs(ul[${bnvars+1}]-c0l),abs(ul[${bnvars+1}]+c0l));
    a = max(abs(ur[${bnvars+1}]-c0r), a);
    a = max(abs(ur[${bnvars+1}]+c0r), a);

    // fpdtype_t nbvl = ${pyfr.dot('n[{i}]', 'vbl[{i}]', i=ndims)};
    // fpdtype_t nbvr = ${pyfr.dot('n[{i}]', 'vbr[{i}]', i=ndims)};

    // fpdtype_t a = max(abs(nbvl)+c0l, abs(nbvr)+c0r);

    ${pyfr.expand('inviscid_flux', 'ul', 'fl', 'pl', 'vl')};
    ${pyfr.expand('inviscid_flux', 'ur', 'fr', 'pr', 'vr')};

    // Sum the left and right velocities and take the normal
    // fpdtype_t nv = ${pyfr.dot('n[{i}]', 'vl[{i}] + vr[{i}]', i=ndims)};

    // Estimate the maximum wave speed / 2
    // 2023/05/08: remove the perturbation variables in estimating the wave speed
    //fpdtype_t a = sqrt(${0.25*c['gamma']}*(pbl + pbr)/(ul[${bnvars}] + ur[${bnvars}]))
    //            + 0.25*fabs(nv);

    // Output
% for i in range(bnvars):
    nf[${i}] = 0.5*(${' + '.join('n[{j}]*(fl[{j}][{i}] + fr[{j}][{i}])'
                                 .format(i=i, j=j) for j in range(ndims))})
             + a*0.5*(ul[${i}] - ur[${i}]);

% endfor
</%pyfr:macro>
