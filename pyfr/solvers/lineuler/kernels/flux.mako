# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='inviscid_flux' params='s, f, p, v'>
% for i in range(ndims):
    v[${i}] = s[${i+1}]/s[${bnvars}];
% endfor

  // Notation for the pressure perturbation
  p = s[${bnvars - 1}];

  // Density and pressure flux
% for i in range(ndims):
    f[${i}][0] = s[0]*s[${i+1+bnvars}] + s[${i+1}];
    f[${i}][${bnvars - 1}] = ${c['gamma']}*s[${nvars - 1}]*v[${i}] + s[${i+1+bnvars}]*p;
% endfor


  // Momentum fluxes
% for i, j in pyfr.ndrange(ndims, ndims):
    f[${i}][${j + 1}] = s[${j+1}]*s[${i+1+bnvars}]${' + p' if i == j else ''};
% endfor
</%pyfr:macro>
