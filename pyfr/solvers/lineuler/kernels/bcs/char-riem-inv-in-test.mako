# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<% gmo = c['gamma'] - 1.0 %>
<% gamma = c['gamma'] %>

<%pyfr:macro name='bc_rsolve_state' params='ul, nl, ur' externs='ploc, t, q1, q2'>
    fpdtype_t inv = 1.0 / ul[${bnvars}];
    // baseflow uses primitive variables
    fpdtype_t V_n0 = ${' + '.join('ul[{0}]*nl[{1}]'.format(i + 1 + bnvars, i)
                                 for i in range(ndims))};
    // perturbation uses conservative variables, Vn = rhob * u dot n
    fpdtype_t V_n = ${' + '.join('ul[{0}]*nl[{1}]'.format(i + 1, i)
                                 for i in range(ndims))};
    // prescribed values use conservative form
    fpdtype_t V_e = ul[${bnvars}] * (${' + '.join('{0}*nl[{1}]'.format(c['uvw'[i]], i)
                                 for i in range(ndims))});

    fpdtype_t c = sqrt(${gamma}*ul[${nvars-1}]*inv);


    fpdtype_t h1 = ${c['rho']} - ${c['p']} / (c * c) ;  //incoming characteristics

    // h4 = (rhob*u + p/cb)/2
    fpdtype_t h4 = V_e/2.0 + ${c['p']}/(2.0*c);

    // h5 = (rhob*u - p/cb)/2
    fpdtype_t h5 = -V_n/2.0 + ul[${bnvars-1}]/(2.0*c);

    ur[0] = h1 + (h4 + h5)/c;
    ur[${bnvars-1}] = c * (h4 + h5);
    ur[1] = (h4 - h5) * nl[0];
% for i in range(ndims-1):
    ur[${i + 1 + 1}] = ul[${bnvars}] * ${c['uvw'[i + 1]]};
% endfor

//copy baseflow (necessary for LEE because ghost.mako only works for LNS)
% for i in range(bnvars):
    ur[${i + bnvars}] = ul[${i + bnvars}];
% endfor
</%pyfr:macro>