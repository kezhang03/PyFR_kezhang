# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<% gmo = c['gamma'] - 1.0 %>
<% gamma = c['gamma'] %>

<%pyfr:macro name='bc_rsolve_state' params='ul, nl, ur' externs='ploc, t, q1, q2'>
    fpdtype_t inv = 1.0 / ul[${bnvars}];
    fpdtype_t V_n0 = inv*${' + '.join('ul[{1}]*nl[{0}]'.format(i + bnvars, i + 1)
                                 for i in range(ndims))};

    fpdtype_t V_n = inv*${' + '.join('ul[{1}]*nl[{0}]'.format(i, i + 1)
                                 for i in range(ndims))};

    fpdtype_t V_e = ${' + '.join('{0}*nl[{1}]'.format(c['uvw'[i]], i)
                                 for i in range(ndims))};

    fpdtype_t c = sqrt(${gamma}*ul[${nvars-1}]*inv);


    fpdtype_t h1 = (V_n0 > 0)
                 ? ul[0] - ul[${bnvars-1}] / (c * c)    //outgoing characteristics
                 : ${c['rho']} - ${c['p']} / (c * c) ;  //incoming characteristics


    fpdtype_t h4 = (V_n0 > c)
                 ? V_n/2.0 - ul[${bnvars-1}]/(2.0*c)
                 : ul[${bnvars}]*V_e/2.0 - ${c['p']}/(2.0*c);

    fpdtype_t h5 = (V_n0 + c >0)
                 ? V_n/2.0 + ul[${bnvars-1}]/(2.0*c)
                 : ul[${bnvars}]*V_e/2.0 + ${c['p']}/(2.0*c);

    ur[0] = h1 + (h5 - h4)/c;

    ur[${bnvars-1}] = c * (h5 - h4);

    ur[1] = (h4 + h5)*nl[0];

% for i in range(ndims - 1):
    ur[${i + 1 + 1}] = ul[${bnvars}]*${c['uvw'[i + 1]]}; // + (h4 + h5)*nl[${i + 1}];
% endfor

% for i in range(bnvars):
    ur[${i + bnvars}] = ul[${i + bnvars}];
% endfor
</%pyfr:macro>
