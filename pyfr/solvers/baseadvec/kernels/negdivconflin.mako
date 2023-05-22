# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='negdivconflin' ndim='2'
              t='scalar fpdtype_t'
              tdivtconf='inout fpdtype_t[${str(nvars)}]'
              ploc='in fpdtype_t[${str(ndims)}]'
              u='in fpdtype_t[${str(nvars)}]'
              rcpdjac='in fpdtype_t'
              cu='in fpdtype_t[${str(bnvars)}]'>

// Compute the C@U term in the formula and combine it as a forcing term
% for i, ex in enumerate(srcex):
    tdivtconf[${i}] = -rcpdjac*tdivtconf[${i}] + ${ex};
     //- cu[${i}];
% endfor
% for j in range(bnvars):
    // set baseflow flux to 0
    tdivtconf[${j + bnvars}] = 0;
% endfor
</%pyfr:kernel>
