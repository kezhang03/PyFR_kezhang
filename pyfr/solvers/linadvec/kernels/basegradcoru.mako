# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='basegradcoru' ndim='2'
              gradu='inout fpdtype_t[${str(ndims)}][${str(nvars)}]'
              smats='in fpdtype_t[${str(ndims)}][${str(ndims)}]'
              rcpdjac='in fpdtype_t'
              >
    // gradu='in fpdtype_t[${str(ndims)}][${str(nvars)}]'
    // gradbaseu='out fpdtype_t[${str(ndims)}][${str(bnvars)}]'
    fpdtype_t tmpgradu[][${nvars}] = ${pyfr.array('gradu[{i}][{j}]', i=ndims, j=nvars)};

% for i, j in pyfr.ndrange(ndims, bnvars):
    gradu[${i}][${j+bnvars}] = rcpdjac*(${' + '.join(f'smats[{k}][{i}]*tmpgradu[{k}][{j+bnvars}]'
                                              for k in range(ndims))});
% endfor
</%pyfr:kernel>
