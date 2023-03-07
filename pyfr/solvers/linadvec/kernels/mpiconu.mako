# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='mpiconu' ndim='1'
              ulin='in view fpdtype_t[${str(nvars)}]'
              urin='in mpi fpdtype_t[${str(nvars)}]'
              ulout='out view fpdtype_t[${str(nvars)}]'>
//% for i in range(nvars):
% for i in range(bnvars):
% if c['ldg-beta'] == -0.5:
    ulout[${i+bnvars}] = ulin[${i+bnvars}];
% elif c['ldg-beta'] == 0.5:
    ulout[${i+bnvars}] = urin[${i+bnvars}];
% else:
    ulout[${i+bnvars}] = urin[${i+bnvars}]*${0.5 + c['ldg-beta']}
                + ulin[${i+bnvars}]*${0.5 - c['ldg-beta']};
% endif
% endfor
</%pyfr:kernel>

