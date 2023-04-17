# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='bc_common_grad_zero' params='ul, nl, gradul, gradur'>
% for i, j in pyfr.ndrange(ndims, bnvars):
    gradur[${i}][${j}] = 0;
% endfor
</%pyfr:macro>

<%pyfr:macro name='bc_common_grad_copy' params='ul, nl, gradul, gradur'>
% for i, j in pyfr.ndrange(ndims, bnvars):
    gradur[${i}][${j}] = gradul[${i}][${j}];
% endfor
</%pyfr:macro>

<%pyfr:macro name='bc_common_base_var_copy' params='ul, ur'>
% for i in range(bnvars):
    ur[${i + bnvars}] = ul[${i + bnvars}];
% endfor
</%pyfr:macro>

<%pyfr:macro name='bc_common_base_grad_copy' params='ul, nl, gradul, gradur'>
% for i, j in pyfr.ndrange(ndims, bnvars):
    gradur[${i}][${bnvars+j}] = gradul[${i}][${bnvars+j}];
% endfor
</%pyfr:macro>

<%pyfr:macro name='bc_common_base_grad_zero' params='ul, ur, gradul, gradur'>
% for i, j in pyfr.ndrange(ndims, bnvars):
    gradur[${i}][${bnvars+j}] = 0;
% endfor
</%pyfr:macro>