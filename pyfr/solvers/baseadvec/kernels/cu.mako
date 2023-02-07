# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='cu' ndim='2'
              t='scalar fpdtype_t'
              u='in fpdtype_t[${str(nvars)}]'
              divub='in fpdtype_t[${str(bnvars)}][${str(bnvars)}]'
              cu='out fpdtype_t[${str(bnvars)}]'>

    // Compute the average quantities
    fpdtype_t p;
    fpdtype_t v[${bnvars}];
    fpdtype_t invrhob = 1.0/u[${bnvars}];
% for i in range(ndims):
    v[${i}] = invrhob*u[${i + 1}];
% endfor
% if ndims == 2:
  // continuity equation
  cu[0] = 0;
  // momentum equation (x)
  cu[1] =  u[0]*(u[5]*divub[1][1] + u[6]*divub[1][2]) + u[4]*(v[1]*divub[1][1] + v[2]*divub[1][2]);
  // momentum equation (y)
  cu[2] =  u[0]*(u[5]*divub[2][1] + u[6]*divub[2][2]) + u[4]*(v[1]*divub[2][1] + v[2]*divub[2][2]);
  // energy equation
  cu[3] =  (${c['gamma'] - 1})*u[3]*(divub[1][1] + divub[2][2]) + cu[3] + (${1 - c['gamma']})*(v[1]*divub[3][1] + v[2]*divub[3][2]);
% elif ndims == 3:
  cu[0] = 0;
  cu[1] =  u[0]*(u[6]*divub[1][1] + u[7]*divub[1][2] + u[8]*divub[1][3]) + u[5]*v[1]*(divub[1][1] + divub[2][2] + divub[3][3]);
  cu[2] =  u[0]*(u[6]*divub[2][1] + u[7]*divub[2][2] + u[8]*divub[2][3]) + u[5]*v[2]*(divub[1][1] + divub[2][2] + divub[3][3]);
  cu[3] =  u[0]*(u[6]*divub[3][1] + u[7]*divub[3][2] + u[8]*divub[3][3]) + u[5]*v[3]*(divub[1][1] + divub[2][2] + divub[3][3]);
  cu[4] =  (${c['gamma'] - 1})*u[4]*(divub[1][1] + divub[2][2] + divub[3][3]) + (${1 - c['gamma']})*(v[1]*divub[4][1] + v[2]*divub[4][2] + v[3]*divub[4][3]);
% endif
</%pyfr:kernel>