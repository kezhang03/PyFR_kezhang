# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='culee' ndim='2'
              t='scalar fpdtype_t'
              u='in fpdtype_t[${str(nvars)}]'
              divub='in fpdtype_t[${str(ndims)}][${str(nvars)}]'
              cu='out fpdtype_t[${str(bnvars)}]'>

    // Compute the average quantities
    fpdtype_t p;
    fpdtype_t v[${ndims}];
    fpdtype_t invrhob = 1.0/u[${bnvars}];
% for i in range(ndims):
    v[${i}] = invrhob*u[${i + 1}];
% endfor
% if ndims == 2:
  // continuity equation
  cu[0] = 0;
  // momentum equation (x)
  // grad: 1st index is the dimension
  // grad: 2nd index is the variable
  cu[1] =  u[0]*(u[5]*divub[0][5] + u[6]*divub[1][5]) + u[4]*(v[0]*divub[0][5] + v[1]*divub[1][5]);
  // momentum equation (y)
  cu[2] =  u[0]*(u[5]*divub[0][6] + u[6]*divub[1][6]) + u[4]*(v[0]*divub[0][6] + v[1]*divub[1][6]);
  // energy equation
  //cu[3] =  (${c['gamma'] - 1})*u[3]*(divub[0][5] + divub[1][6]) + (${1 - c['gamma']})*(v[0]*divub[0][7] + v[1]*divub[1][7]);
  cu[3] =  (1-${c['gamma']}) * (v[0]*divub[0][7] + v[1]*divub[1][7] - u[3]*(divub[0][5] + divub[1][6])) ;
% elif ndims == 3:
  cu[0] = 0;
  cu[1] =  u[0]*(u[5]*divub[0][5] + u[6]*divub[1][5] + u[7]*divub[2][5]) + u[4]*(v[1]*divub[0][5] + v[2]*divub[1][5] + v[3]*divub[2][5]);
  cu[2] =  u[0]*(u[5]*divub[0][6] + u[6]*divub[1][6] + u[7]*divub[2][6]) + u[4]*(v[1]*divub[0][6] + v[2]*divub[1][6] + v[3]*divub[2][6]);
  cu[3] =  u[0]*(u[5]*divub[0][7] + u[6]*divub[1][7] + u[7]*divub[2][7]) + u[4]*(v[1]*divub[0][7] + v[2]*divub[1][7] + v[3]*divub[2][7]);
  cu[4] =  (1-${c['gamma']}) * (v[0]*divub[0][7] + v[1]*divub[1][7] + v[2]*divub[2][7] - u[3]*(divub[0][5] + divub[1][6] + divub[2][7]));
% endif
</%pyfr:kernel>