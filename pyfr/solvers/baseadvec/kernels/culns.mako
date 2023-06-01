# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='culns' ndim='2'
              t='scalar fpdtype_t'
              u='in fpdtype_t[${str(nvars)}]'
              grad_uin='in fpdtype_t[${str(ndims)}][${str(nvars)}]'
              cu='out fpdtype_t[${str(bnvars)}]'>

    // Compute the average quantities
    fpdtype_t v[${ndims}];
    fpdtype_t invrhob = 1.0/u[${bnvars}];
% for i in range(ndims):
    v[${i}] = invrhob*u[${i + 1}];
% endfor
    fpdtype_t pb = u[${nvars-1}];
    fpdtype_t rhob = u[${bnvars}];
    fpdtype_t p = u[${bnvars-1}];
    fpdtype_t rho = u[0];

% if visc_corr == 'sutherland':
    // Compute the temperature and viscosity
    // Use baseflow
    fpdtype_t cpT = ${c['gamma']}/${c['gamma']-1}*(pb/rhob);
    fpdtype_t Trat = ${1/c['cpTref']}*cpT;
    fpdtype_t mu_c = ${c['mu']*(c['cpTref'] + c['cpTs'])}*Trat*sqrt(Trat)
                   / (cpT + ${c['cpTs']});
    // Compute viscosity perturbation
    fpdtype_t dmudcpT = mu_c*sqrt(Trat)/${c['cpTref']*(c['cpTref'] + c['cpTs'])}/(cpT + ${c['cpTs']})*(1.5-cpT/(cpT + ${c['cpTs']}));
    fpdtype_t mu_p = dmudcpT*cpT*(p/pb-rho/rhob);
% else:
    fpdtype_t mu_c = ${c['mu']};
    fpdtype_t mu_p = 0.0;
% endif

% if ndims == 2:
// baseflow derivatives (grad[rhob,ub,vb,pb])
  fpdtype_t rhob_x = grad_uin[0][4];
  fpdtype_t rhob_y = grad_uin[1][4];

  fpdtype_t ub_x = grad_uin[0][5];
  fpdtype_t ub_y = grad_uin[1][5];

  fpdtype_t vb_x = grad_uin[0][6];
  fpdtype_t vb_y = grad_uin[1][6];

  fpdtype_t pb_x = grad_uin[0][7];
  fpdtype_t pb_y = grad_uin[1][7];

  // perturbation derivatives (rho*grad[u,v,w])
  fpdtype_t u_x = grad_uin[0][1] - v[0]*rhob_x;
  fpdtype_t u_y = grad_uin[1][1] - v[0]*rhob_y;

  fpdtype_t v_x = grad_uin[0][2] - v[1]*rhob_x;
  fpdtype_t v_y = grad_uin[1][2] - v[1]*rhob_y;

  // negated stress tensors
  fpdtype_t txx = -2*mu_c*invrhob*(u_x - ${1.0/3.0}*(u_x + v_y));
  fpdtype_t tyy = -2*mu_c*invrhob*(v_y - ${1.0/3.0}*(u_x + v_y));
  fpdtype_t txy = -mu_c*invrhob*(v_x + u_y);

  fpdtype_t t0xx = -2*mu_c*(ub_x - ${1.0/3.0}*(ub_x + vb_y));
  fpdtype_t t0yy = -2*mu_c*(vb_y - ${1.0/3.0}*(ub_x + vb_y));
  fpdtype_t t0xy = -mu_c*(vb_x + ub_y);

  // stress tensor due to viscosity perturbation (use baseflow stress to save time)
  txx += t0xx / mu_c * mu_p;
  tyy += t0yy / mu_c * mu_p;
  txy += t0xy / mu_c * mu_p;


  // continuity equation
  cu[0] = 0;
  // momentum equation (x)
  // grad: 1st index is the dimension
  // grad: 2nd index is the variable
  cu[1] =  u[0]*(u[5]*ub_x + u[6]*ub_y) + u[4]*(v[0]*ub_x + v[1]*ub_y);
  // momentum equation (y)
  cu[2] =  u[0]*(u[5]*vb_x + u[6]*vb_y) + u[4]*(v[0]*vb_x + v[1]*vb_y);
  // energy equation
  cu[3] = (v[0]*pb_x + v[1]*pb_y) - u[3]*(ub_x + vb_y);
  cu[3] += (t0xx * u_x + t0xy * (u_y + v_x) + t0yy * v_y) * invrhob;
  cu[3] += txx * ub_x + txy * (ub_y + vb_x) +  tyy * vb_y;
  cu[3] *= (1 - ${c['gamma']});

% elif ndims == 3:
// baseflow derivatives (grad[u,v,w])
  fpdtype_t rhob_x = grad_uin[0][5];
  fpdtype_t rhob_y = grad_uin[1][5];
  fpdtype_t rhob_z = grad_uin[2][5];

  fpdtype_t ub_x = grad_uin[0][6];
  fpdtype_t ub_y = grad_uin[1][6];
  fpdtype_t ub_z = grad_uin[2][6];

  fpdtype_t vb_x = grad_uin[0][7];
  fpdtype_t vb_y = grad_uin[1][7];
  fpdtype_t vb_z = grad_uin[2][7];

  fpdtype_t wb_x = grad_uin[0][8];
  fpdtype_t wb_y = grad_uin[1][8];
  fpdtype_t wb_z = grad_uin[2][8];

  fpdtype_t pb_x = grad_uin[0][9];
  fpdtype_t pb_y = grad_uin[1][9];
  fpdtype_t pb_z = grad_uin[2][9];

  // perturbation derivatives (rho*grad[u,v,w])
  fpdtype_t u_x = grad_uin[0][1] - v[0]*rhob_x;
  fpdtype_t u_y = grad_uin[1][1] - v[0]*rhob_y;
  fpdtype_t u_z = grad_uin[2][1] - v[0]*rhob_z;

  fpdtype_t v_x = grad_uin[0][2] - v[1]*rhob_x;
  fpdtype_t v_y = grad_uin[1][2] - v[1]*rhob_y;
  fpdtype_t v_z = grad_uin[2][2] - v[1]*rhob_z;

  fpdtype_t w_x = grad_uin[0][3] - v[2]*rhob_x;
  fpdtype_t w_y = grad_uin[1][3] - v[2]*rhob_y;
  fpdtype_t w_z = grad_uin[2][3] - v[2]*rhob_z;

  // negated stress tensors
  fpdtype_t txx = -2*mu_c*invrhob*(u_x - ${1.0/3.0}*(u_x + v_y + w_z));
  fpdtype_t tyy = -2*mu_c*invrhob*(v_y - ${1.0/3.0}*(u_x + v_y + w_z));
  fpdtype_t tzz = -2*mu_c*invrhob*(w_z - ${1.0/3.0}*(u_x + v_y + w_z));
  fpdtype_t txy = -mu_c*invrhob*(v_x + u_y);
  fpdtype_t txz = -mu_c*invrhob*(u_z + w_x);
  fpdtype_t tyz = -mu_c*invrhob*(w_y + v_z);

  fpdtype_t t0xx = -2*mu_c*(ub_x - ${1.0/3.0}*(ub_x + vb_y + wb_z));
  fpdtype_t t0yy = -2*mu_c*(vb_y - ${1.0/3.0}*(ub_x + vb_y + wb_z));
  fpdtype_t t0zz = -2*mu_c*(wb_z - ${1.0/3.0}*(ub_x + vb_y + wb_z));
  fpdtype_t t0xy = -mu_c*(vb_x + ub_y);
  fpdtype_t t0xz = -mu_c*(ub_z + wb_x);
  fpdtype_t t0yz = -mu_c*(wb_y + vb_z);

  // stress tensor due to viscosity perturbation (use baseflow stress to save time)
  txx += txx / mu_c * mu_p;
  tyy += tyy / mu_c * mu_p;
  tzz += tzz / mu_c * mu_p;
  txy += txy / mu_c * mu_p;
  txz += txz / mu_c * mu_p;
  tyz += tyz / mu_c * mu_p;

  cu[0] = 0;
  cu[1] =  u[0]*(u[6]*ub_x + u[7]*ub_y + u[8]*ub_z) + u[5]*(v[0]*ub_x + v[1]*ub_y + v[2]*ub_z);
  cu[2] =  u[0]*(u[6]*vb_x + u[7]*vb_y + u[8]*vb_z) + u[5]*(v[0]*vb_x + v[1]*vb_y + v[2]*vb_z);
  cu[3] =  u[0]*(u[6]*wb_x + u[7]*wb_y + u[8]*wb_z) + u[5]*(v[0]*wb_x + v[1]*wb_y + v[2]*wb_z);
  cu[4] =  (v[0]*pb_x + v[1]*pb_y + v[2]*pb_z) - u[4]*(ub_x + vb_y + wb_z);
  cu[4] += invrhob * (t0xx * u_x + t0xy * (u_y + v_x) + t0yy * v_y + t0xz * (w_x + u_z) + t0yz * (w_y + v_z) + t0zz * w_z);
  cu[4] += (txx * ub_x + txy * (ub_y + vb_x) +  tyy * vb_y + txz * (wb_x + ub_z) + tyz * (wb_y + vb_z) + tzz * wb_z);
  cu[4] *= (1 - ${c['gamma']});
% endif
</%pyfr:kernel>