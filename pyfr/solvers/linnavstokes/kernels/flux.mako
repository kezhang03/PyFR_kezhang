# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

% if ndims == 2:
<%pyfr:macro name='viscous_flux_add' params='u_in, grad_uin, fout'>

    fpdtype_t invrhob = 1.0/u_in[${bnvars}];
    fpdtype_t rho = u_in[0];
    fpdtype_t u = invrhob*u_in[1], v = invrhob*u_in[2];
    fpdtype_t p = u_in[3];
    fpdtype_t rhob = u_in[4];
    fpdtype_t pb = u_in[7];
    fpdtype_t invgmo = ${1/(c['gamma']-1)};
    fpdtype_t gmo = ${c['gamma']-1};

    // baseflow derivatives
    fpdtype_t rhob_x = grad_uin[0][4];
    fpdtype_t rhob_y = grad_uin[1][4];

    fpdtype_t ub_x = grad_uin[0][5];
    fpdtype_t ub_y = grad_uin[1][5];

    fpdtype_t vb_x = grad_uin[0][6];
    fpdtype_t vb_y = grad_uin[1][6];

    fpdtype_t pb_x = grad_uin[0][7];
    fpdtype_t pb_y = grad_uin[1][7];

    // perturbation derivatives (rhob*grad[u,v])
    fpdtype_t rho_x = grad_uin[0][0];
    fpdtype_t rho_y = grad_uin[1][0];

    fpdtype_t u_x = grad_uin[0][1] - u*rhob_x;
    fpdtype_t u_y = grad_uin[1][1] - u*rhob_y;

    fpdtype_t v_x = grad_uin[0][2] - v*rhob_x;
    fpdtype_t v_y = grad_uin[1][2] - v*rhob_y;

    fpdtype_t p_x = grad_uin[0][3];
    fpdtype_t p_y = grad_uin[1][3];

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

    // Compute temperature derivatives (Cv*dT/d[x,y])
    fpdtype_t Tb = invgmo*pb/rhob;
    fpdtype_t Tb_x = invgmo*(invrhob*pb_x-invrhob*invrhob*pb*rhob_x);
    fpdtype_t Tb_y = invgmo*(invrhob*pb_y-invrhob*invrhob*pb*rhob_y);

    fpdtype_t T_x = Tb_x * (p/pb-rho/rhob) + Tb * (p_x/pb - p/pb/pb*pb_x - rho_x*invrhob + rho*invrhob*invrhob*rhob_x);
    fpdtype_t T_y = Tb_y * (p/pb-rho/rhob) + Tb * (p_y/pb - p/pb/pb*pb_y - rho_y*invrhob + rho*invrhob*invrhob*rhob_y);

    // Negated stress tensor elements
    fpdtype_t txx = -2*mu_c*invrhob*(u_x - ${1.0/3.0}*(u_x + v_y));
    fpdtype_t tyy = -2*mu_c*invrhob*(v_y - ${1.0/3.0}*(u_x + v_y));
    fpdtype_t txy = -mu_c*invrhob*(v_x + u_y);

    // stress tensor due to viscosity perturbation
    txx += -2 * mu_p * (ub_x - ${1.0/3.0} * (ub_x + vb_y));
    tyy += -2 * mu_p * (vb_y - ${1.0/3.0} * (ub_x + vb_y));
    txy += -mu_p * (vb_x + ub_y);

    fout[0][1] += txx;
    fout[1][1] += txy;
    fout[0][2] += txy;
    fout[1][2] += tyy;

    // fout[0][3] += -mu_c*${c['gamma']/c['Pr']}*T_x;
    // fout[1][3] += -mu_c*${c['gamma']/c['Pr']}*T_y;

    // add viscosity perturbation to temperature derivatives
    fout[0][3] += gmo*(-mu_c*${c['gamma']/c['Pr']}*T_x - mu_p*${c['gamma']/c['Pr']}*Tb_x);
    fout[1][3] += gmo*(-mu_c*${c['gamma']/c['Pr']}*T_y - mu_p*${c['gamma']/c['Pr']}*Tb_y);

</%pyfr:macro>
% elif ndims == 3:
<%pyfr:macro name='viscous_flux_add' params='uin, grad_uin, fout'>
    fpdtype_t invrhob = 1.0/u_in[${bnvars}];
    fpdtype_t rho = u_in[0];
    fpdtype_t u = invrhob*u_in[1], v = invrhob*u_in[2], w = invrhob*u_in[3];
    fpdtype_t p = u_in[4];
    fpdtype_t rhob = u_in[5];
    fpdtype_t pb = u_in[9];

   // baseflow derivatives (grad[rho,u,v,w,p])
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

      // perturbation derivatives (rhob*grad[u,v,w])
      fpdtype_t u_x = grad_uin[0][1] - u*rhob_x;
      fpdtype_t u_y = grad_uin[1][1] - u*rhob_y;
      fpdtype_t u_z = grad_uin[2][1] - u*rhob_z;

      fpdtype_t v_x = grad_uin[0][2] - v*rhob_x;
      fpdtype_t v_y = grad_uin[1][2] - v*rhob_y;
      fpdtype_t v_z = grad_uin[2][2] - v*rhob_z;

      fpdtype_t w_x = grad_uin[0][3] - w*rhob_x;
      fpdtype_t w_y = grad_uin[1][3] - w*rhob_y;
      fpdtype_t w_z = grad_uin[2][3] - w*rhob_z;

% if visc_corr == 'sutherland':
    // Compute the temperature and viscosity
    // Use baseflow
    fpdtype_t cpT = ${c['gamma']}/${c['gamma']-1}*(pb/rhob);
    fpdtype_t Trat = ${1/c['cpTref']}*cpT;
    fpdtype_t mu_c = ${c['mu']*(c['cpTref'] + c['cpTs'])}*Trat*sqrt(Trat)
                   / (cpT + ${c['cpTs']});
    // Compute viscosity perturbation
    fpdtype_t dmudcpT = mu_c*sqrt(Trat)/c['cpTref']*(c['cpTref'] + c['cpTs'])/(cpT + ${c['cpTs']})*(1.5-cpT/(cpT + ${c['cpTs']}));
    fpdtype_t mu_p = dmudcpT*cpT*(p/pb-rho/rhob);
% else:
    fpdtype_t mu_c = ${c['mu']};
    fpdtype_t mu_p = 0.0;
% endif

    // Compute temperature derivatives (Cv*dT/d[x,y,z])
    fpdtype_t Tb = pb/rhob;
    fpdtype_t Tb_x = (invrhob*pb_x-pb*invrhob*invrhob*rhob_x);
    fpdtype_t Tb_y = (invrhob*pb_y-pb*invrhob*invrhob*rhob_y);
    fpdtype_t Tb_z = (invrhob*pb_z-pb*invrhob*invrhob*rhob_z);

    fpdtype_t T_x = Tb_x * (p/pb-rho/rhob) + Tb * (p_x/pb - p/pb/pb*pb_x - rho_x*invrhob + rho*invrhob*invrhob*rhob_x);
    fpdtype_t T_y = Tb_y * (p/pb-rho/rhob) + Tb * (p_y/pb - p/pb/pb*pb_y - rho_y*invrhob + rho*invrhob*invrhob*rhob_y);
    fpdtype_t T_z = Tb_z * (p/pb-rho/rhob) + Tb * (p_z/pb - p/pb/pb*pb_z - rho_z*invrhob + rho*invrhob*invrhob*rhob_z);

    // Negated stress tensor elements
    fpdtype_t txx = -2*mu_c*invrhob*(u_x - ${1.0/3.0}*(u_x + v_y + w_z));
    fpdtype_t tyy = -2*mu_c*invrhob*(v_y - ${1.0/3.0}*(u_x + v_y + w_z));
    fpdtype_t tzz = -2*mu_c*invrhob*(w_z - ${1.0/3.0}*(u_x + v_y + w_z));
    fpdtype_t txy = -mu_c*invrhob*(v_x + u_y);
    fpdtype_t txz = -mu_c*invrhob*(u_z + w_x);
    fpdtype_t tyz = -mu_c*invrhob*(w_y + v_z);

    // stress tensor due to viscosity perturbation
    txx += -2 * mu_p * (ub_x - ${1.0/3.0} * (ub_x + vb_y + wb_z));
    tyy += -2 * mu_p * (vb_y - ${1.0/3.0} * (ub_x + vb_y + wb_z));
    tzz += -2 * mu_p * (wb_z - ${1.0/3.0} * (ub_x + vb_y + wb_z));
    txy += -mu_p * (vb_x + ub_y);
    txz += -mu_p * (ub_z + wb_x);
    tyz += -mu_p * (wb_y + vb_z);

    fout[0][1] += txx;     fout[1][1] += txy;     fout[2][1] += txz;
    fout[0][2] += txy;     fout[1][2] += tyy;     fout[2][2] += tyz;
    fout[0][3] += txz;     fout[1][3] += tyz;     fout[2][3] += tzz;

    // fout[0][4] += -mu_c*${c['gamma']/c['Pr']}*T_x;
    // fout[1][4] += -mu_c*${c['gamma']/c['Pr']}*T_y;
    // fout[2][4] += -mu_c*${c['gamma']/c['Pr']}*T_z;

    // add viscosity perturbation to temperature derivatives
    fout[0][4] += -mu_c*${c['gamma']/c['Pr']}*T_x - mu_p*${c['gamma']/c['Pr']}*Tb_x;
    fout[1][4] += -mu_c*${c['gamma']/c['Pr']}*T_y - mu_p*${c['gamma']/c['Pr']}*Tb_y;
    fout[2][4] += -mu_c*${c['gamma']/c['Pr']}*T_z - mu_p*${c['gamma']/c['Pr']}*Tb_z;
</%pyfr:macro>
% endif
