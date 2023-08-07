# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.linnavstokes.kernels.bcs.common'/>

<%include file='pyfr.solvers.lineuler.kernels.bcs.char-riem-inv-in-test'/>
<%pyfr:alias name='bc_ldg_state' func='bc_rsolve_state'/>
<%pyfr:alias name='bc_ldg_grad_state' func='bc_common_grad_copy'/>
<%pyfr:alias name='bc_ldg_base_grad_state' func='bc_common_base_grad_copy'/>