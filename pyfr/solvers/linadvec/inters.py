# -*- coding: utf-8 -*-

from pyfr.solvers.baseadvecdiff import (BaseAdvectionDiffusionIntInters,
                                    BaseAdvectionDiffusionMPIInters,
                                    BaseAdvectionDiffusionBCInters)

class TplargsMixin:
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        rsolver = self.cfg.get('solver-interfaces', 'riemann-solver')
        self._tplargs = dict(ndims=self.ndims, nvars=self.nvars,
                             rsolver=rsolver, c=self.c)

        self._tplargs['bnvars'] = self.bnvars


class LinearAdvectionIntInters(TplargsMixin, BaseAdvectionDiffusionIntInters):
    def __init__(self, be, lhs, rhs, elemap, cfg):
        super().__init__(be, lhs, rhs, elemap, cfg)

        be.pointwise.register('pyfr.solvers.linadvec.kernels.baseintconu')

        self.kernels['con_u'] = lambda: be.kernel(
            'baseintconu', tplargs=self._tplargs, dims=[self.ninterfpts],
            ulin=self._scal_lhs, urin=self._scal_rhs,
            ulout=self._vect_lhs, urout=self._vect_rhs
        )

class LinearAdvectionMPIInters(TplargsMixin, BaseAdvectionDiffusionMPIInters):
    def __init__(self, be, lhs, rhsrank, rallocs, elemap, cfg):
        super().__init__(be, lhs, rhsrank, rallocs, elemap, cfg)

        be.pointwise.register('pyfr.solvers.linadvec.kernels.basempiconu')

        self.kernels['con_u'] = lambda: be.kernel(
            'basempiconu', tplargs=self._tplargs, dims=[self.ninterfpts],
            ulin=self._scal_lhs, urin=self._scal_rhs, ulout=self._vect_lhs
        )

class LinearAdvectionBCInters(TplargsMixin, BaseAdvectionDiffusionBCInters):
    def __init__(self, be, lhs, elemap, cfgsect, cfg):
        super().__init__(be, lhs, elemap, cfgsect, cfg)

        # Additional BC specific template arguments
        #self._tplargs['bctype'] = self.type

        be.pointwise.register('pyfr.solvers.linadvec.kernels.basebcconu')

        self.kernels['con_u'] = lambda: be.kernel(
            'basebcconu', tplargs=self._tplargs, dims=[self.ninterfpts],
            extrns=self._external_args, ulin=self._scal_lhs,
            ulout=self._vect_lhs, nlin=self._pnorm_lhs,
            **self._external_vals
        )
