# -*- coding: utf-8 -*-

import numpy as np

# from pyfr.solvers.linadvec import LinearAdvectionElements
from pyfr.solvers.baseadvecdiff import BaseAdvectionDiffusionElements
from pyfr.solvers.lineuler.elements import BaseFluidElements


class LinearNavierStokesElements(BaseFluidElements, BaseAdvectionDiffusionElements):

    @staticmethod
    def grad_con_to_pri(cons, grad_cons, cfg):
        ptr = int(len(cons)/2)
        rho, *rhouvw = cons[:ptr-2]
        # unsure about the dimension of grad_cons

        grad_rho, *grad_rhouvw, grad_p = grad_cons[:][:ptr-1]
        rhob = cons[ptr]

        # Divide momentum components by ρ
        uvw = [rhov / rhob for rhov in rhouvw]

        # Velocity gradients: ∇u⃗ = 1/ρ·[∇(ρu⃗) - u⃗ ⊗ ∇ρ]
        grad_rhob = grad_cons[:][ptr]
        grad_uvw = [(grad_rhov - v*grad_rhob) / rhob
                    for grad_rhov, v in zip(grad_rhouvw, uvw)]

        return [grad_rho] + grad_uvw + [grad_p]

    def set_backend(self, *args, **kwargs):
        super().set_backend(*args, **kwargs)

        # Can elide interior flux calculations at p = 0
        if self.basis.order == 0:
            return

        # Register our flux kernels
        self._be.pointwise.register('pyfr.solvers.linnavstokes.kernels.tflux')
        self._be.pointwise.register('pyfr.solvers.linnavstokes.kernels.tfluxlin')

        # Handle Sutherland's law
        # remove shock capturing
        visc_corr = self.cfg.get('solver', 'viscosity-correction', 'none')
        if visc_corr not in {'sutherland', 'none'}:
            raise ValueError('Invalid viscosity-correction option')

        # Template parameters for the flux kernels
        tplargs = {
            'ndims': self.ndims,
            'nvars': self.nvars,
            'bnvars': self.bnvars,
            'nverts': len(self.basis.linspts),
            'c': self.cfg.items_as('constants', float),
            'jac_exprs': self.basis.jac_exprs,
            'visc_corr': visc_corr
        }

        # Helpers
        c, l = 'curved', 'linear'
        r, s = self._mesh_regions, self._slice_mat

        if c in r and 'flux' not in self.antialias:
            self.kernels['tdisf_curved'] = lambda uin: self._be.kernel(
                'tflux', tplargs=tplargs, dims=[self.nupts, r[c]],
                u=s(self.scal_upts[uin], c), f=s(self._vect_upts, c),
                smats=self.curved_smat_at('upts')
            )
        elif c in r:
            self.kernels['tdisf_curved'] = lambda: self._be.kernel(
                'tflux', tplargs=tplargs, dims=[self.nqpts, r[c]],
                u=s(self._scal_qpts, c), f=s(self._vect_qpts, c),
                smats=self.curved_smat_at('qpts')
            )

        if l in r and 'flux' not in self.antialias:
            self.kernels['tdisf_linear'] = lambda uin: self._be.kernel(
                'tfluxlin', tplargs=tplargs, dims=[self.nupts, r[l]],
                u=s(self.scal_upts[uin], l), f=s(self._vect_upts, l),
                verts=self.ploc_at('linspts', l),
                upts=self.upts
            )
        elif l in r:
            self.kernels['tdisf_linear'] = lambda: self._be.kernel(
                'tfluxlin', tplargs=tplargs, dims=[self.nqpts, r[l]],
                u=s(self._scal_qpts, l), f=s(self._vect_qpts, l),
                verts=self.ploc_at('linspts', l),
                upts=self.qpts
            )
