# -*- coding: utf-8 -*-

import numpy as np
import scipy.io
import os
import scipy.interpolate
import matplotlib.pyplot as plt
def _get_inter_objs(interside, getter, elemap):
    # Map from element type to view mat getter
    emap = {type: getattr(ele, getter) for type, ele in elemap.items()}

    # Get the data from the interface
    return [emap[type](eidx, fidx) for type, eidx, fidx, flags in interside]


class BaseInters:
    def __init__(self, be, lhs, elemap, cfg):
        self._be = be
        self.elemap = elemap
        self.cfg = cfg

        # Get the number of dimensions and variables
        self.ndims = next(iter(elemap.values())).ndims
        self.nvars = next(iter(elemap.values())).nvars
        """
        MODIFICATION FOR LINEAR SOLVER
        """
        if cfg.get('solver','solver-type','None') == 'linear':
            self.bnvars = next(iter(elemap.values())).bnvars
        """
        MODIFICATION FOR LINEAR SOLVER
        """

        # Get the number of interfaces
        self.ninters = len(lhs)

        # Compute the total number of interface flux points
        self.ninterfpts = sum(elemap[etype].nfacefpts[fidx]
                              for etype, eidx, fidx, flags in lhs)

        # By default do not permute any of the interface arrays
        self._perm = Ellipsis

        # Kernel constants
        self.c = cfg.items_as('constants', float)

        # Kernels and MPI requests we provide
        self.kernels = {}
        self.mpireqs = {}

        # Global kernel arguments
        self._external_args = {}
        self._external_vals = {}

    def prepare(self, t):
        pass

    def _set_external(self, name, spec, value=None):
        self._external_args[name] = spec

        if value is not None:
            self._external_vals[name] = value

    def _const_mat(self, inter, meth):
        m = _get_inter_objs(inter, meth, self.elemap)

        # Swizzle the dimensions and permute
        m = np.concatenate(m)
        m = np.atleast_2d(m.T)
        m = m[:, self._perm]

        return self._be.const_matrix(m)

    """
    MODIFICATION FOR LINEAR SOLVER
    """
    def _const_mat_inlet(self, inter, meth):
        coords = _get_inter_objs(inter, meth, self.elemap)
        coords = np.concatenate(coords)
        coords = np.atleast_2d(coords.T)
        coords = coords[:, self._perm]

        # read disturbance data
        dirname = self.cfg.get('soln-UDI', 'UDI-dirname', 'None')
        filenames = self.cfg.get('soln-UDI', 'UDI-filenames', 'None')
        filenames_list = filenames.split(',')
        filenames_list = [filename.strip() for filename in filenames_list]
        fmt = self.cfg.get('soln-UDI', 'UDI-format', 'None')
        # load perturbation
        # !! note that our data is from large to small, so we flipped the data !!
        qs = [np.flip(scipy.io.loadmat(os.path.join(dirname, filename + fmt))[filename]) for filename in
              filenames_list]
        # load coordinates of perturbation
        coord_filename = self.cfg.get('soln-UDI', 'UDI-coords-filename', 'None')
        y = scipy.io.loadmat(os.path.join(dirname, coord_filename+fmt))[coord_filename]
        y = np.flip(np.squeeze(y))
        cs = [scipy.interpolate.CubicSpline(y, q) for q in qs]
        # interpolate value at inlet boundary
        q_tmp = [c(coords[1,:]).flatten() for c in cs]
        # we need to separate the real and imagine part
        q_intp = []
        for q in q_tmp:
            q_intp = np.concatenate((q_intp, np.real(q)))
            q_intp = np.concatenate((q_intp, np.imag(q)))

        q_dim = (self.ndims+2)*2
        q_intp = np.reshape(q_intp, (q_dim, self.ninterfpts))
        # plot all interpolated perturbation
        # for i in range(0, (self.ndims+2)):
        #     plt.scatter(coords[1,:], q_intp[2 * i, :], label='interpolation')
        #     plt.legend(loc='best')
        #     plt.show()

        return self._be.const_matrix(q_intp)

    """
    MODIFICATION FOR LINEAR SOLVER
    """

    def _get_perm_for_view(self, inter, meth):
        vm = _get_inter_objs(inter, meth, self.elemap)
        vm = [np.concatenate(m) for m in zip(*vm)]
        mm = self._be.view(*vm, vshape=()).mapping.get()

        return np.argsort(mm[0])

    def _view(self, inter, meth, vshape=()):
        vm = _get_inter_objs(inter, meth, self.elemap)
        vm = [np.concatenate(m)[self._perm] for m in zip(*vm)]
        return self._be.view(*vm, vshape=vshape)

    def _scal_view(self, inter, meth):
        return self._view(inter, meth, (self.nvars,))

    def _vect_view(self, inter, meth):
        return self._view(inter, meth, (self.ndims, self.nvars))

    def _xchg_view(self, inter, meth, vshape=()):
        vm = _get_inter_objs(inter, meth, self.elemap)
        vm = [np.concatenate(m)[self._perm] for m in zip(*vm)]
        return self._be.xchg_view(*vm, vshape=vshape)

    def _scal_xchg_view(self, inter, meth):
        return self._xchg_view(inter, meth, (self.nvars,))

    def _vect_xchg_view(self, inter, meth):
        return self._xchg_view(inter, meth, (self.ndims, self.nvars))
