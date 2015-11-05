# -*- coding: utf-8 -*-

import os

from pyfr.inifile import Inifile
from pyfr.readers.native import NativeReader
from pyfr.solvers.base import BaseSystem
from pyfr.util import subclass_where


class BaseWriter(object):
    def __init__(self, args):
        self.outf = args.outf

        # Solution file type; usually .pyfrs or .pyfrd
        self.solnextn = os.path.splitext(args.solnf)[1]

        # Load the mesh and solution files
        self.soln = NativeReader(args.solnf)
        self.mesh = NativeReader(args.meshf)

        # Get element types and array shapes
        self.mesh_inf = self.mesh.array_info
        self.soln_inf = self.soln.array_info

        # Dimensions
        self.ndims = next(iter(self.mesh_inf.values()))[1][2]
        self.nvars = next(iter(self.soln_inf.values()))[1][1]

        # Check solution and mesh are compatible
        if self.mesh['mesh_uuid'] != self.soln['mesh_uuid']:
            raise RuntimeError('Solution "%s" was not computed on mesh "%s"' %
                               (args.solnf, args.meshf))

        # Load the configuration and stats files
        self.cfg = Inifile(self.soln['config'])
        self.stats = Inifile(self.soln['stats'])

        # System and elements classs
        self.systemscls = subclass_where(
            BaseSystem, name=self.cfg.get('solver', 'system')
        )
        self.elementscls = self.systemscls.elementscls
