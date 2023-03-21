# -*- coding: utf-8 -*-

# from pyfr.solvers.linadvec import LinearAdvectionSystem
from pyfr.solvers.baseadvecdiff import BaseAdvectionDiffusionSystem
from pyfr.solvers.linnavstokes.elements import LinearNavierStokesElements
from pyfr.solvers.linnavstokes.inters import (LinearNavierStokesBaseBCInters,
                                           LinearNavierStokesIntInters,
                                           LinearNavierStokesMPIInters)


class LinearNavierStokesSystem(BaseAdvectionDiffusionSystem):
    name = 'linear-navier-stokes'

    elementscls = LinearNavierStokesElements
    intinterscls = LinearNavierStokesIntInters
    mpiinterscls = LinearNavierStokesMPIInters
    bbcinterscls = LinearNavierStokesBaseBCInters
