#!/usr/bin/env python
'''Initialization of pyDOG, a module for creating receptive fields
based on the eDOG model proposed in "Einevoll GT and Plesser HE
(2012), Extended difference-of-Gaussians model incorporating cortical
feedback for relay cells in the lateral geniculate nucleus of
cat. Cogn Neurodyn 6:307-324"

Group of Computational Neuroscience (compneuro.umb.no),
Department of Mathematical Sciences and Technology,
Norwegian University of Life Sciences.

Copyright (C) 2012 Computational Neuroscience Group, UMB.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
'''

from tcloop import DOG, full_eDOG, eDOG, eDOG_FastLoop
from stimgen import slidingBar, flashingSpots, patchGrating, driftingGrating
import plot_tools
