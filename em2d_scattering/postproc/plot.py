#! /usr/bin/env python
# -*- coding: utf-8 -*-

from sys import argv
from pylab import *

frame = int(argv[1])
filename = "./output/Ez_%06d.dat" %frame
d = loadtxt(filename)

profile = "./output/Ga_%06d.dat" %0
p = loadtxt(profile)

#interp = 'bilinear';
interp = 'nearest';
title("T = %d" %frame)
imshow(p, origin='lower', interpolation=interp)
#imshow(d, origin='lower', interpolation=interp)
colorbar()
show()

