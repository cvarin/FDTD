#! /usr/bin/env python
# -*- coding: utf-8 -*-

from sys import argv
from pylab import *

frame = int(argv[1])
filename = "./output/Ez_%06d.dat" %frame
d = loadtxt(filename)

profile = "./output/Ga_%06d.dat" %0
p = loadtxt(profile)
rect = Rectangle((409.5, 359.5), 587.5-409.5, 540.5-359.5, facecolor="#aaaaaa", alpha=0.4) # 1000
gca().add_patch(rect)

#interp = 'bilinear';
interp = 'nearest';
title("T = %d" %frame)
#imshow(p, origin='lower', interpolation=interp)
imshow(d, origin='lower', interpolation=interp)
colorbar()
show()

