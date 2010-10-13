#! /usr/bin/env python
# -*- coding: utf-8 -*-

from sys import argv
from pylab import *

frame = int(argv[1])
filename = "./output/Ez_%06d.dat" %frame
d = loadtxt(filename)

interp = 'bilinear';
#interp = 'nearest';
title("T = %d" %frame)
imshow(d, origin='lower', interpolation=interp)
colorbar()
show()

