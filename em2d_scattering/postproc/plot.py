#! /usr/bin/env python
# -*- coding: utf-8 -*-

from sys import argv
from pylab import *

d = loadtxt("./output/Ez.dat")

interp = 'bilinear';
#interp = 'nearest';
imshow(d, origin='lower', interpolation=interp)
colorbar()
show()

