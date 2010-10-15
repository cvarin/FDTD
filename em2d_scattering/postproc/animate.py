#! /usr/bin/env python
# -*- coding: utf-8 -*-

from sys import argv
from pylab import *
import gobject

frame = int(argv[1])
step = int(argv[2])

###########################################################
def load_file(n):
     filename = "./output/Ez_%06d.dat" %n
     return loadtxt(filename)

###########################################################
n = step
d = load_file(n)

fig = figure()
interp = 'bilinear';
title("T = %d" %n)
im = imshow(d, origin='lower', interpolation=interp, vmin = -1.0, vmax = 1.0)
colorbar()

######### Figure update ###################################
def updatefig(*args):
     global n
     title('T = %d' %n)
     
     ##### Update data on graph ###########################
     d = load_file(n)
     im.set_array(d)
     fig.canvas.draw()
     n += step
     if(n <= frame): return True
     else:
          print "Done"
          return False
  
###########################################################   
wait = 100
n = step
gobject.timeout_add(wait,updatefig)
show()

######### End of script ###################################


