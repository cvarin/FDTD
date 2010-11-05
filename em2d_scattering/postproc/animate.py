#! /usr/bin/env python
# -*- coding: utf-8 -*-

from sys import argv
from pylab import *
import gobject

frame = int(argv[1])
step = int(argv[2])

profile = "./output/Ga_%06d.dat" %0
p = loadtxt(profile)

###########################################################
def load_file(n):
     filename = "./output/Ez_%06d.dat" %n
     return loadtxt(filename)

###########################################################
n = step
d = load_file(n)

minmax = 2.0

fig = figure(figsize=(8.5,8.5))
interp = 'bilinear';
title("T = %d" %n)
#subplot(2,1,1)
im = imshow(d, origin='lower', interpolation=interp, vmin = -minmax, vmax = minmax)
grid(True)

#rect = Rectangle((0, 110), 340, 20, facecolor="#aaaaaa", alpha=0.3)
#gca().add_patch(rect)

xlabel("x")
ylabel("y")
#subplot(2,1,2)
#curve = plot(d[170,:])
#axis([0, 340, -1.0, 1.0])
colorbar()

######### Figure update ###################################
def updatefig(*args):
     global n
     title('T = %d' %n)
     
     ##### Update data on graph ###########################
     d = load_file(n)
     im.set_array(d)
     fig.canvas.draw()
     #curve[0].set_ydata(d[170,:])
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


