#!/usr/bin/env python
# -*- coding: utf-8 -*-
from pylab import *
import gobject

nsteps = 1000
fileout = 10
ncell = 400
m2start  = int(ncell/2)

fig = figure(figsize=(16.0,8.0))
grid(True)
xlabel('z')
ylabel('Field')
filename = "../output/out_%.06d.dat"%0
d = loadtxt(filename)
curve1 = plot(d[:,0],d[:,1],label="Ex")
curve2 = plot(d[:,0],d[:,2],label="Hy")
legend()
axvspan(m2start,ncell,facecolor='0.5', alpha=0.5)
axis([0, ncell, -1.5, 1.5])

n = 0
def updatefig(*args):
     global n
     title('step %d' % n)
     
     #### Open files and read their content ###############
     filename = "../output/out_%.06d.dat"%n
     d = loadtxt(filename)
     print filename
     
     #### Update data on graph ############################
     curve1[0].set_ydata(d[:,1])
     curve2[0].set_ydata(d[:,2])
     
     fig.canvas.draw()
     
     #### Jump to next step and check if done #############
     n += fileout
     if(n <= nsteps): return True
     else: return False
     
     #### end of def ######################################

######### Update until updatefig returns false ############
gobject.timeout_add(30,updatefig)
show()

######### End of script ###################################