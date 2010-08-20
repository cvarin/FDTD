#!/usr/bin/env python
# -*- coding: utf-8 -*-
from pylab import *
import gobject
from sys import argv
from matplotlib.widgets import Slider

eta0 = 376.730313461

#########################################################################################
# Get the data directory from the command line and fetch the data file names
if size(argv) == 1: 
  print "Please provide a data directory."
  exit()
else:
  print "Data wil be read from %s" %argv[1]
  def filename(n): return argv[1] + "out_%.06d.dat"%n

#########################################################################################
# Simulation parameters
import os
class read_input_file: 
     def __init__(self, filename):
          if not os.path.exists(filename): 
               print "Can't find " + filename
               exit()
          infile = open(filename)
          for line in infile:
               if len(line) is not 1:
                    array = line.split()
                    if not array[0].find("nsteps"): self.nsteps = int(array[2])
                    if not array[0].find("step"): self.step = int(array[2])
                    if not array[0].find("ncell"): self.ncell = int(array[2])
                    if not array[0].find("m2start"): self.m2start = int(array[2])
                    if not array[0].find("m2stop"): self.m2stop = int(array[2])
                    if not array[0].find("I"): 
                      self.I = float(array[2])
                      self.E0 = sqrt(2.0*eta0*self.I*1.0e4)
          infile.close()

f = read_input_file(argv[1] + "input.txt")
normE = 1.0/f.E0
normH = eta0*normE

#########################################################################################
# The figure is created with frame 0
fig = figure(figsize=(16.0,8.0))
grid(True)
xlabel('z')
ylabel('Field (%.2e V/m)'%f.E0)
d = loadtxt(filename(0))
curve1 = plot(d[:,0],d[:,1]*normE,label=r"$E_x$")
curve2 = plot(d[:,0],d[:,2]*normH,label=r"$\eta_0H_y$")
legend()
axvspan(f.m2start,f.m2stop,facecolor='0.5', alpha=0.5)
axis([0, f.ncell, -1.5, 1.5])

#########################################################################################
# Animation widgets
axcolor = 'lightgoldenrodyellow'
axdef  = axes([0.15, 0.02, 0.7, 0.02], axisbg=axcolor)
progress = Slider(axdef,'', 0, f.nsteps, valinit=0)

resetax = axes([0.01, 0.02, 0.05, 0.1])
playbutton = Button(resetax, '>||', color=axcolor, hovercolor='0.975')
def play_pause(event):
  global run
  run = not run

#########################################################################################
# Function to update the figure
def updatefig(*args):
  global n

  #######################################################################################
  # Open file, read the content, and update data on the graph
  d = loadtxt(filename(n))
  curve1[0].set_ydata(d[:,1]*normE)
  curve2[0].set_ydata(d[:,2]*normH)
  fig.canvas.draw()

  #######################################################################################
  # Update n
  playbutton.on_clicked(play_pause)
  if(run): n += f.step
  else:
    val = int(progress.val)
    n = val - val%10
  if(n <= f.nsteps): 
    progress.set_val(n)
    return True
  else: return False

#########################################################################################
# Do the animation
n = 0
run = True
if(run): gobject.timeout_add(30,updatefig)
show()

######### End of script #################################################################