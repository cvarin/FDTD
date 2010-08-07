#!/usr/bin/env python
# -*- coding: utf-8 -*-
from pylab import *
import gobject
from sys import argv
from matplotlib.widgets import Slider

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
nsteps = 1000
step = 10
ncell = 400
m2start  = int(ncell/2)

#########################################################################################
# The figure is created with frame 0
fig = figure(figsize=(16.0,8.0))
grid(True)
xlabel('z')
ylabel('Field')
d = loadtxt(filename(0))
curve1 = plot(d[:,0],d[:,1],label="Ex")
curve2 = plot(d[:,0],d[:,2],label="Hy")
legend()
axvspan(m2start,ncell,facecolor='0.5', alpha=0.5)
axis([0, ncell, -1.5, 1.5])

#########################################################################################
# Animation widgets
axcolor = 'lightgoldenrodyellow'
axdef  = axes([0.15, 0.02, 0.7, 0.02], axisbg=axcolor)
progress = Slider(axdef,'', 0, nsteps, valinit=0)

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
  curve1[0].set_ydata(d[:,1])
  curve2[0].set_ydata(d[:,2])
  fig.canvas.draw()

  #######################################################################################
  # Update n
  playbutton.on_clicked(play_pause)
  if(run): n += step
  else:
    val = int(progress.val)
    n = val - val%10
  if(n <= nsteps): 
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