#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
     Python script for plotting the results from Plane_Wave
"""
######### Module imports ##################################
from numpy import arange, array, reshape, squeeze, zeros
from pylab import *
import gobject
import re
from mSetFigure import *
from mIO import *
from mField import *

######### Parse the command line arguments ################
cmd = Parse_command_line_input(sys.argv)
start,stop,step = cmd.start,cmd.stop,cmd.step
inputfile = cmd.inputfile

###### Get some information from the simulation inputfile #
basename = Get_io_basename(inputfile)
I        = Get_attribute(inputfile,"EMsource","I","double")
nxyz     = Get_vector_attribute(inputfile,"space","nxyz","int")
nabc     = Get_attribute(inputfile,"boundaries","nabc","int")
if Node_is_present(inputfile,"conducting_plane"):
     cp_position = Get_attribute(inputfile,"conducting_plane","position","int")
     cp_thickness = Get_attribute(inputfile,"conducting_plane","thickness","int")
else:
     cp_position = 0
     cp_thickness = 0

######### Set and evaluate some derived quantities ########
asize = array([nxyz[0]+1, nxyz[1]+1, nxyz[2] + 2*nabc]) # x,y,z
E0 = sqrt(2.0*120.0*pi*1.0e4*I)
Enorm = 1.0/E0

######### Arrays creation #################################
Ex3D = Read_component_in_file_and_convert('x',0,asize,basename)
t = arange(start,stop+step,step)
Et = zeros(size(t))
x,y,z = asize/2.0

######### Graph update function ###########################
print 'Preparing field plot'
fig = set_figure()

######### Subplots ########################################
subd = fig.add_subplot(4,1,4)
grid(True)
xlabel('t')
ylabel('Energy (J)')
em = loadtxt(basename + "_EM_energy.log")
energy = plot(em[:,0],em[:,1])

###########################################################
subc = fig.add_subplot(4,1,3)
grid(True)
ylabel('Ex')
pulse = plot(t,Et*Enorm,label='Field at cell (%d,%d,%d)' %(x,y,z))
legend(loc='upper right')
axis([start,stop, -1, 1])

###########################################################
subb = fig.add_subplot(4,1,2)
grid(True)
xlabel('z')
ylabel('Ex')
curve = plot(Ex3D[x,y,:]*Enorm)
axvspan(0.0,nabc, facecolor='0.5', alpha=0.5)
axvspan(nabc+nxyz[2],asize[2],facecolor='0.5', alpha=0.5)
annotate('PML', xy=(0.05,0.825),  xycoords='axes fraction')
annotate('PML', xy=(0.925,0.825),  xycoords='axes fraction')
axvspan(cp_position,cp_position + cp_thickness, facecolor='0.5', alpha=0.3)
annotate('Medium', xy=(cp_position + 1,0.6),  xycoords='data')
axis([0, asize[2]-1, -1, 1])

###########################################################
suba = fig.add_subplot(4,1,1)
minmax = 1.0e-0;
ylabel('x')
im = imshow(squeeze(Ex3D[:,y,:]*Enorm),origin='lower',\
          vmin=-minmax, vmax=minmax,)
annotate('PML', xy=(0.05,0.825),  xycoords='axes fraction')
annotate('PML', xy=(0.925,0.825),  xycoords='axes fraction')
axvspan(0.0,nabc, facecolor='0.5', alpha=0.5)
axvspan(nabc+nxyz[2],asize[2]-1,facecolor='0.5', alpha=0.5)
axvspan(cp_position,cp_position + cp_thickness, facecolor='0.5', alpha=0.3)
annotate('Medium', xy=(cp_position + 1,8),  xycoords='data')
axis([0, asize[2]-1, 0, asize[0] - 1])

######### Keyboard events #################################
def press(event):
     global n,run,reverse
     
     #### Play, pause #####################################
     if(event.key==' '):
          if(run==False): 
               run = True
          else: run = False
          
     #### Change direction ################################
     if(event.key=='r'):
          if(reverse==False): 
               reverse = True
          else: reverse = False
     
     #### Next, previous ##################################
     if(event.key=='n'):
          run = False
          n += step
          updatefig()
     if(event.key=='p'):
          run = False
          n -= step
          updatefig()
     if event.key=='g':
          run = False
          next_image = raw_input("\nEnter the frame number: ")
          if re.match("^[-+]?[0-9]+$", next_image):
               frame = int(next_image)
               if(frame <= stop and frame >= start):
                    if(frame%step == 0):
                         n = int(next_image)
                         updatefig()
                    else: print "Frame does not exist!"
               else: print "Frame is out of range."

######### Figure update ###################################
def updatefig(*args):
     global n
     #if n%int(stop/10) == 0: print 'step %d of %d' %(n,stop)
     title('step %d' % n)
     
     #### Read field and fetch the pulse shape ############
     Ex3D = Read_component_in_file_and_convert('x',n,asize,basename)
     i = n/step
     Et[i] = Ex3D[x,y,z]
     
     #### Update data on graph ############################
     im.set_array(squeeze(Ex3D[:,y,:]*Enorm))
     #curve[0].set_ydata(1.0 + Ex3D[x,y,:]*Enorm - 1.0)
     curve[0].set_ydata(Ex3D[x,y,:]*Enorm)
     pulse[0].set_ydata(Et*Enorm)
     fig.canvas.draw()

     # Version plus simple sans entrée au clavier
     if(run): n += step
     if(reverse): n -= step
     if(n <= stop or n >= start): return True
     else:
          print "Done"
          return False
  
     #### end of def ######################################
n = start
run = True
reverse = False
updatefig()
fig.canvas.mpl_connect('key_press_event', press)
if(run): gobject.idle_add(updatefig)
#wait = 50
#gobject.timeout_add(wait,updatefig)
show()

######### End of script ###################################
