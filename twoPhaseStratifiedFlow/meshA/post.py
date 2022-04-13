import sys, os, shutil, math, commands, glob 
import os.path, time
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pylab import *
from pygnuplot import gnuplot, pyplot
from subprocess import Popen
from subprocess import call

# input required ===============================================================
muA = 0.0005
muB = 1.85e-05
H = 0.02
L = 0.04
dp = 0.0021
endTime = 2.5
#===============================================================================

# change working directory to current working directory
testcase=os.getcwd()
os.chdir(testcase)  
          
# Numerical solution
## sample  
cmd='sample; '
pipefile = open('output', 'w')
retcode = call(cmd,shell=True,stdout=pipefile)
pipefile.close()
os.remove('output')

## Read samples from sets folder and save as a csv file
Ux_num=pd.read_csv(testcase+'/sets/'+str(endTime)+'/DataFile_U.xy', delimiter='\s+',header=None,names=["y","Ux","Uy","Uz"])

Ux_num.to_csv('Numerical.csv', index=False)

# exact solution
n1 = float(muA)
n2 = float(muB)
d = H/2.0
dpdx = dp/L

Ux_ex = []
for y in Ux_num['y']: #take same Y from sampled data
    denom = (H+d*(n2/n1-1.0))
    num1= ((H**2+(d**2)*(n2/n1-1))*y)
    Ux_ex.append(abs((1/(2*n1))*(dpdx)*(y**2-num1/denom)))

Ux_df=pd.DataFrame({'y':Ux_num['y'], 'Ux':Ux_ex})
Ux_df.to_csv('Analytical.csv', index=False)
        
# plot numerical only
plt.plot(abs(Ux_num['Ux']), Ux_num['y'], 'r--') 
plt.legend(["Ux"])
plt.title('Velocity Profile (Numerical)', fontsize=14)
plt.xlabel('horizontal velocity (m/s)', fontsize=12)
plt.ylabel('channel height (m)', fontsize=12)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.savefig('numProfile.png')
plt.close("all")

# plot Numerical vs Exact
plt.plot(abs(Ux_num['Ux']), Ux_num['y'], 'r--')
plt.plot(Ux_ex, Ux_num['y'],'k')
plt.legend(["Numerical", "Anlytical"])
plt.xlabel('horizontal velocity (m/s)', fontsize=12)
plt.ylabel('channel height (m)', fontsize=12)
plt.title('Velocity Profile', fontsize=14)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.savefig('profile.png')

