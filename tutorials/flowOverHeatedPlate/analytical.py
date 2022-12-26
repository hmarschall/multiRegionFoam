import sys, os, shutil, math, commands, glob 
import os.path, time
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pylab import *
from pygnuplot import gnuplot, pyplot
from subprocess import Popen
from subprocess import call
from os import path
from PyFoam.RunDictionary.SolutionDirectory import SolutionDirectory
from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile

# input required ===============================================================
# none
#===============================================================================

# change working directory to current working directory
testcase=os.getcwd()
os.chdir(testcase)  

# read from dictionaries
case = SolutionDirectory(testcase)
controlDict = ParsedParameterFile(path.join(case.name,"system", "controlDict"))
endTime = controlDict["endTime"]  
          
# Numerical solution
## sample  
cmd='sample -region fluid; '
pipefile = open('output', 'w')
retcode = call(cmd,shell=True,stdout=pipefile)
pipefile.close()
os.remove('output')

## Read samples from sets folder and save as a csv file
T_num=pd.read_csv(testcase+'/postProcessing/sets/fluid/'+str(endTime)+'/DataFile_T.xy', delimiter='\s+',header=None,names=["x","T"])

## calculate non-dimentional temperature (theta) along the interface 
theta_num = []
for t in T_num['T']: 
    theta_num.append((t-300)/(310-300))
    
theta_df=pd.DataFrame({'x':T_num['x'], 'theta':theta_num})
theta_df.to_csv('theta.csv', index=False)

# Analytical solution
theta_anal=pd.read_csv(testcase+'/analyticalTheta.csv', delimiter=',',header=None,names=["x","theta_anal"])
        
# plot numerical only
plt.plot(theta_df['x'], theta_df['theta'], 'r--') 
plt.legend([r'$\theta$'])
plt.title('numerical non-dimensional temperature' + r' $\theta$', fontsize=14)
plt.xlabel('x', fontsize=12)
plt.ylabel(r'$\theta$', fontsize=12)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.savefig('theta_num.png')
plt.close("all")

# plot Numerical vs Exact
plt.plot(theta_df['x'], theta_df['theta'], 's') #'r') 
plt.plot(theta_anal['x'], theta_anal['theta_anal'], '.')#'k')
plt.legend(["Numerical", "Anlytical"])
plt.xlabel('x', fontsize=12)
plt.ylabel(r'$\theta$', fontsize=12)
plt.title('non-dimensional temperature' + r' $\theta$', fontsize=14)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.ylim([0, 1])
plt.savefig('compare_theta.png')

