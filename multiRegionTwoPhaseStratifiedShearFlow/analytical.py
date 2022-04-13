import sys, os, math, glob 
import os.path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
#from pygnuplot import gnuplot, pyplot

# input required ===============================================================
muA = 5e-04
muB = 1.85e-05
H = 0.02
L = 0.04
dp = 0.0021
#===============================================================================

# change working directory to current working directory
testcase=os.getcwd()
os.chdir(testcase)  
          
# Numerical solution
## sample  
script = """sample -latestTime -region fluidA; sample -latestTime -region fluidB;"""
os.system("bash -c '%s'" % script)

latestTimeA = max(glob.glob(testcase+'/postProcessing/sets/fluidA/*') , key=os.path.getctime)
latestTimeB = max(glob.glob(testcase+'/postProcessing/sets/fluidB/*') , key=os.path.getctime)

## combine samples into one dataFrame
Ux1_num=pd.read_csv(latestTimeA+'/DataFile_U.xy', delimiter='\s+',header=None,names=["y","Ux","Uy","Uz"])
Ux2_num=pd.read_csv(latestTimeB+'/DataFile_U.xy', delimiter='\s+',header=None,names=["y","Ux","Uy","Uz"])

# exact solution
n1 = float(muA)
n2 = float(muB)
d = H/2.0
dpdx = dp/L

## meshA
Ux1_ex = []
for y in Ux1_num['y']: #take same Y from sampling data
    denom = (H+d*(n2/n1-1.0))
    num1= ((H**2+(d**2)*(n2/n1-1))*y)
    Ux1_ex.append(abs((1/(2*n1))*(dpdx)*(y**2-num1/denom)))

## meshB
Ux2_ex = []
for y in Ux2_num['y']:
    denom = (H+d*(n2/n1-1))
    num1= ((H**2+(d**2)*(n2/n1-1))*y)    
    num2= ((n2/n1)-1)*(H*(d**2)-(H**2)*d)
    Ux2_ex.append(abs((1/(2*n2))*(dpdx)*((y**2)-(num1/denom)+num2/denom)))
        
Ux1_df=pd.DataFrame({'y':Ux1_num['y'], 'Ux':Ux1_ex})
Ux2_df=pd.DataFrame({'y':Ux2_num['y'], 'Ux':Ux2_ex})
Ux_df= Ux1_df.append(Ux2_df)
        
# plot numerical (both meshes with different colors)
plt.plot(abs(Ux1_num['Ux']), Ux1_num['y'], 'g--') 
plt.plot(abs(Ux2_num['Ux']), Ux2_num['y'],'r--')
plt.legend(["Ux1 (meshA)", "Ux2 (meshB)"])
plt.title('Numerical Velocity Profile', fontsize=14)
plt.xlabel('horizontal velocity (m/s)', fontsize=12)
plt.ylabel('channel height (m)', fontsize=12)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.savefig('numProfile.png')
plt.close("all")

# plot Numerical vs Exact
plt.plot(abs(Ux1_num['Ux']), Ux1_num['y'], 'r--')
plt.plot(Ux1_ex, Ux1_num['y'],'k')
plt.plot(abs(Ux2_num['Ux']), Ux2_num['y'],'r--')
plt.plot(Ux2_ex, Ux2_num['y'], 'k')
plt.legend(["Numerical", "Anlytical"])
plt.xlabel('horizontal velocity (m/s)', fontsize=12)
plt.ylabel('channel height (m)', fontsize=12)
plt.title('Velocity Profile', fontsize=14)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.savefig('numVsAnalProfile.png')
