import matplotlib.pyplot as plt
import pandas as pd

file='log.checkMesh.analyzed/meshvolume'
percent = 1

df = pd.read_csv(file, sep='\t', names=['time','volume'], header=0)

lowBound = (1- percent/100)*df['volume'][0]
upBound = (1+ percent/100)*df['volume'][0]

plt.plot(df['time'],df['volume'])
plt.xlabel('time')
plt.ylabel('bubble volume')
ymin = min(0.999*lowBound,df['volume'].min())
ymax = max(1.001*upBound,df['volume'].max())
plt.axhline(y=lowBound,color='r')
plt.axhline(y=upBound,color='r')
plt.ylim(ymin,ymax)
plt.show()
