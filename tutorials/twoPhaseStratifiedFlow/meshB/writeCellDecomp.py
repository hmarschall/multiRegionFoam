import os

path = os.getcwd()
interfaceCells = open(path+"/constant/polyMesh/sets/interfaceCells", "r") 
i=0
cellNumber = 0
cells = []
for line in interfaceCells: 
    if i==18:
        cellNumber = [int(x) for x in line.split()]
    i+=1
    if (i >20) and (i<=20+cellNumber[0]):
        cell = [int(x) for x in line.split()][0]
        cells.append(cell)
interfaceCells.close()
cells.sort()

cellCo=0
lineCo=0
procList=[]
with open(path+"/constant/cellDecomposition","r") as cellDecomposition:
    for line in cellDecomposition:
        if cellCo < len(cells):
            if (lineCo-20) == cells[cellCo]:
                cellCo+=1
                procList.append('0\n')
            else:
                procList.append(line)
            lineCo +=1
        else:
            procList.append(line)
            
with open(path+"/constant/cellDecomposition2", 'w') as f:
    for item in procList:
        f.write("%s" % item)
