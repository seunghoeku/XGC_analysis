import adios2
import sys
import matplotlib.pyplot as plt
import matplotlib.tri as tri

print(sys.argv)
print(sys.argv[1])

if len(sys.argv) != 4 :
    print('usage: ', sys.argv[0], ' infile outfile skip')
    sys.exit()


inFile = sys.argv[1]
outFile = sys.argv[2]
skip = int(sys.argv[3])

f=adios2.open(inFile, 'r')
tp=f.read('ThetaPsi')
rz=f.read('RZ')
ID=f.read('ID')
TS=f.read('TimeStep')

n = int(len(tp) / 2)

fig = plt.figure(figsize=[8,12])
id0 = -1
x = []
y = []

cnt = 0
for i in range(0, n, skip) :
    id = ID[i]
    if id != id0 and len(x) > 0 :
        plt.scatter(x, y, s=1, marker='x', edgecolor='none')
        x = []
        y = []
        id0 = id
        cnt = cnt+1

    x.append(tp[i*2+0])
    y.append(tp[i*2+1])

plt.scatter(x, y, s=1, marker='x', edgecolor='none')
print('num field lines= %d' % cnt)

plt.title('%s ts=%d' % (inFile, TS[0]))
plt.savefig(outFile + '.' + str(TS[0]) + '.TP.png', format='png')
plt.show()

#Now, plot the R,Z points
fig = plt.figure(figsize=[12,12])
id0 = -1
x = []
y = []

cnt = 0
for i in range(0, n, skip) :
    id = ID[i]
    if id != id0 and len(x) > 0 :
        plt.scatter(x, y, s=1, marker='x', edgecolor='none')
        x = []
        y = []
        id0 = id
        cnt = cnt+1

    x.append(rz[i*2+0])
    y.append(rz[i*2+1])

plt.scatter(x, y, s=1, marker='x', edgecolor='none')
print('num field lines= %d' % cnt)

plt.title('%s ts=%d' % (inFile, TS[0]))
plt.savefig(outFile + '.' + str(TS[0]) + '.RZ.png', format='png')
