import adios2
import sys
import matplotlib.pyplot as plt
import matplotlib.tri as tri

print(sys.argv)
print(sys.argv[1])
print('hi there\n')

if len(sys.argv) != 3 :
    print('usage: ', sys.argv[0], ' infile skip')
    sys.exit()


inFile = sys.argv[1]
skip = int(sys.argv[2])

f=adios2.open(inFile, 'r')
rz=f.read('RZ')
tp=f.read('ThetaPsi')
ID=f.read('ID')

n = int(len(rz) / 2)
print('n= ', n)

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

plt.title('%s' % inFile)
plt.show()
