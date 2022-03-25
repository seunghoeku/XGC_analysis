import adios2
import sys

print(sys.argv)
print(sys.argv[1])
print('hi there\n')

if len(sys.argv) != 4 :
    print('usage: ', sys.argv[0], ' infile outfile skip')
    sys.exit()


inFile = sys.argv[1]
out = sys.argv[2]
skip = int(sys.argv[3])

outRZ = out + '.txt'
outTP = out + '.TP.txt'

f=adios2.open(inFile, 'r')
rz=f.read('RZ')
tp=f.read('ThetaPsi')
ID=f.read('ID')

n = int(len(rz) / 2)
print('n= ', n)
#print(rz)

frz = open(outRZ, 'w')
ftp = open(outTP, 'w')
frz.write('ID, R, Z\n')
ftp.write('ID, THETA, PSI\n')

for i in range(0,n,skip) :
    id = ID[i]
    #print('%d: id=%d\n' % (i, id))
    if id < 0 : continue

    frz.write('%d, %lf, %lf\n' % (id, rz[i*2+0], rz[i*2+1]))
    ftp.write('%d, %lf, %lf\n' % (id, tp[i*2+0], tp[i*2+1]))

frz.close()
ftp.close()
