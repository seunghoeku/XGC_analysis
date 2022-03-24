
def diffMag(a, b) :
    d = (a[0]-b[0], a[1]-b[1], a[2]-b[2])
    return d[0]*d[0] + d[1]*d[1] + d[2]*d[2]

traces = open('traces.txt', 'r').readlines()
tracesXGC = open('../data/sku_8000/POINC/xgc_poinc.txt', 'r').readlines()

N = min(len(traces), len(tracesXGC))

print N

outF = open('diff.txt', 'w')
outF.write('IDX, DIFF\n')

cnt = 0
for i in range(N-1) :
    if i == 0 : continue

    b = traces[i].split(',')
    xgc = tracesXGC[i].split(',')
    bRZP = (float(b[1]), float(b[2]), float(b[3]))
    xgcRZP = (float(xgc[1]), float(xgc[2]), float(xgc[3]))
    diff = diffMag(bRZP, xgcRZP)

    if (diff > 5) : continue

    outF.write('%d, %f\n' % (cnt, diff))
    cnt = cnt+1
