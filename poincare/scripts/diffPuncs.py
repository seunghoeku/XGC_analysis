import math
import sys

which = int(sys.argv[1])
print 'Doing pid= ', which

def parsePt(x) :
    p = x.split(',')
    pt = (float(p[1]), float(p[2]))
    return pt

def diffMagPsi(a, b) :
    d = (a[1]-b[1], a[1]-b[1])
    return math.sqrt(d[0]*d[0] + d[1]*d[1])

def diffMag(a, b) :
    d = (a[0]-b[0], a[1]-b[1])
    return math.sqrt(d[0]*d[0] + d[1]*d[1])
#    d = b[0]-b[1]
#    return math.sqrt(d*d)


puncs = open('bumm.punc.txt', 'r').readlines()
puncs_tp = open('bumm.punc.theta_psi.txt', 'r').readlines()
xgcFile = '../data/sku_8000/xgc_punctures.%d.txt'%which
xgcFile_tp = '../data/sku_8000/xgc_punctures_tp.%d.txt'%which
print 'reading: ',  xgcFile
puncsXGC = open(xgcFile, 'r').readlines()
puncsXGC_tp = open(xgcFile_tp, 'r').readlines()

N = min(len(puncs), len(puncsXGC))
print 'Num samples= ', N

outF = open('diff.txt', 'w')
outF.write('IDX, ERR\n')
outF_tp = open('diff_tp.txt', 'w')
outF_tp.write('IDX, ERR_tp\n')

outF_p = open('diff_p.txt', 'w')
outF_p.write('IDX, ERR_PSI\n')

cnt = 0
errSum = 0
maxErr = 0
for i in range(N-1) :
    if i == 0 : continue

    p = parsePt(puncs[i])
    pxgc = parsePt(puncsXGC[i])
    err = diffMag(p, pxgc)
    if err > maxErr : maxErr = err
    errSum = errSum + err
    err = math.log(err, 2.0)

    outF.write('%d, %f\n' % (cnt, err))

    p = parsePt(puncs_tp[i])
    pxgc = parsePt(puncsXGC_tp[i])
    err = diffMag(p, pxgc)
    err = math.log(err, 2.0)
    outF_tp.write('%d, %f\n' % (cnt, err))

    ##psi only
    p = (p[1],p[1])
    pxgc = (pxgc[1],pxgc[1])
    err = diffMagPsi(p, pxgc)
    err = math.log(err, 2.0)
    outF_p.write('%d, %f\n' % (cnt, err))

    cnt = cnt+1

print 'Total error= ', errSum, 'maxErr=', maxErr
