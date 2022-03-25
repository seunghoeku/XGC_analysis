input = open('RZ.punc.theta_psi.txt', 'r').readlines()
input = open('../data/sku_8000/xgc_theta_psi.txt', 'r').readlines()

tout = open('RZ.trunc.theta_psi.txt', 'w')
tout.write("ID, THETA, PSI, DUMMY\n")

N = 325
N = 190

id0 = -1
id = -1
cnt = 0
for str in input[1:] :
    x = str.split(',')
    id = int(x[0])
    cnt = cnt+1
    if id0 < id : cnt = 0
    id0 = id
    if cnt <= N :
        tout.write('%s' % str)
