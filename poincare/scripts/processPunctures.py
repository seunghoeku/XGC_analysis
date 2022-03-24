import string, math

eq_axis_r = 2.8
eq_axis_z = 0.0

punctures = open('punctures.theta_psi.500.txt', 'r').readlines()

thetaPsi = [[], [], [], [], [], []]
for p in punctures[1:] :
    x = string.split(p, ',')
    id = int(x[0])
    theta = float(x[1])
    psi = float(x[2])
    entry = (theta, psi)
    thetaPsi[id].append(entry)


pout = open('p.theta.psi.txt', 'w')
pout.write('ID, Theta, Psi, Dummy\n')

for id in range(len(thetaPsi)) :
    for pt in thetaPsi[id] :
        theta = pt[0] + math.pi
        psi = pt[1]
        pout.write('%d, %f, %f, 0\n' % (id, theta, psi))

print len(thetaPsi[0])
