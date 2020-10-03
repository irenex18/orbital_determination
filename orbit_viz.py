from vpython import *

a = 2.773017979589484
e = 0.1750074901308245
M = radians(336.0050001501443)
Oprime = radians(108.032597191534)
iprime = radians(16.34548466739393)
wprime = radians(74.95130563682554)

def solvekep(M):
    Eguess = M
    Mguess = Eguess - e*sin(Eguess)
    Mtrue = M
    while abs(Mguess - Mtrue) < 1e-004:
        Mguess = Eguess - e*sin(Eguess)
        Eguess = Eguess - (Eguess - e*sin(Eguess) - Mtrue) / (1 - e*cos(Eguess))
    return Eguess

sqrtmu = 0.01720209895
mu = sqrtmu**2
time = 0
dt = .05
period = sqrt(4*pi**2*a**3/mu)
r1ecliptic = vector(0, 0, 0)
Mtrue = 2*pi/period*(time) + M
Etrue = solvekep(Mtrue)
r1ecliptic.x = (cos(wprime)*cos(Oprime) - sin(wprime)*cos(iprime)*sin(Oprime))*(a*cos(Etrue)-a*e) - (cos(wprime)*cos(iprime)*sin(Oprime) + sin(wprime)*cos(Oprime))*(a*sqrt(1-e**2)*sin(Etrue))
r1ecliptic.y = (cos(wprime)*sin(Oprime) + sin(wprime)*cos(iprime)*cos(Oprime))*(a*cos(Etrue)-a*e)) + (cos(wprime)*cos(iprime)*cos(Oprime) - sin(wprime)*sin(Oprime))*(a*sqrt(1-e**2)*sin(Etrue))
r1ecliptic.z = sin(wprime)*sin(iprime)*(a*cos(Etrue)-a*e) + cos(wprime)*sin(iprime)*(a*sqrt(1-e**2)*sin(Etrue))
asteroid = sphere(pos=r1ecliptic*150, radius=(15), color=color.white)
asteroid.trail = curve(color=color.white)
sun = sphere(pos=vector(0,0,0), radius=(50), color=color.yellow)

while (1==1):
    rate(200)
    time = time + 1
    Mtrue = 2*pi/period*(time) + M
    Etrue = solvekep(Mtrue)
    r1ecliptic.x = (cos(wprime)*cos(Oprime) - sin(wprime)*cos(iprime)*sin(Oprime))*(a*cos(Etrue)-a*e) - (cos(wprime)*cos(iprime)*sin(Oprime) + sin(wprime)*cos(Oprime))*(a*sqrt(1-e**2)*sin(Etrue))
    r1ecliptic.y = ((cos(wprime)*sin(Oprime) + sin(wprime)*cos(iprime)*cos(Oprime))*(a*cos(Etrue)-a*e)) + (cos(wprime)*cos(iprime)*cos(Oprime) - sin(wprime)*sin(Oprime))*(a*sqrt(1-e**2)*sin(Etrue))
    r1ecliptic.z = sin(wprime)*sin(iprime)*(a*cos(Etrue)-a*e) + cos(wprime)*sin(iprime)*(a*sqrt(1-e**2)*sin(Etrue))
asteroid.pos = r1ecliptic*150
asteroid.trail.append(pos=asteroid.pos)  
