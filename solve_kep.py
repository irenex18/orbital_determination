from vpython import *
import numpy as np
import math

a = 2.773017979589484
e = 0.1750074901308245
M = radians(336.0050001501443)
Oprime = radians(108.032597191534)
iprime = radians(16.34548466739393)
wprime = radians(74.95130563682554)

#get E by inverting the Kepler equation
def solvekep(M):
    Mtrue = M
    Eguess = M
    Mguess = Eguess - math.e*math.sin(Eguess)
    while abs(Mguess - Mtrue) > 1e-4:
        Eguess = Eguess - (Mtrue - (Eguess - math.e*math.sin(Eguess)))/(math.e*math.cos(Eguess) - 1)
        Mguess = Eguess - math.e*math.sin(Eguess)
    return Eguess
#call the method above in future code/give it a variable
x = a*math.cos(solvekep(M)) - a*math.e
y = (a*(1-(math.e)**.5)**.5)*math.sin(solvekep(M))
z = 0
r1ecliptic = vector(0,0,0)
arr = np.array([[],[y],[z]])
wArr = np.array([[math.cos(wprime),-1*math.sin(wprime),0][math.sin(wprime),math.cos(wprime),0][0,0,1]])
iArr = np.array([[1,0,0][0,math.cos(iprime),-1*math.sin(iprime)][0,math.sin(iprime),math.cos(iprime)]])
OArr = np.array([[math.cos(Oprime),-1*math.sin(Oprime),0][math.sin(Oprime),math.cos(Oprime),0][0,0,1]])

rotate vector by -wprime
finalV1 = np.array(np.dot(wArr, r1ecliptic))
#rotate vector by iprime
finalV2 =np.array(np.dot(iArr,finalV1))
#rotate vector by Oprime
finalV3 = np.array(np.dot(OArr, finalV2))


sqrtmu = 0.01720209895
mu = sqrtmu**2
time = 0
dt = .05
period = sqrt(4*math.pi**2*a**3/mu)
Mtrue = 2*math.pi/period*(time) + M
Etrue = solvekep(Mtrue)
asteroid = sphere(pos=r1ecliptic*150, radius=(15), color=color.white)
asteroid.trail = curve(color=color.white)
sun = sphere(pos=vector(0,0,0), radius=(50), color=color.yellow)

while (True):
    rate(200)
    time = time + 1
    Mtrue = 2*pi/period*(time) + M
    Etrue = solvekep(Mtrue)
    r1ecliptic.x = (cos(wprime)*cos(Oprime) - sin(wprime)*cos(iprime)*sin(Oprime))*(a*cos(Etrue)-a*e) - (cos(wprime)*cos(iprime)*sin(Oprime) + sin(wprime)*cos(Oprime))*(a*sqrt(1-e**2)*sin(Etrue))
    r1ecliptic.y = ((cos(wprime)*sin(Oprime) + sin(wprime)*cos(iprime)*cos(Oprime))*(a*cos(Etrue)-a*e)) + (cos(wprime)*cos(iprime)*cos(Oprime) - sin(wprime)*sin(Oprime))*(a*sqrt(1-e**2)*sin(Etrue))
    r1ecliptic.z = sin(wprime)*sin(iprime)*(a*cos(Etrue)-a*e) + cos(wprime)*sin(iprime)*(a*sqrt(1-e**2)*sin(Etrue))
    asteroid.pos = r1ecliptic*150
    asteroid.trail.append(pos=asteroid.pos)  
