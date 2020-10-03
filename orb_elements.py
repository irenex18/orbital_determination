from math import *
import numpy as np
#IMPORTANT - CHANGE THE SUN-EARTH VECTOR WHEN RUNNING CODE FOR YOUR ASTEROID
def ephemeris(e, a, I, OPrime, wPrime, M0, t, t0, sunVec):

    mu = 0.01720209895
    #calculate mean motion
    n = mu/((a)**3)**0.5 #mean motion
    print(n)
    #calculate M
    M = M0 + n*(t-t0)
    print(M)
    Eguess = M
    prevEguess = Eguess
    Eguess = M + e*sin(Eguess)
    while abs(Eguess-prevEguess) > 1.0e-4:
        prevEguess = Eguess
        Eguess = M + e*sin(Eguess)
    print(Eguess)
    E = Eguess
    #calculate the x and y coordinates using physics
    xCoord = a*cos(E) - e*a
    yCoord = a*sin(E)*((1-e**2)**0.5)
    print(yCoord)
    #determine r and nu
    r = (xCoord**2 + yCoord**2)**0.5
    nu = acos(xCoord/r)
    print(nu)
    zCoord = 0
    #calculate the ecliptic coordinates
    rcoords = np.array([[xCoord],[yCoord],[zCoord]])
    wArr = np.array([[cos(wPrime),-1*sin(wPrime),0],[sin(wPrime),cos(wPrime),0],[0,0,1]])
    iArr = np.array([[1,0,0],[0,cos(I),-1*sin(I)],[0,sin(I),cos(I)]])
    OArr = np.array([[cos(OPrime),-1*sin(OPrime),0],[sin(OPrime),cos(OPrime),0],[0,0,1]])

    #rotate vector by -wprime
    rcoords = np.array(np.dot(wArr, rcoords))
    #rotate vector by iprime
    rcoords =np.array(np.dot(iArr,rcoords))
    #rotate vector by Oprime
    rcoords = np.array(np.dot(OArr, rcoords))
    print(rcoords)
    #Earth to sun vector from JPL horizons
    Rvector = np.array([[-2.027873566936922e-1],[9.963238789875005e-1],[-4.453100906916791e-5]])
    #Find Earth to asteroid vector
    EtoAst = rcoords + Rvector
    #Equatorial coordinates 
    xEq = (EtoAst[0])
    obliquity = 23.43676*(pi/180)
    yEq = EtoAst[1]*cos(obliquity) - EtoAst[2]*sin(obliquity)
    zEq = EtoAst[1]*sin(obliquity) + EtoAst[2]*cos(obliquity)
    #normalize the vectors
    magnitude = (xEq**2 + yEq**2 + zEq**2)**0.5
    xEq = xEq/magnitude
    yEq = yEq/magnitude
    zEq = zEq/magnitude
    print(xEq)
    print(yEq)
    print(zEq)
    #calculate RA and dec
    decDecimal = asin(zEq)*(180/pi)
    decRad = asin(zEq)
    cosRA = (xEq/(cos(decRad)))
    sinRA = yEq / cos(decRad)
    RADecimal = atan2(sinRA, cosRA)*(180/pi)*(24/360)
    RAmin = (RADecimal - trunc(RADecimal))*60
    RAminStr = str(trunc(RAmin))
    if RAmin< 10.0:
        RAminStr = "0" + RAminStr
    RAsec = (RAmin - trunc(RAmin))*60
    RAsecStr = str(RAsec)
    if RAsec < 10.0:
        RAsecStr = "0"+RAsecStr
    RA = str(trunc(RADecimal)) + ":" + RAminStr + ":" + RAsecStr
    decArcMin = (decDecimal - trunc(decDecimal))*60
    decArcMinStr = str(trunc(decArcMin))
    if decArcMin<10.0:
        decArcMinStr = "0" + decArcMinStr
    decArcSec = (decArcMin - trunc(decArcMin))*60
    decArcSecStr = str(decArcSec)
    if decArcSec < 10.0:
        decArcSecStr = "0" + decArcSecStr
    dec = str(trunc(decDecimal)) + ":" + decArcMinStr + ":" + decArcSecStr
    #returns the dec and RA values
    return "dec: " + str(dec) + " RA: " + str(RA)   

#math and physics pset 5 test case: 46P/Wirtanen
e1 = 0.6587595515873473
a1 = 3.092704185336301 #AU
I1 = 11.74759129647092 * pi/180 #radians
OPrime1 = 82.15763948051409 *pi/180 #radians
wPrime1 = 356.34109239 * pi/180 #radians
M01 = 0.01246738682149958 * pi/180 #radians
t1 = 2458668.5 #JHD
t01 = 2458465.5 #JHD
print (ephemeris(e1, a1, I1, OPrime1, wPrime1, M01, t1, t01))
