from math import*
import numpy as np
import numpy.polynomial.polynomial as poly

def positiveAng(angle):
    if angle < 0:
        angle = angle + 360
    return angle

def orbElements(r2EclV, r2dotEclV, t2, sunVec):
    #get a
    r2SquaredSum = 0
    for i in range(len(r2EclV)):
        r2SquaredSum = r2SquaredSum + r2EclV[i]**2
    rMag = sqrt(r2SquaredSum)
    a = ((2.0/rMag)-np.dot(r2dotEclV, r2dotEclV))**-1
    #get e
    r2rdot2Crossed = np.cross(r2EclV, r2dotEclV)
    r2rdot2CrossedMag = 0
    for i in range(len(r2dotEclV)):
        r2rdot2CrossedMag = r2rdot2CrossedMag + r2rdot2Crossed[i]**2
    r2rdot2CrossedMag = sqrt(r2rdot2CrossedMag)
    e = (1-r2rdot2CrossedMag**2/a)**0.5
    #convert coordinates to ecliptic (not necessary because the coordinates from the Method of Gauss are alreadyin ecliptic
    r2Ecl = 0
    for i in range(len(r2EclV)):
        r2Ecl = r2Ecl + r2EclV[i]**2
    r2Ecl = sqrt(r2Ecl)
    #get I in degrees
    hVector = np.cross(r2EclV,r2dotEclV)
    h = 0
    for i in range(len(hVector)):
        h = h + hVector[i]**2
    h = sqrt(h)
    hz = r2EclV[0]*r2dotEclV[1]-r2EclV[1]*r2dotEclV[0]
    I = acos(1.0*hz/h)
    IDeg = positiveAng(I*(180/pi))
    #get Oprime in degrees
    hx = r2EclV[1]*r2dotEclV[2]-r2EclV[2]*r2dotEclV[1]
    hy = r2EclV[2]*r2dotEclV[0]-r2EclV[0]*r2dotEclV[2]
    sinOprime = hx/(h*sin(I))
    cosOprime = hy/(-h*sin(I))
    Oprime = atan2(sinOprime, cosOprime)
    OprimeDeg = positiveAng(Oprime*(180/pi))
    #get w in degrees (w = U-v)
    #get U
    cosU = (r2EclV[0]*cos(Oprime) + r2EclV[1]*sin(Oprime))/r2Ecl
    sinU = r2EclV[2]/(r2Ecl*sin(I))
    U = atan2(sinU, cosU)
    #get nu
    cosNu = ((a*(1-e**2))/r2Ecl - 1)/e
    rdprdot = np.dot(r2EclV, r2dotEclV)
    sinNu = (((a*(1-e**2))/h)*(rdprdot/r2Ecl))/e
    nu = atan2(sinNu, cosNu)
    #get w
    w = U-nu
    wDeg = positiveAng(w*(180/pi))
    #get M in degrees
    E = acos((1/e)*(1-(r2Ecl/a)))
    if nu < 0:
        nu = nu + 2*pi
    if nu > pi and nu < 2*pi:
        E = 2*pi-E
    M2 = E-e*sin(E)
    M2Deg = M2*(180/pi)
    #get n
    mu = 0.01720209895
    n = mu/((a)**3)**0.5
    #precess the true anomaly
    tOriginal = [2019, 7, 21, 6, 00, 00]
    J0 = 367*tOriginal[0] - int((7*(tOriginal[0]+int((tOriginal[1]+9)/12)))/4) + int((275*tOriginal[1])/9) + tOriginal[2] + 1721013.5
    TDueJD = J0 + tOriginal[3]/24
    M2Precessed = M2 + n*(TDueJD-t2)
    M2PrecessedDeg = positiveAng(M2Precessed*(180/pi))
    #get E
    EDeg = positiveAng(E*(180/pi))
    #get T
    T = t2 - (M2/n)
    #get P
    P = ((2*pi)/n)/365.25
    RV = sunVec
    orbElements = np.array([a, e, IDeg, OprimeDeg, wDeg, M2Deg, t2, M2PrecessedDeg])
    return orbElements

##r2 = np.array([0.13216218, -1.23287955, -0.32258829])
##rdot2 = np.array([0.96835101, 0.13676487, -0.15293932])
##t2 = 2458304.74796
##print(orbElements(r2, rdot2, t2))

def ephemeris(e, a, I, OPrime, wPrime, M0, t, t0, sunVec):
    mu = 0.01720209895
    #calculate mean motion
    n = mu/((a)**3)**0.5 #mean motion
    M = M0 + n*(t-t0)
    #calculate E
    Eguess = M
    prevEguess = Eguess
    Eguess = M + e*sin(Eguess)
    while abs(Eguess-prevEguess) > 1.0e-4:
        prevEguess = Eguess
        Eguess = M + e*sin(Eguess)
    E = Eguess
    #calculate the x and y coordinates using physics
    xCoord = a*cos(E) - e*a
    yCoord = a*sin(E)*((1-e**2)**0.5)
    #determine r and nu
    r = (xCoord**2 + yCoord**2)**0.5
    nu = acos(xCoord/r)
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
    #Earth to sun vector from JPL horizons
    Rvector = sunVec
    #Find Earth to asteroid vector
    EtoAst = np.add(rcoords, Rvector)
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
    #calculate RA and dec
    RAsAndDecs = []
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
    RAsAndDecs.append(decDecimal)
    if RADecimal < 0:
        RADecimal = RADecimal + 24
    RAsAndDecs.append(RADecimal*(360/24))
    return RAsAndDecs
    print( "dec: " + str(dec) + " RA: " + str(RA))  
    
def MoG(dataList):
    #list the constants
    #Gaussian constant
    k = 0.0172020985
    #mu
    mu = 1
    #speed of light
    cLight = 173.145
    #obliquity
    obliquity = 23.43676*(pi/180)
    #constant for the starting value of the Newtonian Rhapson method
    constant = 0.85
    #open and add the contents of the file into an array list
    #[date, time, RA, Dec, Sun Vector]
    #[yyyy, mm, dd, hh, mm, ss, hh, mm, ss.ss, dd, mm, ss.s, sunX, sunY, sunZ]
    data = []
    for i in range(len(dataList)):
        subdata = []
        for j in range(len(dataList[1])):
            subdata.append(dataList[i][j])
        data.append(subdata)
    #create a sun vector
    RVec1 = np.array([data[0][12], data[0][13], data[0][14]])
    RVec2 = np.array([data[1][12], data[1][13], data[1][14]])
    RVec3 = np.array([data[2][12], data[2][13], data[2][14]])
    #find the magnitude of the 2nd sun vector
    R2Mag = 0
    for i in range(len(RVec2)):
        R2Mag = R2Mag + RVec2[i]**2
    R2Mag = sqrt(R2Mag)
    #convert the times to JD
    #create a list of just the dates and times
    dateAndTimes = []
    for i in range(len(data)):
        oneDate = []
        for j in range(6):
            oneDate.append(data[i][j])
        dateAndTimes.append(oneDate)
    #convert the times to fractional days
    for i in range(len(dateAndTimes)):
        dateAndTimes[i][3] = dateAndTimes[i][3] + dateAndTimes[i][4]/60.0 + dateAndTimes[i][5]/3600.0
    #convert the dates and times to JD
    JDTimes = []
    for i in range(len(dateAndTimes)):
        J0 = 367*dateAndTimes[i][0] - int((7*(dateAndTimes[i][0]+int((dateAndTimes[i][1]+9)/12)))/4) + int((275*dateAndTimes[i][1])/9) + dateAndTimes[i][2] + 1721013.5
        JD = J0 + dateAndTimes[i][3]/24
        JDTimes.append(JD)
    t1 = JDTimes[0]
    t2 = JDTimes[1]
    t3 = JDTimes[2]
    t10 = t1
    t20 = t2
    t30 = t3
    #calculate the taos
    tao1 = k*(t1-t2)
    tao3 = k*(t3-t2)
    tao = tao3-tao1
    #convert the inputted RAs and Decs into radians
    for i in range(len(data)):
        RADec = data[i][6] + data[i][7]/60.0 + data[i][8]/3600.0
        RARad = RADec * (360/24.0) * (pi/180)
        data[i][6] = RARad
        decDec = data[i][9] + data[i][10]/60.0 + data[i][11]/3600.0
        decRad = decDec * (pi/180.0)
        data[i][9] = decRad
    #calculate the rho hats
    rhoHats = []
    for i in range(len(data)):
        rhoHats.append([[cos(data[i][6])*cos(data[i][9])], [sin(data[i][6])*cos(data[i][9])], [sin(data[i][9])]])
    rhoHat1 = rhoHats[0]
    rhoHat2 = rhoHats[1]
    rhoHat3 = rhoHats[2]
    #convert the rhohats to ecliptic coordinates
    rotation = np.array([[1,0,0],[0,cos(obliquity),sin(obliquity)],[0, -sin(obliquity), cos(obliquity)]])
    rhoHat1V = np.dot(rotation,rhoHat1)
    rhoHat2V = np.dot(rotation,rhoHat2)
    rhoHat3V = np.dot(rotation,rhoHat3)
    rhoHat1 = []
    rhoHat2 = []
    rhoHat3 = []
    for i in range(len(rhoHat1V)):
        rhoHat1.append(float(rhoHat1V[i]))
        rhoHat2.append(float(rhoHat2V[i]))
        rhoHat3.append(float(rhoHat3V[i]))
    rhoHat1 = np.array(rhoHat1)
    rhoHat2 = np.array(rhoHat2)
    rhoHat3 = np.array(rhoHat3)
    #calculate the rhos
    A1 = tao3/(1.0*tao)
    A3 = -tao1/(1.0*tao)
    B1 = (1/6.0)*A1*((tao**2)-(tao3**2))
    B3 = (1/6.0)*A3*(tao**2-tao1**2)
    #calculate the c's later because you need rMag for that
    rhoHat2Cross3 = np.cross(rhoHat2, rhoHat3)
    D0 = np.dot(rhoHat1, rhoHat2Cross3)
    D11 = np.dot(np.cross(RVec1, rhoHat2), rhoHat3)
    D12 = np.dot(np.cross(RVec2, rhoHat2), rhoHat3)
    D13 = np.dot(np.cross(RVec3, rhoHat2), rhoHat3)
    D21 = np.dot(np.cross(rhoHat1, RVec1), rhoHat3)
    D22 = np.dot(np.cross(rhoHat1, RVec2), rhoHat3)
    D23 = np.dot(np.cross(rhoHat1, RVec3), rhoHat3)
    D31 = np.dot(rhoHat1, np.cross(rhoHat2, RVec1))
    D32 = np.dot(rhoHat1, np.cross(rhoHat2, RVec2))
    D33 = np.dot(rhoHat1, np.cross(rhoHat2, RVec3))
    A = (A1*D21-D22+A3*D23)/(-1*D0)
    B = (B1*D21 + B3*D23)/(-D0)
    E = -2*(np.dot(rhoHat2, RVec2))
    F = R2Mag**2
    #solve for rMag
    a = -1*(A**2 + A*E + F)
    b = -mu*(2*A*B + B*E)
    c = -1*(mu**2)*(B**2)
    r2 = poly.polyroots([c, 0, 0, b, 0, 0, a, 0, 1])
    isReal = np.isreal(r2)
    #lists for the r2Vecs and r2Dots
    r2Vecs = []
    r2DotVecs = []
    ##################################################
    for i in range(len(r2)):
        if isReal[i] == True and r2[i]>0:
            r2Mag = float(np.real(r2[i]))
            #solve for the c's
            c1 = A1 + ((mu*B1)/(r2Mag**3))
            c2 = -1
            c3 = A3 + ((mu*B3)/(r2Mag**3))
            #solve for the rho's
            rho1 = (c1*D11 + c2*D12 + c3*D13)/(1.0*c1*D0)
            rho2 = (c1*D21+c2*D22+c3*D23)/(c2*D0)
            rho3 = (c1*D31 + c2*D32 + c3*D33)/(1.0*c3*D0)
            rhoVec2 = rho2*rhoHat2
            r2Vec = rho2*rhoHat2 - RVec2
            r1Vec = rho1*rhoHat1 - RVec1
            r3Vec = rho3*rhoHat3 - RVec3
            r2Vecs.append(r2Vec)
            #calculate the initial r2DotVec
            #calculate the initial f and g values using the first 2 terms of the taylor series
            fInit1 = 1 - ((mu*(tao1**2))/(2*(r2Mag**3)))
            gInit1 = tao1 - ((mu*(tao1**3))/(6*(r2Mag**3)))
            fInit3 = 1 - ((mu*tao3**2)/(2*(r2Mag**3)))
            gInit3 = tao3 - ((mu*(tao3**3))/(6*(r2Mag**3)))
            f1 = fInit1
            f3 = fInit3
            g1 = gInit1
            g3 = gInit3
            d1 = (-f3)/(f1*g3-f3*g1)
            d3 = (f1)/(f1*g3-f3*g1)
            r2dotV = d1*r1Vec + d3*r3Vec
            r2dotPrev = r2dotV
            #find the new c1, c2, c3 values
            c1 = (g3)/(f1*g3-g1*f3)
            c2 = -1
            c3 = (-g1)/(f1*g3-g1*f3)
            #find the new rho's
            rho1Prev = rho1
            rho1 = (c1*D11 + c2*D12 + c3*D13)/(1.0*c1*D0)
            rho2Prev = rho2
            rho2 = (c1*D21+c2*D22+c3*D23)/(c2*D0)
            rho3Prev = rho3
            rho3 = (c1*D31 + c2*D32 + c3*D33)/(1.0*c3*D0)
            #find the new r2Vecs
            r2Vec = rho2*rhoHat2 - RVec2
            r1Vec = rho1*rhoHat1 - RVec1
            r3Vec = rho3*rhoHat3 - RVec3
            r2MagPrev = r2Mag
            r2Mag = np.linalg.norm(r2Vec)
            #correct for light travel time
            t1 = t10 - (rho1)/cLight
            t2 = t20 - rho2/cLight
            t3 = t30 - rho3/cLight
            #calculate the new taos
            tao1 = k*(t1-t2)
            tao3 = k*(t3-t2)
            tao = tao3-tao1
            #list of the possible values
            results = []
            #Taylor Series continues (iterations)
            count = 0
            while abs(rho2-rho2Prev) > 1e-19 and count < 10000:
                u = (mu)/(r2Mag**3)
                z = (np.dot(r2Vec, r2dotV))/(r2Mag**2)
                q = ((np.dot(r2dotV, r2dotV))/(r2Mag**2))-u
                f1 = 1 - ((mu*tao1**2)/(2*r2Mag**3)) + ((mu*(np.dot(r2Vec, r2dotV)))/(2*(r2Mag**5)))*tao1**3  + (1/24.0)*(tao1**4)*(3*u*q - 15*u*(z**2) + (u**2))
                f3 = 1 - ((mu*tao3**2)/(2*r2Mag**3)) + ((mu*(np.dot(r2Vec, r2dotV)))/(2*r2Mag**5))*tao3**3 + (1/24.0)*(tao3**4)*(3*u*q - 15*u*(z**2) + (u**2))
                g1 = tao1 - ((mu*tao1**3)/(6*r2Mag**3)) + (1/24.0)*(tao1**4)*(6*u*z)
                g3 = tao3 - ((mu*tao3**3)/(6*r2Mag**3)) + (1/24.0)*(tao3**4)*(6*u*z)
                d1 = (-f3)/(f1*g3-f3*g1)
                d3 = (f1)/(f1*g3-f3*g1)
                r2dotV = d1*r1Vec + d3*r3Vec
                #find the new c1, c2, c3 values
                c1 = (g3)/(f1*g3-g1*f3)
                c2 = -1
                c3 = (-g1)/(f1*g3-g1*f3)
                #find the new rho's
                rho1Prev = rho1
                rho1 = (c1*D11 + c2*D12 + c3*D13)/(1.0*c1*D0)
                rho2Prev = rho2
                rho2 = (c1*D21+c2*D22+c3*D23)/(c2*D0)
                rho3Prev = rho3
                rho3 = (c1*D31 + c2*D32 + c3*D33)/(1.0*c3*D0)
                #find the new r2Vecs
                r2Vec = rho2*rhoHat2 - RVec2
                r1Vec = rho1*rhoHat1 - RVec1
                r3Vec = rho3*rhoHat3 - RVec3
                r2MagPrev = r2Mag
                r2Mag = np.linalg.norm(r2Vec)
                #correct for light travel time
                t1 = t10 - (rho1)/cLight
                t2 = t20 - rho2/cLight
                t3 = t30 - rho3/cLight
                #calculate the new taos
                tao1 = k*(t1-t2)
                tao3 = k*(t3-t2)
                tao = tao3-tao1
                count = count + 1
            subResults = []
            subResults.append(r2Vec)
            subResults.append(r2dotV)
            subResults.append(t2)
            results.append(subResults)
    return results

#for loop to read the file and split it into lists with three rows each
file = open("XuInput.txt")
data = []
for line in file:
    subData = []
    for word in line.split(" "):
        subData.append(float(word))
    data.append(subData)
data1 = []
data2 = []
data1.append(data[0])
data1.append(data[1])
data1.append(data[3])
data2.append(data[0])
data2.append(data[2])
data2.append(data[3])

#list that will contain all the results
results = []
results.append(MoG(data1))
results.append(MoG(data2))

#get the sun vectors
RV2 = []
RV3 = []
RV2.append(data[1][12])
RV2.append(data[1][13])
RV2.append(data[1][14])
RV3.append(data[2][12])
RV3.append(data[2][13])
RV3.append(data[2][14])
RVs = []
RVs.append(RV2)
RVs.append(RV3)

#list of the elements from the baby OD
elements = []


#list for the fit RAs and Decs
RAsAndDecs = []

#for loop that runs the results through the baby OD to get the orbital elements
for i in range(len(results)):
    r2 = results[i][0][0]
    rdot2 = results[i][0][1]
    t2 = results[i][0][2]
    elements.append(orbElements(r2, rdot2, t2, RVs[i]))
    
#prints the orbital elements in a nice format
print("******************************************************************************")
print("******************************************************************************")
print("For permutations 1,2,4")
print("a: " + str(elements[0][0]))
print("e: " + str(elements[0][1]))
print("I: " + str(elements[0][2]))
print("O: " + str(elements[0][3]))
print("w: " + str(elements[0][4]))
print("M: " + str(elements[0][7]))
print("******************************************************************************")
print("******************************************************************************")
print("For permutations 1,3,4")
print("a: " + str(elements[1][0]))
print("e: " + str(elements[1][1]))
print("I: " + str(elements[1][2]))
print("O: " + str(elements[1][3]))
print("w: " + str(elements[1][4]))
print("M: " + str(elements[1][7]))
print("******************************************************************************")
print("******************************************************************************")

#Create lists with each set of orbital elements and the average orbital elements of the permutations
elements1 = elements[0]
elements2 = elements[1]
elements = np.add(elements1, elements2)
elements = elements/2.0
#print the averages of the permutations
print("The average of the orbital elements")
print("a: " + str(elements[0]))
print("e: " + str(elements[1]))
print("I: " + str(elements[2]))
print("O: " + str(elements[3]))
print("w: " + str(elements[4]))
print("M: " + str(elements[7]))
print("******************************************************************************")
print("******************************************************************************")

#perform differential correction on the obs RA and dec
def diffCorr(data, elements, r2Vec, r2DotVec):
    #calculate the average of the orbital elements
    aMean = elements[0]
    eMean = elements[1]
    IMean = elements[2]*(pi/180)
    OMean = elements[3]*(pi/180)
    wMean = elements[4]*(pi/180)
    M2Mean = elements[5]*(pi/180)

    #create the fitData list needed for differential correction
    fitData = []
    t2 = 2458660.739861111
    t1 = 2458655.7809375
    t3 = 2458675.7224189816
    t4 = 2458679.758599537
    times = np.array([t1, t2, t3, t4])
    sunVecs = np.array([[[data[0][12]], [data[0][13]], [data[0][14]]], [[data[1][12]], [data[1][13]], [data[1][14]]], [[data[2][12]], [data[2][13]], [data[2][14]]], [[data[3][12]], [data[3][13]], [data[3][14]]]])
    fitData.append(ephemeris(eMean, aMean, IMean, OMean, wMean, M2Mean, t1, t2, sunVecs[0]))
    fitData.append(ephemeris(eMean, aMean, IMean, OMean, wMean, M2Mean, t2, t2, sunVecs[1]))
    fitData.append(ephemeris(eMean, aMean, IMean, OMean, wMean, M2Mean, t3, t2, sunVecs[2]))
    fitData.append(ephemeris(eMean, aMean, IMean, OMean, wMean, M2Mean, t4, t2, sunVecs[3]))

    #calculate the decimal values for the RA and decs in data
    for i in range(len(data)):
        RADec = data[i][6] + data[i][7]/60.0 + data[i][8]/3600.0
        RARad = RADec * (360/24.0) * (pi/180)
        data[i][6] = RARad
        decDec = data[i][9] + data[i][10]/60.0 + data[i][11]/3600.0
        decRad = decDec * (pi/180.0)
        data[i][9] = decRad

    #calculate the delta RAs and decs
    deltaRAs = np.array([data[0][6] - fitData[0][1], data[1][6] - fitData[1][1], data[2][6] - fitData[2][1], data[3][6] - fitData[3][1]])
    deltaDecs = np.array([data[0][9] - fitData[0][0], data[1][9] - fitData[1][0], data[2][9] - fitData[2][0], data[3][9] - fitData[3][0]])

    
    #for loop that creates a list with all the individual partial derivatives for RA and dec that changes the xyz values
    RApartials = []
    decpartials = []
    delta = 1.0e-3

    #so that referencing orbElements inside a double for loop won't throw an error
    def orbElements(r2EclV, r2dotEclV, t2, sunVec):
        #get a
        r2SquaredSum = 0
        for i in range(len(r2EclV)):
            r2SquaredSum = r2SquaredSum + r2EclV[i]**2
        rMag = sqrt(r2SquaredSum)
        a = ((2.0/rMag)-np.dot(r2dotEclV, r2dotEclV))**-1
        #get e
        r2rdot2Crossed = np.cross(r2EclV, r2dotEclV)
        r2rdot2CrossedMag = 0
        for i in range(len(r2dotEclV)):
            r2rdot2CrossedMag = r2rdot2CrossedMag + r2rdot2Crossed[i]**2
        r2rdot2CrossedMag = sqrt(r2rdot2CrossedMag)
        e = (1-r2rdot2CrossedMag**2/a)**0.5
        #convert coordinates to ecliptic (not necessary because the coordinates from the Method of Gauss are alreadyin ecliptic
        r2Ecl = 0
        for i in range(len(r2EclV)):
            r2Ecl = r2Ecl + r2EclV[i]**2
        r2Ecl = sqrt(r2Ecl)
        #get I in degrees
        hVector = np.cross(r2EclV,r2dotEclV)
        h = 0
        for i in range(len(hVector)):
            h = h + hVector[i]**2
        h = sqrt(h)
        hz = r2EclV[0]*r2dotEclV[1]-r2EclV[1]*r2dotEclV[0]
        I = acos(1.0*hz/h)
        IDeg = positiveAng(I*(180/pi))
        #get Oprime in degrees
        hx = r2EclV[1]*r2dotEclV[2]-r2EclV[2]*r2dotEclV[1]
        hy = r2EclV[2]*r2dotEclV[0]-r2EclV[0]*r2dotEclV[2]
        sinOprime = hx/(h*sin(I))
        cosOprime = hy/(-h*sin(I))
        Oprime = atan2(sinOprime, cosOprime)
        OprimeDeg = positiveAng(Oprime*(180/pi))
        #get w in degrees (w = U-v)
        #get U
        cosU = (r2EclV[0]*cos(Oprime) + r2EclV[1]*sin(Oprime))/r2Ecl
        sinU = r2EclV[2]/(r2Ecl*sin(I))
        U = atan2(sinU, cosU)
        #get nu
        cosNu = ((a*(1-e**2))/r2Ecl - 1)/e
        rdprdot = np.dot(r2EclV, r2dotEclV)
        sinNu = (((a*(1-e**2))/h)*(rdprdot/r2Ecl))/e
        nu = atan2(sinNu, cosNu)
        #get w
        w = U-nu
        wDeg = positiveAng(w*(180/pi))
        #get M in degrees
        E = acos((1/e)*(1-(r2Ecl/a)))
        if nu < 0:
            nu = nu + 2*pi
        if nu > pi and nu < 2*pi:
            E = 2*pi-E
        M2 = E-e*sin(E)
        M2Deg = M2*(180/pi)
        #get n
        mu = 0.01720209895
        n = mu/((a)**3)**0.5
        #precess the true anomaly
        tOriginal = [2019, 7, 21, 6, 00, 00]
        J0 = 367*tOriginal[0] - int((7*(tOriginal[0]+int((tOriginal[1]+9)/12)))/4) + int((275*tOriginal[1])/9) + tOriginal[2] + 1721013.5
        TDueJD = J0 + tOriginal[3]/24
        M2Precessed = M2 + n*(TDueJD-t2)
        M2PrecessedDeg = positiveAng(M2Precessed*(180/pi))
        #get E
        EDeg = positiveAng(E*(180/pi))
        #get T
        T = t2 - (M2/n)
        #get P
        P = ((2*pi)/n)/365.25
        RV = sunVec
        orbElements = np.array([a, e, IDeg, OprimeDeg, wDeg, M2Deg, t2])
        return orbElements
    
    for i in range(len(deltaRAs)):
        for j in range(3):
            #copies of the r2Vec to change the xyz values
            r2VecTemp1 = r2Vec.copy()
            r2VecTemp2 = r2Vec.copy()
            r2VecTemp1[j] = r2VecTemp1[j] + delta
            r2VecTemp2[j] = r2VecTemp2[j] - delta
            #lists that will contain the orbital elements form the baby OD
            orbElements1 = []
            orbElements2 = []
            #for the first ra and dec
            orbElements1 = orbElements(r2VecTemp1, r2DotVec, t2, sunVecs[i])
            a1 = orbElements1[0]
            e1 = orbElements1[1]
            I1 = orbElements1[2]
            O1 = orbElements1[3]
            w1 = orbElements1[4]
            M1 = orbElements1[5]
            RAandDec1 = ephemeris(e1, a1, I1, O1, w1, M1, times[i], t2, sunVecs[i])
            #for the second ra and dec
            orbElements2 = orbElements(r2VecTemp2, r2DotVec, t2, sunVecs[i])
            a2 = orbElements2[0]
            e2 = orbElements2[1]
            I2 = orbElements2[2]
            O2 = orbElements2[3]
            w2 = orbElements2[4]
            M2 = orbElements2[5]
            RAandDec2 = ephemeris(e2, a2, I2, O2, w2, M2, times[i], t2, sunVecs[i])
            #find the partial derivatives for RA and dec
            RApartials.append((RAandDec1[1]-RAandDec2[1])/(2.0*delta))
            decpartials.append((RAandDec1[0]-RAandDec2[0])/(2.0*delta))
                
    #for loop that creates a list with all the individual partial derivatives for RA and dec that changes the xyz dot values
    RApartialsdots = []
    decpartialsdots = []
    for i in range(len(deltaRAs)):
        for j in range(3):
            #copies of the r2dotVec to change the xyz dot values
            r2dotVecTemp1 = r2DotVec.copy()
            r2dotVecTemp2 = r2DotVec.copy()
            r2dotVecTemp1[j] = r2dotVecTemp1[j] + delta
            r2dotVecTemp2[j] = r2dotVecTemp2[j] - delta
            #lists that will contain the orbital elements form the baby OD
            orbElements1 = []
            orbElements2 = []
            #for the first ra and dec
            orbElements1 = orbElements(r2Vec, r2dotVecTemp1, t2, sunVecs[i])
            a1 = orbElements1[0]
            e1 = orbElements1[1]
            I1 = orbElements1[2]
            O1 = orbElements1[3]
            w1 = orbElements1[4]
            M1 = orbElements1[5]
            RAandDec1 = ephemeris(e1, a1, I1, O1, w1, M1, times[i], t2, sunVecs[i])
            #for the second ra and dec
            orbElements2 = orbElements(r2Vec, r2dotVecTemp2, t2, sunVecs[i])
            a2 = orbElements2[0]
            e2 = orbElements2[1]
            I2 = orbElements2[2]
            O2 = orbElements2[3]
            w2 = orbElements2[4]
            M2 = orbElements2[5]
            RAandDec2 = ephemeris(e2, a2, I2, O2, w2, M2, times[i], t2, sunVecs[i])
            #find the partial derivaties for RA and dec
            RApartialsdots.append((RAandDec1[1]-RAandDec2[1])/(2.0*delta))
            decpartialsdots.append((RAandDec1[0]-RAandDec2[0])/(2.0*delta))
    #Making the Jacobian Matrix
    #Get lists of the individual components
    daxList = []
    dayList = []
    dazList = []
    daxdotList = []
    daydotList = []
    dazdotList = []
    ddxList = []
    ddyList = []
    ddzList = []
    ddxdotList = []
    ddydotList = []
    ddzdotList = []
    for i in range(len(RApartials)):
        if i % 3 == 0:
            daxList.append(RApartials[i])
            ddxList.append(decpartials[i])
            daxdotList.append(RApartialsdots[i])
            ddxdotList.append(decpartialsdots[i])
        if i % 3 == 1:
            dayList.append(RApartials[i])
            ddyList.append(decpartials[i])
            daydotList.append(RApartialsdots[i])
            ddydotList.append(decpartialsdots[i])
        if i % 3 == 2:
            dazList.append(RApartials[i])
            ddzList.append(decpartials[i])
            dazdotList.append(RApartialsdots[i])
            ddzdotList.append(decpartialsdots[i])
    
    #trying it another way
    #get the sums of the delta RAs and decs
    #Even though they variables are labeled a, they include all positional data
    ddddx = sum(np.multiply(deltaDecs, ddxList))
    ddddy = sum(np.multiply(deltaDecs, ddyList))
    ddddz = sum(np.multiply(deltaDecs, ddzList))
    ddddxdot = sum(np.multiply(deltaDecs, ddxdotList))
    ddddydot = sum(np.multiply(deltaDecs, ddydotList))
    ddddzdot = sum(np.multiply(deltaDecs, ddzdotList))
    dadax = sum(np.multiply(deltaRAs, daxList)) + ddddx
    daday = sum(np.multiply(deltaRAs, dayList)) + ddddy
    dadaz = sum(np.multiply(deltaRAs, dazList)) + ddddz
    dadaxdot = sum(np.multiply(deltaRAs, daxdotList)) + ddddxdot
    dadaydot = sum(np.multiply(deltaRAs, daydotList)) + ddddydot
    dadazdot = sum(np.multiply(deltaRAs, dazdotList)) + ddddzdot
    

    #the jacobian matrix 
    j2 = np.array([[sum(np.multiply(daxList, daxList)), sum(np.multiply(daxList, dayList)), sum(np.multiply(daxList, dazList)), sum(np.multiply(daxList, daxdotList)), sum(np.multiply(daxList, daydotList)), sum(np.multiply(daxList, dazdotList))],
                   [sum(np.multiply(dayList, daxList)), sum(np.multiply(dayList, dayList)), sum(np.multiply(dayList, dazList)), sum(np.multiply(dayList, daxdotList)), sum(np.multiply(dayList, daydotList)), sum(np.multiply(dayList, dazdotList))],
                   [sum(np.multiply(dazList, daxList)), sum(np.multiply(dazList, dayList)), sum(np.multiply(dazList, dazList)), sum(np.multiply(dazList, daxdotList)), sum(np.multiply(dazList, daydotList)), sum(np.multiply(dazList, dazdotList))],
                   [sum(np.multiply(daxdotList, daxList)), sum(np.multiply(daxdotList, dayList)), sum(np.multiply(daxdotList, dazList)), sum(np.multiply(daxdotList, daxdotList)), sum(np.multiply(daxdotList, daydotList)), sum(np.multiply(daxdotList, dazdotList))],
                   [sum(np.multiply(daydotList, daxList)), sum(np.multiply(daydotList, dayList)), sum(np.multiply(daydotList, dazList)), sum(np.multiply(daydotList, daxdotList)), sum(np.multiply(daydotList, daydotList)), sum(np.multiply(daydotList, dazdotList))],
                   [sum(np.multiply(dazdotList, daxList)), sum(np.multiply(dazdotList, dayList)), sum(np.multiply(dazdotList, dazList)), sum(np.multiply(dazdotList, daxdotList)), sum(np.multiply(dazdotList, daydotList)), sum(np.multiply(dazdotList, dazdotList))]])
    deltaMat = np.array([[dadax], [daday], [dadaz], [dadaxdot], [dadaydot], [dadazdot]])


    #calculate the RMS
    #separate the original fit data in lists of just the fit RAs and dec and then calculate the original RMS data
    summation1 = sum(np.multiply(deltaRAs, deltaRAs)) + sum(np.multiply(deltaDecs, deltaDecs))
    RMSOriginal = sqrt(summation1/2.0)

    #find the new fit's RMS
    newFitData = []
    correctedVs = np.dot(np.linalg.pinv(j2), deltaMat)
    r2VCorrected = []
    r2dotVCorrected = []
    for i in range(3):
        r2VCorrected.append(float(correctedVs[i] + results[0][0][0][i]))
    for i in range(3, 6):
        r2dotVCorrected.append(float(correctedVs[i] + results[0][0][1][i-3]))
    t2 = 2458660.739861111
    finalElements = orbElements(r2VCorrected, r2dotVCorrected, t2, RV2)
    
    #recalculate the RAs and decs using the final elements garnered from the differential correction
    a = finalElements[0]
    e = finalElements[1]
    I = finalElements[2]*(pi/180)
    O = finalElements[3]*(pi/180)
    w = finalElements[4]*(pi/180)
    M = finalElements[5]*(pi/180)
    newFitData.append(ephemeris(e, a, I, O, w, M, t1, t2, sunVecs[0]))
    newFitData.append(ephemeris(e, a, I, O, w, M, t2, t2, sunVecs[1]))
    newFitData.append(ephemeris(e, a, I, O, w, M, t3, t2, sunVecs[2]))
    newFitData.append(ephemeris(e, a, I, O, w, M, t4, t2, sunVecs[3]))
    #get the new fit RAs and Decs and calculate the RMS
    
    deltaNewRAs = np.array([data[0][6] - newFitData[0][1], data[1][6] - newFitData[1][1], data[2][6] - newFitData[2][1], data[3][6] - newFitData[3][1]])
    deltaNewDecs = np.array([data[0][9] - newFitData[0][0], data[1][9] - newFitData[1][0], data[2][9] - newFitData[2][0], data[3][9] - newFitData[3][0]])
    summation2 = sum(np.multiply(deltaNewRAs, deltaNewRAs)) + sum(np.multiply(deltaNewDecs, deltaNewDecs))
    RMSNew = sqrt(summation2/2.0)
    return RMSNew, RMSOriginal
  
#print the RMS difference
r2Vec = results[0][0][0]
r2dotVec = results[0][0][1]

#The new and original RMS values
print("The new and original RMS values: ")
print(diffCorr(data, elements1, r2Vec, r2dotVec))
print("Conclusion, the difference between the RMS values is not significant enough to result in better orbital element values.")



