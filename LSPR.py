import numpy as np
from numpy.linalg import inv
from math import *


#function to calculate the 6 plate constants, ra/dec uncertainty. and ra/dec coordinates
#redo in another doc with substring to convert to decimal and to get rid of the space
def lspr(fileName, positionX, positionY):
    file = open(fileName)
    listStars = []
    lineNumber = 0
    for line in file:
        for word in line.split(" "):
            listStars.append(word)
        lineNumber = lineNumber + 1
    starCoords = np.array(listStars)
    #converts ras to decimals and finds N
    raCount = 0
    N=0
    for i in range(len(starCoords)):
        if raCount % 4 == 2:
##            raString = starCoords[i]
##            raDecimal = int(raString[0:2]) + float(raString[3:5])/60. + float(raString[6:11])/3600.
            starCoords[i] = starCoords[i]*(pi/180)
            N = N+1
        raCount = raCount+1
##    #converts decs to decimals
    decCount = 0
    for i in range(len(starCoords)):
        if decCount % 4 == 3:
##            decString = starCoords[i]
##            decDecimal = int(decString[1:3]) + float(decString[4:6])/60. + float(decString[7:11])/3600.
            starCoords[i] = starCoords[i]*(pi/180)
        decCount = decCount + 1
    #create lists for each star component
    xCoord = []
    yCoord = []
    starRA = []
    starDec = []
    count = 0
    for i in range(len(starCoords)):
        if count % 4 == 0:
            xCoord.append(float(starCoords[i]))
        if count % 4 == 1:
            yCoord.append(float(starCoords[i]))
        if count % 4 == 2:
            starRA.append(float(starCoords[i]))
        if count % 4 == 3:
            starDec.append(float(starCoords[i]))
        count = count + 1
    xCoordArr = np.array(xCoord)
    yCoordArr = np.array(yCoord)
    starRAArr = np.array(starRA)
    starDecArr = np.array(starDec)
    
    #determining the plate constants
    raSum = sum(starRA)
    raAndXSum = sum(starRAArr*xCoord)
    raAndYSum = sum(starRAArr*yCoord)
    xSum = sum(xCoord)
    ySum = sum(yCoord)
    xSquaredSum = sum(xCoordArr*xCoordArr)
    xAndYSum = sum(xCoordArr*yCoordArr)
    ySquaredSum = sum(yCoordArr*yCoordArr)
    matrix = np.array([[N, xSum, ySum],[xSum, xSquaredSum, xAndYSum],[ySum, xAndYSum, ySquaredSum]])
    matrixinv = np.array(inv(matrix))
    matrix1 = np.array([raSum, raAndXSum, raAndYSum])
    plates1 = np.array(np.dot(matrixinv, matrix1))
    b1Str = plates1[0]*(360/24)
    a11Str = plates1[1]*(360/24)
    a12Str = plates1[2]*(360/24)
    b1 = plates1[0]
    a11 = plates1[1]
    a12 = plates1[2]
    decSum = sum(starDec)
    decAndXSum = sum(starDecArr*xCoord)
    decAndYSum = sum(starDecArr*yCoord)
    matrix2 = np.array([decSum, decAndXSum, decAndYSum])
    plates2 = np.array(np.dot(matrixinv, matrix2))
    b2 = plates2[0]
    a21 = plates2[1]
    a22 = plates2[2]

    #determining the best fits RA and dec
    RADecimal = b1 + a11*positionX + a12*positionY
    RAmin = (RADecimal - trunc(RADecimal))*60
    RAminStr = str(trunc(RAmin))
    if RAmin< 10.0:
        RAminStr = "0" + RAminStr
    RAsec = (RAmin - trunc(RAmin))*60
    RAsecStr = str(RAsec)
    if RAsec < 10.0:
        RAsecStr = "0"+RAsecStr
    RA = str(trunc(RADecimal)) + ":" + RAminStr + ":" + RAsecStr
    decDecimal = b2 + a21*positionX + a22*positionY
    decArcMin = (decDecimal - trunc(decDecimal))*60
    decArcMinStr = str(trunc(decArcMin))
    if decArcMin<10.0:
        decArcMinStr = "0" + decArcMinStr
    decArcSec = (decArcMin - trunc(decArcMin))*60
    decArcSecStr = str(decArcSec)
    if decArcSec < 10.0:
        decArcSecStr = "0" + decArcSecStr
    dec = str(trunc(decDecimal)) + ":" + decArcMinStr + ":" + decArcSecStr

    #calculate the RA uncertainty
    #starRAArr has the actual RAs
    starRAArr = (360/24)*starRAArr
    fitRAs = []
    for i in range(len(starRAArr)):
        fitRAs.append(b1Str + a11Str*xCoordArr[i] + a12Str*yCoordArr[i])
    RASquaredDiff = []
    for i in range(len(starRAArr)):
        RASquaredDiff.append((starRAArr[i] - fitRAs[i])**2)
    RASquaredDiffArr = np.array(RASquaredDiff)
    RAUncertaintyDeg = ((1/(N-3))*(sum(RASquaredDiffArr)))**0.5
    RAUncert = RAUncertaintyDeg*3600.
    
    #calculate the dec uncertainty
    fitDecs = []
    for i in range(len(starDecArr)):
        fitDecs.append(b2 + a21*xCoordArr[i] + a22*yCoordArr[i])
    DecSquaredDiff = []
    for i in range(len(starDecArr)):
        DecSquaredDiff.append((starDecArr[i] - fitDecs[i])**2)
    DecSquaredDiffArr = np.array(DecSquaredDiff)
    SumDecSquaredDiff = sum(DecSquaredDiffArr)
    DecUncertainty = (((1/(N-3))*(SumDecSquaredDiff))**0.5)*3600

    return "Unflattened: \nb1: " + str(b1Str) + " b2: " + str(b2) + "\na11: " + str(a11Str) + " a12: " + str(a12Str) + "\na21: " + str(a21) + " a22: " + str(a22) \
+ "\nRA: " + RA + " dec: " + dec + "\nRA Uncertainty: " + str(RAUncert) + " Dec Uncertainty: " + str(DecUncertainty)

def lsprFlat(fileName, positionX, positionY):
    file = open(fileName)
    listStars = []
    lineNumber = 0
    for line in file:
        for word in line.split(" "):
            listStars.append(word)
        lineNumber = lineNumber + 1
    starCoords = np.array(listStars)
    print(starCoords)
    #converts ras to decimals and finds N
    raCount = 0
    N=0
    for i in range(len(starCoords)):
        if raCount % 4 == 2:
            raString = starCoords[i]
            raDecimal = int(raString[0:2]) + float(raString[3:5])/60. + float(raString[6:11])/3600.
            starCoords[i] = raDecimal
            N = N+1
        raCount = raCount+1
    #converts decs to decimals
    decCount = 0
    for i in range(len(starCoords)):
        if decCount % 4 == 3:
            decString = starCoords[i]
            decDecimal = int(decString[0:2]) + float(decString[3:5])/60. + float(decString[6:11])/3600.
            starCoords[i] = decDecimal
        decCount = decCount + 1
    #create lists for each star component
    xCoord = []
    yCoord = []
    starRA = []
    starDec = []
    count = 0
    for i in range(len(starCoords)):
        if count % 4 == 0:
            xCoord.append(float(starCoords[i]))
        if count % 4 == 1:
            yCoord.append(float(starCoords[i]))
        if count % 4 == 2:
            starRA.append(float(starCoords[i]))
        if count % 4 == 3:
            starDec.append(float(starCoords[i]))
        count = count + 1
    xCoordArr = np.array(xCoord)
    yCoordArr = np.array(yCoord)
    starRAArr = np.array(starRA)
    starDecArr = np.array(starDec)
    starRARad = []
    starDecRad = []
    for i in range(len(starDecArr)):
        starRARad.append(starRAArr[i]*(360.0/24)*(pi/180))
        starDecRad.append(starDecArr[i]*(pi/180))
    starDecArrRad = np.array(starDecRad)
    starRAArrRad = np.array(starRARad)
    
    #calculate the flat RAs and flat decs
    Dsum = 0
    for i in range(len(starDecArrRad)):
        Dsum = Dsum + starDecArrRad[i]
    D = Dsum/N*1.0
    RAsum = 0
    for i in range(len(starRAArrRad)):
        RAsum = RAsum + starRAArrRad[i]
    A = RAsum/N*1.0
    L = 3911/(9e-3)

    flatRA = []
    flatDec = []
    for i in range(len(starRAArrRad)):
        H = sin(starDecArrRad[i])*sin(D)+cos(starDecArrRad[i])*cos(D)*cos(starRAArrRad[i] - A)
        flatRA.append(((cos(starDecArrRad[i])*sin(starRAArrRad[i]-A))/H)-(xCoordArr[i]/L))
        flatDec.append(((sin(starDecArrRad[i])*cos(D)-cos(starDecArrRad[i])*sin(D)*cos(starRAArrRad[i]-A))/H)-(yCoordArr[i]/L))
    flatRAArr = np.array(flatRA)
    flatDecArr = np.array(flatDec)
    

    #determining the plate constants
    raSum = sum(flatRAArr)
    raAndXSum = sum(flatRAArr*xCoord)
    raAndYSum = sum(flatRAArr*yCoord)
    xSum = sum(xCoord)
    ySum = sum(yCoord)
    xSquaredSum = sum(xCoordArr*xCoordArr)
    xAndYSum = sum(xCoordArr*yCoordArr)
    ySquaredSum = sum(yCoordArr*yCoordArr)
    matrix = np.array([[N, xSum, ySum],[xSum, xSquaredSum, xAndYSum],[ySum, xAndYSum, ySquaredSum]])
    matrixinv = np.array(inv(matrix))
    matrix1 = np.array([raSum, raAndXSum, raAndYSum])
    plates1 = np.array(np.dot(matrixinv, matrix1))
    b1Str = plates1[0]* (360/24)
    a11Str = plates1[1] *(360/24)
    a12Str = plates1[2] * (360/24)
    b1 = plates1[0]
    a11 = plates1[1]
    a12 = plates1[2]
    decSum = sum(flatDec)
    decAndXSum = sum(flatDecArr*xCoord)
    decAndYSum = sum(flatDecArr*yCoord)
    matrix2 = np.array([decSum, decAndXSum, decAndYSum])
    plates2 = np.array(np.dot(matrixinv, matrix2))
    b2 = plates2[0]
    a21 = plates2[1]
    a22 = plates2[2]
    print(matrixinv, matrix1, matrix2)

    #determining the uncertainty of the stars
    flatFitDec = []
    flatFitRA = []
    for i in range(len(flatRAArr)):
        flatFitRA.append(b1 + a11*xCoord[i] + a12*yCoord[i] + xCoord[i]/L)
        flatFitDec.append(b2 + a21*xCoord[i] + a22*yCoord[i] + yCoord[i]/L)
    #determines the fit RAs and Decs for the stars
    fitDecs = []
    fitRAs = []
    for i in range(len(flatRAArr)):
        deltaStar = cos(D)-flatFitDec[i]*sin(D)
        rStar = ((flatFitRA[i])**2 + deltaStar**2)**0.5
        fitRAs.append((A + atan(flatFitRA[i]/deltaStar))*(180/pi))
        fitDecs.append((atan((sin(D) + flatFitDec[i]*cos(D))/rStar))*(180/pi))
    
    #calculate the RA uncertainty
    #starRAArr has the actual RAs
    starRAArr = (360/24)*starRAArr
    RASquaredDiff = []
    for i in range(len(starRAArr)):
        RASquaredDiff.append((starRAArr[i] - fitRAs[i])**2)
    RASquaredDiffArr = np.array(RASquaredDiff)
    RAUncertaintyDeg = ((1/(N-3))*(sum(RASquaredDiffArr)))**0.5
    RAUncert = RAUncertaintyDeg*3600.
    
    #calculate the dec uncertainty
    DecSquaredDiff = []
    for i in range(len(starDecArr)):
        DecSquaredDiff.append((starDecArr[i] - fitDecs[i])**2)
    DecSquaredDiffArr = np.array(DecSquaredDiff)
    SumDecSquaredDiff = sum(DecSquaredDiffArr)
    DecUncertainty = (((1/(N-3))*(SumDecSquaredDiff))**0.5)*3600
    
    #determining the best fits RA and dec
    RADecimalFlatAst = b1 + a11*positionX + a12*positionY + positionX/L
    decDecimalFlatAst = b2 + a21*positionX + a22*positionY + positionY/L
    delta = cos(D) - decDecimalFlatAst*sin(D)
    r = (RADecimalFlatAst**2 + delta**2)**0.5
    RADecimalAst = A + atan(RADecimalFlatAst/delta)
    DecDecimalAst = atan((sin(D) + decDecimalFlatAst*cos(D))/r)
    RADecimal = RADecimalAst*(180/pi)*(24/360)
    decDecimal = DecDecimalAst*(180/pi)
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
    
    return "Flattened & Unflattened: \nb1: " + str(b1Str) + " b2: " + str(b2) + "\na11: " + str(a11Str) + " a12: " + str(a12Str) + "\na21: " + str(a21) + " a22: " + str(a22) \
+ "\nRA: " + RA + " dec: " + dec + "\nRA Uncertainty: " + str(RAUncert) + " Dec Uncertainty: " + str(DecUncertainty)

    
    
#print("input 2: \n" + lspr("LSPRtestinput2.txt", 1403.6, 1585.9))
#print("input 2: \n" + lsprFlat("LSPRtestinput2.txt", 1403.6, 1585.9))

#print(lsprFlat("Obs2Seq1Im5LSPR.txt", 733.9526315982242, 463.95642280204237))
print(lsprFlat("lsprASTRO1.txt", 603, 423))
