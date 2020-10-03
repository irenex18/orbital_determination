#code to calculate the signal, SNR and instrumental magnitude with uncertainty
from math import*
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

def fractional(radius):
    #find the fraction values for the top half of the circle (minus the middle circle)
    pixelValues = []
    for i in range(1, radius+1):
        subValues = []
        for j in range(radius*2+1):
            result = 0
            for N in range(10000):
                if radius**2 - ((N/10000.0)+j-radius-0.5)**2 > 0 and (1/10000.0)*((-1*radius+i-0.5)+sqrt(radius**2 - ((N/10000.0)+j-radius-0.5)**2)) > 0:
                    result = result + (1/10000.0)*((((-1*radius+i-0.5)+sqrt(radius**2 - ((N/10000.0)+j-radius-0.5)**2))+((-1*radius+i-0.5)+sqrt(radius**2 - ((N/10000.0)+j-radius-0.5)**2)))/2)
            subValues.append(result)
        pixelValues.append(subValues)
    for i in range(1,radius):
        for j in range(2*radius):
            pixelValues[radius-i][j] = pixelValues[radius-i][j] - pixelValues[radius-1-i][j]
    for i in range(radius):
        for j in range(radius):
            pixelValues[i][2*radius-j] = pixelValues[i][j]
    pixelArrHalf = np.array(pixelValues)
    #find the fractional values for the middle part of the circle and add it to the list
    middleL = (pixelArrHalf.T)[radius]
    lrFlippedTArr = np.fliplr(pixelArrHalf.T)
    middleR = lrFlippedTArr[radius]
    midSubValues = []
    for i in middleL:
        midSubValues.append(i)
    midSubValues.append(1)
    for i in middleR:
        midSubValues.append(i)
    pixelValues.append(midSubValues)
    #find the fractional values for the bottom part of the circle, add it to the list, and create the final array
    udFlippedArr = np.flipud(pixelArrHalf)
    for i in range(radius):
        subValues = []
        for j in range(2*radius + 1):
            subValues.append(udFlippedArr[i][j])
        pixelValues.append(subValues)
    fractionalArr = np.array(pixelValues)
    return fractionalArr

def De(filename):
    fileContents = fits.getdata(filename)
    slicedData = np.array(fileContents[20:len(fileContents)-20, 20:len(fileContents)-20])
    return np.mean(slicedData)*0.8
    
######################################################################################################################################################################
##the coords still stay the same!! don't go to 0 - need to account for xcoord and ycoord
def apPhotometryRejected(filename, xCoord, yCoord, apR, anInnerR, anOuterR, darkFileName):
    fileContents = fits.getdata(filename)
    apertureSq = fileContents[yCoord-apR-1:yCoord+apR, xCoord-apR-1:xCoord+apR]
    #calculate the signal, SNR, and magnitude if border pixels accepted
    proportionsAp = fractional(apR)
    for i in range(2*apR+1):
        for j in range(2*apR+1):
            proportionsAp[i][j] = round(proportionsAp[i][j],2)
    #get a list of aperture values
    apValues = []
    nap = 0
    for i in range(apR*2+1):
        for j in range(apR*2+1):
            if proportionsAp[i][j] == 1:
                apValues.append(apertureSq[i][j])
                nap = nap + 1
    #create an annulus for the aperture
    annulusSq = fileContents[yCoord-anOuterR-1:yCoord+anOuterR, xCoord-anOuterR-1:xCoord+anOuterR]
    proportionsAnO = fractional(anOuterR)
    annulusInnerSq = fileContents[yCoord-anInnerR-1:yCoord+anInnerR, xCoord-anInnerR-1:xCoord+anInnerR]
    proportionsAnI = fractional(anInnerR)
    annulusInner = annulusInnerSq*proportionsAnI
    for i in range(2*anOuterR+1):
        for j in range(2*anOuterR+1):
            proportionsAnO[i][j]=round(proportionsAnO[i][j],2)
    for i in range(2*anInnerR+1):
        for j in range(2*anInnerR+1):
            proportionsAnI[i][j] = round(proportionsAnI[i][j],2)
    for i in range(2*anOuterR+1):
        for j in range(2*anOuterR+1):
            if proportionsAnO[i][j] < 1:
                annulusSq[i][j] = 0
    for y in range(anOuterR-anInnerR,(anOuterR*2+1)-anInnerR):
        for x in range(anOuterR-anInnerR, (anOuterR*2+1)-anInnerR):
            if proportionsAnI[y-anInnerR][x-anInnerR] > 0:
                annulusSq[y][x] = 0
    #calculate the average Sky ADU Count
    annulusValues = []
    skyCount = 0
    for i in range(2*anOuterR+1):
        for j in range(2*anOuterR+1):
            if abs(annulusSq[i][j] - 0)>1e-4:
                annulusValues.append(annulusSq[i][j])
                skyCount = skyCount + 1
    avSky = sum(annulusValues)/skyCount*1.0
    #calculate the signal
    signal = sum(apValues) - nap*avSky
    #calcualte SNR
    De = 10
    nan = skyCount
    skyE = avSky*.8
    Se = signal*.8
    p2 = 11**2+(.8**2)/12
    SNR = sqrt(Se)/sqrt(1+nap*(1+ nap*1.0/nan)*((skyE+De+p2)/Se))
    #calculate instant magnitude
    instMag = -2.5*log(signal,10)
    uncertainty = 1.087/SNR
    print(nap,nan)
    return "Border Pixels Rejecteded: " + "\nsignal: " + str(signal) + " SNR: " + str(SNR) + "\ninst mag: " + str(instMag) + " +/- " + str(uncertainty)
    

     
def apPhotometryAccepted(filename, xCoord, yCoord, apR, anInnerR, anOuterR):
    fileContents = fits.getdata(filename)
    apertureSq = fileContents[yCoord-apR-1:yCoord+apR, xCoord-apR-1:xCoord+apR]
    #calculate the signal, SNR, and magnitude if border pixels accepted
    proportionsAp = fractional(apR)
    for i in range(2*apR+1):
        for j in range(2*apR+1):
            proportionsAp[i][j] = round(proportionsAp[i][j],2)
    #get a list of aperture values
    apValues = []
    nap = 0
    for i in range(apR*2+1):
        for j in range(apR*2+1):
            if proportionsAp[i][j] > 0:
                apValues.append(apertureSq[i][j])
                nap = nap + 1
    #create an annulus for the aperture
    annulusSq = fileContents[yCoord-anOuterR-1:yCoord+anOuterR, xCoord-anOuterR-1:xCoord+anOuterR]
    proportionsAnO = fractional(anOuterR)
    annulusInnerSq = fileContents[yCoord-anInnerR-1:yCoord+anInnerR, xCoord-anInnerR-1:xCoord+anInnerR]
    proportionsAnI = fractional(anInnerR)
    annulusInner = annulusInnerSq*proportionsAnI
    for i in range(2*anOuterR+1):
        for j in range(2*anOuterR+1):
            proportionsAnO[i][j]=round(proportionsAnO[i][j],2)
    for i in range(2*anInnerR+1):
        for j in range(2*anInnerR+1):
            proportionsAnI[i][j] = round(proportionsAnI[i][j],2)
    for i in range(2*anOuterR+1):
        for j in range(2*anOuterR+1):
            if proportionsAnO[i][j] < 1:
                annulusSq[i][j] = 0
    for y in range(anOuterR-anInnerR,(anOuterR*2+1)-anInnerR):
        for x in range(anOuterR-anInnerR, (anOuterR*2+1)-anInnerR):
            if proportionsAnI[y-anInnerR][x-anInnerR] > 0:
                annulusSq[y][x] = 0
    #calculate the average Sky ADU Count
    annulusValues = []
    skyCount = 0
    for i in range(2*anOuterR+1):
        for j in range(2*anOuterR+1):
            if abs(annulusSq[i][j] - 0)>1e-4:
                annulusValues.append(annulusSq[i][j])
                skyCount = skyCount + 1
    avSky = sum(annulusValues)/skyCount*1.0
    #calculate the signal
    signal = sum(apValues) - nap*avSky
    #calcualte SNR
    De = 10
    nan = skyCount
    skyE = avSky*.8
    Se = signal*.8
    p2 = 11**2+(.8**2)/12
    SNR = sqrt(Se)/sqrt(1+nap*(1+ nap*1.0/nan)*((skyE+De+p2)/Se))
    #calculate instant magnitude
    instMag = -2.5*log(signal,10)
    uncertainty = 1.087/SNR
    print(nap,nan)
    return "Border Pixels Accepted: " + "\nsignal: " + str(signal) + " SNR: " + str(SNR) + "\ninst mag: " + str(instMag) + " +/- " + str(uncertainty)


def apPhotometryFractional(filename, xCoord, yCoord, apR, anInnerR, anOuterR, darkFileName):
    fileContents = fits.getdata(filename)
    apertureSq = np.array(fileContents[yCoord-apR-1:yCoord+apR, xCoord-apR-1:xCoord+apR])
    #calculate the signal, SNR, and magnitude if border pixels accepted
    proportionsAp = np.array(fractional(apR))
    for i in range(2*apR+1):
        for j in range(apR*2+1):
            apertureSq[i][j] = proportionsAp[i][j]*apertureSq[i][j]
    apertureValues = []
    nap = np.sum(proportionsAp)
    for i in range(apR*2+1):
        for j in range(apR*2+1):
            if proportionsAp[i][j] != 0:
                apertureValues.append(round(apertureSq[i][j],3))
    apertureVs = np.array(apertureValues)
    #create an annulus for the aperture
    annulusSq = fileContents[yCoord-anOuterR-1:yCoord+anOuterR, xCoord-anOuterR-1:xCoord+anOuterR]
    proportionsAnO = fractional(anOuterR)
    annulusInnerSq = fileContents[yCoord-anInnerR-1:yCoord+anInnerR, xCoord-anInnerR-1:xCoord+anInnerR]
    proportionsAnI = fractional(anInnerR)
    annulusInner = annulusInnerSq*proportionsAnI
    for i in range(2*anOuterR+1):
        for j in range(2*anOuterR+1):
            if abs(proportionsAnO[i][j] - 1) >1e-4:
                annulusSq[i][j] = 0
    for y in range(anOuterR-anInnerR,(anOuterR*2+1)-anInnerR):
        for x in range(anOuterR-anInnerR, (anOuterR*2+1)-anInnerR):
            if proportionsAnI[y-anInnerR][x-anInnerR] != 0:
                annulusSq[y][x] = 0
    #calculate the average Sky ADU Count
    annulusValues = []
    skyCount = 0
    for i in range(2*anOuterR+1):
        for j in range(2*anOuterR+1):
            if annulusSq[i][j] != 0:
                annulusValues.append(annulusSq[i][j])
                skyCount = skyCount + 1
    avSky = sum(annulusValues)/skyCount*1.0
    
    #calculate signal in ADU
    signal = np.sum(apertureVs) - nap*avSky
    #calcualte SNR
    De = De(darkFileName)
    nan = skyCount
    skyE = avSky*.8
    Se = signal*.8
    p2 = 11**2+(.8**2)/12
    SNR = sqrt(Se)/sqrt(1+nap*(1+ nap*1.0/nan)*((skyE+De+p2)/Se))
    #calculate instant magnitude
    instMag = -2.5*log(signal,10)
    uncertainty = 1.087/SNR
    return "Border Pixels Fractional: " + "\nsignal: " + str(signal) + " SNR: " + str(SNR) + "\ninst mag: " + str(instMag) + " +/- " + str(uncertainty)
    
#print(apPhotometryFractional("aptest.FIT", 490, 293, 5, 8, 13))
#print(apPhotometryRejected("aptest.FIT", 490, 293, 5, 8, 13))
#print(apPhotometryAccepted("aptest.FIT", 490, 293, 5, 8, 13))
print(apPhotometryFractional("1998OH_OBS1_SEQ2.00000001.Entered Coordinates.red.fits", 419, 702, 4, 8, 12, "mfdark1.fits"))
print(apPhotometryFractional("1998OH_OBS1_SEQ2.00000001.Entered Coordinates.red.fits", 244, 298, 3, 7, 12, "mfdark1.fits"))
print(apPhotometryFractional("1998OH_OBS1_SEQ2.00000001.Entered Coordinates.red.fits", 408, 342, 3, 7, 12, "mfdark1.fits"))
print(apPhotometryFractional("1998OH_OBS1_SEQ2.00000001.Entered Coordinates.red.fits", 463, 180, 6, 10, 15, "mfdark1.fits"))
print(apPhotometryFractional("1998OH_OBS1_SEQ2.00000001.Entered Coordinates.red.fits", 594, 219, 5, 9, 13, "mfdark1.fits"))
print(apPhotometryFractional("1998OH_OBS1_SEQ2.00000001.Entered Coordinates.red.fits", 529, 398, 4, 9, 13, "mfdark1.fits"))
print(apPhotometryFractional("1998OH_OBS1_SEQ2.00000001.Entered Coordinates.red.fits", 758, 419, 4, 9, 13, "mfdark1.fits"))
