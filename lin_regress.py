import matplotlib.pyplot as plt
from math import*
import numpy as np
from scipy import stats

def astMag(x, y, xErrors, astInst):
    slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
    astMag = slope*astInst + intercept
    print("Asteroid Magnitude: " + str(round(astMag, 2)))
    print("slope: %.3f    intercept: %.3f" % (slope, intercept))
    print("R-squared: %f" % r_value**2)
    plt.plot(x, y, 'o', label = "Raw Data")
    plt.xlabel("Instrumental Magnitude")
    plt.ylabel("Apparent Visual Magnitude")
    plt.plot(x, intercept + slope*x, 'b', label = "Curve Fit")
    plt.errorbar(x, y, xerr=xErrors, yerr = None, fmt=" ", label = "Error Bars")
    plt.legend()
    plt.show()
    
def standardDeviations(x, y):
    slope, intercept, r, prob2, see = stats.linregress(x,y)
    xMean = x.mean()
    sx2 = ((x-xMean)**2).sum()
    sd_intercept = see*sqrt(1./len(x) + xMean*xMean/sx2)
    sd_slope = see*sqrt(1.0/sx2)
    print("Intercept Standard Deviation: " + str(round(sd_intercept,3)) + "\nSlope Standard Deviation: " + str(round(sd_slope, 2)))
    
x1 = np.array([-11.866, -10.478, -14.508, -12.931, -11.593, -12.215])
y1 = np.array([13.22, 14.574, 10.62, 12.237, 13.566, 12.924])
xUnc = np.array([.006, .014, .002, .004, .008, .005])
print("obst: " + str(standardDeviations(x1,y1)))
print("obs1: " + str(astMag(x1,y1,xUnc, -8.546)))
print("*********************************************************")
##x2 = np.array([-10.112, -9.0967, -10.261, -9.376, -11.555, -9.683])
##y2 = np.array([15.023, 15.972, 14.816, 15.717, 13.58, 15.395])
##print("obst: " + str(standardDeviations(x2,y2)))
##print("obs1: " + str(astMag(x2,y2,-8.645)))
print("*********************************************************")
x3 = np.array([-8.267, -9.883, -10.734, -10.866, -10.334, -10.904, -11.173])
y3 = np.array([16.639, 15.066, 14.234, 14.122, 14.626, 14.077, 13.838])
xUnc = np.array([.084, .022, .012, .012, .016, .012, .01])
print("obs2: " + str(standardDeviations(x3,y3)))
print("obs2: " + str(astMag(x3,y3,xUnc, -8.358)))

