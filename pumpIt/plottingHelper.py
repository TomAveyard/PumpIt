import numpy as np
from math import sin, cos, radians, degrees, sqrt, atan2, tan, atan, acos, atan2
import matplotlib.pyplot as plt
import sys
from scipy.interpolate import interp1d

def plotStraightLine(startCoords: tuple, endCoords: tuple, numberOfPoints: int):

    xCoords = np.linspace(startCoords[0], endCoords[0], numberOfPoints)
    yCoords = np.linspace(startCoords[1], endCoords[1], numberOfPoints)

    coords = np.array([xCoords, yCoords])

    return coords

def plotArc(startCoords: tuple, originCoords: tuple, angle: float, numberOfPoints: int):
    
    deltax =  originCoords[0] - startCoords[0]
    deltay =  originCoords[1] - startCoords[1]
    radiusOfCurvature = sqrt((deltax ** 2) + (deltay ** 2))

    startAngle = -90 - degrees(atan2(deltax, deltay))
    endAngle = startAngle + angle

    angles = np.linspace(startAngle, endAngle, numberOfPoints)
    xCoords = np.zeros(numberOfPoints)
    yCoords = np.zeros(numberOfPoints)

    for i in range(numberOfPoints):

        xCoords[i] = originCoords[0] + radiusOfCurvature * cos(radians(angles[i]))
        yCoords[i] = originCoords[1] + radiusOfCurvature * sin(radians(angles[i]))

    coords = np.array([xCoords, yCoords])

    return coords

def findIntersection(point1, angle1, point2, angle2):

    m1 = tan(radians(angle1))
    m2 = tan(radians(angle2))

    x = (point2[1] - point1[1] + (m1 * point1[0]) - (m2 * point2[0])) / (m1 - m2)
    y = m1 * x + (point1[1] - m1 * point1[0])

    return [x,y]

def simplePlot(xCoords, yCoords, show=True, equalAxes=True, plotPoints=False):

    ax = plt.axes()

    if type(xCoords) is list:
        for i in range(len(xCoords)):
            ax.plot(xCoords[i], yCoords[i])
    else:
        ax.plot(xCoords, yCoords)
    if plotPoints:
        ax.scatter(xCoords, yCoords)
    if equalAxes:
        ax.axis("equal")
    if show:
        plt.show()

def lengthBetween2Points(point1, point2):

    return sqrt(((point1[0] - point2[0]) ** 2) + ((point1[1] - point2[1]) ** 2))

def cosineRuleAngle(a, b, c):

    return degrees(acos((a ** 2 + b ** 2 - c ** 2) / (2 * a * b)))

def findIntersectionOfCoords(line1, line2, method="interpolate", side="left", returnIndex=False, n=10000):

    if type(line1) == list:
        line1 = np.array(line1)
    if type(line2) == list:
        line2 = np.array(line2)

    f1 = interp1d(line1[0], line1[1], kind='linear')
    f2 = interp1d(line2[0], line2[1], kind='linear')

    xx = np.linspace(max(line1[0][0],  line2[0][0]), min(line1[0][-1], line2[0][-1]), n)

    line1Interp = f1(xx)
    line2Interp = f2(xx)
    
    idxInterp = np.argwhere(np.diff(np.sign(line1Interp - line2Interp))).flatten()
    
    if method.lower() == "interpolate" or method.lower() == "interp":

        if returnIndex:
            return (line1Interp[idxInterp][0], xx[idxInterp][0]), idxInterp
        else:
            return (line1Interp[idxInterp][0], xx[idxInterp][0])

    elif method.lower() == "nearest" or method.lower() == "nearest point":

        idx = np.searchsorted(line1[1], line1Interp[idxInterp][0])

        if side.lower() != "right":

            idx -= 1

        if returnIndex:
            return (line1[0][idx], line1[1][idx]), idx
        else:
            return (line1[0][idx], line1[1][idx])

def polarToCartesian(r, theta):

    x = r * cos(radians(theta))
    y = r * sin(radians(theta))

    return [x, y]

def polarToCartesianLines(rs, thetas):
    xs = []
    ys = []
    for i in range(len(rs)):
        cartesianCoord = polarToCartesian(rs[i], thetas[i])
        xs.append(cartesianCoord[0])
        ys.append(cartesianCoord[1])

    return xs, ys

def cartesianToPolar(x, y):

    r = sqrt(x**2 + y**2)
    theta = degrees(atan2(y, x))

    return [r, theta]

def cartesianToPolarLines(xs, ys):

    rs = []
    thetas = []
    for i in range(len(xs)):
        polarCoord = cartesianToPolar(xs[i], ys[i])
        rs.append(polarCoord[0])
        thetas.append(polarCoord[1])

    return rs, thetas

def rotateCartesianCoord(x, y, rotation, origin=[0,0]):

    polarCoords = cartesianToPolar(x-origin[0], y-origin[1])
    rotatedCartesianCoords = polarToCartesian(polarCoords[0], polarCoords[1] + rotation)
    rotatedCartesianCoords[0] += origin[0]
    rotatedCartesianCoords[1] += origin[1]

    return rotatedCartesianCoords

# Generates a bezier curve given a list of control point coordinates
class Bezier:

    def __init__(self, controlPoints: list):
        
        self.controlPoints = controlPoints
        self.n = len(controlPoints) - 1
        self.length = None

    # Binomial lookup table, more efficient than generating each time
    def binomial(self, k: float) -> float:

        self.pascalsTriangle = [[1],
                                [1,1],
                                [1,2,1],
                                [1,3,3,1],
                                [1,4,6,4,1],
                                [1,5,10,10,5,1],
                                [1,6,15,20,15,6,1]]

        # Generates new lines if it's not already in the above table
        while self.n > len(self.pascalsTriangle):
            nextRow = []
            nextRow.append(1)
            for i in range(len(self.pascalsTriangle[-1]) - 1):
                nextRow.append(self.pascalsTriangle[-1][i] + self.pascalsTriangle[-1][i+1])
            nextRow.append(1)
            self.pascalsTriangle.append(nextRow)

        return self.pascalsTriangle[self.n][k]

    # Finds the x/y coordinate for a given t, intended as a helper function for the getPoint function 
    def bezier(self, t: float, coord: str) -> float:

        sum = 0

        for i in range(self.n + 1):

            if coord == 'x':
                weight = self.controlPoints[i][0]
            elif coord == 'y':
                weight = self.controlPoints[i][1]
            else:
                sys.exit("Invalid coordinate type given")

            sum = sum + weight * self.binomial(i) * ((1-t) ** (self.n - i)) * (t ** i)

        return sum

    # Returns the x,y coordinates for a given t
    def getPoint(self, t: float) -> list:

        if t > 1 or t < 0:
            sys.exit("t has to be 0 <= t <= 1")

        x = self.bezier(t, 'x')
        y = self.bezier(t, 'y')

        return [x,y]

    # Estimates length by adding up discrete sections
    def estimateLength(self, numberOfPoints: int = 100) -> float:

        ts = np.linspace(0, 1, numberOfPoints)
        sum = 0
        i = 1
        prevPoint = self.getPoint(0)
        while i < numberOfPoints:
            currentPoint = self.getPoint(ts[i])
            sum += sqrt(((currentPoint[0] - prevPoint[0]) ** 2) + ((currentPoint[1] - prevPoint[1]) ** 2))

            prevPoint = currentPoint
            i += 1
        
        self.length = sum
        return sum
