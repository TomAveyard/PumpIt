import numpy as np
from math import sin, cos, radians, degrees, sqrt, atan2
import matplotlib.pyplot as plt

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

def simplePlot(xCoords, yCoords, show=True, equalAxes=True):

    ax = plt.axes()

    ax.plot(xCoords, yCoords)
    if equalAxes:
        ax.axis("equal")
    if show:
        plt.show()
