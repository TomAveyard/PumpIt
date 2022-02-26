import sys
from meridional import Meridional
from math import cos, atan2, degrees, radians, sqrt, tan, pi
from plottingHelper import Bezier
import matplotlib.pyplot as plt

class Blade:

    def __init__(self, meridionalSection: Meridional, numberOfStreamlines: int = 3, bladeDevelopmentControlPoints: list = [[0,0], [1,1]]) -> None:
        
        self.meridionalSection = meridionalSection
        self.numberOfStreamlines = numberOfStreamlines
        self.bladeDevelopmentControlPoints = bladeDevelopmentControlPoints

        outerStreamlineXCoords = self.meridionalSection.outerStreamlineXCoords
        outerStreamlineYCoords = self.meridionalSection.outerStreamlineYCoords
        innerStreamlineXCoords = self.meridionalSection.innerStreamlineXCoords
        innerStreamlineYCoords = self.meridionalSection.innerStreamlineYCoords
        
        if self.bladeDevelopmentControlPoints[0] != [0,0] or bladeDevelopmentControlPoints[-1] != [1,1]:

            sys.exit("Error: First and last blade development control points are required to be [0,0] and [1,1] respectively")

        # Initialise lists to hold coords for each streamline
        self.streamlinesXCoords = []
        self.streamlinesYCoords = []
        # Append on the number of lists required
        for i in range(self.numberOfStreamlines):
            self.streamlinesXCoords.append([])
            self.streamlinesYCoords.append([])

        i = 0
        while i < self.meridionalSection.numberOfPoints:
            
            # Calculate inner to outer streamline geometry
            deltax = abs(innerStreamlineXCoords[i] - outerStreamlineXCoords[i])
            deltay = abs(outerStreamlineYCoords[i] - innerStreamlineYCoords[i])
            theta = degrees(atan2(deltax, deltay))
            b = sqrt((deltax ** 2) + (deltay ** 2))
            
            streamlinesY = [outerStreamlineYCoords[i]]
            streamlinesX = [outerStreamlineXCoords[i]]
            streamTubeWidths = []

            for j in range(1, self.numberOfStreamlines + 1):
                
                # Calculate coords for each streamline so that each streamtube formed by the streamlines has an equal amount of flow through it
                # Meridional velocity cm is assumed constant
                streamlinesY.append(sqrt((streamlinesY[j-1] ** 2) - ((b * cos(radians(theta)) *  (outerStreamlineYCoords[i] + innerStreamlineYCoords[i])) / (self.numberOfStreamlines + 1))))
                streamlinesX.append(streamlinesX[j-1] + (abs(streamlinesY[j-1] - streamlinesY[j])) * tan(radians(theta)))
            
            streamlinesY.append(innerStreamlineYCoords[i])

            for j in range(1, len(streamlinesY)-1):

                width = sqrt(((streamlinesX[j] - streamlinesX[j-1]) ** 2) + ((streamlinesY[j-1] - streamlinesY[j]) ** 2))
                streamTubeWidths.append(width)

                self.streamlinesXCoords[j-1].append(streamlinesX[j])
                self.streamlinesYCoords[j-1].append(streamlinesY[j])

            i+=1

        self.streamlinesXCoords.insert(0, outerStreamlineXCoords)
        self.streamlinesYCoords.insert(0, outerStreamlineYCoords)
        self.streamlinesXCoords.append(innerStreamlineXCoords)
        self.streamlinesYCoords.append(innerStreamlineYCoords)

        self.bladeDevelopmentBezier = Bezier(self.bladeDevelopmentControlPoints)
        self.streamlineMeridionalLengths = []
        self.streamlineCircumferentialLengths = []
        self.streamlineTotalLengths = []
        self.streamlinesBladeAngles = []
        self.streamlinesDeltaMs = []
        self.streamlinesDeltaUs = []
        self.streamlinesDeltaLs = []
        self.streamlinesDeltaZs = []
        self.streamlinesDeltaRs = []

        tIncrement = 1 / self.meridionalSection.numberOfPoints
        t = 0
        ts = []
        while t <= 1:
            ts.append(t)
            t += tIncrement

        for i in range(len(self.streamlinesXCoords)):
            
            xCoords = self.streamlinesXCoords[i]
            yCoords = self.streamlinesYCoords[i]

            bladeAngles = []
            deltaMs = []
            deltaUs = []
            deltaLs = []
            deltaZs = []
            deltaRs = []
            streamlineMeridionalLength = 0
            streamlineCircumferentialLength = 0
            streamlineTotalLength = 0

            for j in range(self.meridionalSection.numberOfPoints-2, 0, -1):
                
                bladeAnglej = self.meridionalSection.impeller.outletBladeAngle - ts[j] * (self.meridionalSection.impeller.outletBladeAngle - self.meridionalSection.impeller.inletBladeAngle)
                bladeAngles.append(bladeAnglej)

                deltaMj = sqrt(((xCoords[j] - xCoords[j+1]) ** 2) + ((yCoords[j] - yCoords[j+1]) ** 2))
                deltaUj = -deltaMj / tan(radians(bladeAnglej))
                deltaLj = sqrt(deltaMj ** 2 + deltaUj ** 2)
                deltaZj = xCoords[j] - xCoords[j+1]
                deltaRj = sqrt(deltaMj ** 2 - deltaZj ** 2)

                deltaMs.append(deltaMj)
                deltaUs.append(deltaUj)
                deltaLs.append(deltaLj)
                deltaZs.append(deltaZj)
                deltaRs.append(deltaRj)

                streamlineMeridionalLength += deltaMj
                streamlineCircumferentialLength += deltaUj
                streamlineTotalLength += deltaLj

            self.streamlineMeridionalLengths.append(streamlineMeridionalLength)
            self.streamlineCircumferentialLengths.append(streamlineCircumferentialLength)
            self.streamlineTotalLengths.append(streamlineTotalLength)

            self.streamlinesBladeAngles.append(bladeAngles)
            self.streamlinesDeltaMs.append(deltaMs)
            self.streamlinesDeltaUs.append(deltaUs)
            self.streamlinesDeltaLs.append(deltaLs)
            self.streamlinesDeltaZs.append(deltaZs)
            self.streamlinesDeltaRs.append(deltaRs)

        self.streamlinesEpsilonSchs = []
        self.streamlinesEpsilonSchsRadians = []
        self.streamlinesRadiuses = []

        for i in range(len(self.streamlinesBladeAngles)):

            epsilonsch = 0

            r = self.streamlinesYCoords[i][-1]
            rs = [r]
            epsilonschs = [epsilonsch]
            epsilonschsRadians = [radians(epsilonsch)]

            for j in range(len(self.streamlinesDeltaRs[i])):

                r -= self.streamlinesDeltaRs[i][j]
                rs.append(r)
                deltaEpsilonsch = 360 * self.streamlinesDeltaUs[i][j] / (2 * pi * r)
                epsilonsch += deltaEpsilonsch
                epsilonschs.append(epsilonsch)
                epsilonschsRadians.append(radians(epsilonsch))

            self.streamlinesEpsilonSchs.append(epsilonschs)
            self.streamlinesEpsilonSchsRadians.append(epsilonschsRadians)
            self.streamlinesRadiuses.append(rs)
