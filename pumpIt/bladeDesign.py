import sys
from meridional import Meridional
from math import cos, atan2, degrees, radians, sqrt, tan, pi
from plottingHelper import Bezier, findIntersectionOfCoords, polarToCartesian

class Blade:

    def __init__(self, meridionalSection: Meridional, numberOfStreamlines: int = 3, bladeDevelopmentControlPoints: list = [[0,0], [1,1]]) -> None:
        
        self.meridionalSection = meridionalSection
        self.numberOfStreamlines = numberOfStreamlines
        self.bladeDevelopmentControlPoints = bladeDevelopmentControlPoints

        outerStreamlineMeridionalXCoords = self.meridionalSection.outerStreamlineXCoords
        outerStreamlineMeridionalYCoords = self.meridionalSection.outerStreamlineYCoords
        innerStreamlineMeridionalXCoords = self.meridionalSection.innerStreamlineXCoords
        innerStreamlineMeridionalYCoords = self.meridionalSection.innerStreamlineYCoords
        
        if self.bladeDevelopmentControlPoints[0] != [0,0] or bladeDevelopmentControlPoints[-1] != [1,1]:

            sys.exit("Error: First and last blade development control points are required to be [0,0] and [1,1] respectively")

        # Initialise lists to hold coords for each streamline
        self.streamlinesMeridionalXCoords = []
        self.streamlinesMeridionalYCoords = []
        # Append on the number of lists required
        for i in range(self.numberOfStreamlines):
            self.streamlinesMeridionalXCoords.append([])
            self.streamlinesMeridionalYCoords.append([])

        i = 0
        while i < self.meridionalSection.numberOfPoints:
            
            # Calculate inner to outer streamline geometry
            deltax = abs(innerStreamlineMeridionalXCoords[i] - outerStreamlineMeridionalXCoords[i])
            deltay = abs(outerStreamlineMeridionalYCoords[i] - innerStreamlineMeridionalYCoords[i])
            theta = degrees(atan2(deltax, deltay))
            b = sqrt((deltax ** 2) + (deltay ** 2))
            
            streamlinesY = [outerStreamlineMeridionalYCoords[i]]
            streamlinesX = [outerStreamlineMeridionalXCoords[i]]
            streamTubeWidths = []

            for j in range(1, self.numberOfStreamlines + 1):
                
                # Calculate coords for each streamline so that each streamtube formed by the streamlines has an equal amount of flow through it
                # Meridional velocity cm is assumed constant
                streamlinesY.append(sqrt((streamlinesY[j-1] ** 2) - ((b * cos(radians(theta)) *  (outerStreamlineMeridionalYCoords[i] + innerStreamlineMeridionalYCoords[i])) / (self.numberOfStreamlines + 1))))
                streamlinesX.append(streamlinesX[j-1] + (abs(streamlinesY[j-1] - streamlinesY[j])) * tan(radians(theta)))
            
            streamlinesY.append(innerStreamlineMeridionalYCoords[i])

            for j in range(1, len(streamlinesY)-1):

                width = sqrt(((streamlinesX[j] - streamlinesX[j-1]) ** 2) + ((streamlinesY[j-1] - streamlinesY[j]) ** 2))
                streamTubeWidths.append(width)

                self.streamlinesMeridionalXCoords[j-1].append(streamlinesX[j])
                self.streamlinesMeridionalYCoords[j-1].append(streamlinesY[j])

            i+=1

        self.streamlinesMeridionalXCoords.insert(0, outerStreamlineMeridionalXCoords)
        self.streamlinesMeridionalYCoords.insert(0, outerStreamlineMeridionalYCoords)
        self.streamlinesMeridionalXCoords.append(innerStreamlineMeridionalXCoords)
        self.streamlinesMeridionalYCoords.append(innerStreamlineMeridionalYCoords)

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

        # Finds various deltas to achieve blade angle development
        for i in range(len(self.streamlinesMeridionalXCoords)):
            
            xCoords = self.streamlinesMeridionalXCoords[i]
            yCoords = self.streamlinesMeridionalYCoords[i]

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

        # Find the coords of the streamlines in the plan view for a single blade using the deltas
        # Done for one blade, coords for other blades can be found by adding on a degree rotation to the coords for the one blade
        for i in range(len(self.streamlinesBladeAngles)):

            epsilonsch = 0

            r = self.streamlinesMeridionalYCoords[i][-1]
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

        # Find indexes for each streamline where the streamline intersects the LE of the blade

        self.streamlinesLEIntersectionIdxs = []
        self.streamlineLEIntersectionCoords = []

        LECoords = [
            [self.meridionalSection.bladeLEShroudCoords[0], self.meridionalSection.bladeLEHubCoords[0]],
            [self.meridionalSection.bladeLEShroudCoords[1], self.meridionalSection.bladeLEHubCoords[1]]
        ]

        for i in range(len(self.streamlinesMeridionalXCoords)):

            streamlineLine = [
                self.streamlinesMeridionalXCoords[i],
                self.streamlinesMeridionalYCoords[i]
            ]

            # Finds index at which blade LE intersect each streamline
            intersectionCoords, intersectionIdx = findIntersectionOfCoords(streamlineLine, LECoords, method="nearest", side="right", n=10000, returnIndex=True)
            self.streamlinesLEIntersectionIdxs.append(intersectionIdx)
            self.streamlineLEIntersectionCoords.append(intersectionCoords)

        # New coords lists for the blade coords
        self.bladesEpsilonSchs = []
        self.bladesEpsilonSchsRadians = []
        self.bladesRadiuses = []

        for i in range(len(self.streamlinesEpsilonSchs)):

            self.bladesEpsilonSchs.append(self.streamlinesEpsilonSchs[i][0:-intersectionIdx])
            self.bladesEpsilonSchsRadians.append(self.streamlinesEpsilonSchsRadians[i][0:-intersectionIdx])
            self.bladesRadiuses.append(self.streamlinesRadiuses[i][0:-intersectionIdx])

        # Gets coords of LE in plan view
        self.bladeLEEpsilons = []
        self.bladeLEEpsilonsRadians = []
        self.bladeLERadiuses = []
        for i in range(len(self.bladesEpsilonSchs)):
            self.bladeLEEpsilons.append(self.bladesEpsilonSchs[i][-1])
            self.bladeLEEpsilonsRadians.append(self.bladesEpsilonSchsRadians[i][-1])
            self.bladeLERadiuses.append(self.bladesRadiuses[i][-1])

        # Get streamline coords in cartesian coordinates
        self.streamlinesXCoords = []
        self.streamlinesYCoords = []
        for i in range(len(self.streamlinesEpsilonSchs)):
            self.streamlinesXCoords.append([])
            self.streamlinesYCoords.append([])
            for j in range(len(self.streamlinesEpsilonSchs[i])):
                coords = polarToCartesian(self.streamlinesRadiuses[i][j], self.streamlinesEpsilonSchs[i][j])
                self.streamlinesXCoords[i].append(coords[0])
                self.streamlinesYCoords[i].append(coords[1])

        # Get blade coords in cartesian coordinates
        self.bladesXCoords = []
        self.bladesYCoords = []
        for i in range(len(self.bladesEpsilonSchs)):
            self.bladesXCoords.append([])
            self.bladesYCoords.append([])
            for j in range(len(self.bladesEpsilonSchs[i])):
                coords = polarToCartesian(self.bladesRadiuses[i][j], self.bladesEpsilonSchs[i][j])
                self.bladesXCoords[i].append(coords[0])
                self.bladesYCoords[i].append(coords[1])

        # Get blade LE coords in cartesian coordinates
        self.bladeLEXCoords = []
        self.bladeLEYCoords = []
        for i in range(len(self.bladeLEEpsilons)):
            coords = polarToCartesian(self.bladeLERadiuses[i], self.bladeLEEpsilons[i])
            self.bladeLEXCoords.append(coords[0])
            self.bladeLEYCoords.append(coords[1])
