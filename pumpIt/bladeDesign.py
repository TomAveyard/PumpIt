import imp
import sys
from meridional import Meridional
from math import cos, atan2, degrees, radians, sqrt, tan, pi, asin, sin
from plottingHelper import Bezier, findIntersectionOfCoords, polarToCartesian

class Blade:

    def __init__(self, meridionalSection: Meridional, numberOfStreamlines: int = 3, bladeDevelopmentControlPoints: list = [[0,0], [1,1]], bladeLoadingBandwidth=15) -> None:
        
        self.meridionalSection = meridionalSection
        self.numberOfStreamlines = numberOfStreamlines
        self.bladeDevelopmentControlPoints = bladeDevelopmentControlPoints
        self.bladeLoadingBandwidth = bladeLoadingBandwidth

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
                
                bladeAnglej = self.meridionalSection.impeller.inletBladeAngle - (ts[j] * (self.meridionalSection.impeller.inletBladeAngle - self.meridionalSection.impeller.outletBladeAngle))
                bladeAngles.append(bladeAnglej)

                deltaMj = sqrt(((xCoords[j] - xCoords[j+1]) ** 2) + ((yCoords[j] - yCoords[j+1]) ** 2))
                deltaUj = deltaMj / tan(radians(bladeAnglej))
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

        self.outletBladeInclination = self.meridionalSection.impeller.outletBladeInclination
        self.outletBladeTwist = 90 - self.outletBladeInclination

        for i in range(len(self.streamlinesBladeAngles)):

            if self.outletBladeTwist != 0:
                width = abs(self.streamlinesMeridionalXCoords[0][-1] - self.streamlinesMeridionalXCoords[i][-1])
                deltaU = width * tan(radians(self.outletBladeTwist))
                epsilonsch = (360 * deltaU) / (2 * pi * (self.meridionalSection.impeller.impellerOutletDiameter/2))
            else:
                epsilonsch = 0

            r = self.streamlinesMeridionalYCoords[i][-1]
            rs = [r]
            epsilonschs = [epsilonsch]
            epsilonschsRadians = [radians(epsilonsch)]

            for j in range(len(self.streamlinesDeltaRs[i])):

                r -= self.streamlinesDeltaRs[i][j]
                rs.append(r)
                deltaEpsilonsch = (360 * self.streamlinesDeltaUs[i][j]) / (2 * pi * r)
                epsilonsch -= deltaEpsilonsch
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
        self.bladesAxialCoords = []

        for i in range(len(self.streamlinesEpsilonSchs)):
            intersectionIdx = self.streamlinesLEIntersectionIdxs[i] + 1
            self.bladesEpsilonSchs.append(self.streamlinesEpsilonSchs[i][0:-intersectionIdx])
            self.bladesEpsilonSchsRadians.append(self.streamlinesEpsilonSchsRadians[i][0:-intersectionIdx])
            self.bladesRadiuses.append(self.streamlinesRadiuses[i][0:-intersectionIdx])
            self.bladesAxialCoords.append(self.streamlinesMeridionalXCoords[i][0:-intersectionIdx])
            self.bladesAxialCoords[i] = self.bladesAxialCoords[i][::-1]

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

        # Calculate blade lengths

        self.bladeLengths = []
        for i in range(len(self.bladesAxialCoords)):
            length = 0
            for j in range(1, len(self.bladesAxialCoords[i])-1):
                deltaR = self.bladesRadiuses[i][j] - self.bladesRadiuses[i][j-1]
                deltaZ = self.bladesAxialCoords[i][j] - self.bladesAxialCoords[i][j-1]
                deltaM = sqrt(deltaR ** 2 + deltaZ ** 2)
                deltaEpsilon = self.bladesEpsilonSchs[i][j] - self.bladesEpsilonSchs[i][j-1]
                averageR = (self.bladesRadiuses[i][j] + self.bladesRadiuses[i][j-1]) / 2
                deltaU = deltaEpsilon * 2 * pi * averageR / 360
                deltaL = sqrt(deltaM ** 2 + deltaU ** 2)
                length += deltaL
            self.bladeLengths.append(length)

        self.wrapAngles = []
        for i in range(len(self.bladesEpsilonSchs)):
            self.wrapAngles.append(abs(self.bladesEpsilonSchs[i][-1]))

        self.LELengthFractions = []
        for i in [0, len(self.streamlineMeridionalLengths)-1]:
            length = 0
            for j in range(1, self.streamlinesLEIntersectionIdxs[i]):
                length += sqrt((self.streamlinesMeridionalXCoords[i][j] - self.streamlinesMeridionalXCoords[i][j-1]) ** 2 + (self.streamlinesMeridionalYCoords[i][j] - self.streamlinesMeridionalYCoords[i][j-1]) ** 2)
            self.LELengthFractions.append(length / self.streamlineMeridionalLengths[i])

        # Calculate effective blade loading

        impeller = self.meridionalSection.impeller
        minimumNonDimensionalBladeLength = min(self.bladeLengths) / impeller.d2
        w1NonDimensional = impeller.w1 / impeller.u2
        w2NonDimensional = impeller.w2 / impeller.u2

        self.effectiveBladeLoading = (2 * pi * impeller.headCoefficient) \
            / (impeller.hydraulicEfficiency * impeller.numberOfBlades * minimumNonDimensionalBladeLength \
                * (w1NonDimensional + w2NonDimensional))
            
        self.allowableBladeLoading = (40 / impeller.nq) ** 0.77
        bladeLoadingDelta = self.allowableBladeLoading * (self.bladeLoadingBandwidth / 100)
        self.allowableBladeLoadingBand = [self.allowableBladeLoading - bladeLoadingDelta, self.allowableBladeLoading + bladeLoadingDelta]
        self.recommendedBladeLoading = self.allowableBladeLoading - (self.allowableBladeLoading * 0.1)

        if self.effectiveBladeLoading < self.allowableBladeLoadingBand[0] or self.effectiveBladeLoading > self.allowableBladeLoadingBand[1]:

            print("Warning: Effective blade loading not within recommended allowable band")
            print("Effective blade loading: " + str(round(self.effectiveBladeLoading, 5)))
            print("Allowable blade loading: " + str(round(self.allowableBladeLoadingBand[0], 5)) + " to " + str(round(self.allowableBladeLoadingBand[1], 5)))
            print("It is recommended to aim for a blade loading of " + str(round(self.recommendedBladeLoading, 2)) + "\n")

        # Calculate blade distance at outlet

        rotation = 360 / impeller.numberOfBlades
        secondBladeEpsilonschs = [epsilon + rotation for epsilon in self.bladesEpsilonSchs[0]]
        minDistance = 10000000
        for i in range(len(secondBladeEpsilonschs)):
            distance = sqrt((self.bladesRadiuses[0][-1] ** 2 + self.bladesRadiuses[0][i] ** 2) - (2 * self.bladesRadiuses[0][-1] * self.bladesRadiuses[0][i] * cos(radians(secondBladeEpsilonschs[i]))))
            if distance < minDistance:
                minDistance = distance
        
        self.outletBladeDistance = minDistance
        self.a2 = self.outletBladeDistance # Alias
        
        self.betaa2 = degrees(asin(self.a2 / impeller.t2))
        self.outletAngleRatio = sin(radians(self.betaa2)) / sin(radians(impeller.beta2B))
        if self.outletAngleRatio < 0.7 or self.outletAngleRatio > 0.9:
            print("Warning: Outlet angle ratio should be between 0.7 and 0.9")
            print("Outlet angle ratio: " + str(round(self.outletAngleRatio, 2)) + "\n")
