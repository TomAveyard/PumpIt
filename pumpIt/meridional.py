from impeller import Impeller
import fluid as fl
import plottingHelper as ph
from math import atan2, floor, pi, degrees, cos, sin, radians, atan
import numpy as np
import matplotlib.pyplot as plt

class Meridional:

    def __init__(self, 
    impeller: Impeller,
    outerStreamlineRadiusOfCurvatureFactor: float = 0.7,
    outerStreamlineOutletAngle: float = 0,
    innerStreamlineOutletAngle: float = 0,
    inletAxialSectionFactor: float = 0.25,
    outerStreamlineCircularSectionArcLength: float = 45,
    numberOuterStreamlineControlPoints: int = 1,
    crossSectionVariationType: str = "linear",
    numberOfPoints: int = 400
    ) -> None:
        
        self.impeller = impeller
        self.outerStreamlineRadiusOfCurvatureFactor = outerStreamlineRadiusOfCurvatureFactor
        self.outerStreamlineOutletAngle = outerStreamlineOutletAngle
        self.epsilonDS = self.outerStreamlineOutletAngle
        self.innerStreamlineOutletAngle = innerStreamlineOutletAngle
        self.epsilonTS = self.innerStreamlineOutletAngle
        self.inletAxialSectionFactor = inletAxialSectionFactor
        self.outerstreamlineCircularSectionArcAngle = outerStreamlineCircularSectionArcLength
        self.numberOuterStreamlineControlPoints = numberOuterStreamlineControlPoints
        self.controlPoints = []
        self.crossSectionVariationType = crossSectionVariationType.lower()
        self.numberOfPoints = numberOfPoints

        self.outerStreamlineHeight = (self.impeller.impellerOutletDiameter - self.impeller.impellerInletDiameter) * ((self.impeller.nq / 74) ** 1.07)
        self.ZE = self.outerStreamlineHeight # Alias
        self.outerStreamlineRadiusOfCurvature = self.outerStreamlineRadiusOfCurvatureFactor * self.impeller.impellerInletWidth
        self.RDs = self.outerStreamlineRadiusOfCurvature # Alias
        self.inletAxialSection = self.inletAxialSectionFactor * self.impeller.impellerInletWidth
        self.g1 = self.inletAxialSection # Alias

        # /// Coords for outer streamline ///

        # Coords for circular section
        arcLength = (self.outerstreamlineCircularSectionArcAngle / 360) * 2 * pi * self.outerStreamlineRadiusOfCurvature
        arcEndCoords = [self.g1 + self.outerStreamlineRadiusOfCurvature * cos(radians(self.outerstreamlineCircularSectionArcAngle)), (self.impeller.impellerInletDiameter / 2) + (self.outerStreamlineRadiusOfCurvature * (1 - sin(radians(self.outerstreamlineCircularSectionArcAngle))))]
        outerStreamlineEndCoords = [self.ZE, self.impeller.impellerOutletDiameter]
        # Generate radial section bezier
        self.controlPoints.append(arcEndCoords)
        if numberOuterStreamlineControlPoints == 1:
            self.controlPoints.append(ph.findIntersection(arcEndCoords, self.outerstreamlineCircularSectionArcAngle, outerStreamlineEndCoords, 270 - self.outerStreamlineOutletAngle))
        self.controlPoints.append(outerStreamlineEndCoords)
        radialSectionBezier = ph.Bezier(self.controlPoints)
        radialSectionLength = radialSectionBezier.estimateLength()
        # Calculate a constant delta m
        outerStreamlineLength = self.g1 + arcLength + radialSectionLength
        deltam = outerStreamlineLength / self.numberOfPoints

        #Initialise for coord finding
        self.outerStreamlineXCoords = np.zeros(self.numberOfPoints)
        self.outerStreamlineYCoords = np.zeros(self.numberOfPoints)

        self.outerStreamlineXCoords[0] = 0
        self.outerStreamlineYCoords[0] = self.impeller.d1 / 2
        i = 1

        # Coords for initial axial section
        while self.outerStreamlineXCoords[i-1] < self.g1:
            self.outerStreamlineXCoords[i] = self.outerStreamlineXCoords[i-1] + deltam
            self.outerStreamlineYCoords[i] = self.impeller.d1 / 2
            i += 1

        # Coords for circular arc
        transitionIndex = i

        extraArcAngle = degrees(atan2(deltam - (self.g1 - self.outerStreamlineXCoords[i-1]), self.outerStreamlineRadiusOfCurvature))
        startCoords = (self.g1 + ((deltam - (self.g1 - self.outerStreamlineXCoords[i-1])) * cos(radians(extraArcAngle))), (self.impeller.d1 / 2) + ((deltam - (self.g1 - self.outerStreamlineXCoords[i-1])) * sin(radians(extraArcAngle))))
        arcNumberOfPoints = floor((arcLength / deltam))
        arcCoords = ph.plotArc(startCoords, (self.g1, (self.impeller.d1 / 2) + self.outerStreamlineRadiusOfCurvature), self.outerstreamlineCircularSectionArcAngle - extraArcAngle, arcNumberOfPoints)
        
        while i < transitionIndex + arcNumberOfPoints:
            self.outerStreamlineXCoords[i] = arcCoords[0][i - transitionIndex]
            self.outerStreamlineYCoords[i] = arcCoords[1][i - transitionIndex]
            i += 1

        # Control point for arc end needs to be shifted to the drawn arc end rather than the exact calculated end as it causes issues drawing the inner streamline
        self.controlPoints = []
        arcEndCoords = [self.outerStreamlineXCoords[i-1], self.outerStreamlineYCoords[i-1]]
        self.controlPoints.append(arcEndCoords)
        if numberOuterStreamlineControlPoints == 1:
            self.controlPoints.append(ph.findIntersection(arcEndCoords, self.outerstreamlineCircularSectionArcAngle, outerStreamlineEndCoords, 270 - self.outerStreamlineOutletAngle))
        self.controlPoints.append(outerStreamlineEndCoords)
        radialSectionBezier = ph.Bezier(self.controlPoints)

        # Coords for radial section
        bezierNumberOfPoints = self.numberOfPoints - i + 1
        tIncrement = 1 / bezierNumberOfPoints
        t = tIncrement
        while i < self.numberOfPoints:
            coord = radialSectionBezier.getPoint(t)
            self.outerStreamlineXCoords[i] = coord[0]
            self.outerStreamlineYCoords[i] = coord[1]
            t += tIncrement
            i += 1
        
        # Impeller width variation
        startb = self.impeller.impellerInletWidth
        endb = self.impeller.impellerOutletWidth
        if crossSectionVariationType == "linear":
            bs = np.linspace(startb, endb, self.numberOfPoints)

        # Initialising loop
        self.innerStreamlineXCoords = np.zeros(self.numberOfPoints)
        self.innerStreamlineYCoords = np.zeros(self.numberOfPoints)
        self.innerStreamlineXCoords[0] = self.outerStreamlineXCoords[0]
        self.innerStreamlineYCoords[0] = self.outerStreamlineYCoords[0] - startb
        i = 1
        # Draws the inner streamline from impeller width variation
        while i < self.numberOfPoints:
            
            b = bs[i]
            angle = degrees(atan2((self.outerStreamlineYCoords[i] - self.outerStreamlineYCoords[i-1]), (self.outerStreamlineXCoords[i] - self.outerStreamlineXCoords[i-1])))
            deltax = b * sin(radians(angle))
            deltay = b * cos(radians(angle))
            self.innerStreamlineXCoords[i] = self.outerStreamlineXCoords[i] + deltax
            self.innerStreamlineYCoords[i] = self.outerStreamlineYCoords[i] - deltay
            i += 1
        
        # Draw blade leading edge
        i = 0
        bladeLEy = self.impeller.d1i / 2
        while self.innerStreamlineYCoords[i] < bladeLEy:
            i += 1
        x2, x1 = self.innerStreamlineXCoords[i], self.innerStreamlineXCoords[i-1]
        y2, y1 = self.innerStreamlineYCoords[i], self.innerStreamlineYCoords[i-1]

        bladeLEx = x2 - (((x2 - x1) * (y2- bladeLEy)) / (y2 - y1))

        self.bladeLEHubCoords = [bladeLEx, bladeLEy]
        self.bladeLEShroudCoords = [self.outerStreamlineXCoords[0], self.outerStreamlineYCoords[0]]
        self.bladeLEAngle = 180 - degrees(atan2((self.bladeLEShroudCoords[1] - self.bladeLEHubCoords[1]), (self.bladeLEShroudCoords[0] - self.bladeLEHubCoords[0])))
        self.epsilonEK = self.bladeLEAngle # Alias

        if self.epsilonEK < 30 or self.epsilonEK > 40:
            print("Warning: Blade angle is " + str(round(self.epsilonEK, 2)) + " degrees. Recommended value is between 30 and 40 degrees")

        self.bladeTEHubCoords = [self.innerStreamlineXCoords[-1], self.innerStreamlineYCoords[-1]]
        self.bladeTEShroudCoords = [self.outerStreamlineXCoords[-1], self.outerStreamlineYCoords[-1]]
        
    def plot(self, show=True, full=False, shaft=True):

        ax = plt.axes()
        ax.plot(self.outerStreamlineXCoords, self.outerStreamlineYCoords, color="black")
        ax.plot(self.innerStreamlineXCoords, self.innerStreamlineYCoords, color="black")
        ax.plot([self.bladeLEShroudCoords[0], self.bladeLEHubCoords[0]], [self.bladeLEShroudCoords[1], self.bladeLEHubCoords[1]], color="grey")
        ax.plot([self.bladeTEShroudCoords[0], self.bladeTEHubCoords[0]], [self.bladeTEShroudCoords[1], self.bladeTEHubCoords[1]], color="grey")

        if shaft:
            ax.plot([self.outerStreamlineXCoords[0], self.innerStreamlineXCoords[-1]], [self.impeller.shaftDiameter/2, self.impeller.shaftDiameter/2], color="grey", ls="--")

        ax.axis("equal")
        if full:
            ax.plot(self.outerStreamlineXCoords, -self.outerStreamlineYCoords, color="black")
            ax.plot(self.innerStreamlineXCoords, -self.innerStreamlineYCoords, color="black")
            ax.plot([self.bladeLEShroudCoords[0], self.bladeLEHubCoords[0]], [-self.bladeLEShroudCoords[1], -self.bladeLEHubCoords[1]], color="grey")
            ax.plot([self.bladeTEShroudCoords[0], self.bladeTEHubCoords[0]], [-self.bladeTEShroudCoords[1], -self.bladeTEHubCoords[1]], color="grey")
            if shaft:
                ax.plot([self.outerStreamlineXCoords[0], self.innerStreamlineXCoords[-1]], [-self.impeller.shaftDiameter/2, -self.impeller.shaftDiameter/2], color="grey", ls="--")
        else:
            ax.set_ylim(ymin=0)

        if show:
            plt.show()
            


fluid = fl.Fluid(density=787, viscosity=2.86e-3, vapourPressure=4100)
impeller = Impeller(suctionSpecificSpeedEU=650,
            suctionSidePressure=3e5,
            kgPerSec=1.32,
            headRise=295,
            fluid=fluid,
            shaftAllowableShearStress=8e7,
            numberOfBlades=6,
            approachFlowAngle=90,
            outletBladeAngle=22.5,
            inletBladeInnerDiameterRatio=1.15
)
meridionalSection = Meridional(impeller)
meridionalSection.plot()
