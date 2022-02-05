from impeller import Impeller
import fluid as fl
import plottingHelper as ph
from math import atan2, floor, pi, degrees, cos, sin, radians
import numpy as np

class Meridional:

    def __init__(self, 
    impeller: Impeller,
    outerStreamlineRadiusOfCurvatureFactor: float = 0.7,
    outerStreamlineOutletAngle: float = 0,
    innerStreamlineOutletAngle: float = 0,
    inletAxialSectionFactor: float = 0.25,
    crossSectionVariationType: str = "linear",
    numberOfPoints: int = 100
    ) -> None:
        
        self.impeller = impeller
        self.outerStreamlineRadiusOfCurvatureFactor = outerStreamlineRadiusOfCurvatureFactor
        self.outerStreamlineOutletAngle = outerStreamlineOutletAngle
        self.epsilonDS = self.outerStreamlineOutletAngle
        self.innerStreamlineOutletAngle = innerStreamlineOutletAngle
        self.epsilonTS = self.innerStreamlineOutletAngle
        self.inletAxialSectionFactor = inletAxialSectionFactor
        self.crossSectionVariationType = crossSectionVariationType.lower()
        self.numberOfPoints = numberOfPoints

        self.outerStreamlineHeight = (self.impeller.impellerOutletDiameter - self.impeller.impellerInletDiameter) * (self.impeller.nq / 74) ** 1.07
        self.ZE = self.outerStreamlineHeight # Alias
        self.outerStreamlineRadiusOfCurvature = self.outerStreamlineRadiusOfCurvatureFactor * self.impeller.impellerInletWidth
        self.RDs = self.outerStreamlineRadiusOfCurvature # Alias
        self.inletAxialSection = self.inletAxialSectionFactor * self.impeller.impellerInletWidth
        self.g1 = self.inletAxialSection # Alias

        # Coords for outer streamline

        arcAngle = 90 - self.outerStreamlineOutletAngle
        arcLength = (arcAngle / 360) * 2 * pi * self.outerStreamlineRadiusOfCurvature
        arcEndY = self.outerStreamlineRadiusOfCurvature - self.outerStreamlineRadiusOfCurvature * cos(radians((arcAngle)))
        radialSectionLength = (self.impeller.d2 / 2) - arcEndY
        outerStreamlineLength = self.g1 + arcLength + radialSectionLength
        deltam = outerStreamlineLength / self.numberOfPoints

        self.outerStreamlineXCoords = np.zeros(self.numberOfPoints)
        self.outerStreamlineYCoords = np.zeros(self.numberOfPoints)

        self.outerStreamlineXCoords[0] = 0
        self.outerStreamlineYCoords[0] = self.impeller.d1 / 2
        i = 1

        while self.outerStreamlineXCoords[i-1] < self.g1:
            self.outerStreamlineXCoords[i] = self.outerStreamlineXCoords[i-1] + deltam
            self.outerStreamlineYCoords[i] = self.impeller.d1 / 2
            i += 1

        transitionIndex = i

        extraArcAngle = degrees(atan2(deltam - (self.g1 - self.outerStreamlineXCoords[i-1]), self.outerStreamlineRadiusOfCurvature))
        startCoords = (self.g1 + ((deltam - (self.g1 - self.outerStreamlineXCoords[i-1])) * cos(radians(extraArcAngle))), (self.impeller.d1 / 2) + ((deltam - (self.g1 - self.outerStreamlineXCoords[i-1])) * sin(radians(extraArcAngle))))
        arcNumberOfPoints = floor((arcLength / deltam))
        arcCoords = ph.plotArc(startCoords, (self.g1, (self.impeller.d1 / 2) + self.outerStreamlineRadiusOfCurvature), arcAngle - extraArcAngle, arcNumberOfPoints)
        
        while i < transitionIndex + arcNumberOfPoints:
            self.outerStreamlineXCoords[i] = arcCoords[0][i - transitionIndex]
            self.outerStreamlineYCoords[i] = arcCoords[1][i - transitionIndex]
            i += 1

        while i < self.numberOfPoints:
            self.outerStreamlineXCoords[i] = self.outerStreamlineXCoords[i-1]
            self.outerStreamlineYCoords[i] = self.outerStreamlineYCoords[i-1] + deltam
            i += 1
        
        ph.simplePlot(self.outerStreamlineXCoords, self.outerStreamlineYCoords)

        if crossSectionVariationType == "linear":

            startArea = 2 * pi

fluid = fl.Fluid(density=787, viscosity=2.86e-3, vapourPressure=4100)
impeller = Impeller(suctionSpecificSpeedEU=650,
            suctionSidePressure=3e5,
            kgPerSec=1.32,
            headRise=295,
            fluid=fluid,
            shaftAllowableShearStress=8e7,
            numberOfBlades=6,
            approachFlowAngle=90,
            outletBladeAngle=22.5
)
meridionalSection = Meridional(impeller)
print(impeller.d2)
print(meridionalSection.ZE)