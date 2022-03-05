from math import atan2, degrees
from impeller import Impeller
from plottingHelper import Bezier, findIntersection
import matplotlib.pyplot as plt
import numpy as np

class CrossSection:

    def __init__(self, numberOfSections: int = 100):

        self.type = None
        self.numberOfSections = numberOfSections
        self.rCoords = []
        self.aCoords = []

    def generateCoords(self):

        pass

    def calculateSummation(self) -> float:

        i = 0
        sum = 0

        while i < self.numberOfSections - 1:

            r = self.rCoords[i]
            deltar = self.rCoords[i+1] - r
            b = 2 * abs(self.aCoords[i])

            sum += (r / b) * deltar
            i += 1

        return sum

    def plotShape(self):

        ax = plt.axes()
        ax.plot(self.aCoords, self.rCoords)
        ax.set_title("Volute Cross Section")
        ax.set_xlabel("Axial Distance [m]")
        ax.set_ylabel("Radial Distance [m]")
        ax.axis("equal")

        plt.show()

class TrapezoidalCrossSection(CrossSection):

    def __init__(self, angle, numberOfSections: int = 100):
        super().__init__(numberOfSections)
        self.angle = angle

    def generateCoords(self, rBase, width: float, height: float) -> float:

        self.width = width
        self.height = height
        self.rBase = rBase
        base = [-self.width / 2, self.rBase]
        top = [0, self.height]

        control = findIntersection(base, 90 + self.angle, top, 180)

        bzr = Bezier([base, control, top])
        tIncrement = 1 / self.numberOfSections
        t = 0
        for i in range(self.numberOfSections):
            coord = bzr.getPoint(t)
            self.aCoords.append(coord[0])
            self.rCoords.append(coord[1])
            t += tIncrement
        aCoordsReflected = [-x for x in self.aCoords]
        aCoordsReflected.reverse()
        rCoordsReflected = self.rCoords.copy()
        rCoordsReflected.reverse()
        self.aCoords += aCoordsReflected
        self.rCoords += rCoordsReflected

        return self.calculateSummation()

class RectangularCrossSection(CrossSection):

    def __init__(self, numberOfSections: int = 100):
        super().__init__(numberOfSections)

    def generateCoords(self, rBase: float, width: float, height: float) -> float:
        
        self.width = width
        self.height = height
        self.rBase = rBase

        self.aCoords = [-self.width/2 for x in range(self.numberOfSections)]
        self.rCoords = list(np.linspace(rBase, self.height + rBase, self.numberOfSections))
        aCoordsReflected = [self.width/2 for x in range(self.numberOfSections)]
        aCoordsReflected.reverse()
        rCoordsReflected = self.rCoords.copy()
        rCoordsReflected.reverse()

        self.aCoords += aCoordsReflected
        self.rCoords += rCoordsReflected

        return self.calculateSummation()
   
test = RectangularCrossSection(numberOfSections=1000)
print(test.generateCoords(1.1, 1.4, 1.2))
test.plotShape()

class Volute:

    def __init__(
        self, 
        impeller: Impeller, 
        partialVolutes: int = 1, 
        wrapAngleIncrement: float = 12.5, 
        diffuser = None, 
        dzd2RatioOverride: float = None,
        b3b2RatioOverride: float = None,
        voluteCrossSection: TrapezoidalCrossSection | RectangularCrossSection = None
        ) -> None:

        self.partialVolutes = partialVolutes
        self.ZLe = self.partialVolutes # Alias
        self.impeller = impeller
        self.wrapAngleIncrement = wrapAngleIncrement
        self.diffuser = diffuser
        self.dzd2RatioOverride = dzd2RatioOverride
        self.b3b2RatioOverride = b3b2RatioOverride

        # Determine wrap angles
        self.wrapAngles = [None for x in range(self.partialVolutes)]
        self.wrapAngles[0] = 360 / self.partialVolutes
        if self.partialVolutes > 1:
            for i in range(1, len(self.wrapAngles)):
                if self.impeller.numberOfBlades % 2 == 0:
                    self.wrapAngleIncrement = 0
                self.wrapAngles[i] = self.wrapAngles[i-1] - self.wrapAngleIncrement

        # Determine the flow rate through volute
        self.casingDesignFlowRate = self.impeller.impellerVolumeFlowRate
        self.QLe = self.casingDesignFlowRate # Alias

        # If no diffuser, inlet velocity = circumferntial component of absolute velocity of impeller outlet
        if self.diffuser == None:
            self.inletVelocity = self.impeller.c2u
            self.c2u = self.inletVelocity # ALias
        
        # Determine cutwater diameter
        # Currently unused
        if self.diffuser != None:
            if self.impeller.nq <= 40:
                self.d3d2Ratio = 1.015 + 0.08 * (((self.impeller.fluid.density * self.impeller.headRise) / (1000 * 1000)) - 0.1) ** 0.8
            else:
                self.d3d2Ratio = 1.04 + 0.001 * (self.impeller.nq - 40)
        
        # Determine cutwater diameter
        if self.dzd2RatioOverride == None:
            self.dzd2Ratio = 1.03 + 0.1 * (self.impeller.nq / 40) + 0.07 * ((self.impeller.fluid.density * self.impeller.headRise) / (1000 * 1000))
        else:
            self.dzd2Ratio = self.dzd2RatioOverride

        self.cutwaterDiameter = self.dzd2Ratio * self.impeller.d2
        self.dz = self.cutwaterDiameter

        # Determine inlet width
        # TODO: Research better ways of calculating this
        if self.b3b2RatioOverride == None:
            if self.impeller.nq >= 40:
                self.b3b2Ratio = 4 - (4 - 2) * (self.impeller.nq / 40)
            else:
                self.b3b2Ratio = 1.2 - (1.2 - 1.05) * (self.impeller.nq / 160)
        else:
            self.b3b2Ratio = self.b3b2RatioOverride

        self.b3 = self.b3b2Ratio * self.impeller.b2
        self.inletWidth = self.b3 # Alias

        # Calculate cutwater LE

        self.cutwaterLEThickness = 0.2 * self.impeller.d2
        self.e3 = self.cutwaterLEThickness # Alias

        self.cutwaterIncidence = 3
        self.i3 = self.cutwaterIncidence # ALias
        self.inletFlowAngle = degrees(atan2(self.impeller.c2mDash, self.impeller.c2u))
        self.alpha3 = self.inletFlowAngle # Alias
        self.cutwaterCamberAngle = self.alpha3 + self.cutwaterIncidence
        self.alpha3B = self.cutwaterCamberAngle # Alias

        # Calculate volute throat area
         
        
        

        
