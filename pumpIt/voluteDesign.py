from cmath import pi
from math import atan2, degrees, sqrt, tan, radians, asin, sin, cos, acos
from impeller import Impeller
from plottingHelper import Bezier, findIntersection, findIntersectionOfCoords, polarToCartesian, cartesianToPolar, polarToCartesianLines, cartesianToPolarLines
import matplotlib.pyplot as plt
import numpy as np
import csv
import os

class CrossSection:

    def __init__(self, numberOfSections: int = 100):

        self.type = None
        self.numberOfSections = numberOfSections
        self.rCoords = []
        self.aCoords = []

        self.width = None
        self.height = None
        self.diameter = None
        self.rBase = None

    def generateCoords(self):

        pass

    def calculateArea(self) -> float:

        i = 0
        sum = 0

        while i < self.numberOfSections - 1:

            deltar = self.rCoords[i+1] - self.rCoords[i]
            b = 2 * abs(self.aCoords[i])

            sum += b * deltar
            i += 1

        return sum

    def calculateSummation(self) -> float:

        i = 0
        sum = 0

        while i < self.numberOfSections - 1:

            r = self.rCoords[i]
            deltar = self.rCoords[i+1] - r
            b = 2 * abs(self.aCoords[i])

            sum += (b / r) * deltar
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

    def __init__(self, angle: float = 45, numberOfSections: int = 1000):
        super().__init__(numberOfSections)
        self.angle = angle

    def generateCoords(self, rBase, width: float, height: float) -> float:

        self.aCoords = []
        self.rCoords = []

        self.width = width
        self.height = height
        self.rBase = rBase
        base = [-self.width / 2, self.rBase]
        top = [0, self.rBase + self.height]

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

    def __init__(self, numberOfSections: int = 1000):
        super().__init__(numberOfSections)

    def generateCoords(self, rBase: float, width: float, height: float) -> float:
        
        self.aCoords = []
        self.rCoords = []

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

class CircularCrossSection(CrossSection):

    def __init__(self, numberOfSections: int = 100):
        super().__init__(numberOfSections)

    def generateCoords(self, rBase: float, width: float, height: float) -> float:

        self.aCoords = []
        self.rCoords = []

        self.rBase = rBase
        self.width = width
        self.height = height

        self.diameter = 2 * (((self.width / 2) ** 2) + (self.height ** 2)) / (2 * self.height)

        centre = (0, self.rBase + sqrt((self.diameter / 2) ** 2 - (self.width/2) ** 2))
        if self.height > self.diameter / 2:
            startAngle = degrees(asin(self.width / self.diameter))
        else:
            startAngle = degrees(acos(self.width / self.diameter)) + 90

        angles = np.linspace(startAngle, 180, self.numberOfSections)

        self.aCoords = [((self.diameter / 2) * sin(radians(theta))) for theta in angles]
        self.rCoords = [((-self.diameter / 2) * cos(radians(theta)) + centre[1]) for theta in angles]
        aCoordsReflected = [-((self.diameter / 2) * sin(radians(theta))) for theta in angles]
        aCoordsReflected.reverse()
        rCoordsReflected = self.rCoords.copy()
        rCoordsReflected.reverse()

        self.aCoords += aCoordsReflected
        self.rCoords += rCoordsReflected

        return self.calculateSummation()

class Volute:

    def __init__(
        self, 
        impeller: Impeller,
        voluteCrossSection: TrapezoidalCrossSection | RectangularCrossSection,
        diffuser = None, 
        dzd2RatioOverride: float = None,
        b3b2RatioOverride: float = None,
        rAIncrementFactor: float = 1000,
        applyCutwaterCorrection: bool = True,
        dischargeAreaRatio: float = None,
        dischargeLength: float = None,
        dischargeExitDiameter: float = None,
        numberOfPointsCutwater: int = 100,
        dischargeInnerAngle: float = 2.5

        ) -> None:

        self.impeller = impeller
        self.diffuser = diffuser
        self.dzd2RatioOverride = dzd2RatioOverride
        self.b3b2RatioOverride = b3b2RatioOverride
        self.voluteCrossSection = voluteCrossSection
        self.rAIncrementFactor = rAIncrementFactor
        self.dischargeAreaRatio = dischargeAreaRatio
        self.dischargeLength = dischargeLength
        self.dischargeExitDiameter = dischargeExitDiameter
        self.numberOfPointsCutwater = numberOfPointsCutwater
        self.dischargeInnerAngle = dischargeInnerAngle

        if self.dischargeExitDiameter != None:
            self.dischargeArea = pi * (self.dischargeExitDiameter/2) ** 2

        # TODO: Support partial volutes
        self.ZLe = 1
        self.wrapAngle = 360

        # Determine the flow rate through volute
        self.casingDesignFlowRate = self.impeller.metreCubedPerSec
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
        self.rz = self.dz / 2

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

        self.cutwaterLEThickness = 0.02 * self.impeller.d2
        self.e3 = self.cutwaterLEThickness # Alias

        self.cutwaterIncidence = 3
        self.i3 = self.cutwaterIncidence # ALias
        self.inletFlowAngle = degrees(atan2(self.impeller.c2mDash, self.impeller.c2u))
        self.alpha3 = self.inletFlowAngle # Alias
        self.cutwaterCamberAngle = self.alpha3 + self.cutwaterIncidence
        self.alpha3B = self.cutwaterCamberAngle # Alias

        # Calculate volute
        
        self.rzDash = self.rz + (self.e3 / 2)
        self.rAs = []
        self.epsilons = []
        self.areas = []
        rAIncrement = self.rz / self.rAIncrementFactor
        epsilon = 0
        area = 0
        rA = self.rzDash
        while epsilon <= self.wrapAngle:
            
            self.rAs.append(rA)
            self.epsilons.append(epsilon)
            self.areas.append(area)

            rA += rAIncrement

            summation = voluteCrossSection.generateCoords(self.rz, self.b3, rA - self.rzDash)
            area = voluteCrossSection.calculateArea()
            epsilon = ((360 * self.c2u * (self.impeller.d2 / 2)) / self.QLe) * summation
            
        self.a3 = self.rAs[-1] - self.rzDash
        self.throatWidth = self.a3 # Alias
        self.A3q = self.areas[-1]
        self.throatArea = self.A3q # Alias

        # Cutwater correction
        # Applies linear correction from 0 to 360 degrees
        self.deltaa3 = self.a3 * 0.2 * ((self.e3 / self.a3) ** 2)
        correction = np.linspace(0, self.deltaa3, len(self.rAs))
        self.rAsCorrected = [None for x in range(len(self.rAs))]
        for i in range(len(self.rAs)):
            self.rAsCorrected[i] = self.rAs[i] + correction[i]
        
        if applyCutwaterCorrection:
            self.rAs = self.rAsCorrected

        # Discharge nozzle

        if type(voluteCrossSection) == RectangularCrossSection:
            self.hd = self.voluteCrossSection.height
        else:
            self.Rdeq = sqrt(self.throatArea / pi)

        if self.dischargeExitDiameter != None:
            self.dischargeAreaRatio = self.dischargeArea / self.throatArea

        # Returns either area ratio or length depending on which was specified, from the specified file
        def readCpData(file, length=None, areaRatio=None):
            
            if type(self.voluteCrossSection) == RectangularCrossSection:
                dimension = self.hd
            else:
                dimension = self.Rdeq

            reader = csv.reader(file)
            cpData = [[], []]
            for row in reader:
                cpData[0].append(float(row[0]))
                cpData[1].append(float(row[1]))
            if length != None:
                x = length / dimension
                line2 = [[x, x], [0, 10]]
                point = findIntersectionOfCoords(cpData, line2)
                return point[0] + 1
            elif areaRatio != None:
                y = areaRatio - 1
                line2 = [[1, 40], [y, y]]
                point = findIntersectionOfCoords(cpData, line2)
                return point[1] * dimension

        # Uses above function with the correct cp data file
        scriptDir = os.path.dirname(__file__)
        if self.dischargeAreaRatio != None:
            if type(voluteCrossSection) == RectangularCrossSection:
                relPath = "Data\planarcp__.csv"
                absFilePath = os.path.join(scriptDir, relPath)
                with open(absFilePath, 'r') as file:
                    self.dischargeLength = readCpData(file, areaRatio=self.dischargeAreaRatio)
            else:
                relPath = "Data\conicalcp__.csv"
                absFilePath = os.path.join(scriptDir, relPath)
                with open(absFilePath, 'r') as file:
                    self.dischargeLength = readCpData(file, areaRatio=self.dischargeAreaRatio)
        elif self.dischargeLength != None:
            if type(voluteCrossSection) == RectangularCrossSection:
                relPath = "Data\planarcp_.csv"
                absFilePath = os.path.join(scriptDir, relPath)
                with open(absFilePath, 'r') as file:
                    self.dischargeAreaRatio = readCpData(file, length=self.dischargeLength)
            else:
                relPath = "Data\conicalcp_.csv"
                absFilePath = os.path.join(scriptDir, relPath)
                with open(absFilePath, 'r') as file:
                    self.dischargeAreaRatio = readCpData(file, length=self.dischargeLength)

        self.dischargeArea = self.throatArea * self.dischargeAreaRatio

        if self.dischargeExitDiameter == None:

            self.dischargeExitDiameter = sqrt(self.dischargeArea / pi) * 2

        # Convert to cartesian coords for easier construction of discarge nozzle

        self.xCoords = []
        self.yCoords = []
        self.rzDashXCoords = []
        self.rzDashYCoords = []
        for i in range(len(self.rAsCorrected)):
            coords = polarToCartesian(self.rAsCorrected[i], self.epsilons[i])
            self.xCoords.append(coords[0])
            self.yCoords.append(coords[1])
            coords = polarToCartesian(self.rzDash, self.epsilons[i])
            self.rzDashXCoords.append(coords[0])
            self.rzDashYCoords.append(coords[1])

        # Calculate cutwater coords
        self.cutwaterXCoords = []
        self.cutwaterYCoords = []
        degreeIncrement = 180 / self.numberOfPointsCutwater
        degree = 180
        for i in range(self.numberOfPointsCutwater):
            cutwaterCentre = [self.xCoords[0] + (self.e3 / 2), self.yCoords[0]]
            cartesianCoords = polarToCartesian(self.e3 / 2, degree)
            self.cutwaterXCoords.append(cutwaterCentre[0] + cartesianCoords[0])
            self.cutwaterYCoords.append(cutwaterCentre[1] + cartesianCoords[1])
            degree += degreeIncrement
        
        # Convert cutwater coords to polar for polar plotting
        self.cutwaterEpsilons = []
        self.cutwaterRs = []
        for i in range(len(self.cutwaterXCoords)):
            polarCoords = cartesianToPolar(self.cutwaterXCoords[i], self.cutwaterYCoords[i])
            self.cutwaterRs.append(polarCoords[0])
            self.cutwaterEpsilons.append(polarCoords[1])

        # Calcaulte coordinates of discharge outlet
        xModifier = self.dischargeLength * tan(radians(self.dischargeInnerAngle))
        self.dischargeOutletXCoords = [self.cutwaterXCoords[-1] + self.dischargeExitDiameter - xModifier, self.cutwaterXCoords[-1] - xModifier]
        self.dischargeOutletYCoords = [self.cutwaterYCoords[-1] + self.dischargeLength, self.cutwaterYCoords[-1] + self.dischargeLength]

        # Join all coords into another coordinate set for easy plotting
        self.cutwaterXCoords.reverse()
        self.cutwaterYCoords.reverse()

        self.totalXCoords = self.xCoords + self.dischargeOutletXCoords + self.cutwaterXCoords
        self.totalYCoords = self.yCoords + self.dischargeOutletYCoords + self.cutwaterYCoords
        self.totalRCoords, self.totalEpsilons = cartesianToPolarLines(self.totalXCoords, self.totalYCoords)
        