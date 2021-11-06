import pumpIt.conversions as convert
from pumpIt.constants import G
from pumpIt.fluid import Fluid
import sys
from math import pi, sqrt

class Pump:

    def __init__(self, massFlowRate: float = None, volumeFlowRate: float = None, inletPressure: float = None, outletPressure: float = None, inletHead: float = None, outletHead: float = None, fluid: Fluid = None):
        
        if massFlowRate != None and volumeFlowRate == None:

            self.massFlowRate = massFlowRate
            self.volumeFlowRate = convert.massToVolumeFlowRate(massFlowRate, fluid.density)

        elif massFlowRate == None and volumeFlowRate != None:

            self.massFlowRate = convert.volumeToMassFlowRate(volumeFlowRate, fluid.density)
            self.volumeFlowRate = volumeFlowRate

        else:

            sys.exit("Error [Overspecified]: Please only specify one of mass flow rate or volume flow rate")

        # Following allows specification of either pressure or head at the inlet/outlet, with error catching so that pressure and head can't both be specified at the same time
        
        if inletPressure != None:

            if inletHead != None:

                sys.exit("Error [Overspecified]: Please only specify one of inlet pressure or head")

            self.inletPressure = inletPressure
            self.inletHead = convert.pressureToHead(self.inletPressure, fluid.density)

        elif inletHead != None:

            if inletPressure != None:

                sys.exit("Error [Overspecified]: Please only specify one of inlet pressure or head")

            self.inletHead = inletHead
            self.inletPressure = convert.headToPressure(inletHead, fluid.density)

        if outletPressure != None:

            if outletHead != None:

                sys.exit("Error [Overspecified]: Please only specify one of outlet pressure or head")

            self.outletPressure = outletPressure
            self.outletHead = convert.pressureToHead(self.outletPressure, fluid.density)

        elif outletHead != None:

            if outletPressure != None:

                sys.exit("Error [Overspecified]: Please only specify one of outlet pressure or head")

            self.outletHead = outletHead
            self.outletPressure = convert.headToPressure(outletHead, fluid.density)
        
        self.pressureRise = self.outletPressure - self.inletPressure
        self.headRise = self.outletHead - self.inletHead

        self.fluid = fluid
        self.npsh = self.inletHead - self.fluid.vapourHead

        # Variables to be calculated in methods below
        self.angularSpeed = None
        self.rpm = None
        self.suctionSpecificSpeed = None
        self.specificSpeed = None
        self.headCoefficient = None
        self.blockageEffect = None
        self.impellerOuterRadius = None


    def setSpeed(self, radPerSec: float = None, rpm: float = None, revPerSec: float = None):

        if radPerSec != None and rpm == None:

            self.angularSpeed = radPerSec
            self.rpm = convert.angularSpeedToRPM(radPerSec)

        elif rpm != None and radPerSec == None:

            self.angularSpeed = convert.rpmToAngularSpeed(rpm)
            self.rpm = rpm

        self.suctionSpecificSpeed = self.angularSpeed / (((G * self.npsh) / (self.volumeFlowRate ** (2/3))) ** (3/4))

    def setSuctionSpecificSpeed(self, suctionSpecificSpeed: float = None):

        self.suctionSpecificSpeed = suctionSpecificSpeed

        self.angularSpeed = (((G * self.npsh) / (self.volumeFlowRate ** (2/3))) ** (3/4)) * self.suctionSpecificSpeed
        self.rpm = convert.angularSpeedToRPM(self.angularSpeed)

    def setBlockageEffect(self, blockageEffect: float = 1):

        self.blockageEffect = blockageEffect

    def calculate(self):

        # Calculate specific speed
        self.specificSpeed = self.angularSpeed * sqrt(self.volumeFlowRate) / ((G * self.headRise) ** (3/4))

        # Calculate head coefficient (based on equations found in Pump Handbook 2.28)
        if self.specificSpeed < 1:

            self.headCoefficient = 0.4 / (self.specificSpeed ** (1/4))

        else:

            self.headCoefficient = 0.4 / sqrt(self.specificSpeed)

        self.impellerOuterRadius = (1 / self.angularSpeed) * sqrt((G * self.headRise) / self.headCoefficient)
        self.impellerOuterDiameter = self.impellerOuterRadius * 2
        self.exitFlowCoefficient = 0.1715 * sqrt(self.specificSpeed)
        self.outerWidth = self.volumeFlowRate / 2 * pi * self.angularSpeed * self.impellerOuterRadius ** 2 * self.exitFlowCoefficient

        