import conversions
import constants
import fluid as fl
from math import radians, degrees, sqrt, pi, log10, exp, tan, atan2
from scipy.interpolate import interp1d
from matplotlib import pyplot as plt

class Pump:

    def __init__(self, 
    rpm: float = None, 
    metreCubedPerSec: float = None, 
    headRise: float = None,
    fluid: fl.Fluid = None,
    approachFlowAngle: float = 90,
    axialThrustBalanceHoles: bool = False,
    shaftAllowableShearStress: float = None,
    shaftDiameterSafetyFactor: float = 1.1,
    headCoefficientCorrelation: str = "Karassik",
    numberOfBlades: int = 6,
    hubShaftDiameterRatio: float = 1.5,
    headCoefficientOverride: bool = False,
    isSuctionImpeller: bool = False,
    incidenceAngle: float = 2,
    bladeThickness: float = None,
    convergenceCriteria: float = 0.001
    ):

        self.rpm = rpm
        self.metreCubedPerSec = metreCubedPerSec
        self.headRise = headRise
        self.fluid = fluid
        self.approachFlowAngle = approachFlowAngle
        self.alpha1 = self.approachFlowAngle # Alternate name
        self.axialThrustBalanceHoles = axialThrustBalanceHoles
        self.shaftAllowableShearStress = shaftAllowableShearStress
        self.shaftDiameterSafetyFactor = shaftDiameterSafetyFactor
        self.headCoefficientCorrelation = headCoefficientCorrelation
        self.headCoefficientOverride = headCoefficientOverride
        self.hubShaftDiameterRatio = hubShaftDiameterRatio
        self.numberOfBlades = numberOfBlades
        self.convergenceCriteria = convergenceCriteria
        self.isSuctionImpeller = isSuctionImpeller
        self.incidenceAngle = incidenceAngle
        self.bladeThickness = bladeThickness
        
        # Unit conversions

        self.radPerSec = conversions.RPMToAngularSpeed(self.rpm)
        self.revPerSec = conversions.RPMToRevPerSec(self.rpm)
        self.metreCubedPerMin = self.metreCubedPerSec * 60
        self.metreCubedPerHour = self.metreCubedPerMin * 60
        self.kgPerSec = self.metreCubedPerSec * self.fluid.density
        self.gallonPerMin = conversions.metreCubedPerSecToGallonPerMin(self.metreCubedPerSec)
        self.headRiseFeet = conversions.metreToFeet(self.headRise)
        self.pressureRise = conversions.headToPressure(self.headRise, self.fluid.density)

        # Calculate specific speed

        self.specificSpeedDimensionless = (self.radPerSec * sqrt(self.metreCubedPerSec)) / ((constants.G * self.headRise) ** 0.75)
        self.specificSpeedEU = (self.rpm * sqrt(self.metreCubedPerSec)) / (self.headRise ** 0.75)
        self.specificSpeedUS = (self.rpm * sqrt(self.gallonPerMin)) / (self.headRiseFeet ** 0.75)

        # Estimate efficiencies
        # Table 3.9 Gulich

        if metreCubedPerSec < 0.005:

            print("Warning: hydraulic efficiency estimate not valid for Q < 0.005 [m^3 s^-1] so value may be unreliable")

        Qref = 1

        if self.metreCubedPerSec >= 1:
            a = 1
        else:
            a = 0.5

        m = 0.1 * a * ((Qref / metreCubedPerSec) ** 0.15) * ((45 / self.specificSpeedEU) ** 0.06)

        self.overallEfficiency = 1 - 0.095 * ((Qref / self.metreCubedPerSec) ** m) - 0.3 * ((0.35 - log10(self.specificSpeedEU / 23)) ** 2) * ((Qref / self.metreCubedPerSec) ** 0.05)

        m = 0.08 * a * ((Qref / metreCubedPerSec) ** 0.15) * ((45 / self.specificSpeedEU) ** 0.06)

        self.hydraulicEfficiency = 1 - 0.055 * ((Qref / self.metreCubedPerSec) ** m) - 0.2 * ((0.26 - log10(self.specificSpeedEU / 25)) ** 2) * ((Qref / self.metreCubedPerSec) ** 0.1)

        # Estimate volumetric losses
        # Table 3.5 Gulich

        if self.specificSpeedEU >= 27:
            a = 0.15
            m = 0.6
        else:
            a = 4.1
            m = 1.6

        if self.axialThrustBalanceHoles:
            z_H = 2
        else:
            z_H = 1

        self.impellerInletSealLeakageFlowRate = self.metreCubedPerSec * (a * z_H) / (self.specificSpeedEU ** m)

        self.volumetricEfficiency = self.metreCubedPerSec / (self.metreCubedPerSec + self.impellerInletSealLeakageFlowRate)

        # Calculate shaft power

        self.shaftPower = (self.kgPerSec * self.headRise * constants.G) / self.overallEfficiency

        # Calculate shaft diameter
        # Eq T7.1.2 Gulich

        self.minShaftDiameter = 3.65 * ((self.shaftPower / (self.rpm * self.shaftAllowableShearStress)) ** (1/3))
        self.shaftDiameter = self.minShaftDiameter * self.shaftDiameterSafetyFactor

        # Estimation of head coefficient

        if self.headCoefficientOverride:

            self.headCoefficient = headCoefficientOverride

        else:

            if self.headCoefficientCorrelation.lower() == "gulich":

                # Eq. 3.26 Gulich
                self.headCoefficient = 1.21 * exp(-0.77 * self.specificSpeedEU / 100) 

            elif self.headCoefficientCorrelation.lower() == "karassik":

                # Fig 12 Karassik
                if self.specificSpeedDimensionless >= 1:
                    self.headCoefficient = 0.4 / (self.specificSpeedDimensionless ** 0.25)
                else:
                    self.headCoefficient = 0.4 / (self.specificSpeedDimensionless ** 0.5)

        # Calculate impeller outer diameter

        self.impellerOuterDiameter = (2 / self.radPerSec) * sqrt((constants.G * self.headRise) / self.headCoefficient)
        self.u2 = (self.impellerOuterDiameter / 2) * self.radPerSec

        # Calculate impeller inlet diameter

        self.impellerHubDiameter = self.shaftDiameter * self.hubShaftDiameterRatio
        self.impellerHubDiameterDimensionless = self.impellerHubDiameter / self.impellerOuterDiameter

        self.impellerInletDiameter = self.impellerHubDiameter * 1.5
        impellerInletDiameterPrev = 10

        if self.isSuctionImpeller:

            self.fd1 = 1.20

        else:

            if self.specificSpeedEU <= 15:

                self.fd1 = 1.15

            elif self.specificSpeedEU >= 40:

                self.fd1 = 1.05

            else:

                self.fd1 = interp1d([40, 15], [1.15, 1.05])(self.specificSpeedEU)

        while abs(abs(self.impellerInletDiameter) - abs(impellerInletDiameterPrev)) > self.convergenceCriteria:

            impellerInletDiameterPrev = self.impellerInletDiameter

            self.u1 = (self.impellerInletDiameter / 2) * self.radPerSec
            self.c1m = (4 * self.metreCubedPerSec) / (pi * ((self.impellerInletDiameter ** 2) - (self.impellerHubDiameter ** 2)))
            self.swirlNumber = 1 - (self.c1m / (self.u1 * tan(radians(self.approachFlowAngle))))

            self.impellerInletDiameterDimensionless = self.fd1 * sqrt((self.impellerHubDiameterDimensionless ** 2) + (1.5e-3 * self.headCoefficient * ((self.specificSpeedEU ** 1.33) / (self.swirlNumber ** 0.67))))
            self.impellerInletDiameter = self.impellerInletDiameterDimensionless * self.impellerOuterDiameter

        # Calculate inlet velocity triangle without blockage

        self.c1u = self.c1m / tan(radians(self.approachFlowAngle))
        self.c1 = sqrt((self.c1m ** 2) + (self.c1u ** 2))
        self.w1 = sqrt((self.c1m ** 2) + ((self.u1 - self.c1u) ** 2))
        self.inletFlowCoefficient = self.c1m / self.u1
        self.beta1 = degrees(atan2(self.c1m, self.u1 - self.c1u))

        # Calculate required blade thickness if not overridden

        if self.bladeThickness == None:

            upperHead = 600
            lowerHead = 10
            self.bladeThicknessDimensionless = interp1d([upperHead, lowerHead], [0.022, 0.016])(self.headRise)
            self.bladeThickness = self.bladeThicknessDimensionless * self.impellerOuterDiameter


    def printResults(self, 
    inputs: bool = True,
    performance: bool = True,
    sizes: bool = True,
    velocityTriangles: bool = True,
    decimalPlaces: int = 4
    ):

        print("")

        def formatNumber(num):

            return str(round(num, decimalPlaces))

        printSeparator = "---------- \n"

        if inputs:
            
            inputsOutput = [
                printSeparator, "Inputs \n", printSeparator,
                "\nPump Speed: \n",
                formatNumber(self.rpm), " [RPM] \n",
                formatNumber(self.revPerSec), " [s^-1] \n",
                formatNumber(self.radPerSec), " [rad s^-1] \n",
                "\nVolume Flow Rate: \n",
                formatNumber(self.metreCubedPerSec), " [m^3 s^-1] \n",
                formatNumber(self.metreCubedPerHour), " [m^3 hr^-1] \n",
                formatNumber(self.gallonPerMin), " [gpm] \n",
                "\nMass Flow Rate: \n",
                formatNumber(self.kgPerSec), " [kg s^-1] \n",
                "\nHead Rise: \n",
                formatNumber(self.headRise), " [m] \n",
                formatNumber(self.headRiseFeet), " [ft] \n",
                "\nPressure Rise: \n",
                formatNumber(self.pressureRise/1e5), " [Bar] \n"
            ]

            inputsString = "".join(inputsOutput)

            print(inputsString)

        if performance:

            specificSpeedOutput = [
                printSeparator, "Specific Speed \n", printSeparator,
                "\nDimensionless Specific Speed: \n",
                formatNumber(self.specificSpeedDimensionless), " \n",
                "\nEU Specific Speed: \n",
                formatNumber(self.specificSpeedEU), " \n",
                "\nUS Specific Speed: \n",
                formatNumber(self.specificSpeedUS), " \n"
            ]

            specificSpeedString = "".join(specificSpeedOutput)

            print(specificSpeedString)

            efficienciesOutput = [
                printSeparator, "Efficiencies \n", printSeparator,
                "\nOverall Efficiency: \n",
                formatNumber(self.overallEfficiency), " \n",
                "\nHydraulic Efficiency: \n",
                formatNumber(self.hydraulicEfficiency), " \n",
                "\nVolumetric Efficiency: \n",
                formatNumber(self.volumetricEfficiency), " \n"
            ]

            efficienciesString = "".join(efficienciesOutput)

            print(efficienciesString)

            coefficientsOutput = [
                printSeparator, "Coefficients \n", printSeparator,
                "\nHead Coefficient: \n",
                formatNumber(self.headCoefficient), " \n"
            ]

            coefficientsString = "".join(coefficientsOutput)

            print(coefficientsString)

        if sizes:

            sizesOutput = [
                printSeparator, "Sizes \n", printSeparator,
                "\nShaft Diameter: \n",
                formatNumber(self.shaftDiameter*1e3), " [mm]\n",
                "\nHub Diameter: \n",
                formatNumber(self.impellerHubDiameter*1e3), " [mm]\n",
                "\nImpeller Inlet Diameter: \n",
                formatNumber(self.impellerInletDiameter*1e3), " [mm]\n",
                "\nImpeller Outer Diameter: \n",
                formatNumber(self.impellerOuterDiameter*1e3), " [mm]\n",
                "\nBlade Thickness: \n",
                formatNumber(self.bladeThickness*1e3), " [mm]\n"
            ]

            sizesString = "".join(sizesOutput)

            print(sizesString)

        if velocityTriangles:

            velocityTriangleOutput = [
                printSeparator, "Inlet Velocity Triangle Without Blockage\n", printSeparator,
                "\nTip Speed u1: \n",
                formatNumber(self.u1), " [m s^-1]\n",
                "\nAbsolute Speed c1: \n",
                formatNumber(self.c1), " [m s^-1]\n", "(Meridional: ", formatNumber(self.c1m), " [m s^-1])\n", "(Circumferential: ", formatNumber(self.c1u), " [m s^-1])\n",
                "\nRelative Speed w1: \n", 
                formatNumber(self.w1), " [m s^-1]\n",
                "\nApproach Flow Angle alpha1: \n", 
                formatNumber(self.approachFlowAngle), " [deg]\n",
                "\nRelative Flow Angle beta1: \n", 
                formatNumber(self.beta1), " [deg]\n",
                "\nFlow Coefficient: \n", 
                formatNumber(self.inletFlowCoefficient), " \n",
                printSeparator, "Oulet Velocity Triangle \n", printSeparator,
                "\nTip Speed u2: \n",
                formatNumber(self.u2), " [m s^-1]\n"
            ]

            velocityTriangleString = "".join(velocityTriangleOutput)

            print(velocityTriangleString)

        print("")

    def plotVelocityTriangle(self, area: str):

        pass


fluid = fl.Fluid(density=787, viscosity=2.86e-3, vapourPressure=3000)
pump = Pump(rpm=25000,
            metreCubedPerSec=0.00167,
            headRise=291.4,
            fluid=fluid,
            shaftAllowableShearStress=8e7
)

"""

fluid = fl.Fluid(density=1000, viscosity=2.86e-3, vapourPressure=3000)
pump = Pump(rpm=1450,
            metreCubedPerSec=0.0778,
            headRise=20,
            fluid=fluid,
            shaftAllowableShearStress=8e7
)

"""

pump.printResults()


