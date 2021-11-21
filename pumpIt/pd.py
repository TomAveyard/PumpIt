import conversions
import constants
import fluid as fl
from math import sqrt, pi, log10, exp

class Pump:

    def __init__(self, 
    rpm: float = None, 
    metreCubedPerSec: float = None, 
    headRise: float = None,
    fluid: fl.Fluid = None,
    inletFlowAngle: float = 90,
    axialThrustBalanceHoles: bool = False,
    allowableShaftShearStress: float = None,
    shaftDiameterSafetyFactor: float = 1.1,
    headCoefficientCorrelation: str = "Karassik",
    headCoefficientOverride: bool = False
    ):

        self.rpm = rpm
        self.metreCubedPerSec = metreCubedPerSec
        self.headRise = headRise
        self.fluid = fluid
        self.inletFlowAngle = inletFlowAngle
        self.axialThrustBalanceHoles = axialThrustBalanceHoles
        self.allowableShaftShearStress = allowableShaftShearStress
        self.shaftDiameterSafetyFactor = shaftDiameterSafetyFactor
        self.headCoefficientOverride = headCoefficientOverride
        self.headCoefficientCorrelation = headCoefficientCorrelation

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

        # Estimate pump hydraulic efficiency

        if metreCubedPerSec < 0.005:

            print("Warning: hydraulic efficiency estimate not valid for Q < 0.005 [m^3 s^-1] so value may be unreliable")

        Qref = 1

        if self.metreCubedPerSec >= 1:
            a = 1
        else:
            a = 0.5

        m = 0.08 * a * ((Qref / metreCubedPerSec) ** 0.15) * ((45 / self.specificSpeedEU) ** 0.06)

        self.hydraulicEfficiency = 1 - 0.055 * ((Qref / self.metreCubedPerSec) ** m) - 0.2 * ((0.26 - log10(self.specificSpeedEU / 25)) ** 2) * ((Qref / self.metreCubedPerSec) ** 0.1)

        # Estimate volumetric losses

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

        # Calculate shaft diameter


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


    def printResults(self, 
    inputs: bool = True,
    performance: bool = True,
    sizes: bool = True,
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
                "\nHydraulic Efficiency: \n",
                formatNumber(self.hydraulicEfficiency), " \n"
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
                "\nImpeller Outer Diameter: \n",
                formatNumber(self.impellerOuterDiameter*1e3), " mm\n"
            ]

            sizesString = "".join(sizesOutput)

            print(sizesString)

        print("")

"""
fluid = fl.Fluid(density=787, viscosity=2.86e-3, vapourPressure=3000)
pump = Pump(rpm=30000,
            metreCubedPerSec=0.00167,
            headRise=291.4,
            fluid=fluid,
)
"""

fluid = fl.Fluid(density=1000, viscosity=2.86e-3, vapourPressure=3000)
pump = Pump(rpm=1450,
            metreCubedPerSec=0.0778,
            headRise=20,
            fluid=fluid,
)

pump.printResults()


