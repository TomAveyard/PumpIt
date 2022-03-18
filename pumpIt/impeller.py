import conversions
import constants
import fluid as fl
from math import atan, radians, degrees, sqrt, pi, log10, exp, tan, atan2, sin, cos
from scipy.interpolate import interp1d
from matplotlib import pyplot as plt

class Impeller:

    def __init__(self, 
    rpm: float = None, 
    suctionSpecificSpeedEU: float = None,
    suctionSidePressure: float = None,
    NPSHA: float = None,
    specificSpeedEU: float = None,
    metreCubedPerSec: float = None,
    kgPerSec: float = None,
    headRise: float = None,
    fluid: fl.Fluid = None,
    approachFlowAngle: float = 90,
    axialThrustBalanceHoles: bool = False,
    shaftAllowableShearStress: float = 42e8,
    shaftDiameterSafetyFactor: float = 1.1,
    shaftDiameterOverride: float = None,
    headCoefficientCorrelation: str = "gulich",
    numberOfBlades: int = 6,
    hubShaftDiameterRatio: float = 1.5,
    inletBladeInnerDiameterRatio: float = 1.15,
    headCoefficientOverride: bool = False,
    isSuctionImpeller: bool = False,
    incidenceAngle: float = 2,
    bladeThickness: float = None,
    leadingEdgeThicknessLengthRatio: float = 0.2,
    inletBladeInclination: float = 90,
    outletBladeAngle: float = 22.5,
    outletBladeInclination: float = 90,
    outletWidthOverride: float = None,
    convergenceCriteria: float = 0.001
    ):

        self.rpm = rpm
        self.suctionSpecificSpeedEU = suctionSpecificSpeedEU
        self.nss = self.suctionSpecificSpeedEU # Alias
        self.suctionSidePressure = suctionSidePressure
        self.NPSHA = NPSHA
        self.specificSpeedEU = specificSpeedEU
        self.metreCubedPerSec = metreCubedPerSec
        self.kgPerSec = kgPerSec
        self.headRise = headRise
        self.deltaH = self.headRise # Alias
        self.fluid = fluid
        self.approachFlowAngle = approachFlowAngle
        self.alpha1 = self.approachFlowAngle # Alias
        self.axialThrustBalanceHoles = axialThrustBalanceHoles
        self.shaftAllowableShearStress = shaftAllowableShearStress
        self.shaftDiameterSafetyFactor = shaftDiameterSafetyFactor
        self.shaftDiameterOverride = shaftDiameterOverride
        self.headCoefficientCorrelation = headCoefficientCorrelation
        self.headCoefficientOverride = headCoefficientOverride
        self.hubShaftDiameterRatio = hubShaftDiameterRatio
        self.inletBladeInnerDiameterRatio = inletBladeInnerDiameterRatio
        self.d1iRatio = self.inletBladeInnerDiameterRatio #Alias
        self.numberOfBlades = numberOfBlades
        self.ZLa = self.numberOfBlades #Alias
        self.convergenceCriteria = convergenceCriteria
        self.isSuctionImpeller = isSuctionImpeller
        self.incidenceAngle = incidenceAngle
        self.bladeThickness = bladeThickness
        self.leadingEdgeThicknessLengthRatio = leadingEdgeThicknessLengthRatio
        self.inletBladeInclination = inletBladeInclination
        self.outletBladeAngle = outletBladeAngle
        self.outletWidthOverride = outletWidthOverride
        self.outletBladeInclination = outletBladeInclination

        # Convert kg/s to m3/s or vice versa

        if self.kgPerSec != None and self.metreCubedPerSec == None:

            self.metreCubedPerSec = self.kgPerSec / self.fluid.density

        elif self.metreCubedPerSec != None and self.kgPerSec == None:

            self.kgPerSec = self.metreCubedPerSec * self.fluid.density

        else:

            exit("Error: both m3PerSec and kgPerSec have been specified - please specify only one")

        # Convert suction side pressure to NPSHA or vice versa if either is specified

        if self.suctionSidePressure != None:

            self.NPSHA = conversions.pressureToHead(self.suctionSidePressure - self.fluid.vapourPressure, self.fluid.density)

        elif self.NPSHA != None:

            self.suctionSidePressure = conversions.headToPressure(self.fluid.vapourHead + self.NPSHA, self.fluid.density)

        # Calculate suction performance based on specified parameters

        if self.suctionSpecificSpeedEU != None and self.suctionSidePressure != None:

            self.NPSH3 = 0.878 * (self.NPSHA ** (0.14))
            self.rpm = self.suctionSpecificSpeedEU * (self.NPSH3 ** 0.75) / sqrt(self.metreCubedPerSec)

        elif self.suctionSidePressure != None and self.rpm != None:

            self.NPSH3 = 0.878 * (self.NPSHA ** (0.14))
            self.suctionSpecificSpeedEU = self.rpm * sqrt(self.metreCubedPerSec) / ((self.NPSH3) ** 0.75)

        elif self.rpm != None and self.suctionSpecificSpeedEU != None:

            self.NPSH3 = ((self.rpm / self.suctionSpecificSpeedEU) * sqrt(self.metreCubedPerSec)) ** (1.333)

        elif self.specificSpeedEU != None and self.suctionSpecificSpeedEU != None:

            self.rpm = self.specificSpeedEU * (self.headRise ** 0.75) / sqrt(self.metreCubedPerSec)
            self.NPSH3 = ((self.rpm / self.specificSpeedEU) * sqrt(self.metreCubedPerSec)) ** (1.333)
            self.NPSHA = 1.16 * (self.NPSH3 ** 1.14)

        else:

            exit("Error: invalid pump specification")

        # Unit conversions

        self.radPerSec = conversions.RPMToAngularSpeed(self.rpm)
        self.revPerSec = conversions.RPMToRevPerSec(self.rpm)
        self.metreCubedPerMin = self.metreCubedPerSec * 60
        self.metreCubedPerHour = self.metreCubedPerMin * 60
        self.gallonPerMin = conversions.metreCubedPerSecToGallonPerMin(self.metreCubedPerSec)
        self.headRiseFeet = conversions.metreToFeet(self.headRise)
        self.pressureRise = conversions.headToPressure(self.headRise, self.fluid.density)

        # Calculate specific speed

        self.specificSpeedDimensionless = (self.radPerSec * sqrt(self.metreCubedPerSec)) / ((constants.G * self.headRise) ** 0.75)
        self.specificSpeedEU = (self.rpm * sqrt(self.metreCubedPerSec)) / (self.headRise ** 0.75)
        self.nq = self.specificSpeedEU # ALias
        self.specificSpeedUS = (self.rpm * sqrt(self.gallonPerMin)) / (self.headRiseFeet ** 0.75)

        # Estimate efficiencies
        # Table 3.9 Gulich

        if self.metreCubedPerSec < 0.005:

            print("Warning: hydraulic efficiency estimate not valid for Q < 0.005 [m^3 s^-1] so value may be unreliable")

        Qref = 1

        if self.metreCubedPerSec >= 1:
            a = 1
        else:
            a = 0.5

        m = 0.1 * a * ((Qref / self.metreCubedPerSec) ** 0.15) * ((45 / self.specificSpeedEU) ** 0.06)

        self.overallEfficiency = 1 - 0.095 * ((Qref / self.metreCubedPerSec) ** m) - 0.3 * ((0.35 - log10(self.specificSpeedEU / 23)) ** 2) * ((Qref / self.metreCubedPerSec) ** 0.05)

        m = 0.08 * a * ((Qref / self.metreCubedPerSec) ** 0.15) * ((45 / self.specificSpeedEU) ** 0.06)

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

        self.impellerVolumeFlowRate = self.metreCubedPerSec / self.volumetricEfficiency

        # Calculate shaft power

        self.shaftPower = (self.kgPerSec * self.headRise * constants.G) / self.overallEfficiency

        # Calculate shaft diameter
        # Eq T7.1.2 Gulich

        self.minShaftDiameter = 3.65 * ((self.shaftPower / (self.rpm * self.shaftAllowableShearStress)) ** (1/3))
        if self.shaftDiameterOverride == None:
            self.shaftDiameter = self.minShaftDiameter * self.shaftDiameterSafetyFactor
        else:
            self.shaftDiameter = self.shaftDiameterOverride
        self.shaftTorque = self.shaftPower / self.radPerSec

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

        self.impellerOutletDiameter = (60 / (pi * self.rpm)) * sqrt((2 * constants.G * self.headRise) / self.headCoefficient)
        self.u2 = (self.impellerOutletDiameter / 2) * self.radPerSec

        # Calculate impeller inlet diameter

        self.impellerHubDiameter = self.shaftDiameter * self.hubShaftDiameterRatio
        self.impellerHubDiameterDimensionless = self.impellerHubDiameter / self.impellerOutletDiameter

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

            self.u1 = (self.impellerInletDiameter * pi * self.rpm) / 60
            self.c1m = (4 * self.impellerVolumeFlowRate) / (pi * ((self.impellerInletDiameter ** 2) - (self.impellerHubDiameter ** 2)))
            self.swirlNumber = 1 - (self.c1m / (self.u1 * tan(radians(self.approachFlowAngle))))

            self.impellerInletDiameterDimensionless = self.fd1 * sqrt((self.impellerHubDiameterDimensionless ** 2) + (1.5e-3 * self.headCoefficient * ((self.specificSpeedEU ** 1.33) / (self.swirlNumber ** 0.67))))
            self.impellerInletDiameter = self.impellerInletDiameterDimensionless * self.impellerOutletDiameter

        self.impellerInletWidth = 0.5 * (self.impellerInletDiameter - self.impellerHubDiameter)
        self.b1 = self.impellerInletWidth

        # Calculate inlet blade inner diameter

        self.inletBladeInnerDiameter = self.impellerHubDiameter * self.inletBladeInnerDiameterRatio
        self.d1i = self.inletBladeInnerDiameter # Alias

        # Calculate inlet blade outer diameter

        self.inletBladeOuterDiameter = self.impellerInletDiameter

        # Calculate geometric average of diameters at impeller inlet

        self.inletDiametersGeometricAverage = sqrt(0.5 * ((self.inletBladeInnerDiameter ** 2) + (self.inletBladeOuterDiameter ** 2)))

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
            self.bladeThickness = self.bladeThicknessDimensionless * self.impellerOutletDiameter

        # Calculate leading edge tip radius
        profileLength = self.bladeThickness / self.leadingEdgeThicknessLengthRatio
        ellipseLocus = sqrt((profileLength ** 2) - ((self.bladeThickness ** 2) / 4))
        self.leadingEdgeTipRadius = profileLength - ellipseLocus

        # Calculate inlet blade angle and inlet blockage

        self.beta1B = self.beta1 + self.incidenceAngle
        beta1Bprev = self.beta1B + 10

        while abs(abs(self.beta1B) - abs(beta1Bprev)) > convergenceCriteria:

            beta1Bprev = self.beta1B

            self.inletBladeBlockage = (1 - ((self.numberOfBlades * self.leadingEdgeTipRadius) / (pi * self.impellerInletDiameter * sin(radians(self.beta1B)) * sin(radians(self.inletBladeInclination))))) ** (-1)
            self.beta1dash = degrees(atan2(self.c1m * self.inletBladeBlockage, self.u1 - self.c1u))
            self.beta1B = self.beta1dash + self.incidenceAngle

        self.inletBladeAngle = self.beta1B

        # Calculate inlet velocity triangle with blockage

        self.c1mdash = (self.u1 - self.c1u) * tan(radians(self.beta1dash))
        self.w1dash = (self.u1 - self.c1u) / cos(radians(self.beta1dash))
        self.c1dash = sqrt((self.c1mdash ** 2) + (self.c1u ** 2))
        self.inletFlowCoefficientDash = self.c1mdash / self.u1

        # Estimate outlet width
        # Eq 7.1 Gulich

        if self.outletWidthOverride == None:
            specificSpeedEURef = 100
            self.impellerOutletWidthDimensionless = 0.017 + (0.262 * (self.specificSpeedEU / specificSpeedEURef)) - (0.08 * ((self.specificSpeedEU / specificSpeedEURef) ** 2)) + (0.0093 * ((self.specificSpeedEU / specificSpeedEURef) ** 3))
            self.impellerOutletWidth = self.impellerOutletWidthDimensionless * self.impellerOutletDiameter
        else:
            self.impellerOutletWidth = self.outletWidthOverride

        # Calculate outlet blade blockage

        self.outletBladeBlockage = (1 - ((self.bladeThickness * self.numberOfBlades) / (pi * self.impellerOutletDiameter * sin(radians(self.outletBladeAngle)) * sin(radians(self.outletBladeInclination))))) ** (-1)

        # Calculate slip factor

        self.etalim = exp((-8.16 * sin(radians(self.outletBladeAngle))) / self.numberOfBlades)
        self.inletDiametersGeometricAverageDimensionless = self.inletDiametersGeometricAverage / self.impellerOutletDiameter

        if self.inletDiametersGeometricAverageDimensionless <= self.etalim:
            self.kw = 1
        else:
            self.kw = 1 - ((self.inletDiametersGeometricAverageDimensionless - self.etalim) / (1 - self.etalim))

        self.slipFactor = 0.98 * (1 - (sqrt(sin(radians(self.outletBladeAngle))) / (self.numberOfBlades ** 0.7))) * self.kw

        # Calculate areas

        self.inletArea = pi * (((self.impellerInletDiameter / 2) ** 2) - ((self.impellerHubDiameter / 2) ** 2))
        self.outletArea = pi * self.impellerOutletDiameter * self.impellerOutletWidth

        # Calculate head

        self.calculatedHead = ((self.hydraulicEfficiency * (self.u2 ** 2) / constants.G) * 
                                (self.slipFactor - ((self.impellerVolumeFlowRate / (self.outletArea * self.u2 * tan(radians(self.outletBladeAngle)))) *
                                (self.outletBladeBlockage + ((self.outletArea * self.inletDiametersGeometricAverageDimensionless * tan(radians(self.outletBladeAngle))) / (self.inletArea * tan(radians(self.approachFlowAngle)))))))
                                )

        # Calculate outlet velocity triangle 

        self.beta2B = self.outletBladeAngle

        self.c2m = self.impellerVolumeFlowRate / self.outletArea
        self.c2mDash = self.c2m * self.outletBladeBlockage
        self.outletFlowCoefficient = self.c2m / self.u2
        self.c2u = self.u2 * (self.slipFactor - ((self.c2m * self.outletBladeBlockage) / (self.u2 * tan(radians(self.beta2B)))))
        self.c2 = sqrt((self.c2m ** 2) + (self.c2u ** 2))
        self.c2Dash = sqrt((self.c2mDash ** 2) + (self.c2u ** 2))
        self.w2u = self.u2 - self.c2u
        self.w2 = sqrt((self.c2m ** 2) + (self.w2u ** 2))
        self.w2dash = sqrt((self.c2mDash ** 2) + (self.w2u ** 2))
        self.alpha2 = degrees(atan2(self.c2m, self.c2u))
        self.alpha2Dash = degrees(atan2(self.c2mDash, self.c2u))
        self.beta2 = degrees(atan2(self.c2m, self.w2u))
        self.beta2Dash = degrees(atan2(self.c2mDash, self.w2u))
        self.deviationAngle = self.beta2B - self.beta2
        self.deviationAngleDash = self.beta2B - self.beta2Dash

        # Outlet pitch

        self.outletPitch = pi * self.impellerOutletDiameter / self.numberOfBlades
        self.t2 = self.outletPitch # Alias

        # Creating aliases

        self.dn = self.impellerHubDiameter
        self.d1 = self.impellerInletDiameter
        self.d2 = self.impellerOutletDiameter
        self.b2 = self.impellerOutletWidth
