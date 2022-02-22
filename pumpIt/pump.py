from numpy import number
from impeller import Impeller
from meridional import Meridional
from bladeDesign import Blade
import matplotlib.pyplot as plt
from math import tan, radians, pi

class Pump:

    def __init__(self, impeller: Impeller, meridional: Meridional, blade: Blade):

        self.impeller = impeller
        self.meridional = meridional
        self.blade = blade

    def plotMeridional(self, show=True, full=False, shaft=True, internalStreamlines=True):

        ax = plt.axes()
        ax.plot(self.meridional.outerStreamlineXCoords, self.meridional.outerStreamlineYCoords, color="black")
        ax.plot(self.meridional.innerStreamlineXCoords, self.meridional.innerStreamlineYCoords, color="black")
        ax.plot([self.meridional.bladeLEShroudCoords[0], self.meridional.bladeLEHubCoords[0]], [self.meridional.bladeLEShroudCoords[1], self.meridional.bladeLEHubCoords[1]], color="grey")
        ax.plot([self.meridional.bladeTEShroudCoords[0], self.meridional.bladeTEHubCoords[0]], [self.meridional.bladeTEShroudCoords[1], self.meridional.bladeTEHubCoords[1]], color="grey")

        if internalStreamlines:
            for i in range(1, len(self.blade.streamlinesXCoords)-1):
                ax.plot(self.blade.streamlinesXCoords[i], self.blade.streamlinesYCoords[i], color='grey', ls='--', alpha=0.5, lw=0.5)

        if shaft:
            ax.plot([self.meridional.outerStreamlineXCoords[0], self.meridional.innerStreamlineXCoords[-1]], [self.impeller.shaftDiameter/2, self.impeller.shaftDiameter/2], color="grey", ls="dashdot")

        ax.axis("equal")
        if full:
            ax.plot(self.meridional.outerStreamlineXCoords, -self.meridional.outerStreamlineYCoords, color="black")
            ax.plot(self.meridional.innerStreamlineXCoords, -self.meridional.innerStreamlineYCoords, color="black")
            ax.plot([self.meridional.bladeLEShroudCoords[0], self.meridional.bladeLEHubCoords[0]], [-self.meridional.bladeLEShroudCoords[1], -self.meridional.bladeLEHubCoords[1]], color="grey")
            ax.plot([self.meridional.bladeTEShroudCoords[0], self.meridional.bladeTEHubCoords[0]], [-self.meridional.bladeTEShroudCoords[1], -self.meridional.bladeTEHubCoords[1]], color="grey")
            if shaft:
                ax.plot([self.meridional.outerStreamlineXCoords[0], self.meridional.innerStreamlineXCoords[-1]], [-self.impeller.shaftDiameter/2, -self.impeller.shaftDiameter/2], color="grey", ls="dashdot")
            if internalStreamlines:
                for i in range(1, len(self.blade.streamlinesXCoords)-1):
                    ax.plot(self.blade.streamlinesXCoords[i], -self.blade.streamlinesYCoords[i], color='grey', ls='--', alpha=0.5, lw=0.5)
        else:
            ax.set_ylim(ymin=0)

        if show:
            plt.show()

    def plotVelocityTriangle(self, area: str, withBlockage: bool = True, withoutBlockage: bool = True, decimalPoints: int = 2):

        if area.lower() == "inlet":
            title = "Inlet"
            u = self.impeller.u1
            c = self.impeller.c1
            cu = self.impeller.c1u
            cm = self.impeller.c1m
            w = self.impeller.w1
            beta = self.impeller.beta1
            betaB = self.impeller.beta1B
            alpha = self.impeller.alpha1
            cdash= self.impeller.c1dash
            cmdash = self.impeller.c1mdash
            wdash = self.impeller.w1dash
            betadash = self.impeller.beta1dash
            i = self.impeller.incidenceAngle
        if area.lower() == "outlet":
            title = "Outlet"
            u = self.impeller.u2
            c = self.impeller.c2
            cu = self.impeller.c2u
            cm = self.impeller.c2m
            w = self.impeller.w2
            beta = self.impeller.beta2
            betaB = self.impeller.beta2B
            alpha = self.impeller.alpha2
            cdash= self.impeller.c2Dash
            cmdash = self.impeller.c2mDash
            wdash = self.impeller.w2dash
            betadash = self.impeller.beta2Dash
            i = self.impeller.deviationAngle

        hw = 0.25
        hl = 0.5
        standardOffset = 2
        lineClashOffset = 0.1
        padding = 1

        ax = plt.axes()

        # Impeller velocity
        ax.arrow(u, 0, -u, 0, head_width=hw, head_length=hl, length_includes_head=True, edgecolor="orange", facecolor="orange", label="U: " + str(round(u, 2)) + " [m/s]")
        #ax.annotate("U: \n" + str(round(u, 2)) + " [m/s]", xy=(u/2, 0 - standardOffset))

        # Absolute velocity
        ax.arrow(u, 0, -cu, cm, head_width=hw, head_length=hl, length_includes_head=True, edgecolor="red", facecolor="red", label="C: " + str(round(c, 2)) + " [m/s]")
        #ax.annotate("C: \n" + str(round(c, 2)) + " [m/s]", xy=(u-2*standardOffset, (0+cm)/2 - standardOffset))

        # Absolute velocity components
        if alpha != 90:
            ax.arrow(u, lineClashOffset, -cu, 0, head_width=hw, head_length=hl, length_includes_head=True, edgecolor="lightgrey", facecolor="lightgrey", label="Cu: " + str(round(cu, 2)) + " [m/s]")
            ax.arrow(u-cu, 0, 0, cm, head_width=hw, head_length=hl, length_includes_head=True, edgecolor="grey", facecolor="grey", label="Cm: " + str(round(cm, 2)) + " [m/s]")

        # Relative velocity
        ax.arrow(0, 0, u-cu, cm, head_width=hw, head_length=hl, length_includes_head=True, edgecolor="blue", facecolor="blue", label="W: " + str(round(w, 2)) + " [m/s]")
        #ax.annotate("W: \n" + str(round(w, 2)) + " [m/s]", xy=((u+cu) / 2, (0+cm)/2 - standardOffset))

        # Absolute velocity with blockage
        ax.arrow(u, 0, -cu, cmdash, head_width=hw, head_length=hl, length_includes_head=True, edgecolor="red", facecolor="red", linestyle="dotted", label="C': " + str(round(cdash, 2)) + " [m/s]")
        # Relative velocity with blockage
        ax.arrow(0, 0, u-cu, cmdash, head_width=hw, head_length=hl, length_includes_head=True, edgecolor="blue", facecolor="blue", linestyle="dotted", label="W': " + str(round(wdash, 2)) + " [m/s]")
        # Blade angle line
        ax.arrow(0, 0, u-cu, (u-cu) * tan(radians(betaB)), head_width=0, head_length=0, length_includes_head=True, edgecolor="black", facecolor="black", linestyle="dotted", label="BetaB: " + str(round(betaB, 2)) + " [deg]")

        ax.legend()
        ax.set_title(title + " Velocity Triangle")
        ax.set_xlabel("Speed [m/s]")
        ax.set_ylabel("Speed [m/s]")

        plt.gca().set_ylim(bottom=-standardOffset - padding)

        plt.show()
    
    def plotBladesPlanView(self, numberOfBlades = 1):

        if type(numberOfBlades) == str:
            if numberOfBlades.lower() == "all":
                numberOfBlades = self.impeller.numberOfBlades
        
        ax = plt.axes(projection='polar')
        colors = ['black', 'pink', 'blue', 'orange', 'green', 'grey', 'red', 'purple']

        for blade in range(numberOfBlades):

            for i in range(len(self.blade.streamlinesBladeAngles)):

                epsilonsch = blade * (360 / self.impeller.numberOfBlades)

                r = self.blade.streamlinesYCoords[i][-1]
                rs = [r]
                epsilonschs = [epsilonsch]
                epsilonschsRadians = [radians(epsilonsch)]

                for j in range(len(self.blade.streamlinesDeltaRs[i])):
                    r -= self.blade.streamlinesDeltaRs[i][j]
                    rs.append(r)
                    deltaEpsilonsch = 360 * self.blade.streamlinesDeltaUs[i][j] / (2 * pi * r)
                    epsilonsch += deltaEpsilonsch
                    epsilonschs.append(epsilonsch)
                    epsilonschsRadians.append(radians(epsilonsch))

                ax.plot(epsilonschsRadians, rs, color=colors[blade])
            
        plt.show()
    
    def printImpellerResults(self, 
    inputs: bool = True,
    performance: bool = True,
    sizes: bool = True,
    velocityTriangles: bool = True,
    suction: bool = True,
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
                formatNumber(self.impeller.rpm), " [RPM] \n",
                formatNumber(self.impeller.revPerSec), " [s^-1] \n",
                formatNumber(self.impeller.radPerSec), " [rad s^-1] \n",
                "\nVolume Flow Rate: \n",
                formatNumber(self.impeller.metreCubedPerSec), " [m^3 s^-1] \n",
                formatNumber(self.impeller.metreCubedPerHour), " [m^3 hr^-1] \n",
                formatNumber(self.impeller.gallonPerMin), " [gpm] \n",
                "\nMass Flow Rate: \n",
                formatNumber(self.impeller.kgPerSec), " [kg s^-1] \n",
                "\nHead Rise: \n",
                formatNumber(self.impeller.headRise), " [m] \n",
                formatNumber(self.impeller.headRiseFeet), " [ft] \n",
                "\nPressure Rise: \n",
                formatNumber(self.impeller.pressureRise/1e5), " [Bar] \n"
            ]

            inputsString = "".join(inputsOutput)

            print(inputsString)

        if performance:

            specificSpeedOutput = [
                printSeparator, "Specific Speed \n", printSeparator,
                "\nDimensionless Specific Speed: \n",
                formatNumber(self.impeller.specificSpeedDimensionless), " \n",
                "\nEU Specific Speed: \n",
                formatNumber(self.impeller.specificSpeedEU), " \n",
                "\nUS Specific Speed: \n",
                formatNumber(self.impeller.specificSpeedUS), " \n"
            ]

            specificSpeedString = "".join(specificSpeedOutput)

            print(specificSpeedString)

            efficienciesOutput = [
                printSeparator, "Efficiencies \n", printSeparator,
                "\nOverall Efficiency: \n",
                formatNumber(self.impeller.overallEfficiency), " \n",
                "\nHydraulic Efficiency: \n",
                formatNumber(self.impeller.hydraulicEfficiency), " \n",
                "\nVolumetric Efficiency: \n",
                formatNumber(self.impeller.volumetricEfficiency), " \n"
            ]

            efficienciesString = "".join(efficienciesOutput)

            print(efficienciesString)

            coefficientsOutput = [
                printSeparator, "Coefficients \n", printSeparator,
                "\nHead Coefficient: \n",
                formatNumber(self.impeller.headCoefficient), " \n"
            ]

            coefficientsString = "".join(coefficientsOutput)

            print(coefficientsString)

        if sizes:

            sizesOutput = [
                printSeparator, "Sizes \n", printSeparator,
                "\nShaft Diameter: \n",
                formatNumber(self.impeller.shaftDiameter*1e3), " [mm]\n",
                "\nHub Diameter: \n",
                formatNumber(self.impeller.impellerHubDiameter*1e3), " [mm]\n",
                "\nImpeller Inlet Diameter: \n",
                formatNumber(self.impeller.impellerInletDiameter*1e3), " [mm]\n",
                "\nImpeller Outer Diameter: \n",
                formatNumber(self.impeller.impellerOutletDiameter*1e3), " [mm]\n",
                "\nImpeller Outlet Width: \n",
                formatNumber(self.impeller.impellerOutletWidth*1e3), " [mm]\n"
                "\nBlade Thickness: \n",
                formatNumber(self.impeller.bladeThickness*1e3), " [mm]\n"
            ]

            sizesString = "".join(sizesOutput)

            print(sizesString)

        if velocityTriangles:

            velocityTriangleOutput = [
                printSeparator, "Inlet Velocity Triangle Without Blockage\n", printSeparator,
                "\nTip Speed u1: \n",
                formatNumber(self.impeller.u1), " [m s^-1]\n",
                "\nAbsolute Speed c1: \n",
                formatNumber(self.impeller.c1), " [m s^-1]\n", "(Meridional c1m: ", formatNumber(self.impeller.c1m), " [m s^-1])\n", "(Circumferential c1u: ", formatNumber(self.impeller.c1u), " [m s^-1])\n",
                "\nRelative Speed w1: \n", 
                formatNumber(self.impeller.w1), " [m s^-1]\n",
                "\nApproach Flow Angle alpha1: \n", 
                formatNumber(self.impeller.approachFlowAngle), " [deg]\n",
                "\nRelative Flow Angle beta1: \n", 
                formatNumber(self.impeller.beta1), " [deg]\n",
                "\nFlow Coefficient: \n", 
                formatNumber(self.impeller.inletFlowCoefficient), " \n",
                printSeparator, "Inlet Velocity Triangle With Blockage\n", printSeparator,
                "\nBlade Blockage: \n",
                formatNumber(self.impeller.inletBladeBlockage), " \n",
                "\nTip Speed u1: \n",
                formatNumber(self.impeller.u1), " [m s^-1]\n",
                "\nAbsolute Speed c1': \n",
                formatNumber(self.impeller.c1dash), " [m s^-1]\n", "(Meridional c1m': ", formatNumber(self.impeller.c1mdash), " [m s^-1])\n", "(Circumferential c1u: ", formatNumber(self.impeller.c1u), " [m s^-1])\n",
                "\nRelative Speed w1': \n", 
                formatNumber(self.impeller.w1dash), " [m s^-1]\n",
                "\nApproach Flow Angle alpha1: \n", 
                formatNumber(self.impeller.approachFlowAngle), " [deg]\n",
                "\nRelative Flow Angle beta1': \n", 
                formatNumber(self.impeller.beta1dash), " [deg]\n",
                "\nBlade Angle beta1B: \n", 
                formatNumber(self.impeller.beta1B), " [deg]\n",
                "\nIncidence Angle: \n", 
                formatNumber(self.impeller.incidenceAngle), " [deg]\n",
                "\nFlow Coefficient: \n", 
                formatNumber(self.impeller.inletFlowCoefficientDash), " \n",
                printSeparator, "Outlet Velocity Triangle Without Blockage\n", printSeparator,
                "\nTip Speed u2: \n",
                formatNumber(self.impeller.u2), " [m s^-1]\n",
                "\nAbsolute Speed c2: \n",
                formatNumber(self.impeller.c2), " [m s^-1]\n", "(Meridional c2m: ", formatNumber(self.impeller.c2m), " [m s^-1])\n", "(Circumferential c2u: ", formatNumber(self.impeller.c2u), " [m s^-1])\n",
                "\nRelative Speed w2: \n", 
                formatNumber(self.impeller.w2), " [m s^-1]\n",
                "\nAbsolute Outlet Flow Angle alpha2: \n", 
                formatNumber(self.impeller.alpha2), " [deg]\n",
                "\nRelative Flow Angle beta2: \n", 
                formatNumber(self.impeller.beta2), " [deg]\n",
                "\nDeviation Angle delta: \n", 
                formatNumber(self.impeller.deviationAngle), " [deg]\n",
                "\nFlow Coefficient: \n", 
                formatNumber(self.impeller.outletFlowCoefficient), " \n",
                printSeparator, "Outlet Velocity Triangle With Blockage\n", printSeparator,
                "\nBlade Blockage: \n",
                formatNumber(self.impeller.outletBladeBlockage), " \n",
                "\nTip Speed u2: \n",
                formatNumber(self.impeller.u2), " [m s^-1]\n",
                "\nAbsolute Speed c2': \n",
                formatNumber(self.impeller.c2Dash), " [m s^-1]\n", "(Meridional c2m': ", formatNumber(self.impeller.c2mDash), " [m s^-1])\n", "(Circumferential c2u: ", formatNumber(self.impeller.c2u), " [m s^-1])\n",
                "\nRelative Speed w2': \n", 
                formatNumber(self.impeller.w2dash), " [m s^-1]\n",
                "\nApproach Flow Angle alpha2': \n", 
                formatNumber(self.impeller.alpha2Dash), " [deg]\n",
                "\nRelative Flow Angle beta2': \n", 
                formatNumber(self.impeller.beta2Dash), " [deg]\n",
                "\nBlade Angle beta2B: \n", 
                formatNumber(self.impeller.beta2B), " [deg]\n",
                "\nDeviation Angle delta': \n", 
                formatNumber(self.impeller.deviationAngleDash), " [deg]\n",
                "\nFlow Coefficient: \n", 
                formatNumber(self.impeller.outletFlowCoefficient), " \n"
            ]

            velocityTriangleString = "".join(velocityTriangleOutput)

            print(velocityTriangleString)

        if suction:

            suctionOutput = [
                printSeparator, "Suction Performance\n", printSeparator,
                "\nSuction Specific Speed: \n",
                formatNumber(self.impeller.suctionSpecificSpeedEU), " \n",
                "\nAvailable NPSH: \n",
                formatNumber(self.impeller.NPSHA), " [m]\n",
                "\nRequired (3% Cavitation) NPSH: \n",
            formatNumber(self.impeller.NPSH3), " [m]\n"
            ]

            suctionString = "".join(suctionOutput)

            print(suctionString)

        print("")