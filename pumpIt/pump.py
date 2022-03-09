from impeller import Impeller
from meridional import Meridional
from bladeDesign import Blade
from plottingHelper import polarToCartesian
from voluteDesign import Volute
import matplotlib.pyplot as plt
import matplotlib
from math import tan, radians, pi, degrees
from plottingHelper import rotateCartesianCoord

class Pump:

    def __init__(self, impeller: Impeller, meridional: Meridional, blade: Blade, volute: Volute):

        self.impeller = impeller
        self.meridional = meridional
        self.blade = blade
        self.volute = volute

    def plotMeridional(self, show=True, full=False, shaft=True, internalStreamlines=True, ax=None):
        
        if not ax:
            ax = plt.axes()
        
        ax.plot(self.meridional.outerStreamlineXCoords, self.meridional.outerStreamlineYCoords, color="black")
        ax.plot(self.meridional.innerStreamlineXCoords, self.meridional.innerStreamlineYCoords, color="black")
        ax.plot([self.meridional.bladeLEShroudCoords[0], self.meridional.bladeLEHubCoords[0]], [self.meridional.bladeLEShroudCoords[1], self.meridional.bladeLEHubCoords[1]], color="grey")
        ax.plot([self.meridional.bladeTEShroudCoords[0], self.meridional.bladeTEHubCoords[0]], [self.meridional.bladeTEShroudCoords[1], self.meridional.bladeTEHubCoords[1]], color="grey")

        if internalStreamlines:
            for i in range(1, len(self.blade.streamlinesXCoords)-1):
                ax.plot(self.blade.streamlinesMeridionalXCoords[i], self.blade.streamlinesMeridionalYCoords[i], color='grey', ls='--', alpha=0.5, lw=0.5)

        if shaft:
            ax.plot([self.meridional.outerStreamlineXCoords[0], self.meridional.innerStreamlineXCoords[-1]], [self.impeller.shaftDiameter/2, self.impeller.shaftDiameter/2], color="grey", ls="dashdot")

        ax.axis("equal")
        ax.set_title("Impeller Meridional View")
        ax.set_xlabel("Axial Distance [m]")
        ax.set_ylabel("Radial Distance [m]")

        if full:
            ax.plot(self.meridional.outerStreamlineXCoords, -self.meridional.outerStreamlineYCoords, color="black")
            ax.plot(self.meridional.innerStreamlineXCoords, -self.meridional.innerStreamlineYCoords, color="black")
            ax.plot([self.meridional.bladeLEShroudCoords[0], self.meridional.bladeLEHubCoords[0]], [-self.meridional.bladeLEShroudCoords[1], -self.meridional.bladeLEHubCoords[1]], color="grey")
            ax.plot([self.meridional.bladeTEShroudCoords[0], self.meridional.bladeTEHubCoords[0]], [-self.meridional.bladeTEShroudCoords[1], -self.meridional.bladeTEHubCoords[1]], color="grey")
            if shaft:
                ax.plot([self.meridional.outerStreamlineXCoords[0], self.meridional.innerStreamlineXCoords[-1]], [-self.impeller.shaftDiameter/2, -self.impeller.shaftDiameter/2], color="grey", ls="dashdot")
            if internalStreamlines:
                for i in range(1, len(self.blade.streamlinesMeridionalXCoords)-1):
                    ax.plot(self.blade.streamlinesMeridionalXCoords[i], [-y for y in self.blade.streamlinesMeridionalYCoords[i]], color='grey', ls='--', alpha=0.5, lw=0.5)
        else:
            ax.set_ylim(ymin=0)

        if show:
            plt.show()

    def plotVelocityTriangle(self, area: str, withBlockage: bool = True, withoutBlockage: bool = True, decimalPoints: int = 2, show=True, ax=None, equalAxes=False):

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

        if not ax:
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
        ax.set_title("Impeller " + title + " Velocity Triangle")
        ax.set_xlabel("Speed [m/s]")
        ax.set_ylabel("Speed [m/s]")

        plt.gca().set_ylim(bottom=-standardOffset - padding)

        if equalAxes:
            ax.axis("equal")

        if show:
            plt.show()
    
    def plotPlanView(self, plotType="polar", numberOfBlades=1, bladesOrStreamlines="blades", color="black", plotLE=True, show=True, ax=None):
        
        plotType = plotType.lower()
        bladesOrStreamlines = bladesOrStreamlines.lower()
        color = color.lower()

        if type(numberOfBlades) == str:
            if numberOfBlades.lower() == "all":
                numberOfBlades = self.impeller.numberOfBlades
        
        if not ax:
            if plotType == "polar":
                ax = plt.axes(projection='polar')
            elif plotType == "cartesian":
                ax = plt.axes()
        
        ax.set_title("Impeller Plan View")

        if plotType == "cartesian":
            ax.axis("equal")
            ax.set_xlabel("x [m]")
            ax.set_ylabel("y [m]")
        if plotType == "polar":
            ax.set_xlabel("Radius [m]")
            ax.set_ylabel("Theta [Degrees]")
        
        colors = ['black', 'pink', 'blue', 'orange', 'green', 'grey', 'red', 'purple']

        for blade in range(numberOfBlades):

            if color == "multi":
                c = colors[blade]
            else:
                c = color

            bladeRotation = blade * (2 * pi / self.impeller.numberOfBlades)

            if plotType == "cartesian":

                bladeRotation = degrees(bladeRotation)

                if bladesOrStreamlines == "streamlines" or bladesOrStreamlines == "streamline":

                    for i in range(len(self.blade.streamlinesXCoords)):

                        rotatedXCoords = []
                        rotatedYCoords = []

                        for j in range(len(self.blade.streamlinesXCoords[i])):

                            rotatedCoords = rotateCartesianCoord(self.blade.streamlinesXCoords[i][j], self.blade.streamlinesYCoords[i][j], bladeRotation)
                            rotatedXCoords.append(rotatedCoords[0])
                            rotatedYCoords.append(rotatedCoords[1])
                        
                        ax.plot(rotatedXCoords, rotatedYCoords, color=c)
                
                if bladesOrStreamlines == "blades" or bladesOrStreamlines == "blade":

                    for i in range(len(self.blade.bladesXCoords)):

                        rotatedXCoords = []
                        rotatedYCoords = []

                        for j in range(len(self.blade.bladesXCoords[i])):
                            
                            rotatedCoords = rotateCartesianCoord(self.blade.bladesXCoords[i][j], self.blade.bladesYCoords[i][j], bladeRotation)
                            rotatedXCoords.append(rotatedCoords[0])
                            rotatedYCoords.append(rotatedCoords[1])
                        
                        ax.plot(rotatedXCoords, rotatedYCoords, color=c)

            if plotType == "polar":

                if bladesOrStreamlines == "streamlines" or bladesOrStreamlines == "streamline":

                    for i in range(len(self.blade.streamlinesEpsilonSchsRadians)):

                        epsilonschsRadians = []

                        for j in range(len(self.blade.streamlinesEpsilonSchsRadians[i])):

                            epsilonschsRadians.append(self.blade.streamlinesEpsilonSchsRadians[i][j] + bladeRotation)

                        ax.plot(epsilonschsRadians, self.blade.streamlinesRadiuses[i], color=c)

                elif bladesOrStreamlines == "blades" or bladesOrStreamlines == "blade":

                    for i in range(len(self.blade.bladesEpsilonSchsRadians)):

                        epsilonschsRadians = []

                        for j in range(len(self.blade.bladesEpsilonSchsRadians[i])):

                            epsilonschsRadians.append(self.blade.bladesEpsilonSchsRadians[i][j] + bladeRotation)

                        ax.plot(epsilonschsRadians, self.blade.bladesRadiuses[i], color=c)

            if plotLE:

                if plotType == "cartesian":

                    rotatedXCoords = []
                    rotatedYCoords = []

                    for i in range(len(self.blade.bladeLEXCoords)):

                        rotatedCoords = rotateCartesianCoord(self.blade.bladeLEXCoords[i], self.blade.bladeLEYCoords[i], bladeRotation)
                        rotatedXCoords.append(rotatedCoords[0])
                        rotatedYCoords.append(rotatedCoords[1])
                    
                    ax.plot(rotatedXCoords, rotatedYCoords, color=c)
                
                if plotType == "polar":

                    LEEpsilonRadians = []

                    for i in range(len(self.blade.bladeLEEpsilonsRadians)):

                        LEEpsilonRadians.append(self.blade.bladeLEEpsilonsRadians[i] + bladeRotation)

                    ax.plot(LEEpsilonRadians, self.blade.bladeLERadiuses, color=c)
        
        if show:
            plt.show()

    def plotVoluteDevelopmentPlan(self, polar=True, show=True, ax=None):
        
        if not ax:
            if polar:
                ax = plt.axes(projection='polar')
            else:
                ax = plt.axes()
                ax.axis("equal")
        
        ax.set_title("Volute Development Plan View")
        if polar:
            ax.plot([radians(x) for x in range(360)], [self.impeller.d2/2 for y in range(360)], color="grey")
            ax.plot([radians(x) for x in range(360)], [self.volute.rzDash for y in range(360)], color="black", ls = "--")
            ax.plot([radians(x) for x in self.volute.totalEpsilons], self.volute.totalRCoords, color="black")
        else:
            xCoordsImpeller = []
            yCoordsImpeller = []
            xCoordsrzDash = []
            yCoordsrzDash = []
            for i in range(len(self.volute.rAs)):
                coords = polarToCartesian(self.impeller.d2/2, self.volute.epsilons[i])
                xCoordsImpeller.append(coords[0])
                yCoordsImpeller.append(coords[1])
                coords = polarToCartesian(self.volute.rzDash, self.volute.epsilons[i])
                xCoordsrzDash.append(coords[0])
                yCoordsrzDash.append(coords[1])
            ax.plot(self.volute.totalXCoords, self.volute.totalYCoords, color="black")
            ax.plot(xCoordsImpeller, yCoordsImpeller, color = "grey")
            ax.plot(xCoordsrzDash, yCoordsrzDash, color="black", ls="--")
            ax.set_xlabel("x Distance [m]")
            ax.set_ylabel("y Distance [m]")

        if show:
            plt.show()

    def plotVoluteDevelopmentCrossSection(self, corrected=True, show=True, ax=None, numberOfIntermediateSections=4):

        if not ax:
            ax = ax = plt.axes()
        if corrected:
            rAs = self.volute.rAsCorrected
        else:
            rAs = self.volute.rAs

        ax.set_title("Volute Development Cross Section View")
        ax.axis("equal")
        ax.set_xlabel("Axial Distance [m]")
        ax.set_ylabel("Radial Distance [m]")
        self.volute.voluteCrossSection.generateCoords(self.volute.rzDash, self.volute.b3, rAs[0] - self.volute.rzDash)
        ax.plot(self.volute.voluteCrossSection.aCoords, self.volute.voluteCrossSection.rCoords, color="black")
        self.volute.voluteCrossSection.generateCoords(self.volute.rzDash, self.volute.b3, rAs[-1] - self.volute.rzDash)
        ax.plot(self.volute.voluteCrossSection.aCoords, self.volute.voluteCrossSection.rCoords, color="black")

        idxIncrement = round(len(rAs) / (numberOfIntermediateSections + 1))
        i = idxIncrement
        while i < len(rAs):
            self.volute.voluteCrossSection.generateCoords(self.volute.rzDash, self.volute.b3, rAs[i] - self.volute.rzDash)
            ax.plot(self.volute.voluteCrossSection.aCoords, self.volute.voluteCrossSection.rCoords, color="grey", ls="--", alpha=0.5)
            i += idxIncrement
        if show:
            plt.show()

    def plotResult(self, fullMeridional=False):
        
        fig = plt.figure(constrained_layout=True)
        gs = fig.add_gridspec(ncols=4, nrows=2, height_ratios=[1, 1], width_ratios=[1, 2, 2, 2])
        ax1 = fig.add_subplot(gs[0,0])
        ax2 = fig.add_subplot(gs[1,0])
        ax3 = fig.add_subplot(gs[0:,1])
        ax4 = fig.add_subplot(gs[0:,2])
        ax5 = fig.add_subplot((gs[0,3]))
        ax6 = fig.add_subplot((gs[1,3]))

        self.plotVelocityTriangle(area="inlet", ax=ax1, show=False, equalAxes=True)
        self.plotVelocityTriangle(area="outlet", ax=ax2, show=False, equalAxes=True)
        self.plotMeridional(show=False, ax=ax3, full=fullMeridional)
        self.plotPlanView(plotType="cartesian", numberOfBlades="all", show=False, ax=ax4)
        self.plotVoluteDevelopmentPlan(polar=False, show=False, ax=ax5)
        self.plotVoluteDevelopmentCrossSection(show=False, ax=ax6)

        fig.suptitle("Pump Design Result", fontsize=20)
        plot_backend = matplotlib.get_backend()

        def moveFigure(f, x, y, fullScreen=False):
            """Move figure's upper left corner to pixel (x, y)"""
            backend = matplotlib.get_backend()
            mng = plt.get_current_fig_manager()
            if backend == 'TkAgg':
                f.canvas.manager.window.wm_geometry("+%d+%d" % (x, y))
                if fullScreen:
                    mng.resize(*mng.window.maxsize())
            elif backend == 'wXAgg':
                f.canvas.manager.window.SetPosition((x, y))
                if fullScreen:
                    mng.frame.Maximize(True)
            elif backend == "Qt4Agg":
                # This works for QT and GTK
                # You can also use window.setGeometry
                f.canvas.manager.window.move(x, y)
                if fullScreen:
                    mng.window.showMaximized()

        moveFigure(fig, 0, 0, fullScreen=True)

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