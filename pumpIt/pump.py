from impeller import Impeller
from meridional import Meridional
from bladeDesign import Blade
import matplotlib.pyplot as plt

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