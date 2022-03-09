import fluid as fl
from impeller import Impeller
from meridional import Meridional
from bladeDesign import Blade
from pump import Pump
from voluteDesign import Volute, TrapezoidalCrossSection, RectangularCrossSection

fluid = fl.Fluid(density=787, viscosity=2.86e-3, vapourPressure=4100)
impeller = Impeller(suctionSpecificSpeedEU=650,
            suctionSidePressure=3e5,
            kgPerSec=1.32,
            headRise=295,
            fluid=fluid,
            shaftAllowableShearStress=8e7,
            numberOfBlades=5,
            approachFlowAngle=90,
            outletBladeAngle=22.5,
            inletBladeInnerDiameterRatio=1.15
)
meridionalSection = Meridional(impeller)
bladeDesign = Blade(meridionalSection, 3)
voluteCrossSection = TrapezoidalCrossSection()
volute = Volute(impeller, voluteCrossSection, dischargeAreaRatio=5)

pump = Pump(impeller, meridionalSection, bladeDesign, volute)

#pump.plotMeridional()
pump.plotResult(fullMeridional=True)
pump.plotVoluteDevelopmentPlan(polar=False)