import fluid as fl
from impeller import Impeller
from meridional import Meridional
from bladeDesign import Blade
from pump import Pump
from voluteDesign import Volute, TrapezoidalCrossSection, RectangularCrossSection, CircularCrossSection

fluid = fl.Fluid(density=787, viscosity=2.86e-3, vapourPressure=4100)
impeller = Impeller(suctionSpecificSpeedEU=650,
            suctionSidePressure=3e5,
            kgPerSec=1.32,
            headRise=295,
            fluid=fluid,
            shaftAllowableShearStress=42e8,
            shaftDiameterSafetyFactor=4,
            shaftDiameterOverride=8e-3,
            numberOfBlades=5,
            approachFlowAngle=90,
            outletBladeAngle=18.5,
            inletBladeInnerDiameterRatio=1.2,
            isSuctionImpeller=True
)
meridionalSection = Meridional(impeller, numberOfPoints=1000)
bladeDesign = Blade(meridionalSection, 3)
voluteCrossSection = CircularCrossSection()
volute = Volute(impeller, voluteCrossSection, dischargeExitDiameter=0.0127)

pump = Pump(impeller, meridionalSection, bladeDesign, volute)

#pump.plotVelocityTriangle("inlet")
#pump.plotMeridional()
#pump.plotPlanView(bladesOrStreamlines="blades", numberOfBlades=1)
pump.plotResult(fullMeridional=True)
#pump.plotVoluteDevelopmentPlan(polar=False)
pump.outputBladeGenFiles()
