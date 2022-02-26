import fluid as fl
from impeller import Impeller
from meridional import Meridional
from bladeDesign import Blade
from pump import Pump

fluid = fl.Fluid(density=1000, viscosity=1e-3, vapourPressure=3173)
impeller = Impeller(rpm=1450,
            suctionSpecificSpeedEU=165,
            metreCubedPerSec=0.07777,
            headRise=20,
            fluid=fluid,
            shaftAllowableShearStress=8e7,
            numberOfBlades=6,
            approachFlowAngle=90,
            outletBladeAngle=22.5,
            inletBladeInnerDiameterRatio=1.7,
            hubShaftDiameterRatio=1.8,
            isSuctionImpeller=True,
            headCoefficientCorrelation="gulich"
)
meridionalSection = Meridional(impeller, outerStreamlineCircularSectionArcLength=35, outerStreamlineOutletAngle=20, innerStreamlineOutletAngle=10)
bladeDesign = Blade(meridionalSection, 3)

pump = Pump(impeller, meridionalSection, bladeDesign)
pump.plotMeridional()
pump.plotPlanView()
