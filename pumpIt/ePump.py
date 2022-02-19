import fluid as fl
from impeller import Impeller
from meridional import Meridional
from bladeDesign import Blade

fluid = fl.Fluid(density=787, viscosity=2.86e-3, vapourPressure=4100)
impeller = Impeller(suctionSpecificSpeedEU=650,
            suctionSidePressure=3e5,
            kgPerSec=1.32,
            headRise=295,
            fluid=fluid,
            shaftAllowableShearStress=8e7,
            numberOfBlades=6,
            approachFlowAngle=90,
            outletBladeAngle=22.5,
            inletBladeInnerDiameterRatio=1.15
)
meridionalSection = Meridional(impeller)
#meridionalSection.plot()
bladeDesign = Blade(meridionalSection, 3)