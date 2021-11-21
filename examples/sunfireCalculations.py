import pumpIt.pumpDesign as pd

fluid = pd.Fluid(density=786, viscosity=0.00177, vapourPressure=5000)
pump = pd.Pump(volumeFlowRate=1.5, inletPressure=1.01325e5, outletPressure=25e5, fluid=fluid, shaftRadius=0)

#pump.setSuctionSpecificSpeed(suctionSpecificSpeed=4.5)
pump.setSpeed(rpm=1780)
pump.setBlockageEffect(blockageEffect=0.9)
pump.calculate()

print(pump.specificSpeed)
print(pump.suctionSpecificSpeed)
print(pump.rpm)
print(pump.npsh)
print(pump.impellerOuterRadius*1e2)
print(pump.outerWidth*1e2)
print(pump.outerTipSpeed)
print(pump.eyeCavitationCoefficient)
print(pump.eyeFlowCoefficient)
print(pump.eyeRadius*1e2)