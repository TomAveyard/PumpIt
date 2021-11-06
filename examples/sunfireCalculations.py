import pumpIt.pumpDesign as pd

fluid = pd.Fluid(density=786, viscosity=0.00177, vapourPressure=5332.9)
pump = pd.Pump(massFlowRate=1.3, inletPressure=3e5, outletPressure=90e5, fluid=fluid)

#pump.setSuctionSpecificSpeed(suctionSpecificSpeed=4)
pump.setSpeed(rpm=25000)
pump.setBlockageEffect(blockageEffect=0.85)
pump.calculate()

print(pump.specificSpeed)
print(pump.suctionSpecificSpeed)
print(pump.rpm)
print(pump.npsh)
print(pump.impellerOuterRadius*1e3)
print(pump.outerWidth*1e3)