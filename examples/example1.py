import pumpIt.pumpDesign as pd

fluid = pd.Fluid(density=1000, viscosity=1.0016e-3, vapourPressure=2300)
pump = pd.Pump(volumeFlowRate=0.1577, inletHead=4.5, outletHead=36.2, fluid=fluid, shaftRadius=0)

#pump.setSuctionSpecificSpeed(suctionSpecificSpeed=4.5)
pump.setSpeed(rpm=1780)
pump.setBlockageEffect(blockageEffect=0.9)
pump.calculate(convergenceCriteria=1e-2)

print("Specific Speed: " + str(pump.specificSpeed))
print("Suction Specific Speed: " + str(pump.suctionSpecificSpeed))
print("RPM: " + str(pump.rpm))
print("NPSH: " + str(pump.npsh))
print("Head Coefficient: " + str(pump.headCoefficient))
print("Impeller Outer Radius: " + str(pump.impellerOuterRadius*1e2))
print("Outer Width: " + str(pump.outerWidth*1e2))
print("Outer Tip Speed: " + str(pump.outerTipSpeed))
print("Eye Cavitation Coefficient: " + str(pump.eyeCavitationCoefficient))
print("Eye Flow Coefficient: " + str(pump.eyeFlowCoefficient))
print("Eye Radius: " + str(pump.eyeRadius*1e2))
print("Shaft Radius: " + str(pump.shaftRadius))
print("Shaft Eye Ratio: " + str(pump.shaftEyeRatio))