from math import pi
from pumpIt.constants import G

def pressureToHead(pressure: float, density: float):

    return pressure / (density * G)

def headToPressure(head: float, density: float):

    return head * density * G

def rpmToAngularSpeed(rpm: float):

    return (rpm / 60) * (2 * pi)

def angularSpeedToRPM(angularSpeed: float):

    return (angularSpeed / (2 * pi)) * 60

def volumeToMassFlowRate(volumeFlowRate: float, density: float):

    return volumeFlowRate * density

def massToVolumeFlowRate(massFlowRate: float, density: float):

    return massFlowRate / density