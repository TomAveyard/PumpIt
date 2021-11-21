from math import pi
from pumpIt.constants import G

def pressureToHead(pressure: float, density: float):

    return pressure / (density * G)

def headToPressure(head: float, density: float):

    return head * density * G

def RPMToAngularSpeed(RPM: float):

    return (RPM / 60) * (2 * pi)

def angularSpeedToRPM(angularSpeed: float):

    return (angularSpeed / (2 * pi)) * 60

def RPMToRevPerSec(RPM: float):

    return RPM / 60

def revPerSecToRPM(revPerSec: float):

    return revPerSec * 60

def volumeToMassFlowRate(volumeFlowRate: float, density: float):

    return volumeFlowRate * density

def massToVolumeFlowRate(massFlowRate: float, density: float):

    return massFlowRate / density

def metreToFeet(m: float):

    return m * 3.281

def feetToMetre(ft: float):

    return ft / 3.281

def metreCubedPerSecToGallonPerMin(m3PerSec: float):

    return m3PerSec * 15850.323141

def GallonPerMinToMetreCubedPerSec(GPM: float):

    return GPM / 15850.323141