from pumpIt.conversions import pressureToHead

class Fluid:

    def __init__(self, density:float = None, viscosity:float = None, vapourPressure: float = None):

        self.density = density
        self.viscosity = viscosity
        self.vapourPressure = vapourPressure
        self.vapourHead = pressureToHead(vapourPressure, self.density)
        