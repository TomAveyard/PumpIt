from meridional import Meridional

class Blade:

    def __init__(self, meridionalSection: Meridional, numberOfStreamlines: int = 3) -> None:
        
        self.meridionalSection = meridionalSection
        self.outerStreamlineXCoords = self.meridionalSection.outerStreamlineXCoords
        self.outerStreamlineYCoords = self.meridionalSection.outerStreamlineYCoords
        self.innerStreamlineXCoords = self.meridionalSection.innerStreamlineXCoords
        self.innerStreamlineYCoords = self.meridionalSection.innerStreamlineYCoords

