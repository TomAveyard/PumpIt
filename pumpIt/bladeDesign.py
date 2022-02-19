from meridional import Meridional
from math import sin, cos, atan2, degrees, radians, sqrt, tan
import matplotlib.pyplot as plt

class Blade:

    def __init__(self, meridionalSection: Meridional, numberOfStreamlines: int = 3) -> None:
        
        self.meridionalSection = meridionalSection
        self.numberOfStreamlines = numberOfStreamlines
        self.outerStreamlineXCoords = self.meridionalSection.outerStreamlineXCoords
        self.outerStreamlineYCoords = self.meridionalSection.outerStreamlineYCoords
        self.innerStreamlineXCoords = self.meridionalSection.innerStreamlineXCoords
        self.innerStreamlineYCoords = self.meridionalSection.innerStreamlineYCoords

        # Initialise lists to hold coords for each streamline
        self.streamlinesXCoords = []
        self.streamlinesYCoords = []
        # Append on the number of lists required
        for i in range(self.numberOfStreamlines):
            self.streamlinesXCoords.append([])
            self.streamlinesYCoords.append([])

        i = 0
        while i < self.meridionalSection.numberOfPoints:
            
            # Calculate inner to outer streamline geometry
            deltax = abs(self.innerStreamlineXCoords[i] - self.outerStreamlineXCoords[i])
            deltay = abs(self.outerStreamlineYCoords[i] - self.innerStreamlineYCoords[i])
            theta = degrees(atan2(deltax, deltay))
            b = sqrt((deltax ** 2) + (deltay ** 2))
            
            streamlinesY = [self.outerStreamlineYCoords[i]]
            streamlinesX = [self.outerStreamlineXCoords[i]]
            streamTubeWidths = []

            for j in range(1, self.numberOfStreamlines + 1):
                
                # Calculate coords for each streamline so that each streamtube formed by the streamlines has an equal amount of flow through it
                # Meridional velocity cm is assumed constant
                streamlinesY.append(sqrt((streamlinesY[j-1] ** 2) - ((b * cos(radians(theta)) *  (self.outerStreamlineYCoords[i] + self.innerStreamlineYCoords[i])) / (self.numberOfStreamlines + 1))))
                streamlinesX.append(streamlinesX[j-1] + (abs(streamlinesY[j-1] - streamlinesY[j])) * tan(radians(theta)))
            
            streamlinesY.append(self.innerStreamlineYCoords[i])

            for j in range(1, len(streamlinesY)-1):

                width = sqrt(((streamlinesX[j] - streamlinesX[j-1]) ** 2) + ((streamlinesY[j-1] - streamlinesY[j]) ** 2))
                streamTubeWidths.append(width)

                self.streamlinesXCoords[j-1].append(streamlinesX[j])
                self.streamlinesYCoords[j-1].append(streamlinesY[j])

            i+=1
