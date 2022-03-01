from impeller import Impeller

class Volute:

    def __init__(self, impeller: Impeller, type: str = "single") -> None:

        self.type = type.lower()
        self.impeller = impeller


        self.casingDesignFlowRate = self.impeller.impellerVolumeFlowRate
        self.QLe = self.casingDesignFlowRate # Alias
        
