import os
import multiprocessing as mp

class SimulationSetup:
    def __init__(self, vSymbol, vGenName, vGenCount, hGenName, hGenCount, refGenDir, vhRatio, rMin, rMax, eStatus, pError, tNum, sNum):
        self.virusSymbol = vSymbol 
        self.virusRefGenName = vGenName 
        self.virusRefGenCount = int(vGenCount)
        self.hostRefGenName = hGenName 
        self.hostRefGenCount = int(hGenCount) 
        self.refGenomeDirectory = refGenDir 
        self.virusHostRatio = vhRatio 
        self.minReadLen = int(rMin) 
        self.maxReadLen = int(rMax) 
        self.errorStatus = eStatus.upper()
        self.percentError= int(pError) 
        self.threadsNumber = int(tNum) 
        self.simulationNumber = sNum 

        print(f"\nNumber of CPU logical cores (threads) in the system: {mp.cpu_count()}") 
        print(f"Number of threads used: {self.threadsNumber}")

    def readRefGenFileFunc(self, name):
        path = os.path.join(self.refGenomeDirectory, name)
        genome = ''
        with open(path, 'r') as f:
            for line in f:
                if not line[0] == '>':
                    genome += line.rstrip()
        return genome
