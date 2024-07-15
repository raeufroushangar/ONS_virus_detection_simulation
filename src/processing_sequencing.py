import random
from operator import itemgetter

class ProcessingSequencing:
    def processingFunc(self, virusRefGenName, hostRefGenName, hostRefGenCount, minReadLen, maxReadLen, percentError, errorStatus, refGenomeDirectory):
        """
        Process the reference genomes and generate read fragments.

        Args:
            virusRefGenName (str): Virus reference genome file name.
            hostRefGenName (str): Host reference genome file name.
            hostRefGenCount (int): Number of host genomes to use.
            minReadLen (int): Minimum read length.
            maxReadLen (int): Maximum read length.
            percentError (int): Percentage of sequencing read errors.
            errorStatus (str): Error status ('YES' or 'NO').
            refGenomeDirectory (str): Path to the reference genome directory.

        Returns:
            tuple: Virus genome size, host genome size, list of generated reads.
        """
        reads = []
        vReadGenome = self.readRefGenFileFunc(virusRefGenName) 
        vGenomeSize = len(vReadGenome)
        
        hReadGenome = self.readRefGenFileFunc(hostRefGenName) 
        hGenomeSize = len(hReadGenome)

        if maxReadLen > vGenomeSize:
            print(f"Can not enter desired read length (max) > {vGenomeSize} bases. Try again!")
        else:
            hGenomes = hostRefGenCount * [hReadGenome]        
            for i in hGenomes:
                reads.extend(self.fragmentsGeneratorFunc(i, minReadLen, maxReadLen, percentError, errorStatus))

            vGenomes = self.virusRefGenCount * [vReadGenome]        
            for j in vGenomes:
                reads.extend(self.fragmentsGeneratorFunc(j, minReadLen, maxReadLen, percentError, errorStatus))

        return (vGenomeSize, hGenomeSize, reads)

    def readsPickerThreadingFunc(self, reads, threadsNumber):
        """
        Randomly select reads from a list of reads to mimic the sequencing process.

        Args:
            reads (list): List of reads.
            threadsNumber (int): Number of threads to use.

        Returns:
            list: List of selected reads.
        """
        try:
            if threadsNumber == 1:
                randomIndex = random.randrange(0, len(reads))
                result = [reads[randomIndex]]
                del reads[randomIndex]
                return result

            if threadsNumber > 1:
                if len(reads) >= threadsNumber:
                    randomIndexes = random.sample(range(len(reads)), threadsNumber)
                    result = list(itemgetter(*randomIndexes)(reads))
                    for index in sorted(randomIndexes, reverse=True):
                        del reads[index]
                    return result

                if len(reads) < threadsNumber and len(reads) > 1:
                    randomIndexes = random.sample(range(len(reads)), len(reads))
                    result = list(itemgetter(*randomIndexes)(reads))
                    for index in sorted(randomIndexes, reverse=True):
                        del reads[index]
                    return result

                if len(reads) == 1:
                    result = [reads[0]]
                    del reads[0]
                    return result
        except Exception as e:
            print(f"Error in readsPickerThreadingFunc: {e}")
            return []
