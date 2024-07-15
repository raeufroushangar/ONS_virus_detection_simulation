import random

class ErrorHandlingFragmentGeneration:
    def addingErrorsFunc(self, seqFragment, percentError):
        errorNum = round(len(seqFragment)*(percentError/100))
        alphabetNoA = ['T', 'C', 'G']
        alphabetNoT = ['A', 'C', 'G']
        alphabetNoC = ['A', 'G', 'T']
        alphabetNoG = ['A', 'C', 'T']

        fragments = list(seqFragment)
        cps = random.sample(range(0, len(seqFragment)), errorNum)
        for i in cps:
            randomACGT = random.randint(0, 3)
            if randomACGT == 3:
                fragments[i] = ''
            else:
                if fragments[i] == 'A':
                    fragments[i] = alphabetNoA[randomACGT]
                elif fragments[i] == 'T':
                    fragments[i] = alphabetNoT[randomACGT]
                elif fragments[i] == 'C':
                    fragments[i] = alphabetNoC[randomACGT]
                else:
                    fragments[i] = alphabetNoG[randomACGT]
        return "".join(fragments)

    def fragmentsGeneratorFunc(self, genome, readLenMin, readLenMax, percentError, errorStatus):
        reads = []
        readLen = random.randint(readLenMin, readLenMax)
        startPos = random.randint(0, len(genome) - readLen)
        read1, read2, read3 = genome[:startPos], genome[startPos:startPos+readLen], genome[startPos+readLen:]
        if len(read2) != 0 and len(read2) >= readLenMin:
            if errorStatus == 'YES':
                read2WithError = self.addingErrorsFunc(read2, percentError)
                reads.append(read2WithError)
            else:
                reads.append(read2)
        for read in [read1, read3]:
            if len(read) != 0 and len(read) >= readLenMin:
                if len(read) > readLenMax:
                    reads.extend(self.fragmentsGeneratorFunc(read, readLenMin, readLenMax, percentError, errorStatus))
                else:
                    if errorStatus == 'YES':
                        readWithError = self.addingErrorsFunc(read, percentError)
                        reads.append(readWithError)
                    else:
                        reads.append(read)
        return reads
