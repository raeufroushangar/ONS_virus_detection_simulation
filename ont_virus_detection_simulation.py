import os
import csv
import sys
import time
import random
import pandas as pd
import concurrent.futures
from datetime import date
import multiprocessing as mp
from operator import itemgetter
from ont_virus_detection_blast_search import blastsearch as bs

# Global
reads= []
COUNT = 0

class simulation:
    """
       virusSymbol: virus symbol (e.g. FeLV, MVM, PCV)
       virusRefGenName: virus fasta file name  
       virusRefGenCount: number of virus genome to use
       hostRefGenName: host fasta file name
       hostRefGenCount: number of host genome to use
       refGenomeDirectory: reference genome file path
       virusHostRatio: string to add to output file name
       minReadLen: smallest read fragmente size
       maxReadLen: largest read fragmente size
       errorStatus: Yes/No
       percentError: sequencing read error percentage, 10% and higher
       threadsNumber: number of logical threads to use 
       simulationNumber: simulation number
    """
    #-------------------------
    def __init__(self, vSymbol, vGenName, vGenCount, hGenName, hGenCount, refGenDir, vhRatio, rMin,
                 rMax, eStatus, pError, tNum, sNum):

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

        # number of physical cores * 2 threads per core = number of logical cores
        print(f"\nNumber of CPU logical cores (threads) in the system: {mp.cpu_count()}") 
        print(f"Number of threads used: {self.threadsNumber}") 
    #-------------------------
    def readRefGenFileFunc(self, name):
        """
        -Input: 
            1: name: refernce genome file name, type: fasta
        -Functionality: open and read reference genome file
        -Output:
            1: single string, type: string
        """
        path = os.path.join(self.refGenomeDirectory, name)
        genome= ''
        with open(path, 'r') as f:
            for line in f:
                # ignore header line with genome information
                if not line[0] == '>':
                    genome += line.rstrip()
        return genome
    #-------------------------
    def addingErrorsFunc(self, seqFragment):
        """
        -Input: 
            1: seqFragment: sequence fragment, type: string
            2: percentError: error number, type: integer
        -Functionality: multiple user input error by length of read fragment to generate number of 
                        errors (errorNum) to use. Generate a list of individual nucleotides from
                        fragment. Loop through the nucleotides and add errors by replacing them 
                        using their index.
        -Output: 
            1: seqFragmentNew: same input genome fragment but with errors, type: string
        """
        errorNum = round(len(seqFragment)*(self.percentError/100))
        alphabetNoA = ['T', 'C', 'G']
        alphabetNoT = ['A', 'C', 'G']
        alphabetNoC = ['A', 'G', 'T']
        alphabetNoG = ['A', 'C', 'T']

        fragments= list (seqFragment)
        # generate unique (without replacement) random numbers with sample size equal to errorNum 
        # (e.g. 10% of the length of the input fragment). The generated random numbers are within 
        # the range 0 and the length of the input fragment.
        cps = random.sample(range(0, len(seqFragment)), errorNum)
        for i  in cps:
            randomACGT= random.randint(0, 3)
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
        seqFragmentNew= "".join(fragments)
        return seqFragmentNew
    #-------------------------
    def fragmentsGeneratorFunc(self, genome, readLenMin, readLenMax, errorRate, errorStatus):
        """ 
        -Input: 
            1: genome: refernce genome string, type: string
            2: readLenMin: smallest read length size, type: integer
            3: readLenMax: largest read length size, type: integer
            4: errorStatus: fragment with or without error, string
        -Functionality: generate read fragments with/without errors from a refernce genome
        -Output: 
            1: list of read fragments type: list
        """
        global reads
        readLen = random.randint(readLenMin, readLenMax)
        startPos = random.randint(0,len(genome)-readLen)
        read1, read2, read3 = (genome[:startPos],
                            genome[startPos:startPos+readLen],
                            genome[startPos+readLen:])
        if len(read2) !=0 and len(read2) >=readLenMin:
            if errorStatus == 'YES':
                read2WithError = self.addingErrorsFunc(read2)
                reads.append(read2WithError)
            else:
                reads.append(read2)
        for read in [read1, read3]:
            if len(read) !=0 and len(read) >=readLenMin:
                if len(read) > readLenMax:
                    self.fragmentsGeneratorFunc(read, readLenMin, readLenMax, errorRate, errorStatus)
                else:
                    if errorStatus == 'YES':
                        readWithError = self.addingErrorsFunc(read)
                        reads.append(readWithError)
                    else:
                        reads.append(read)
        return reads
    #-------------------------
    def processingFunc(self):
        """ 
        -Input: 
            1: virusRefGenName: virus refernce genome file name, type: string
            2: hostRefGenName: host refernce geneome file name, type: string
            3: hostRefGenCount: number of host genome to use, type: integer
            4: minReadLen: size of smallest read fragmente, type: integer
            5: maxReadLen: size of largets read fragmente, type: integer
            6: percentError: sequencing read error percentage, type: integer
            7: errorStatus: error status (Yes/No), type: string
        -Functionality: call 'reading refernce geneome' and 'read fragmentation' functions
        -Output: 
            1: genemoe size type: string
        """
        # read target refrence genome file (fasta)
        vReadGenome = self.readRefGenFileFunc(self.virusRefGenName) 
        vGenomeSize = len(vReadGenome)
        
        # read host refrence genome file (fasta)
        hReadGenome = self.readRefGenFileFunc(self.hostRefGenName) 
        hGenomeSize = len(hReadGenome)

        if self.maxReadLen > vGenomeSize:
            print(f"Can not enter desired read length (max) > {vGenomeSize} bases. Try again!")
        else:
            hGenomes = self.hostRefGenCount * [hReadGenome]        
            [self.fragmentsGeneratorFunc (i, self.minReadLen, self.maxReadLen, self.percentError, 
            self.errorStatus) for i in hGenomes]

            vGenomes = self.virusRefGenCount * [vReadGenome]        
            [self.fragmentsGeneratorFunc(j, self.minReadLen, self.maxReadLen, self.percentError, 
            self.errorStatus) for j in vGenomes]

            return (vGenomeSize, hGenomeSize)
    #-------------------------
    def readsPickerThreadingFunc(self):
        """
        -Input: 
            1: threadsNumber: number of  threads to use, type: integer
        -Functionality: use number threads to randomly select reads from a list of reads. This to 
                        mimic Oxford Nanopore Sequenicng process
        -Output: 
            1: list of sequencing reads, type: list
        """
        try:
            if self.threadsNumber == 1:
                randomIndexes= random.randrange(0, len(reads))
                result= [reads[randomIndexes]]
                del reads[randomIndexes]
                return result

            if self.threadsNumber > 1:
                if len(reads) >= self.threadsNumber:
                    randomIndexes= random.sample(range(len(reads)), self.threadsNumber)
                    result= list(itemgetter(*randomIndexes)(reads))
                    for index in sorted(randomIndexes, reverse=True):
                        del reads[index]
                    return result

                if len(reads) < self.threadsNumber and len(reads) > 1:
                    randomIndexes= random.sample(range(len(reads)), len(reads))
                    result= list(itemgetter(*randomIndexes)(reads))
                    for index in sorted(randomIndexes, reverse=True):
                        del reads[index]
                    return result

                if len(reads) == 1:
                    result = [reads[0]]
                    del reads[0]
                    return result
        except:
            pass
    #-------------------------
    @staticmethod
    def generateFragmentsFilesFunc(fragments, fileName):
        """
        -Input:
             1: fragments: genome fragments, type: list
             2: fileName: name of file to write fragment size, type: string
        -Functionality: Loop through a list of fragments, create file, and write the fragment 
                        order number and sie the size of each 
        -Output:
              1: two column file with fragment size and number, type: csv file
        """
        global COUNT
        counter=0
        try:
            for fragment in fragments:
                COUNT+=1
                counter+=1
                readFileName= 'fragmentsFile_'+str(counter)+'_Threading.fa'
                readFile = open(readFileName,"w") 
                readFile.write(f">read{str(COUNT)}\n") 
                readFile.write(fragment+"\n\n") 
                readFile.close()
                readCountLen = [[COUNT,len(fragment)]]
                my_df = pd.DataFrame(readCountLen)
                my_df.to_csv(fileName+'.csv', mode='a', index=False, header=False)
        except:
            pass
    #-------------------------
    def sequenceAlignmentFunc(self):
        """ write simulation metadata and matched read from blast seach into file"""
        start_time= time.time()

        global COUNT
        processingFuncResult = self.processingFunc()

        resultDirectory= f"{self.virusSymbol}_{str(self.threadsNumber)}Threads_{self.errorStatus}{self.percentError}%ReadPercentError_{self.virusHostRatio}VirusHostRatio_Simulation{self.simulationNumber}" # result directory

        os.mkdir(resultDirectory) # create directory
        os.chdir(resultDirectory) # set directory for result

        fragmentsFileName= f"{self.virusSymbol}_FragmentsCountLength_{str(self.threadsNumber)}Threads_{self.virusHostRatio}VirusHostRatio_Simulation{self.simulationNumber}"

        fields_csv= ['alignment_num','sbjct_title', 'sbjct_id', 'sbjct_accession', 'sbjct_length',
                     'query_id', 'query_length', 'score', 'evalue', 'identities', 'positives',
                      'gaps', 'align_length', 'strand', 'query_start', 'query_end', 'query', 
                      'match', 'sbjct', 'sbjct_start', 'sbjct_end', 'detect_time_min']

        finalResultFile_csv= fragmentsFileName+'_RVDB_result.csv'
        with open(finalResultFile_csv, 'a') as csvfile:  
            csvwriter = csv.writer(csvfile)  
            csvwriter.writerow(fields_csv)  

        today = date.today().strftime("%b-%d-%Y") # mm/dd/y
        finalResultFile_txt= open(fragmentsFileName+'_RVDB_result.txt', 'w')
        finalResultFile_txt.write(f'\nDate: {today}')
        finalResultFile_txt.write(f'\nVirus: {self.virusRefGenName}\n')
        finalResultFile_txt.write(f'Genome size (bp): {processingFuncResult[0]}\n')
        finalResultFile_txt.write(f'Count: {self.virusRefGenCount}\n')
        finalResultFile_txt.write(f'Host name: {self.hostRefGenName}\n')
        finalResultFile_txt.write(f'Genome size (bp): {processingFuncResult[1]}\n')
        finalResultFile_txt.write(f'Count: {self.hostRefGenCount}\n')
        finalResultFile_txt.write(f'\nVirus to host ratio: {self.virusHostRatio}')
        finalResultFile_txt.write(f'\nNumber of total fragments: {len(reads)}')
        finalResultFile_txt.write(f'\n{self.percentError}% Error per read: {self.errorStatus}\n\n')

        for _ in range(len(reads)):
            COUNT+=0
            if reads:
                pickedReads= self.readsPickerThreadingFunc()
                self.generateFragmentsFilesFunc(pickedReads, fragmentsFileName)
                readFilesName= [fname for fname in os.listdir('.') if '_Threading.fa' in fname]
                with concurrent.futures.ThreadPoolExecutor() as executor:
                    blastn_RVDB_results= executor.map(bs.rvdbdbBlastSearchFunc, readFilesName)
                    for result in blastn_RVDB_results:
                        if result:
                            alignment_num= 0
                            for RVDB_result in result:
                                alignment_num+=1
                                bs.searchResultWriterFunc(finalResultFile_txt,
                                                          RVDB_result,
                                                          alignment_num)
                                with open(finalResultFile_csv, 'a') as csvfile:  
                                    csvwriter = csv.writer(csvfile)  
                                    detection_time = (time.time() - start_time)/60
                                    RVDB_result_alignment_num_detection_time= [alignment_num]+RVDB_result+[detection_time]
                                    csvwriter.writerow(RVDB_result_alignment_num_detection_time)  
                            break
                    else:
                        time.sleep(0.2)
                        continue
                    break
        end_time= time.time() - start_time
        print(f"Processing took {round(end_time/60,2)} min using {str(self.threadsNumber)} parallel threads.\n")
        finalResultFile_txt.write(f"Processing took {round(end_time/60,2)} min using {str(self.threadsNumber)} parallel threads.")
        finalResultFile_txt.close()


#-------------------------
#-------------------------
if __name__ == "__main__":

    simul = simulation(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], 
                       sys.argv[6], sys.argv[7], sys.argv[8], sys.argv[9], sys.argv[10], 
                       sys.argv[11], sys.argv[12], sys.argv[13])
    simul.sequenceAlignmentFunc()
   
    # command line arguments:
    # python3 -B ont_virus_detection_simulation.py MVM 'Minute virus of mice.fasta' 50 'GCF_000223135.1_CriGri_1.0_genomic.fna' 1 '/home/raeuf/raeuf_notebook/bioinformatics_tools/apps/blast/db/RefSeq' 10power2 50 500 yes 10 24 3

    # python3 -B ont_virus_detection_simulation.py PCV1 'Porcine circovirus 1.fasta' 50 'GCF_000223135.1_CriGri_1.0_genomic.fna' 1 '/home/raeuf/raeuf_notebook/bioinformatics_tools/apps/blast/db/RefSeq' 10power2 50 500 yes 10 24 3