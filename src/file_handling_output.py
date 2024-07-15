import os
import time
import pandas as pd
import concurrent.futures
import csv
from datetime import date
from src.error_handling_fragment_generation import ErrorHandlingFragmentGeneration
from src.blast_search import BlastSearch

class FileHandlingOutput(ErrorHandlingFragmentGeneration):
    COUNT = 0

    @staticmethod
    def generateFragmentsFilesFunc(fragments, fileName):
        """
        Write genome fragments to files and record their lengths.

        Args:
            fragments (list): List of genome fragments.
            fileName (str): The name of the file to write fragment size information.

        Returns:
            None
        """
        counter = 0
        try:
            for fragment in fragments:
                FileHandlingOutput.COUNT += 1
                counter += 1
                readFileName = 'fragmentsFile_' + str(counter) + '_Threading.fa'
                with open(readFileName, "w") as readFile:
                    readFile.write(f">read{str(FileHandlingOutput.COUNT)}\n") 
                    readFile.write(fragment + "\n\n")
                readCountLen = [[FileHandlingOutput.COUNT, len(fragment)]]
                my_df = pd.DataFrame(readCountLen)
                my_df.to_csv(fileName + '.csv', mode='a', index=False, header=False)
        except Exception as e:
            print(f"Error generating fragment files: {e}")

    def sequenceAlignmentFunc(self, simulation):
        """
        Manage the sequence alignment process and write results to files.

        Args:
            simulation (SimulationSetup): The simulation setup object.

        Returns:
            None
        """
        start_time = time.time()
        processingFuncResult, reads = simulation.processingFunc()
        resultDirectory = f"{simulation.virusSymbol}_{str(simulation.threadsNumber)}Threads_{simulation.errorStatus}{simulation.percentError}%ReadPercentError_{simulation.virusHostRatio}VirusHostRatio_Simulation{simulation.simulationNumber}"
        
        # Create and set result directory
        os.mkdir(resultDirectory)
        os.chdir(resultDirectory)
        
        fragmentsFileName = f"{simulation.virusSymbol}_FragmentsCountLength_{str(simulation.threadsNumber)}Threads_{simulation.virusHostRatio}VirusHostRatio_Simulation{simulation.simulationNumber}"
        
        # Prepare CSV file for results
        fields_csv = ['alignment_num','sbjct_title', 'sbjct_id', 'sbjct_accession', 'sbjct_length', 'query_id', 'query_length', 'score', 'evalue', 'identities', 'positives', 'gaps', 'align_length', 'strand', 'query_start', 'query_end', 'query', 'match', 'sbjct', 'sbjct_start', 'sbjct_end', 'detect_time_min']
        finalResultFile_csv = fragmentsFileName + '_RVDB_result.csv'
        with open(finalResultFile_csv, 'a') as csvfile:  
            csvwriter = csv.writer(csvfile)  
            csvwriter.writerow(fields_csv)
        
        # Prepare TXT file for results
        today = date.today().strftime("%b-%d-%Y")
        finalResultFile_txt = open(fragmentsFileName + '_RVDB_result.txt', 'w')
        finalResultFile_txt.write(f'\nDate: {today}')
        finalResultFile_txt.write(f'\nVirus: {simulation.virusRefGenName}\n')
        finalResultFile_txt.write(f'Genome size (bp): {processingFuncResult[0]}\n')
        finalResultFile_txt.write(f'Count: {simulation.virusRefGenCount}\n')
        finalResultFile_txt.write(f'Host name: {simulation.hostRefGenName}\n')
        finalResultFile_txt.write(f'Genome size (bp): {processingFuncResult[1]}\n')
        finalResultFile_txt.write(f'Count: {simulation.hostRefGenCount}\n')
        finalResultFile_txt.write(f'\nVirus to host ratio: {simulation.virusHostRatio}')
        finalResultFile_txt.write(f'\nNumber of total fragments: {len(reads)}')
        finalResultFile_txt.write(f'\n{simulation.percentError}% Error per read: {simulation.errorStatus}\n\n')
        
        for _ in range(len(reads)):
            if reads:
                pickedReads = simulation.readsPickerThreadingFunc()
                self.generateFragmentsFilesFunc(pickedReads, fragmentsFileName)
                readFilesName = [fname for fname in os.listdir('.') if '_Threading.fa' in fname]
                
                with concurrent.futures.ThreadPoolExecutor() as executor:
                    blastn_RVDB_results = executor.map(BlastSearch.rvdbdbBlastSearchFunc, readFilesName)
                    for result in blastn_RVDB_results:
                        if result:
                            alignment_num = 0
                            for RVDB_result in result:
                                alignment_num += 1
                                BlastSearch.searchResultWriterFunc(finalResultFile_txt, RVDB_result, alignment_num)
                                with open(finalResultFile_csv, 'a') as csvfile:  
                                    csvwriter = csv.writer(csvfile)  
                                    detection_time = (time.time() - start_time) / 60
                                    RVDB_result_alignment_num_detection_time = [alignment_num] + RVDB_result + [detection_time]
                                    csvwriter.writerow(RVDB_result_alignment_num_detection_time)  
                            break
                    else:
                        time.sleep(0.2)
                        continue
                    break
        
        end_time = time.time() - start_time
        print(f"Processing took {round(end_time/60,2)} min using {str(simulation.threadsNumber)} parallel threads.\n")
        finalResultFile_txt.write(f"Processing took {round(end_time/60,2)} min using {str(simulation.threadsNumber)} parallel threads.")
        finalResultFile_txt.close()
