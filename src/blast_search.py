from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML

E_VALUE_THRESH = 1e-10

class BlastSearch:
    @staticmethod
    def rvdbdbBlastSearchFunc(tempFastaFile=''):
        """
        Perform a BLAST search against the RVDB database.

        Args:
            tempFastaFile (str): Path to the sequencing read file in fasta format.

        Returns:
            list: List of sequence alignment results.
        """
        BLAST_xml_file = tempFastaFile.split('_Threading.fa')[0] + '_BlastResult.xml'
        blastn_cline = NcbiblastnCommandline(query=tempFastaFile, 
                                             db="/home/raeuf/raeuf_notebook/bioinformatics_tools/apps/blast/db/RVDB_db/C-RVDBv24.1", 
                                             evalue=E_VALUE_THRESH,
                                             outfmt=5,
                                             num_descriptions=10,
                                             num_alignments=10,
                                             out=BLAST_xml_file,
                                             num_threads=1)
        blastn_cline()
        result_handle = open(BLAST_xml_file, 'r')
        blast_records = NCBIXML.parse(result_handle)
        blast_record = next(blast_records)
        resultList = []
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                if hsp.expect < E_VALUE_THRESH:
                    target_list = ['porcine', 'circovirus']  # Customize target list as needed
                    if any(target in alignment.hit_def.lower() for target in target_list):
                        result = [alignment.hit_def, alignment.hit_id, alignment.accession, 
                                  alignment.length, blast_record.query, blast_record.query_length,
                                  hsp.score, hsp.expect, hsp.identities, hsp.positives, hsp.gaps, 
                                  hsp.align_length, hsp.strand, hsp.query_start, hsp.query_end, 
                                  hsp.query, hsp.match, hsp.sbjct, hsp.sbjct_start, hsp.sbjct_end]
                        resultList.append(result)
        return resultList

    @staticmethod
    def searchResultWriterFunc(resultFile, blastResult, alignmentNum):
        """
        Write BLAST search results to a file.

        Args:
            resultFile (file object): The file to write the results to.
            blastResult (list): List of alignment attributes.
            alignmentNum (int): The alignment number.

        Returns:
            None
        """
        resultFile.write(f'\n****Alignment:{alignmentNum}****\n')
        resultFile.write(f'sbjct_title: {blastResult[0]}\n')
        resultFile.write(f'sbjct_id: {blastResult[1]}\n')
        resultFile.write(f'sbjct_accession: {blastResult[2]}\n')
        resultFile.write(f'sbjct_length: {blastResult[3]}\n')
        resultFile.write(f'query_id: {blastResult[4]}\n')
        resultFile.write(f'query_length: {blastResult[5]}\n')
        resultFile.write(f'score: {blastResult[6]}\n')
        resultFile.write(f'evalue: {blastResult[7]}\n')
        resultFile.write(f'identities: {blastResult[8]}\n')
        resultFile.write(f'positives: {blastResult[9]}\n')
        resultFile.write(f'gaps: {blastResult[10]}\n')
        resultFile.write(f'align_length: {blastResult[11]}\n')
        resultFile.write(f'strand: {blastResult[12]}\n')
        resultFile.write(f'query_start: {blastResult[13]}\n')
        resultFile.write(f'query_end: {blastResult[14]}\n')
        resultFile.write(f'query: {blastResult[15]}\n')
        resultFile.write(f'match: {blastResult[16]}\n')
        resultFile.write(f'sbjct: {blastResult[17]}\n')
        resultFile.write(f'sbjct_start: {blastResult[18]}\n')
        resultFile.write(f'sbjct_end: {blastResult[19]}\n')
        resultFile.write('*****************\n\n')
