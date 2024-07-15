
# Virus Detection in CHO Cells Using Oxford Nanopore Sequencing Reads, BLAST, and the Hypergeometric Probability Distribution

## Project Description

This repository implements the research paper by Roushangar et al. that covers the work for simulating virus detection using Oxford Nanopore sequencing and perform sequence alignments using BLAST searches.


## Directory Structure

- `src`: Contains source code files and a detailed `README.md` describing each script and its functions.
- `FeLV_24Threads_YES10%ReadPercentError_10power2VirusHostRatio_Simulation3`: Simulation data directory.
- `MVM_24Threads_YES10%ReadPercentError_10power2VirusHostRatio_Simulation3`: Simulation data directory.
- `PCV1_24Threads_YES10%ReadPercentError_10power2VirusHostRatio_Simulation3`: Simulation data directory.
- `ons_simulation_main.py`: Main script to run the simulation.
- `submit.sh`: Script for submitting jobs.
- `requirements.txt`: List of required Python packages.

## System Requirements

- macOS 10.15 or higher / Windows 10 or higher / Linux
- Python 3.10.6
- BLAST+ 2.13.0
- *May work with similar versions.

## Dependencies

- pandas >= 1.5.1
- biopython >= 1.79

## Installation Instructions

1. **Clone the repository:**
   ````bash
   git clone https://github.com/raeufroushangar/ONS_virus_detection_simulation.git
   cd ONS_virus_detection_simulation
   ```

2. **Create a virtual environment inside the `ONS_virus_detection_simulation` directory:**
   ```bash
   python3 -m venv venv

3. **Activate the virtual environment:**

   - On macOS and Linux:
     ```bash
     source venv/bin/activate

   - On Windows:
    ```bash
     .\venv\Scripts\activate

4. **Install required packages:**
   ```bash
   pip install -r requirements.txt

5. **Run the analysis script:**
   ```bash
   python3 ons_simulation_main.py MVM 'Minute virus of mice.fasta' 50 'GCF_000223135.1_CriGri_1.0_genomic.fna' 1 '/path/to/refSeq' 10power2 50 500 yes 10 24 3

   ### Explanation of Inputs:
   - `MVM`: Virus symbol (e.g., MVM for Minute Virus of Mice)
   - `'Minute virus of mice.fasta'`: Virus reference genome file name
   - `50`: Number of virus genomes to use
   - `'GCF_000223135.1_CriGri_1.0_genomic.fna'`: Host reference genome file name
   - `1`: Number of host genomes to use
   - `'/path/to/refSeq'`: Reference genome directory path
   - `10power2`: Virus to host ratio string to add to the output file name
   - `50`: Minimum read length
   - `500`: Maximum read length
   - `yes`: Error status (Yes/No)
   - `10`: Sequencing read error percentage
   - `24`: Number of logical threads to use
   - `3`: Simulation number

   Replace the paths and parameters with those suitable for your environment and data.
