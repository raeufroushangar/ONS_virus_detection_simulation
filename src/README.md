# Source Code Directory

This directory contains the source code files for the virus detection simulation project. Below is a summary of the functions contained in each file.

## blast_search.py

- `rvdbdbBlastSearchFunc`: Perform a BLAST search against the RVDB database.
- `searchResultWriterFunc`: Write BLAST search results to a file.

## error_handling_fragment_generation.py

- `addingErrorsFunc`: Adds errors to a sequence fragment based on the specified error percentage.
- `fragmentsGeneratorFunc`: Generates read fragments with or without errors from a reference genome.

## file_handling_output.py

- `generateFragmentsFilesFunc`: Writes the fragments to files and records their lengths.
- `sequenceAlignmentFunc`: Manages the sequence alignment process and writes the results to files.

## processing_sequencing.py

- `processingFunc`: Processes the reference genomes and generates fragments.
- `readsPickerThreadingFunc`: Randomly selects reads from a list of reads to mimic the sequencing process.

## simulation_setup.py

- `__init__`: Initializes the simulation setup with the given parameters.
- `readRefGenFileFunc`: Reads a reference genome file and returns its contents as a string.
