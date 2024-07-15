import sys
from src.simulation_setup import SimulationSetup
from src.error_handling_fragment_generation import ErrorHandlingFragmentGeneration
from src.processing_sequencing import ProcessingSequencing
from src.file_handling_output import FileHandlingOutput

if __name__ == "__main__":
    simulation = SimulationSetup(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7], sys.argv[8], sys.argv[9], sys.argv[10], sys.argv[11], sys.argv[12], sys.argv[13])
    output_handler = FileHandlingOutput()
    output_handler.sequenceAlignmentFunc(simulation)
