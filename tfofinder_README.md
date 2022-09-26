# tfofinder
The tfofinder program identifies purine-only double stranded regions without bulges on the 3' strand.
The input file is a ct file in the RNAstructure or mfold format.
The output is a csv file containing parallel TFO (triple helix forming oligonucleotides) for a given RNA target.
The program first converts the ct file into the sscount format, and a base is considered double stranded if it is predicted to be paired in at least 
one of the suboptimal structures listed in the ct file. Alternatively, the input file should contain only the MFE (minimum free energy structure) if 
it is the only predicted structure to be considered.
!!!!This is a beta version that it is currently undergoing review to remove several checkpoints.!!!!
