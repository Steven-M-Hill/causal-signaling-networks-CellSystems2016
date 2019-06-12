### Code to reproduce the causal signaling network analysis in [Hill, Nesser et al. Context Specificity in Causal Signaling Networks Revealed by Phosphoprotein Profiling, Cell Systems 4, 73-83.e10 (2017)](http://www.cell.com/cell-systems/abstract/S2405-4712(16)30408-2)

### Overview

The [code](code) directory contains Matlab functions to read in the data and perform the following analysis (see paper for details):
- identification of changes under kinase inhibition
- network learning
- assessment of network learning performance (code to follow)

A script [`run.m`](code/run.m) is provided in the [code](code) directory that calls each of the functions. The output from running this script will be saved in the [data](data) and [results](results) directories. For convenience, this output is itself provided as part of this repository.

Before running the script, the supplemental zip file [`DataS1.zip`](data/DataS1.zip) must be in the data directory. This zip file is available to download via the link to the paper above, but for convenience it is already made available in the [data](data) directory.

The following functions are provided (see header comments in individual functions for more details of function inputs and outputs):

- [`readDataCore.m`](code/readDataCore.m) <br />
Read in data (CSV format) from supplemental zip file [DataS1](data/DataS1.zip) and save as MAT files in the [data](data) directory. Reads in "core" data - this is the data focussed on in the analyses in the paper.
- [`readDataComplete.m`](code/readDataComplete.m) <br />
Read in data (CSV format) from supplemental zip file [DataS1](data/DataS1.zip) and save as MAT files in the [data](data) directory. Reads in "complete" data - this is all data, including additional later time points and antibodies not focussed on in the analyses in the paper. The additional time points are used within the function [`causalDescendancyMatrices.m`](code/causalDescendancyMatrices.m) and in a preprocessing step in the function [`networkLearning.m`](code/networkLearning.m).
- [`causalDescendancyMatrices.m`](code/causalDescendancyMatrices.m) <br />
Generates causal descendancy matrices (CDMs) for each inhibitor regime - see Figures 4A and S3 in paper. Also related to Figure S2.
For details see STAR Methods subsection titled 'Identification of Changes Under Kinase Inhibition'.
CDMs are saved in the [results](results) directory in a MAT file and also as CSV files (one file for each inhibitor regime).
- [`networkLearning.m`](code/networkLearning.m) <br />
Learns networks for each *(cell line, stimulus)* context - see Figure 5 and Table S3 in the paper.
For details see STAR Methods subsection titled 'Network Learning' and references therein.
Networks are learned using a modified version of the Joint Network Inference approach of [Oates et al. Annals of Applied Statistics 8, 1892–1919 (2014)](https://doi.org/10.1214/14-AOAS761). This modified version is provided as function [`JNI_serial_mod.m`](code/JNI_serial_mod.m). The original version can be found at https://doi.org/10.1214/14-AOAS761SUPPB.
Network posterior edge probabilities are saved in the [results](results) directory in a MAT file and also as CSV files (one file for each context-specific network).
- [`networkLearningAssessment.m`](code/networkLearningAssessment.m) <br />
Code to follow.
