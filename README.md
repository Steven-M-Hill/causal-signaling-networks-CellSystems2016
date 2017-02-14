### Code to reproduce the causal signaling network analysis in [Hill, Nesser et al. Context Specificity in Causal Signaling Networks Revealed by Phosphoprotein Profiling, Cell Systems 4, 73-83.e10 (2017)](http://www.cell.com/cell-systems/abstract/S2405-4712(16)30408-2)

### Overview

Matlab functions to read in the data and perform the following analysis (see paper for details):
- identification of changes under kinase inhibition
- network learning
- assessment of network learning performance (code to follow)

A script `run.m` is provided that calls each of the functions.

Before running the script, the supplemental zip file `DataS1.zip` (available to download via the link to the paper above) must be in the present working directory (together with the script and the other functions).

The following functions are provided (see header comments in individual functions for more details of function inputs and outputs):

- `readDataCore.m` <br />
Read in data (CSV format) from supplemental zip file DataS1 and save as MAT files. Reads in "core" data - this is the data focussed on in the analyses in the paper.
- `readDataComplete.m` <br />
Read in data (CSV format) from supplemental zip file DataS1 and save as MAT files. Reads in "complete" data - this is all data, including additional later time points and antibodies not focussed on in the analyses in the paper. The additional time points are used within the function `causalDescendancyMatrices.m`and in a preprocessing step in the function `networkLearning.m`.
- `causalDescendancyMatrices.m` <br />
Generates causal descendancy matrices (CDMs) for each inhibitor regime - see Figures 4A and S3 in paper. Also related to Figure S2.
For details see STAR Methods subsection titled 'Identification of Changes Under Kinase Inhibition'.
- `networkLearning.m` <br />
Learns networks for each *(cell line, stimulus)* context - see Figure 5 and Table S3 in the paper.
For details see STAR Methods subsection titled 'Network Learning' and references therein.
Networks are learned using a modified version of the Joint Network Inference approach of [Oates et al. Annals of Applied Statistics 8, 1892–1919 (2014)](https://doi.org/10.1214/14-AOAS761). This modified version is provided as function `JNI_serial_mod`. The original version can be found at https://doi.org/10.1214/14-AOAS761SUPPB.
- `networkLearningAssessment.m` <br />
Code to follow.
