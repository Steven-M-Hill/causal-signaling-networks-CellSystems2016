%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Code to reproduce analysis in
% Hill, Nesser et al. Cell Systems 4, 73-83 (2017), DOI: 10.1016/j.cels.2016.11.013.
%
% See header comments in individual functions for more details of function inputs and outputs
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% read and save data in MAT format
% Read in data (CSV format) from supplemental zip file DataS1 (available at the DOI above) and save as MAT files.
% The file DataS1.zip needs to be in the present working directory.

% Read in "core" data - the data focussed on in the analyses in the manuscript.
excludeSample = 1;
readDataCore(excludeSample)
% Reads in "complete" data - all data, including additional later time points and antibodies not focussed on in the analyses in the manuscript.
readDataComplete

%% Causal descendancy matrices
% Generate causal descendancy matrices (CDMs) for each inhibitor regime - see Figures 4A and S3. Also related to Figure S2.
% For details see STAR Methods subsection titled 'Identification of Changes Under Kinase Inhibition'.

[CDMs,contextLabels,pvals_fdr,effectSizeRatios] = causalDescendancyMatrices;
