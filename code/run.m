%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Code to reproduce analysis in
% Hill, Nesser et al. Cell Systems 4, 73-83.e10 (2017), DOI: 10.1016/j.cels.2016.11.013.
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
% The additional time points are used within the function causalDescendancyMatrices.m.
readDataComplete

%% Causal descendancy matrices
% Generate causal descendancy matrices (CDMs) for each inhibitor regime - see Figures 4A and S3. Also related to Figure S2.
% For details see STAR Methods subsection titled 'Identification of Changes Under Kinase Inhibition'.

[CDMs,contextLabels,pvals_fdr,effectSizeRatios] = causalDescendancyMatrices;

%% Network learning
% Performs network learning - see Figure 5 and Table S3.
% For details see STAR Methods subsection titled 'Network Learning' and references therein.
% Networks are learned using a modified version of the Joint Network Inference approach of Oates et al. Annals of Applied Statistics 8, 1892–1919 (2014). See function JNI_serial_mod.

cellLine = {'UACC812','MCF7','BT20','BT549'};
lambda = 3;
eta = 15;
networks = networkLearning(lambda,eta);
