function readDataComplete
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Code to reproduce analysis in
% Hill, Nesser et al. Cell Systems 4, 73-83 (2017), DOI: 10.1016/j.cels.2016.11.013.
%
% Read in data (CSV format) from supplemental zip file DataS1 (available at the DOI above) and save as MAT files.
% Reads in "complete" data - all data, including additional later time points and antibodies not focussed on in the analyses in the manuscript.
% The file DataS1.zip needs to be in the present working directory.
%
% Input variables:
% 
% Output:
% Function generates 4 MAT files with filenames <cellLine>_log2_complete.mat. Each file contains a 5-dimensional array with dimensions phosphoprotein x timepoint x stimulus x inhibitor x replicate. 
%   Phosphoprotein abundance values are on a log2 scale.
%   Datapoints for the first timepoint (t=0min) are identical across the stimuli, since there was no stimulus here (stimulus was added at t=0min).
%   Each MAT file also contains some metadata: names of proteins, inhibitors, stimuli and time points; number of replicates for each sample; binary vectors indicating which antibodies and samples were used in the analyses.
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% unzip DataS1 zip file    
unzip('DataS1.zip')

% for each cell line, read in data from CSV file and save as MAT file

cellLine = {'UACC812','MCF7','BT20','BT549'};
nProt = [160,150,183,160]; % number of antibodies for each cell line

for c=1:length(cellLine)
    
    % import data
    fileID = fopen(['DataS1/complete/',cellLine{c},'.csv']);
    textscan(fileID,'%s',6,'Delimiter',',');
    proteinNames = textscan(fileID,'%s',nProt(c),'Delimiter',',');
    proteinNames = proteinNames{1}';
    textscan(fileID,'%s',6,'Delimiter',',','HeaderLines',2);
    antibodyUsedInd = textscan(fileID,'%d',nProt(c),'Delimiter',','); % antibody used in analyses indicator
    antibodyUsedInd = antibodyUsedInd{1}'; %#ok<NASGU>
    formatSpec = ['%s %s %s %s', repmat(' %20.15f',1,nProt(c)+2)];
    data = textscan(fileID,formatSpec,'Delimiter',',','HeaderLines',3);
    fclose(fileID);
    
    inhibText = data{3}; % inhibitor metadata
    stimText = data{2}; % stimulus metadata
    timepointText = data{4}; % timepoint metadata
    sampleUsedInd = data{6}; % sample used in analyses indicator
    data = cell2mat(data(7:end)); % data points
    
    % convert into 5-dimensional array with dimensions phosphoprotein x timepoint x stimulus x inhibitor x replicate
    n = size(data,1); % number of samples
    
    if strcmp(cellLine{c},'BT549')
        timeStr = {'0min','5min','15min','30min','60min','2hr','4hr','12hr','24hr','48hr','72hr'};
        timeNum = [0 5 15 30 60 120 240 60*12 60*24 60*48 60*72]; %#ok<NASGU>
    else
        timeStr = {'0min','5min','15min','30min','60min','2hr','4hr','12hr','24hr','48hr'};
        timeNum = [0 5 15 30 60 120 240 60*12 60*24 60*48]; %#ok<NASGU>
    end
    stimulus = {'Serum','PBS','EGF','Insulin','FGF1','HGF','NRG1','IGF1'};
    inhibitor = {'DMSO','BEZ235','PD173074','AZD8055','GSK690693','GSK690693_GSK112021'};
    dimLabel = {'protein','time','stimulus','inhibitor','replicates'}; %#ok<NASGU>
    
    dataNew = NaN*ones(length(proteinNames),length(timeStr),length(stimulus),length(inhibitor),1); % initialise 5-dimensional array
    nReps = zeros(length(timeStr),length(stimulus),length(inhibitor)); % number of replicates
    sampleUsedInd2 = NaN*ones(length(timeStr),length(stimulus),length(inhibitor),32); % sample used in analyses indicator
    
    for i=1:n
                
        if ~isempty(strfind(inhibText{i},'DMSO'))
            inhibIdx = 1;
        elseif ~isempty(strfind(inhibText{i},'BEZ'))
            inhibIdx = 2;
        elseif ~isempty(strfind(inhibText{i},'PD17'))
            inhibIdx = 3;
        elseif ~isempty(strfind(inhibText{i},'AZD'))
            inhibIdx = 4;
        elseif ~isempty(strfind(inhibText{i},'GSK1'))
            inhibIdx = 6;
        elseif ~isempty(strfind(inhibText{i},'GSK6'))
            inhibIdx = 5;
        else
            error('no inhibitor match')
        end
        
        if ~isempty(strfind(stimText{i},'Serum'))
            stimIdx = 1;
        elseif ~isempty(strfind(stimText{i},'PBS'))
            stimIdx = 2;
        elseif ~isempty(strfind(stimText{i},'EGF'))
            stimIdx = 3;
        elseif ~isempty(strfind(stimText{i},'Insulin'))
            stimIdx = 4;
        elseif ~isempty(strfind(stimText{i},'FGF'))
            stimIdx = 5;
        elseif ~isempty(strfind(stimText{i},'HGF'))
            stimIdx = 6;
        elseif ~isempty(strfind(stimText{i},'NRG'))
            stimIdx = 7;
        elseif ~isempty(strfind(stimText{i},'IGF'))
            stimIdx = 8;
        else
            stimIdx = 0; % time 0 points - no stimulus
        end
        
        if ~isempty(strfind(timepointText{i},'15min'))
            timeIdx = 3;
        elseif ~isempty(strfind(timepointText{i},'30min'))
            timeIdx = 4;
        elseif ~isempty(strfind(timepointText{i},'60min'))
            timeIdx = 5;
        elseif ~isempty(strfind(timepointText{i},'12hr'))
            timeIdx = 8;
        elseif ~isempty(strfind(timepointText{i},'24hr'))
            timeIdx = 9;
        elseif ~isempty(strfind(timepointText{i},'48hr'))
            timeIdx = 10;
        elseif ~isempty(strfind(timepointText{i},'72hr'))
            timeIdx = 11;
        elseif ~isempty(strfind(timepointText{i},'0min'))
            timeIdx = 1;
        elseif ~isempty(strfind(timepointText{i},'5min'))
            timeIdx = 2;
        elseif ~isempty(strfind(timepointText{i},'2hr'))
            timeIdx = 6;
        elseif ~isempty(strfind(timepointText{i},'4hr'))
            timeIdx = 7;
        else
            error('no time match')
        end
       
        if stimIdx==0 % no stimulus, entries dataNew(:,:,x,:,:) are set to be the same for all x=1,..,8
            for j=1:length(stimulus)
                nReps(timeIdx,j,inhibIdx) = nReps(timeIdx,j,inhibIdx)+1;
                dataNew(:,timeIdx,j,inhibIdx,nReps(timeIdx,j,inhibIdx)) = data(i,:);
                sampleUsedInd2(timeIdx,j,inhibIdx,nReps(timeIdx,j,inhibIdx)) = sampleUsedInd(i);
            end
        else
            nReps(timeIdx,stimIdx,inhibIdx) = nReps(timeIdx,stimIdx,inhibIdx)+1;
            dataNew(:,timeIdx,stimIdx,inhibIdx,nReps(timeIdx,stimIdx,inhibIdx)) = data(i,:);
            sampleUsedInd2(timeIdx,stimIdx,inhibIdx,nReps(timeIdx,stimIdx,inhibIdx)) = sampleUsedInd(i);
        end
    end
    data = log2(dataNew);
    data(isinf(data)) = NaN; %#ok<NASGU>
    sampleUsedInd = sampleUsedInd2;
    sampleUsedInd(:,:,:,max(nReps(:))+1:end) = []; %#ok<NASGU>
    
    save([cellLine{c},'_log2_complete'],'data','proteinNames','timeStr','timeNum','stimulus','inhibitor','dimLabel','nReps','sampleUsedInd','antibodyUsedInd')
    
end