function networks = networkLearning(lambda,eta)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Code to reproduce analysis in
% Hill, Nesser et al. Cell Systems 4, 73-83 (2017), DOI: 10.1016/j.cels.2016.11.013.
%
% Performs network learning - see Figure 5 and Table S3.
% For details see STAR Methods subsection titled 'Network Learning' and references therein.
% Networks are learned using a modified version of the Joint Network Inference approach of Oates et al. Annals of Applied Statistics 8, 1892–1919 (2014). See function JNI_serial_mod.
% The functions readDataCore and readDataComplete need to be run first and the generated data files must be in the data directory.
%
% networks = networkLearning(lambda,eta)
%
% Input:
% lambda - individual regularisation hyperparameter in Joint Network Inference algorithm
% eta - latent network regularisation hyperparameter in Joint Network Inference algorithm
%
% Outputs:
% networks - 4-dimensional array of size 35x35x8x4; networks(i,j,s,c) is the posterior edge probability for
%               the edge from proteinNames{i} to proteinNames{j} in the network for cellLine{c}, stimulus{s}
%               where the variables proteinNames, cellLine and stimulus can be found in the data files generated 
%               by the function readDataCore.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cellLine = {'UACC812','MCF7','BT20','BT549'}; 
nCellLine = length(cellLine);

% load data
for c=1:nCellLine
    load(['../data/',cellLine{c},'_log2_core'])
    dataCore = data;
    
    % use additional time points in complete data for the interpolation in the preprocessing step
    load(['../data/',cellLine{c},'_log2_complete'],'data','timeNum','antibodyUsedInd') 
    dataComplete = data;
    dataComplete = dataComplete(antibodyUsedInd==1,:,:,:,:);
    dataCoreExtended = dataCore;
    dataCoreExtended(:,end+1:length(timeNum),:,:,:) = dataComplete(:,size(dataCore,2)+1:end,:,:,:);
    dataAll{c} = preprocessData(dataCoreExtended,timeNum); % average replicates and interpolate missing data
    
    % remove additional time points
    dataAll{c} = dataAll{c}(:,1:size(dataCore,2),:,:);
end

[nProt,nTimepoint,nStim,nInhib,~] = size(dataAll{1});

% load prior network
priorGraph = loadPrior;

% create matrix giving targets of inhibitors
inhibitorAction = zeros(nInhib,nProt);
inhibitorAction(1,[]) = 1;       % DMSO
inhibitorAction(2,19) = 1;       % BEZ235 (PI3K/mTOR - mTOR_pS2448)
inhibitorAction(3,[]) = 1;       % PD173074 (FGRFi)
inhibitorAction(4,19) = 1;  % AZD8055 (mTOR - mTOR_pS2448)
inhibitorAction(5,[3 4]) = 1;       % GSK690693 (pan-AKT - AKTpS473, AKTpT308)
inhibitorAction(6,[3 4 18]) = 1;  % GSK690693_GSK112021 (pan-AKT+MEKi - AKTpS473, AKTpT308, MEK1_pS217_S221)

% set parameters for Joint Network Inference procedure and put data into correct format - see JNI_serial function for details
dmax = 3;
G0 = priorGraph;
var_sel = 1;
zeroIntercept = 0;
multiplicityPrior = 1;
selfEdges = 1;
fixedEffectRegime = 2;

% put data and inhibitor target information into correct format
D = cell(nStim*nCellLine,1);
I = cell(nStim*nCellLine,1);
for c=1:nCellLine
    for s=1:nStim
        contextCount = (c-1)*nStim+s;
        includedInhibs = [];
        for i=1:nInhib
            if nnz(isnan(dataAll{c}(:,:,s,i)))>0
                % skips any time courses with missing data
                continue
            else
                D{contextCount}{i} = dataAll{c}(:,:,s,i);
                includedInhibs(end+1) = i;
            end
        end
        I{contextCount} = logical(inhibitorAction(includedInhibs,:))';
    end
end

[~,nets] = JNI_serial_mod(D,I,G0,lambda,eta,dmax,var_sel,zeroIntercept,multiplicityPrior,selfEdges,fixedEffectRegime);

networks = cell2mat(nets');
networks = reshape(networks,[35 35 8 4]);

end

%%
function data = preprocessData(data,time)
% Linearly interpolates missing datapoints and averages replicates. Can handle up to two adjacent missing timepoints
%
% Inputs:
% data - in format proteins x timepoint x stimulus x inhibitor x replicates
% time - timepoints. If not input, assumed to be equally spaced
%
% Output:
% data - processed data

data = nanmean(data,5); % average replicates

nTimepoint = size(data,2);
if nargin<2 || isempty(time)
    time = 1:nTimepoint;
end

% linearly interpolate missing values if in middle of time course - can handle up to two consecutive missing points, otherwise leave as NaN
I = find(isnan(data));
[I J K] = ind2sub(size(data),I);
for i=1:length(I)
    if J(i)>1 && J(i)<nTimepoint
        if ~isnan(data(I(i),J(i)-1,K(i))) && ~isnan(data(I(i),J(i)+1,K(i)))
            td1 = time(J(i))-time(J(i)-1);
            td2 = time(J(i)+1)-time(J(i));
            frac1 = td1/(td1+td2);
            frac2 = td2/(td1+td2);
            data(I(i),J(i),K(i)) = frac2*data(I(i),J(i)-1,K(i))+frac1*data(I(i),J(i)+1,K(i));
        elseif isnan(data(I(i),J(i)-1,K(i))) && ~isnan(data(I(i),J(i)+1,K(i))) && J(i)>2 && ~isnan((data(I(i),J(i)-2,K(i))))
            td1 = time(J(i))-time(J(i)-2);
            td2 = time(J(i)+1)-time(J(i));
            frac1 = td1/(td1+td2);
            frac2 = td2/(td1+td2);
            data(I(i),J(i),K(i)) = frac2*data(I(i),J(i)-2,K(i))+frac1*data(I(i),J(i)+1,K(i));
            td1 = time(J(i)-1)-time(J(i)-2);
            td2 = time(J(i))-time(J(i)-1);
            frac1 = td1/(td1+td2);
            frac2 = td2/(td1+td2);
            data(I(i),J(i)-1,K(i)) = frac2*data(I(i),J(i)-2,K(i))+frac1*data(I(i),J(i),K(i));
        end
    end
end

end

%%
function priorGraph = loadPrior
% defines the prior network
priorGraph = zeros(35);

EBP1 = 1;
ACC = 2;
AKTpS=3; AKTpT=4; AKT = [AKTpS AKTpT];
AMPK = 5;
BAD = 6;
cMET = 7;
cRAF = 8;
CHK1 = 9;
CHK2 = 10;
EGFRpY1068 = 11; EGFRpY1173 = 12; EGFR = [EGFRpY1068 EGFRpY1173];
ERa = 13;
GSK3 = 14;
HER2 = 15;
JNK = 16;
MAPK = 17;
MEK = 18;
mTOR = 19;
NFkB = 20;
p27pT157 = 21; p27pT198 = 22; p27 = [p27pT157 p27pT198];
p38 = 23;
p70 = 24;
p90 = 25;
PDK1 = 26;
PKCa = 27;
PRAS40 = 28;
RB = 29;
S6 = 30;
SRCpY416 = 31; SRCpY527 = 32; SRC = [SRCpY416 SRCpY527];
STAT3 = 33;
YAP = 34;
YB1 = 35;

priorGraph([EGFR HER2 cMET],PKCa)=1;
priorGraph([EGFR HER2 cMET],PDK1)=1;
priorGraph([EGFR HER2 cMET],cRAF)=1;
priorGraph([EGFR HER2 cMET],SRC)=1;
priorGraph([EGFR HER2 cMET],p38)=1;
priorGraph([EGFR HER2 cMET],JNK)=1;
priorGraph(AKT,PRAS40)=1;
priorGraph(AKT,mTOR)=1;
priorGraph(AKT,GSK3)=1;
priorGraph(AKT,BAD)=1;
priorGraph(AKT,YAP)=1;
priorGraph(AMPK,ACC)=1;
priorGraph(AMPK,mTOR)=1;
priorGraph(cRAF,MEK)=1;
priorGraph(GSK3,mTOR)=1;
priorGraph(MAPK,p90)=1;
priorGraph(MAPK,cRAF)=1;
priorGraph(MAPK,MEK)=1;
priorGraph(MAPK,mTOR)=1;
priorGraph(MEK,MAPK)=1;
priorGraph(mTOR,p70)=1;
priorGraph(mTOR,EBP1)=1;
priorGraph(mTOR,AKT)=1;
priorGraph(p70,S6)=1;
priorGraph(p70,PDK1)=1;
priorGraph(p70,mTOR)=1;
priorGraph(p90,S6)=1;
priorGraph(p90,mTOR)=1;
priorGraph(p90,GSK3)=1;
priorGraph(p90,BAD)=1;
priorGraph(PDK1,AKT)=1;
priorGraph(PDK1,PKCa)=1;
priorGraph(PDK1,p70)=1;
priorGraph(PRAS40,mTOR)=1;
priorGraph(SRC,STAT3)=1;

end
