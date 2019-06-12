function [CDMs,contextLabels,pvals_fdr,effectSizeRatios] = causalDescendancyMatrices(pval_threshold,effectSize_threshold)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Code to reproduce analysis in
% Hill, Nesser et al. Cell Systems 4, 73-83 (2017), DOI: 10.1016/j.cels.2016.11.013.
%
% Generates causal descendancy matrices (CDMs) for each inhibitor regime - see Figures 4A and S3. Also related to Figure S2.
% For details see STAR Methods subsection titled 'Identification of Changes Under Kinase Inhibition'.
% The functions readDataCore and readDataComplete need to be run first and the generated data files must be in the present working directory.
%
% [CDMs,contextLabels,pvals_fdr,effectSizeRatios] = causalDescendancyMatrices(pval_threshold,effectSize_threshold)%
% 
% Input:
% pval_threshold - threshold to apply to fdr-corrected p-values to create binary CDMs (defaults to value used in manuscript: 0.05)
% effectSize_threshold - threshold to apply to (effect size)/(replicate standard deviation) ratio (defaults to value used in manuscript: 1)
%
% Outputs:
% CDMs - a 3-dimensional array; proteins x (cell line,stimulus) contexts x inhibitor regime; CDMs(:,:,i) is the CDM for inhibitor regime inhibitor{i+1}, where inhibitor = {'DMSO','BEZ235','PD173074','AZD8055','GSK690693','GSK690693_GSK1120212'};.
% contextLabels - string giving labels for the CDM columns (contexts)
% pvals_fdr - a 4-dimensional array of signed fdr-corrected p-values - proteins x stimuli x cell lines x inhibitor regime.
%               Positive/negative indicates mean DMSO value is larger/smaller than mean inhibitor regime value.
% effectSizeRatios -  a 4-dimensional array of effect size ratios: proteins x stimuli x cell lines x inhibitor regime.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set defaults
if nargin<2
    effectSize_threshold = 1;
    if nargin<1
        pval_threshold = 0.05;
    end
end

cellLine = {'UACC812','MCF7','BT20','BT549'};  
stimulus = {'Serum','PBS','EGF','Insulin','FGF1','HGF','NRG1','IGF1'};
inhibitor = {'DMSO','BEZ235','PD173074','AZD8055','GSK690693','GSK690693_GSK1120212'};
    

%% Calculate paired t-test p-values
% For each phosphoprotein, inhibitor regime and (cell line, stimulus) context, perform a paired t test to assess whether 
%       mean phosphoprotein abundance under DMSO control is significantly different to mean abundance under the inhibitor regime
pvals1 = pvalsPairedT(0); % a 4-dimensional array of signed p-values - proteins x stimuli x cell lines x inhibitor regime
pvals2 = pvalsPairedT(1); % as previous line, but looks for "peak" shapes in time-courses and, when found, performs paired t-test within the peak region

[nProt,nStim,nCellLine,nInhib] = size(pvals1);

% if p-value within peak region is smaller than the one calculated across all 7 time points, use the smaller p-value
for i=1:nProt
    for j=1:nStim
        for k=1:nCellLine
            for l=1:nInhib
                [~,tmp] = min(abs([pvals1(i,j,k,l) pvals2(i,j,k,l)]));
                if tmp==2
                    pvals1(i,j,k,l) = pvals2(i,j,k,l);
                end
            end
        end
    end
end
pvals = pvals1;

% fdr correction
pvals_fdr = NaN*pvals;
for j=1:nStim
    for k=1:nCellLine
        for l=1:nInhib
            pvals_fdr(:,j,k,l) = fdr_medlsu(abs(pvals(:,j,k,l))); % m0 estimated using BKY two-stage method
        end
    end
end
pvals_fdr = min(pvals_fdr,1); % some fdr values can end up being >1; set these to 1
pvals_fdr = pvals_fdr.*sign(pvals); % give fdr values the same sign as the p-values

%% calculate ratio of effect size to replicate standard deviation

effectSizeRatios = getEffectSizeRatios;

%% generate CDMs

pvals_fdr2 = reshape(pvals_fdr,[nProt,nCellLine*nStim,nInhib]);
effectSizeRatios2 = reshape(effectSizeRatios,[nProt,nCellLine*nStim,nInhib]);
CDMs = abs(pvals_fdr2)<=pval_threshold & effectSizeRatios2>=effectSize_threshold;
CDMs = CDMs*1;
nanIdx = find(isnan(pvals_fdr2));
CDMs(nanIdx) = NaN;

contextLabels = {};
for i=1:length(cellLine)
    for j=1:length(stimulus)
        contextLabels{end+1} = [cellLine{i},',',stimulus{j}];
    end
end

end


%%
function pvals = pvalsPairedT(peakShape)
% For each phosphoprotein, inhibitor regime and (cell line, stimulus) context, perform a paired t test to assess whether 
%       mean phosphoprotein abundance under DMSO control is significantly different to mean abundance under the inhibitor regime
%
% inputs:
% peakShape - binary; if set to 1, looks for "peak" shapes in time courses and, if found, applies paired t-test within the peak region only
%
% outputs:
% pvals -  a 4-dimensional array of p-values from paired t-tests: proteins x stimuli x cell lines x inhibitor regime.
%          p-values are signed. Positive/negative indicates mean DMSO value is larger/smaller than mean inhibitor regime value.

cellLine = {'UACC812','MCF7','BT20','BT549'};  
for c=1:length(cellLine)
    
    load([cellLine{c},'_log2_core'])

    [nProt,nTimepoint,nStim,nInhib,~] = size(data);
    
    %average replicates
    data = nanmean(data,5);
    
    % initialise p-value matrix
    pvalsTmp = NaN*ones(nProt,nStim,nInhib-1);
    
    if peakShape
        % determine if there is a "peak" shape in each DMSO time course
        [peak peaking_region] = peakFinder(cellLine{c});
    end
    
    for p=1:nProt
        for d=2:nInhib
            for s=1:nStim
                if peakShape && peak(p,s)
                    t1 = peaking_region(p,s,1); t2 = min(peaking_region(p,s,2),nTimepoint);
                    if t2-t1+1<4
                        % if peak region is small (less than 4 time points), do not use peak region for t-test (use all 7 time points)
                        t1 = 1; t2 = nTimepoint;
                    elseif nnz(~isnan(data(p,t1:t2,s,1))&~isnan(data(p,t1:t2,s,d)))<4
                        % if there are less than 4 non-missing data points in peak region, do not use peak region for t-test
                        t1 = 1; t2 = nTimepoint;
                    end
                else
                    t1 = 1; t2 = nTimepoint;
                end
                
                dmso = squeeze(data(p,t1:t2,s,1)); inhib = squeeze(data(p,t1:t2,s,d));
                [~, pvalsTmp(p,s,d-1) , ~, testStat] = ttest(dmso,inhib);
                pvalsTmp(p,s,d-1) = sign(testStat.tstat)*(pvalsTmp(p,s,d-1)); % make pvalue +ve or -ve to indicate direction of effect
            end
        end
    end
    pvals(:,:,c,:) = pvalsTmp;
end
end

%%
function effectSizeRatios = getEffectSizeRatios
% For each phosphoprotein, inhibitor regime and (cell line, stimulus) context, calculate ratio of effect size (mean phosphoprotein abundance 
%   under DMSO control - mean abundance under the inhibitor regime) to replicate variation
%
% Output:
% effectSizeRatios -  a 4-dimensional array of effect size ratios: proteins x stimuli x cell lines x inhibitor regime.

cellLine = {'UACC812','MCF7','BT20','BT549'};

for c=1:length(cellLine)
    
    load([cellLine{c},'_log2_core'])
    
    [nProt,nTimepoint,nStim,nInhib,~] = size(data);
    
    % calculate replicate standard deviations for each time course
    for i=1:nProt
        for j=1:nStim
            for k=1:nInhib
                d = squeeze(data(i,:,j,k,:));
                s(i,j,k,:) = nanstd(d,[],2);
                n(i,j,k,:) = sum(~isnan(d),2);
            end
        end
    end
    
    % pooled standard deviations for (DMSO, inhib) pairs
    for i=1:nProt
        for j=1:nStim
            for k=2:nInhib
                for t=1:nTimepoint
                    if n(i,j,1,t)>0 && n(i,j,k,t)>0
                        ps(i,j,k-1,t) = sqrt( ( (n(i,j,1,t)-1).*s(i,j,1,t)^2 + (n(i,j,k,t)-1).*s(i,j,k,t)^2 )./(n(i,j,1,t)+n(i,j,k,t)-2) );
                    else
                        ps(i,j,k-1,t) = NaN;
                    end
                end
            end
        end
    end
    
    % mean pooled std deviations across time
    for i=1:nProt
        for j=1:nStim
            for k=2:nInhib
                mps(i,j,k-1) = nanmean(ps(i,j,k-1,:),4);
            end
        end
    end
    
    % effect sizes
    data2 = nanmean(data,5);
    for k=2:nInhib
        es(:,:,k-1) = squeeze(abs(nanmean(data2(:,:,:,1),2)-nanmean(data2(:,:,:,k),2)));
    end
    % if only t=0 available, set effect size to 0
    allNaNidx = find(squeeze(all(isnan(data2(:,2:end,:,2:end)),2))); 
    es(allNaNidx) = NaN;
    
    % effect size ratio
    effectSizeRatios(:,:,:,c) = es./mps;
end

effectSizeRatios = permute(effectSizeRatios,[1,2,4,3]);

end

%%
function [peak peaking_region] = peakFinder(cellLine)
% Heuristic approach to determine if "peak" shapes exist in DMSO time courses and finds 'peaking region' to constrain t-test to within this region.
% 
% inputs:
% cellLine - string specifying cell line being considered
%
% outputs:
% peak - binary matrix. peak(p,s)=1 if protein p, stimulus s has a peak in DMSO time course
% peaking_region - 3-dimensional array: if peak(p,s)==1, peaking_region(p,s,1:2) gives start and end points of peak for protein p, stimulus s. Set to 0 if peak(p,s)==0.


% load complete data set - later time points (> 4hrs) in this file are used within the heuristic approach to determine existence of peaks
load([cellLine,'_log2_complete'])
data = data(antibodyUsedInd==1,:,:,:,:); %#ok<NODEF> % keep only the 35 phosphoproteins contained in core data files

if ndims(data)==5
    data = nanmean(data,5); % average replicates
end
[P,T,S,D] = size(data);

% linearly interpolate missing values if in middle of a time course
I = find(isnan(data));
[I J K] = ind2sub(size(data),I);
for i=1:length(I)
    if J(i)>1 && J(i)<T
        if ~isnan(data(I(i),J(i)-1,K(i))) && ~isnan(data(I(i),J(i)+1,K(i)))
            % interpolate single missing value
            data(I(i),J(i),K(i)) = ( data(I(i),J(i)-1,K(i))+data(I(i),J(i)+1,K(i)) )/2;
        elseif ~isnan(data(I(i),J(i)-1,K(i))) && isnan(data(I(i),J(i)+1,K(i))) && J(i)+1~=T && ~isnan(data(I(i),J(i)+2,K(i)))
            % interpolate two consecutive missing values
            data(I(i),J(i),K(i)) = ( 2*data(I(i),J(i)-1,K(i))+data(I(i),J(i)+2,K(i)) )/3;
            data(I(i),J(i)+1,K(i)) = ( data(I(i),J(i)-1,K(i))+2*data(I(i),J(i)+2,K(i)) )/3;
        end
    end
end

% initialise outputs
peak = false(P,S,D);
peaking_region = zeros(P,S,2,D);

for p=1:P
    for d=1:D
        
        dat = squeeze(data(p,:,:,d)); % data for single protein and inhibitor regime
       
        % to simplify code below, fill in missing last timepoints with last available timepoint
        for s=1:S
            idx = T+1;
            while isnan(dat(idx-1,s))
                idx = idx-1;
            end
            if idx<=T
                dat(idx:end,s) = dat(idx-1,s);
            end
        end
        
        % metrics for capturing amplitude (height) of peaks
        amp1 = max(dat)-dat(1,:); % maximum amplitude of time course for each stimulus (relative to first timepoint)
        amp2 = max(dat)-dat(end,:); % maximum amplitude of time course for each stimulus (relative to last timepoint)
        amp3 = range(dat); % range of time course for each stimulus
         
        % filter out stims where time course amplitude (relative to first time point) is less than a third of the maximum amplitude seen across all stims,
        %   or where time course amplitude (relative to last time point) is less than a third of the maximum amplitude seen across all stims and last time point is below first time point,
        %   or where less than a third of the range of the time course is above the first time point  (i.e. majority of values are below the t=0 value).
        stimIdx1 = find(amp1>max(amp1)/3); % stims with at least a third of maximum amplitude (relative to first time point)
        stimIdx2 = find(amp2>max(amp2)/3 | dat(end,:)<=dat(1,:)); % stims with at least a third of maximum amplitude (relative to last time point) or last time point is below first time point
        stimIdx3 = find(amp1>amp3/3); % stims which have at least a third of range above first timepoint
        stimIdx = intersect(stimIdx1,stimIdx2); stimIdx = intersect(stimIdx,stimIdx3);
        
        for s=stimIdx
            
            dat = squeeze(data(p,:,s,d)); dat = dat(~isnan(dat)); % single time course for a protein, inhibitor regime and one of the stimuli not filtered out. Any trailing NaNs are removed.
            
            [~,maxIdx] = max(dat(:)); % position of max

            if maxIdx>1 && maxIdx<=6 % no longer consider stimulus if peak is not at first time point (t=0) or at the last time point considered in subsequent analyses (7th time point; 4hrs) or later
                
                % remove peak maximum data point + check all amplitude conditions above still hold (check of robustness of peak)
                tmp = dat([1:maxIdx-1 maxIdx+1:end]); 
                if (max(tmp)-dat(1))>max(amp1)/3 && ((max(tmp)-dat(end))>max(amp2)/3 || dat(end)<=dat(1)) && (max(tmp)-dat(1))>range(tmp)/3
                                    
                    % find peaking region - move outwards from peak, to the left and stop before getting within 10% of first timepoint
                    %                       move outwards from peak, to the right and stop before getting below first timepoint OR if go below last timepoint
                    tmp = find((dat(1:maxIdx-1)-dat(1))<0.1*amp1(s),1,'last')+1;
                    if ~isempty(tmp)
                        peaking_region(p,s,1,d) = tmp;
                    else
                        peaking_region(p,s,1,d) = 1;
                    end
                    tmp = find( ( (dat(maxIdx+1:end)-dat(1))<0  ) | (dat(maxIdx+1:end)-dat(end))<0 ,1)+maxIdx-1;
                    if ~isempty(tmp)
                        peaking_region(p,s,2,d) = tmp;
                    else
                        peaking_region(p,s,2,d) = length(dat);
                    end                    
                    
                    % remove peak if there is another peak before peaking region, or a peak after peaking region that is very large
                    if any( (dat(1:peaking_region(p,s,1,d)-1)-dat(1))>amp1(s)/3 )|| any( (dat(peaking_region(p,s,2,d)+1:end)-dat(1))>amp1(s)*0.75 )
                        peaking_region(p,s,1,d) = 0; peaking_region(p,s,2,d) = 0;
                        continue
                    end
                    
                    % check gradient between consecutive time points to check mostly increasing before the peak maximum and decreasing after the peak maximum.
                    % Only one gradient error (decreasing before peak max, increasing after peak max) allowed per 4 consecutive time points - ignore ones where gradient is pretty much flat. 
                    %                                                  Can't have 3 in a row (even flatish ones).
                    %                                                  Any single incorrect gradient can't be steep
                    % before peak max
                    grad = dat(max(2,peaking_region(p,s,1,d)):maxIdx)-dat(max(1,peaking_region(p,s,1,d)-1):maxIdx-1);
                    tmp1 = find(grad<-0.1*amp1(s));
                    tmp1 = diff(tmp1);
                    tmp2 = find(tmp1==1);
                    if any(tmp1<=2) || any(diff(tmp2)==1) || any(grad<-0.5*amp1(s))
                        peaking_region(p,s,1,d) = 0; peaking_region(p,s,2,d) = 0;
                        continue
                    end
                    % after peak max
                    grad = dat(maxIdx+1:min(peaking_region(p,s,2,d)+1,length(dat)))-dat(maxIdx:min(peaking_region(p,s,2,d),length(dat)-1));
                    tmp1 = find(grad>0.1*amp2(s));
                    tmp1 = diff(tmp1);
                    tmp2 = find(tmp1==1);
                    if any(tmp1<=2) || any(diff(tmp2)==1) || any(grad>0.5*amp2(s))
                        peaking_region(p,s,1,d) = 0; peaking_region(p,s,2,d) = 0;
                        continue
                    end
                    
                    peak(p,s,d) = 1;
                    
                end
            end
        end
    end
    
    % For each DMSO peak, check for another inhibitor regime that has a similar peak (peaking region of DMSO is similar to that of a inhibitor regime, and maxima are not very different). Robustness check - peaks typically present for more than one inhibitor regime.
    for s=1:S
        if peak(p,s,1)
            match = 0;
            d = 2;
            while match==0 && d<=D
                if peak(p,s,d) && peaking_region(p,s,1,1)+1>=peaking_region(p,s,1,d) && peaking_region(p,s,2,1)-1<=peaking_region(p,s,2,d) && abs(max(data(p,:,s,d))-max(data(p,:,s,1)))<=(max(data(p,:,s,1))-data(p,1,s,1))/2
                    match = 1;
                end
                d = d+1;
            end
            if match==0
                % if no peak in another inhibitor regime, robustness test failed -> specify no peak in DMSO
                peak(p,s,1) = 0; peaking_region(p,s,1,1) = 0; peaking_region(p,s,2,1) = 0;
            end
        end
    end
end

% return DMSO only
peak = peak(:,:,1); peaking_region = peaking_region(:,:,:,1);
end

%%
function fdr_adj = fdr_medlsu(p_values)
% Executes the adaptive linear step-up FDR procedure with m0 (number of true null hypotheses) estimated using median p-value.
% Reference: Benjamini, Y., Krieger, A.M., & Yekutieli, D. (2006) Adaptive linear 
%     step-up procedures that control the false discovery rate. Biometrika. 93(3), 491-507.
%
% Inputs:
%   p_values - A column vector of p-values (values between 0 and 1)
%
% Output:
%   fdr_adj - FDR adjusted p-values

m = length(p_values); % number of hyptheses
fdr_adj = p_values*0; % initialise fdr adjusted values

[p_values order_pval] = sort(p_values); % sort p-values

m0hat = (m-m/2) ./ (1-median(p_values)); % estimate for number of true null hypotheses

% get fdr-corrected p-values
[~,fdr_cor_values] = fdr(p_values);
fdr_cor_values = fdr_cor_values*m0hat/m; % modify values to take m0 into account

% get fdr-adjusted p-values
for i=1:m
    fdr_adj(i) = min(fdr_cor_values(i:m));
end

[~,order_pval] = sort(order_pval);
fdr_adj = fdr_adj(order_pval);
end

%% FDR code below is taken from http://dl.dropbox.com/u/2785709/brainder/2011/fdr/fdr.m (linked to from https://brainder.org/2011/09/05/fdr-corrected-fdr-adjusted-p-values/)
function varargout = fdr(varargin)
% Computes the FDR-threshold for a vector of p-values.
%
% Usage:
% [pthr,pcor,padj] = fdr(pvals)
%                    fdr(pval,q)
%                    fdr(pval,q,cV)
%
% Inputs:
% pvals  = Vector of p-values.
% q      = Allowed proportion of false positives (q-value).
%          Default = 0.05.
% cV     = If set to anything but 1, uses an harmonic sum for c(V).
%          Default = 1.
%
% Outputs:
% pthr   = FDR threshold.
% pcor   = FDR corrected p-values.
% padj   = FDR adjusted p-values.
%
% Note that the corrected and adjusted p-values do **not** depend
% on the supplied q-value, but they do depend on the choice of c(V).
%
% References:
% * Benjamini & Hochberg. Controlling the false discovery
%   rate: a practical and powerful approach to multiple testing.
%   J. R. Statist. Soc. B (1995) 57(1):289-300.
% * Yekutieli & Benjamini. Resampling-based false discovery rate
%   controlling multiple test procedures for multiple testing
%   procedures. J. Stat. Plan. Inf. (1999) 82:171-96.
%
% ________________________________
% Anderson M. Winkler
% Research Imaging Center/UTHSCSA
% Dec/2007 (first version)
% Aug/2011 (this version)
% http://brainder.org

% Accept arguments
switch nargin,
    case 0,
        error('Error: Not enough arguments.');
    case 1,
        pval = varargin{1};
        qval = 0.05;
        cV   = 1;
    case 2,
        pval = varargin{1};
        qval = varargin{2};
        cV   = 1;
    case 3,
        pval = varargin{1};
        qval = varargin{2};
        if varargin{3}==1, cV = 1;
        else cV = sum(1./(1:numel(pval))) ;
        end
    otherwise
        error('Error: Too many arguments.')
end

% Check if pval is a vector
if numel(pval) ~= length(pval),
    error('p-values should be a row or column vector, not an array.')
end

% Check if pvals are within the interval
if min(pval) < 0 || max(pval) > 1,
    error('Values out of range (0-1).')
end

% Check if qval is within the interval
if qval < 0 || qval > 1,
    error('q-value out of range (0-1).')
end

% ========[PART 1: FDR THRESHOLD]========================================

% Sort p-values
[pval,oidx] = sort(pval);

% Number of observations
V = numel(pval);

% Order (indices), in the same size as the pvalues
idx = reshape(1:V,size(pval));

% Line to be used as cutoff
thrline = idx*qval/(V*cV);

% Find the largest pval, still under the line
thr = max(pval(pval<=thrline));

% Deal with the case when all the points under the line
% are equal to zero, and other points are above the line
if thr == 0,
    thr = max(thrline(pval<=thrline));
end

% Case when it does not cross
if isempty(thr), thr = 0; end

% Returns the result
varargout{1} = thr;

% ========[PART 2: FDR CORRECTED]========================================

if nargout == 2 || nargout == 3,
    
    % p-corrected
    pcor = pval.*V.*cV./idx;

    % Sort back to the original order and output
    [ignore,oidxR] = sort(oidx);
    varargout{2} = pcor(oidxR);
end

% ========[PART 3: FDR ADJUSTED ]========================================

if nargout == 3,

    % Loop over each sorted original p-value
    padj = zeros(size(pval));
    for i = 1:V,
        % The p-adjusted for the current p-value is the smallest slope among
        % all the slopes of each of the p-values larger than the current one
        % Yekutieli & Benjamini (1999), equation #3.
        padj(i) = min(pval(idx>=i).*V.*cV./idx(idx>=i));
    end
    varargout{3} = padj(oidxR);
end

% That's it!
end