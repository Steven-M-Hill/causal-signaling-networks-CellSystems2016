function [G,Gj] = JNI_serial_mod(D,I,G0,lambda,eta,dmax,var_sel,zeroIntercept,multiplicityPrior,selfEdges,fixedEffectRegime)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Code to reproduce analysis in
% Hill, Nesser et al. Cell Systems 4, 73-83 (2017), DOI: 10.1016/j.cels.2016.11.013.
% 
% This is a modified version of the Joint Network Inference approach of Oates et al. Annals of Applied Statistics 8, 1892–1919 (2014).
% Original version of the code is available at DOI: 10.1214/14-AOAS761SUPPB, © Chris Oates 2013.
% Modifications made: added multiplicityPrior, selfEdges and fixedEffectRegime options.
%
% [G,Gj] = JNI_serial(D,I,G0,lambda,eta,dmax,var_sel,zeroIntercept,multiplicityPrior,selfEdges,fixedEffectRegime)
%
% Inputs:
% D = cell(J,1); This is going to hold the time course data.
% D{j} = cell(Ej,1); D{j} holds data on individual j.
% D{j}{e} = Pxn data matrix; Time course e for individual j.
% I = cell(J,1); This is going to hold the experimental design.
% I{j} = PxEj LOGICAL indicator of drug treatment; Entries 1 correspond to perfect inhibition, entries 0 correspond to no inhibition.
% dmax = maximum in degree;
% G0 = prior network (binary matrix);
% lambda = individual regularisation hyperparameter;
% eta = latent network regularisation hyperparameter;
% var_sel = choice of local likelihood function; Options are 1 = classical Zellner, 2 = empirical Bayes Zellner (used in Oates and Mukherjee 2013), 3 = Deltell et al. AOS 2012
% zeroIntercept - binary, include zero intercept in model and include 1st timepoint in response vector?
% multiplicityPrior - binary, include soft multiplicity model prior?
% selfEdges - binary, allow self-edges to be inferred?
% fixedEffectRegime - integer: either 1 or 2. 
%                              1 = inhibitors with same target have one fixed effect (FE) parameter between them, and combinations of inhibitors have additive FEs (sum of individual inhibitor FEs)
%                              2 = inhibitors with same target each have their own FE, as does a combination of inhibitors
% Outputs:
% G = PxP posterior edge probabilities for G;
% Gj = cell(J,1);
% Gj{j} = PxP posterior edge probabilities for Gj;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%
% Initialise %
%%%%%%%%%%%%%%
P = size(D{1}{1},1);
n = size(D{1}{1},2);
J = length(D);
G = zeros(P);
Gj = cell(J,1);
Ej = zeros(J,1);
for j = 1:J
    Gj{j} = zeros(P);
    Ej(j) = length(D{j});
end

%%%%%%%%%%%%%%%%%%%%
% Enumerate models %
%%%%%%%%%%%%%%%%%%%%
% models = PxM binary matrix;
models = enumerate_models(P,dmax);
M = size(models,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cache marginal likelihoods %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MLs = zeros(J,P,M); 
for j = 1:J
    for m = 1:M
        % Build design matrix for line j, model m.
        % All experimental data enters the design matrix.
        num_prot = sum(models(:,m));
        DT = sum(I{j}(models(:,m),:),2)>0; % drug targets over all experiments
        X0 = []; Xm = []; y = []; 
        for e = 1:Ej(j)
            if zeroIntercept
                X0e = [1 zeros(1,n-1); ones(1,n)]'; % or [1 zeros(1,n-1); 0 ones(1,n-1)]' but improper prior so the same
                Xme = [zeros(1,num_prot); D{j}{e}(models(:,m),1:end-1)'];
                ye = D{j}{e};
                % Fixed effect, perfect out interventions:
                DT_e = I{j}(models(:,m),e); % drug targets in experiment e
                n_DT_e = sum(DT_e); % num. drug targets in experiment e
                Xme(:,DT_e) = nan(n,n_DT_e);
                if fixedEffectRegime==1
                    Xme = [Xme repmat(DT_e(DT)',n,1)];
                elseif fixedEffectRegime==2
                    Xme = [Xme zeros(n,Ej(j))];
                    if n_DT_e
                        Xme(:,num_prot+e) = 1;
                    end
                end
            else
                X0e = ones(n-1,1);
                Xme = D{j}{e}(models(:,m),1:end-1)';
                ye = D{j}{e}(:,2:end);
                % Fixed effect, perfect out interventions:
                DT_e = I{j}(models(:,m),e); % drug targets in experiment e
                n_DT_e = sum(DT_e); % num. drug targets in experiment e
                Xme(:,DT_e) = nan(n-1,n_DT_e);
                if fixedEffectRegime==1
                    Xme = [Xme repmat(DT_e(DT)',n-1,1)];
                elseif fixedEffectRegime==2
                    Xme = [Xme zeros(n-1,Ej(j))];
                    if n_DT_e
                        Xme(:,num_prot+e) = 1;
                    end
                end
            end
            X0 = [X0; X0e];
            Xm = [Xm; Xme];
            y = [y ye];
        end
        if fixedEffectRegime==2
            Xm(:,all(Xm==0)) = [];
        end
        
        P0 = X0*pinv(X0); %P0 = X0*((X0'*X0)\(X0'));
        %%%%% model complexity
        [N,a] = size(X0);
        %b = num_prot; % rather than b = size(Xm,2); % model complexity                
        b = size(Xm,2);
        %%%%% Degenerate case where e.g. a predictor is inhibited in all experiments
        ix = or((std(Xm,0,1)==0),(sum(~isnan(Xm),1)==0)); % indices of columns containing no variation
        if sum(ix)>0
            Xm = Xm(:,~ix); % remove problem column(s)         
        end        
        %%%%% Orthogonalisation
        % An extension of Xm = (eye(N)-P0)*Xm; which preserves zeros:
        for p = 1:size(Xm,2)
            wh = (~isnan(Xm(:,p)));
            Xm(wh,p) = (eye(sum(wh)) - X0(wh,:)*pinv(X0(wh,:)))*Xm(wh,p);
        end
        Xm(isnan(Xm)) = 0;
        %%%%%
        Pm = Xm*pinv(Xm); %Pm = Xm*((Xm'*Xm)\(Xm'));
        for p = 1:P 
            if var_sel == 1
                % classical Zellner
                g = N;
                err = y(p,:)*((eye(N)-P0-g*Pm/(g+1))*(y(p,:)'));
                MLs(j,p,m) = (1/2) * gamma((N-a)/2) * pi^(-(N-a)/2) ...
                              * det(X0'*X0)^(-1) * (g+1)^(-b/2) ...
                              * err^(-(N-a)/2);                          
            elseif var_sel == 2
                % empirical Bayes Zellner
                g_w = [N 19 9 4 2 33 1.5 1]; % prior weight = 100/N,5,10,20,30,40,50% respectively
                for i = 1:length(g_w)
                    g = g_w(i);
                    err = y(p,:)*((eye(N)-P0-g*Pm/(g+1))*(y(p,:)'));
                    EBML = (1/2) * gamma((N-a)/2) * pi^(-(N-a)/2) ...
                                 * det(X0'*X0)^(-1) * (g+1)^(-b/2) ...
                                 * err^(-(N-a)/2);                      
                   if EBML > MLs(j,p,m)
                       MLs(j,p,m) = EBML;
                   end
                end
            elseif var_sel == 3
                % Deltell (AOS, 2012)
                SSE0 = y(p,:)*((eye(N)-P0)*(y(p,:)'));
                SSEi = y(p,:)*((eye(N)-Pm)*(y(p,:)'));
                MLs(j,p,m) = (SSEi/SSE0)^(-(N-a)/2) * (1/(b+1)) * ((1+N)/(a+b+1))^(-b/2) ...
                              * hypergeom([(1+b)/2,(N-a)/2],(3+b)/2,(1-SSE0/SSEi)*(a+b+1)/(1+N));
            end
        end 
    end   
    % numerical regularisation:
    for p = 1:P
        MLs(j,p,:) = MLs(j,p,:) / max(MLs(j,p,:));
    end
end

%%%%%%%%%%%%%%%%%%%%
% Cache p(G_p|G^0) %
%%%%%%%%%%%%%%%%%%%%
TP = cell(P,1);
M2 = zeros(1,dmax+1);
for i=0:dmax
    M2(i+1) = nchoosek(p,i);
end
for p = 1:P
    TP{p} = typical_prior(p,G0,eta,models,M,M2,dmax,multiplicityPrior,selfEdges);
end

%%%%%%%%%%%%%%%%%%%%%%
% Cache p(G_p^j|G_p) %
%%%%%%%%%%%%%%%%%%%%%%
[CLP normCnsts] = cell_line_prior(lambda,models,M2,dmax,multiplicityPrior,selfEdges);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cache \mathfrak{P}(\bm{y}_p^j|G_p) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C1 = zeros(P,J,M);
for p = 1:P
    for j = 1:J
        mlsjp = squeeze(MLs(j,p,:));
        for m = 1:M
            if selfEdges
                C1(p,j,m) = sum( mlsjp .* (CLP(:,m)/normCnsts{p}(m)) );
            else
                C1(p,j,m) = sum( mlsjp(~models(p,:)) .* (CLP(~models(p,:),m)/normCnsts{p}(m)) );
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cache p(G_p^j|\bm{y}_p,G_p^0) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C2 = zeros(P,J,M);
for p = 1:P
    tpp = TP{p};
    c1p = prod(squeeze(C1(p,:,:)),1)';
    for j = 1:J
        c1pj = squeeze(C1(p,j,:));
        for m = 1:M
            if selfEdges
                C2(p,j,m) = MLs(j,p,m) * sum( tpp .* (CLP(m,:)./normCnsts{p})' .* (c1p ./ c1pj) );
            else
                C2(p,j,m) = MLs(j,p,m) * sum( tpp(~models(p,:)) .* (CLP(m,(~models(p,:)))./normCnsts{p}(~models(p,:)))' .* (c1p(~models(p,:)) ./ c1pj(~models(p,:))) );
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variable selection decomposition %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for p = 1:P 
    for i = 1:P
        modelsi = models(i,:)';
        if selfEdges
            G(i,p) = sum( modelsi .* TP{p} .* prod(squeeze(C1(p,:,:)),1)' ) ...
                / sum( TP{p} .* prod(squeeze(C1(p,:,:)),1)' );
            for j = 1:J
                Gj{j}(i,p) = sum( modelsi .* squeeze(C2(p,j,:)) ) ...
                    / sum( squeeze(C2(p,j,:)) );
            end
        else
            G(i,p) = sum( modelsi(~models(p,:)) .* TP{p}(~models(p,:)) .* prod(squeeze(C1(p,:,~models(p,:))),1)' ) ...
                / sum( TP{p}(~models(p,:)) .* prod(squeeze(C1(p,:,~models(p,:))),1)' );
            for j = 1:J
                Gj{j}(i,p) = sum( modelsi(~models(p,:)) .* squeeze(C2(p,j,~models(p,:))) ) ...
                    / sum( squeeze(C2(p,j,~models(p,:))) );
            end
        end
    end
end

end

function models = enumerate_models(P,dmax)
% number of possible parent sets
n_models = 0;
for i=0:dmax
%     if selfEdges
        n_models = n_models + nchoosek(P,i);
%     else
%         n_models = n_models + nchoosek(P,i);
%     end
end
% calculate all parent sets
models = false(P,n_models);
counter = 0;
for d = 0:dmax
    nck = nchoosek(1:P,d);
    for i = 1:size(nck,1);
        counter = counter + 1;
        for j = 1:size(nck,2)
            models(nck(i,j),counter) = true;
        end
    end
end
end

% Returns z = p(G_p|G^0), an Mx1 vector.
function z = typical_prior(p,G0,eta,models,M,M2,dmax,multiplicityPrior,selfEdges)
if sum(sum(~ismember(G0,[0 1]))) == 0 % if G0 is a prior network
    tmp = sum(abs(models - repmat(G0(:,p),1,M)));
    tmp = tmp - min(tmp); % numerical regularisation
    z = exp(-eta*tmp);
else % if G0 contains prior edge probabilities
    z = prod((repmat(G0(:,p),1,M).^models).*(repmat(1-G0(:,p),1,M).^(1-models)))';
    %     z = z/sum(z);
end
if multiplicityPrior
    z = z./(M2(sum(models)+1)*(dmax+1));
end
if selfEdges
    z = z'/sum(z);
else
    z = z'/sum(z(~models(p,:)));
end
end

% Returns z = p(G_p^j|G_p), an MxM matrix.
function [z zNormCnsts] = cell_line_prior(lambda,models,M2,dmax,multiplicityPrior,selfEdges)
tmp = pdist(+models','cityblock');
tmp = squareform(tmp);
z = exp(-lambda*tmp);
if multiplicityPrior
    z = bsxfun(@rdivide,z,M2(sum(models)+1)'*(dmax+1));
end
if selfEdges
    for p=1:size(models,1)
        zNormCnsts{p} = sum(z);
    end
%     z = bsxfun(@rdivide,z,sum(z));
else
    for p=1:size(models,1)
        zNormCnsts{p} = sum(z(~models(p,:),:));
    end
end
end
