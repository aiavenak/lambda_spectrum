clear all
close all
clc

path2code = fullfile(pwd);
addpath('~/mat_tools/BCT_2017'); % add BCT to path
load('../data/LAU2_scale3.mat','SC') % load SC from .mat file; SC is NxNxM - N nodes, M subjects

SCall = SC; clear SC;

[N,~,M] = size(SCall);

lambda_vals = logspace(1.45,-2.2,30);
lambda_vals = [0,lambda_vals(end:-1:1)];

model_ops = {'SPL_W_log','SPL_W_inv'}; % weight transforms 

model = model_ops{1};
rho = [.125]; % optional density parameter to ensure same density across all subjects

path2results = fullfile(pwd,'..','results',num2str(rho),model);
if ~exist(path2results,'dir')
    mkdir(path2results);
end

mask = find(~eye(N));

for s = 1:M
    
    SC = SCall(:,:,s);
    
    SC(eye(N)>0) = 0;
    
%     thr = prctile(SC(mask),(1-rho)*100); % threshold SC to desired dens
%     SC(SC <= thr) = 0;
%     if length(unique(get_components(SC))) > 1
%         error('\n Network %d is disconnected \n',s)
%     end   
    
    k = sum(SC,1);
    k = k(ones(N,1),:).*eye(N,N);
    Pref = k\SC;
    
    % calculate Local and Global matrices
    [G,D,B,steps] = f_local_global(SC,[],model);
    
    KLref = zeros(N,N,length(lambda_vals));
    Cinfo = zeros(N,N,length(lambda_vals));
    Ctrans = zeros(N,N,length(lambda_vals));
    visits = zeros(N,N,length(lambda_vals));
    
    tic;
    for lind = 1:length(lambda_vals)
        fprintf('lambda = %f \n',lambda_vals(lind));
        for t = 1:N
            % define G(t) - global term
            Gt = G(:,t)';
            Gt = Gt(ones(N,1),:) + D;
            
            % prob of going from i to j, given a target t
            FG = exp(-((lambda_vals(lind)*Gt) + D)).*(SC>0);
            Z_lamb_t = sum(FG,2);
            Z_lamb_t = Z_lamb_t(:,ones(1,N)); % Normalize to ensure sum(P,2) = ones(N,1)
            P = FG./Z_lamb_t;
            if sum(sum(P,2),1) ~= N  % sanity check 
                error(sprintf('\n ERROR -- subj %d \t lind %d \t target %d',s,lind,t));
            end
            % KL divergence
            if lind > 1
                kk = P.*(log2(P) - log2(Pref));
                kk(isnan(kk)) = 0; kk(isinf(kk)) = 0;
                KLref(:,t,lind) = sum(kk,2);
            end
            % calculate mean absorbtion times
            Q = P;
            Q(t,:) = [];  % make target an absorbing state
            Q(:,t) = [];
            FM = inv(eye(N-1) - Q);    % Fundamental Matrix
            ind = true(N,1); ind(t) = false;
            Cinfo(ind,t,lind) = (bsxfun(@(a,b) a./b,FM,sum(FM,2)))*KLref(ind,t,lind);
            Ctrans(ind,t,lind) = FM*(sum(P(ind,:).*(D(ind,:)),2));
            visits(ind,t,lind) = mean(FM,1)';
            
        end
    end
    fname = fullfile(path2results,sprintf('BRW_%d.mat',s));    
    save(fname,'lambda_vals','KLref','Cinfo','Ctrans','visits')
    toc
end


