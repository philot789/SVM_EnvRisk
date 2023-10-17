function results = sdem_panel_FE_g(y,x,W,T,ndraw,nomit,prior)
% PURPOSE: MCMC SDEM estimates for static spatial panels 
%          (N regions*T time periods) with spatial fixed effects (sfe) 
%          and/or time period fixed effects (tfe)
%          y = X*b + WX*g + sfe(optional) + tfe(optional) + u,  u = p*W*u + e, 
%            e = N(0,sige*V), V = diag(v_1,v_2,...v_N*T) 
%            r/vi = ID chi(r)/r, r = 5 default
%            [b g]' = N(c,C), 
%            1/sige = Gamma(nu,d0), 
% Supply data sorted first by time and then by spatial units, so first region 1,
% region 2, et cetera, in the first year, then region 1, region 2, et
% cetera in the second year, and so on
% sem_panel_FE_g computes y and x in deviation of the spatial and/or time means
% ---------------------------------------------------
%  USAGE: results = sdem_panel_FE_g(y,x,W,T,ndraw,nomit,prior)
%  where:  y = dependent variable vector
%          x = independent variables matrix
%          The function adds W*x for you
%          W = spatial weights matrix (standardized) (N x N)
%          T = number of points in time
%       prior = a structure variable with input options:
%       prior.novi_flag = 1, for e = N(0,sige*I), homoscedastic model
%                       = 0, for e = N(0,sige*V), heteroscedastic model
%                            sets V = diag(v_1,v_2,...v_N*T), r/vi = ID chi(r)/r, r = 5 (default)
%       prior.rval = rval, r prior hyperparameter, default=5
%       prior.model = 0 pooled model without fixed effects (default, x may contain an intercept)
%                  = 1 spatial fixed effects (x may not contain an intercept)
%                  = 2 time period fixed effects (x may not contain an intercept)
%                  = 3 spatial and time period fixed effects (x may not contain an intercept)
%       prior.fe    = report fixed effects and their t-values in prt_sp (default=0=not reported; info.fe=1=report) 
%       prior.beta, prior means for beta,   b (default (k x 1) vector = 0)
%       priov.bcov, prior beta covariance, C above (default eye(k)*1e+12)
%       prior.rval, rval prior hyperparameter, default=5
%       prior.nu,   informative Gamma(nu,d0) prior on sige
%       prior.d0    informative Gamma(nu,d0) prior on sige
%                   default for above: nu=0,d0=0 (diffuse prior)
%       prior.rmin  = (optional) minimum value of rho to use in search  
%       prior.rmax  = (optional) maximum value of rho to use in search    
%       prior.lflag = 0 for full lndet computation (default = 1, fastest)
%                  = 1 for MC lndet approximation (fast for very large problems)
%                  = 2 for Spline lndet approximation (medium speed)
%       prior.order = order to use with info.lflag = 1 option (default = 50)
%       prior.iter  = iterations to use with info.lflag = 1 option (default = 30)  
%       prior.lndet = a matrix returned containing log-determinant information to save time
% ---------------------------------------------------
%  RETURNS: a structure
%         results.meth  = 'sdemsfe_g' if info.model=1
%                       = 'sdemtfe_g' if info.model=2
%                       = 'sdemstfe_g' if info.model=3
%         results.beta  = bhat (includes coefficients on W*X)
%         results.rho   = rho (p above)
%         results.cov   = asymptotic variance-covariance matrix of the parameters b(eta) and rho
%         results.tstat = asymp t-stat (last entry is rho=spatial autoregressive coefficient)
%         results.direct   = nvar x 5 matrix with direct effect, t-stat, t-prob, lower01, upper99
%         results.indirect = nvar x 5 matrix with indirect effect, t-stat, t-prob, lower01, upper99
%         results.total    = nvar x 5 matrix with total effect, t-stat, t-prob, lower01, upper99
%         results.direct_draws   = ndraw x nvar matrix of direct effect draws
%         results.indirect_draws = ndraw x nvar matrix of indirect effect draws
%         results.total_draws    = ndraw x nvar matrix of total effect draws
%         results.bmean = b prior means (prior.beta from input)
%         results.bstd  = b prior std deviation, sqrt(diag(prior.bcov))
%         results.nu    = prior nu-value for sige prior (default = 0)
%         results.d0    = prior d0-value for sige prior (default = 0)
%         results.iprior = 1 for informative prior on beta, 
%                        = 0 for default no prior on beta
%         results.yhat  = [inv(y-p*W)]*[x*b+fixed effects] (according to prediction formula)
%         results.resid = y-x*b
%         results.sige  = e'(I-p*W)'*(I-p*W)*e/nobs
%         results.rsqr  = rsquared
%         results.corr2 = goodness-of-fit between actual and fitted values
%         results.sfe   = spatial fixed effects (if info.model=1 or 3)
%         results.tfe   = time period fixed effects (if info.model=2 or 3)
%         results.tsfe  = t-values spatial fixed effects (if info.model=1 or 3)
%         results.ttfe  = t-values time period fixed effects (if info.model=2 or 3)
%         results.con   = intercept 
%         results.con   = t-value intercept
%         results.lik   = log likelihood
%         results.nobs  = # of observations
%         results.nvar  = # of explanatory variables in x 
%         results.tnvar = # fixed effects
%         results.iter  = # of iterations taken
%         results.rmax  = 1/max eigenvalue of W (or rmax if input)
%         results.rmin  = 1/min eigenvalue of W (or rmin if input)
%         results.lflag = lflag from input
%         results.fe    = fe from input
%         results.liter = info.iter option from input
%         results.order = info.order option from input
%         results.limit = matrix of [rho lower95,logdet approx, upper95] intervals
%                         for the case of lflag = 1
%         results.time1 = time for log determinant calcluation
%         results.time2 = time for eigenvalue calculation
%         results.time4 = time for MCMC sampling
%         results.time  = total time taken      
%         results.lndet = a matrix containing log-determinant information
%                          (for use in later function calls to save time)
% --------------------------------------------------
%  NOTES: if you use lflag = 1 or 2, info.rmin will be set = -1 
%                                    info.rmax will be set = 1
%         For number of spatial units < 500 you should use lflag = 0 to get
%         exact results, 
%         Fixed effects and their t-values are calculated as the deviation
%         from the mean intercept
% ---------------------------------------------------

% written by James P LeSage (Texas State University)
% last updated 11/2020

[nobs, nvar] = size(x);
N = nobs/T;

    [n1,n2] = size(W);
    
    if n1 == N*T
        Wlarge = 1;
    else
        Wlarge = 0;
        W = kron(eye(T),W);
    end

    

xraw = x;

time1 = 0; 
time2 = 0;
time3 = 0;
time4 = 0;

timet = clock; % start the clock for overall timing

fe=0;
model=0;
novi_flag = 0;
rval = 5; % rval = 5 is default
results.rval = rval;
fields = fieldnames(prior);
nf = length(fields);
if nf > 0
    for i=1:nf
        if strcmp(fields{i},'model'), model = prior.model;
        elseif strcmp(fields{i},'fe'), fe = prior.fe;
        elseif strcmp(fields{i},'rval')
                    rval = prior.rval;
        results.rval = rval;
        novi_flag = 0;
    elseif strcmp(fields{i},'novi_flag')
        novi_flag = prior.novi_flag;
        if novi_flag == 1
            results.rval = 0;
        elseif novi_flag == 0
            results.rval = rval;
        end
        results.novi_flag = novi_flag;
    end
end

% check if the user handled the intercept term okay
    n = length(y);
    if sum(x(:,1)) ~= n
        tst = sum(x); % we may have no intercept term
        ind = find(tst == n); % we do have an intercept term
        if length(ind) > 0
            error('sdm_conv_panel_g: intercept term must be in first column of the x-matrix');
        elseif length(ind) == 0 % case of no intercept term
            cflag = 0;
            pp = size(x,2);
        end
    elseif sum(x(:,1)) == n % we have an intercept in the right place
        cflag = 1;
        pp = size(x,2)-1;
    end
      
    if cflag == 0
            xmat = [x W*x];
    elseif cflag == 1
        xmat = [x W*x(:,2:end)];
    end
    
nvarw = size(xmat,2);
    
prior.novi_flag = novi_flag;

results = sem_panel_FE_g(y,xmat,W,T,ndraw,nomit,prior);

results.cflag = cflag;
results.p = pp;
results.fe = fe;
results.nvarw = nvarw;
results.nvar = nvar;

if cflag == 0
results.direct_draws = results.bdraw(:,1:nvar);
results.indirect_draws = results.bdraw(:,nvar+1:end);
results.total_draws =  results.direct_draws + results.indirect_draws;
elseif cflag == 1
 results.direct_draws = results.bdraw(:,2:nvar);
results.indirect_draws = results.bdraw(:,nvar+1:end);
results.total_draws =  results.direct_draws + results.indirect_draws;
end   

% Compute means, std deviation and upper and lower 0.95 intervals
p = results.p;
total_out = zeros(p,5);
direct_out = zeros(p,5);
indirect_out = zeros(p,5);

for i=1:p
    tmp = squeeze(results.total_draws(:,i)); % an ndraw by 1  matrix
    total_mean = mean(tmp);
    total_std = std(tmp);
    % Bayesian 0.95 credible intervals
    % for the cumulative total effects
    bounds = cr_interval(tmp,0.95);
    ubounds = bounds(1,1);
    lbounds = bounds(1,2);
    total_out(i,:) = [total_mean total_mean./total_std tdis_prb(total_mean./total_std,nobs) lbounds ubounds];
    
    tmp = squeeze(results.direct_draws(:,i)); % an ndraw by 1  matrix
    direct_mean = mean(tmp);
    direct_std = std(tmp);
    % Bayesian 0.95 credible intervals
    % for the cumulative total effects
    bounds = cr_interval(tmp,0.95);
    ubounds = bounds(1,1);
    lbounds = bounds(1,2);
    direct_out(i,:) = [direct_mean direct_mean./direct_std tdis_prb(direct_mean./direct_std,nobs) lbounds ubounds];

    tmp = squeeze(results.indirect_draws(:,i)); % an ndraw by 1  matrix
    indirect_mean = mean(tmp);
    indirect_std = std(tmp);
    % Bayesian 0.95 credible intervals
    % for the cumulative total effects
    bounds = cr_interval(tmp,0.95);
    ubounds = bounds(1,1);
    lbounds = bounds(1,2);
    indirect_out(i,:) = [indirect_mean indirect_mean./indirect_std tdis_prb(indirect_mean./indirect_std,nobs) lbounds ubounds];

end

results.total = total_out;
results.direct = direct_out;
results.indirect = indirect_out;


if model==0
    results.meth='psdem_g';
elseif model==1
    results.meth='sdemsfe_g';
elseif model==2
    results.meth='sdemtfe_g';
elseif model==3
    results.meth='sdemstfe_g';
else
    error('sdem_panel_FE_g: wrong input number of prior.model');
end

fields = fieldnames(prior);
nf = length(fields);
% novi_flag = 0;
% rval = 4; % rval = 4 is default
% if results.novi == 0
%     results.rval = rval;
% else
%     results.rval = 0;
% end
nu = 0; d0 = 0; % default to a diffuse prior on sige
c = zeros(nvar,1); 
Tj = eye(nvar)*1e+12;
Q = inv(Tj);
Qpc = Q*c;
iprior = 0;
results.iprior = 0;
for i=1:nf
%     if strcmp(fields{i},'rval')
%         rval = prior.rval;
%         results.rval = rval;
    if strcmp(fields{i},'nu')
        nu = prior.nu;
        results.nu = nu;
    elseif strcmp(fields{i},'d0')
        d0 = prior.d0;
        results.d0 = d0;
%     elseif strcmp(fields{i},'novi_flag')
%         novi_flag = prior.novi_flag;
%         if novi_flag == 1
%             results.rval = 0;
%         elseif novi_flag == 0
%             results.rval = rval;
%         end
%         results.novi_flag = novi_flag;
    elseif strcmp(fields{i},'beta')
        c = prior.beta;
        if size(c,1) ~= nvarw
            error('sdem_panel_FE_g: wrong size prior means, must be 2*k x 1 vector');
        end
        results.iprior = 1;
        results.bmean = c;
    elseif strcmp(fields{i},'bcov')
        TI = prior.bcov;
        if size(TI,2) ~= nvarw
            error('sdem_panel_FE_g: wrong size prior variance-covariance, must be 2*k x 2*k matrix');
        end
        results.iprior = 1;
        results.bstd = diag(sqrt(TI));
        
        Q = inv(TI);
        Qpc = Q*c;
        
    end
end


if novi_flag == 1
    results.homo = 1;
    results.hetero = 0;
elseif novi_flag == 0
    results.hetero = 1;
    results.homo = 0;
end

en=ones(T,1);
et=ones(N,1);
ent=ones(nobs,1);

xraw = xmat;
[ywith,xwith,meanny,meannx,meanty,meantx]=demean(y,xmat,N,T,model);
x = xmat;
bmean = results.beta;
sige = results.sige;
% step 4) find fixed effects and their t-values
if model==1
    intercept=mean(y)-mean(xraw)*bmean;
    results.con=intercept;
    results.sfe=meanny-meannx*bmean-kron(et,intercept);
    xhat=x*bmean+kron(en,results.sfe)+kron(ent,intercept);
    results.tsfe=results.sfe./sqrt(sige/T*ones(N,1)+diag(sige*meannx*(xwith'*xwith)*meannx'));
    results.tcon=results.con/sqrt(sige/nobs+sige*mean(xraw)*(xwith'*xwith)*mean(xraw)');
    tnvar=N;
elseif model==2
    intercept=mean(y)-mean(xraw)*bmean;
    results.con=intercept;
    results.tfe=meanty-meantx*bmean-kron(en,intercept); 
    xhat=x*bmean+kron(results.tfe,et)+kron(ent,intercept);
    results.ttfe=results.tfe./sqrt(sige/N*ones(T,1)+diag(sige*meantx*(xwith'*xwith)*meantx'));
    results.tcon=results.con/sqrt(sige/nobs+sige*mean(xraw)*(xwith'*xwith)*mean(xraw)');
    tnvar=T;
elseif model==3
    intercept=mean(y)-mean(xraw)*bmean; 
    results.con=intercept;
    results.sfe=meanny-meannx*bmean-kron(et,intercept);
    results.tfe=meanty-meantx*bmean-kron(en,intercept);
    results.tsfe=results.sfe./sqrt(sige/T*ones(N,1)+diag(sige*meannx*(xwith'*xwith)*meannx'));
    results.ttfe=results.tfe./sqrt(sige/N*ones(T,1)+diag(sige*meantx*(xwith'*xwith)*meantx'));
    results.tcon=results.con/sqrt(sige/nobs+sige*mean(xraw)*(xwith'*xwith)*mean(xraw)');
    xhat=x*bmean+kron(en,results.sfe)+kron(results.tfe,et)+kron(ent,intercept);
    tnvar=N+T;
else
    xhat=xraw*bmean;
    tnvar=0;
end
results.tnvar=tnvar;
results.resid = y - xhat; 
yme=y-mean(y);
rsqr2=yme'*yme;
rsqr1 = results.resid'*results.resid;
results.rsqr=1.0-rsqr1/rsqr2; %rsquared


 % r-squared and corr-squared between actual and fitted values
results.tnvar=tnvar;
results.resid = y - xhat; 
yme=y-mean(y);
rsqr2=yme'*yme;
rsqr1 = results.resid'*results.resid;
results.rsqr=1.0-rsqr1/rsqr2; %rsquared

yhat=xhat;

ywithhat=xwith*results.beta;
res1=ywith-mean(ywith);
res2=ywithhat-mean(ywith);
rsq1=res1'*res2;
rsq2=res1'*res1;
rsq3=res2'*res2;
results.corr2=rsq1^2/(rsq2*rsq3); %corr2
results.yhat=yhat;

% % calculate t-stats
bsave = results.bdraw;
psave = results.pdraw;
bmean = mean(bsave);
rmean = mean(psave);
bstd = std(bsave);
pstd = std(psave);

mean_parms = [bmean'
              rmean];
std_parms = [bstd'
             pstd];
         
tstat = mean_parms./std_parms;

results.tstat = tstat;

% parms = [bsave psave ssave];
% 
% results.cov = cov(parms);


% return stuff
% results.nobs  = nobs; 
% results.nvar  = nvar;
% results.rmax  = rmax;      
% results.rmin  = rmin;
% results.lflag = ldetflag;
% results.order = order;
% results.miter = miter;
results.fe    = fe;
results.time  = etime(clock,timet);
% results.time1 = time1;
% results.time2 = time2;
% results.time3 = time3;
% results.time4 = time4;
% results.lndet = detval;
results.N     = N;
results.T     = T;
results.model = model;

end

function [ywith,xwith,meanny,meannx,meanty,meantx]=demean(y,x,N,T,model)
% demeaning of the y and x variables, depending on (info.)model
[nobs nvar]=size(x);
meanny=zeros(N,1);
meannx=zeros(N,nvar);
meanty=zeros(T,1);
meantx=zeros(T,nvar);

if (model==1 | model==3);
for i=1:N
    ym=zeros(T,1);
    xm=zeros(T,nvar);
    for t=1:T
        ym(t)=y(i+(t-1)*N,1);
        xm(t,:)=x(i+(t-1)*N,:);
    end
    meanny(i)=mean(ym);
    meannx(i,:)=mean(xm);
end
clear ym wym xm;
end % if statement

if ( model==2 | model==3)
for i=1:T
    t1=1+(i-1)*N;t2=i*N;
    ym=y([t1:t2],1);
    xm=x([t1:t2],:);
    meanty(i)=mean(ym);
    meantx(i,:)=mean(xm);
end
clear ym wym xm;
end % if statement
    
en=ones(T,1);
et=ones(N,1);
ent=ones(nobs,1);

if model==1
    ywith=y-kron(en,meanny);
    xwith=x-kron(en,meannx);
elseif model==2
    ywith=y-kron(meanty,et);
    xwith=x-kron(meantx,et);
elseif model==3
    ywith=y-kron(en,meanny)-kron(meanty,et)+kron(ent,mean(y));
    xwith=x-kron(en,meannx)-kron(meantx,et)+kron(ent,mean(x));
else
    ywith=y;
    xwith=x;
end % if statement

