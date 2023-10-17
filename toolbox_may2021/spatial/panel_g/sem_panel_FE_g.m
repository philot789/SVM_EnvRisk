function results = sem_panel_FE_g(y,x,W,T,ndraw,nomit,prior)
% PURPOSE: computes MCMC spatial error model estimates for static spatial panels 
%          (N regions*T time periods) with spatial fixed effects (sfe) 
%          and/or time period fixed effects (tfe)
%          y = XB + sfe(optional) + tfe(optional) + u,  u = p*W*u + e, 
%            e = N(0,sige*V), V = diag(v_1,v_2,...v_N*T) 
%            r/vi = ID chi(r)/r, r = Gamma(m,k)
%            B = N(c,C), 
%            1/sige = Gamma(nu,d0), 
% Supply data sorted first by time and then by spatial units, so first region 1,
% region 2, et cetera, in the first year, then region 1, region 2, et
% cetera in the second year, and so on
% sem_panel_FE_g computes y and x in deviation of the spatial and/or time means
% ---------------------------------------------------
%  USAGE: results = sem_panel_FE_g(y,x,W,T,ndraw,nomit,prior)
%  where:  y = dependent variable vector
%          x = independent variables matrix
%          W = spatial weights matrix (standardized) (N x N)
%          T = number of points in time
%       prior = a structure variable with input options:
%       prior.novi_flag = 1, for e = N(0,sige*I), homoscedastic model
%                       = 0, for e = N(0,sige*V), heteroscedastic model
%                            sets V = diag(v_1,v_2,...v_N*T), rval/vi = ID chi(rval)/rval, rval = 5 (default)
%       prior.rval = rval, r prior hyperparameter, default=5
%       prior.model = 0 pooled model without fixed effects (default, x may contain an intercept)
%                  = 1 spatial fixed effects (x may not contain an intercept)
%                  = 2 time period fixed effects (x may not contain an intercept)
%                  = 3 spatial and time period fixed effects (x may not contain an intercept)
%       prior.fe    = report fixed effects and their t-values in prt_sp (default=0=not reported; info.fe=1=report) 
%       prior.beta, prior means for beta,   b (default (k x 1) vector = 0)
%       priov.bcov, prior beta covariance, C above (default eye(k)*1e+12)
%       prior.rval, rval prior hyperparameter, default=4
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
%       prior.lndet = a matrix returned by sar containing log-determinant information to save time
% ---------------------------------------------------
%  RETURNS: a structure
%         results.meth  = 'semsfe_g' if info.model=1
%                       = 'semtfe_g' if info.model=2
%                       = 'semstfe_g' if info.model=3
%         results.beta  = bhat
%         results.rho   = rho (p above)
%         results.bdraw = (ndraw-nomit) x k matrix of MCMC draws
%         results.pdraw = (ndraw-nomit) x 1 vector of MCMC draws
%         results.cov   = asymptotic variance-covariance matrix of the parameters b(eta) and rho
%         results.tstat = asymp t-stat (last entry is rho=spatial autoregressive coefficient)
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


time1 = 0; 
time2 = 0;
time3 = 0;
time4 = 0;

timet = clock; % start the clock for overall timing

fe=0;
model=0;
novi_flag = 1;
fields = fieldnames(prior);
nf = length(fields);
if nf > 0
    for i=1:nf
        if strcmp(fields{i},'model'), model = prior.model;
        elseif strcmp(fields{i},'fe'), fe = prior.fe;
        elseif strcmp(fields{i},'novi'), novi_flag = prior.novi_flag;
        elseif strcmp(fields{i},'rval'), novi_flag = 0;
        end
    end
end

[ywith,xwith,meanny,meannx,meanty,meantx]=demean(y,x,N,T,model);

prior.novi = novi_flag;

results = sem_g(ywith,xwith,W,ndraw,nomit,prior);

results.fe = fe;

if model==0
    results.meth='psem_g';
elseif model==1
    results.meth='semsfe_g';
elseif model==2
    results.meth='semtfe_g';
elseif model==3
    results.meth='semstfe_g';
else
    error('sem_panel_FE_g: wrong input number of prior.model');
end

fields = fieldnames(prior);
nf = length(fields);
novi_flag = 0;
rval = 5; % rval = 5 is default
if results.novi == 0
    results.rval = rval;
else
    results.rval = 0;
end
nu = 0; d0 = 0; % default to a diffuse prior on sige
c = zeros(nvar,1); 
Tj = eye(nvar)*1e+12;
Q = inv(Tj);
Qpc = Q*c;
iprior = 0;
results.iprior = 0;
for i=1:nf
    if strcmp(fields{i},'rval')
        rval = prior.rval;
        results.rval = rval;
    elseif strcmp(fields{i},'nu')
        nu = prior.nu;
        results.nu = nu;
    elseif strcmp(fields{i},'d0')
        d0 = prior.d0;
        results.d0 = d0;
    elseif strcmp(fields{i},'novi_flag')
        novi_flag = prior.novi_flag;
        if novi_flag == 1
            results.rval = 0;
        elseif novi_flag == 0
            results.rval = rval;
        end
        results.novi_flag = novi_flag;
    elseif strcmp(fields{i},'beta')
        c = prior.beta;
        if size(c,1) ~= nvar
            error('sem_panel_FE_g: wrong size prior means, must be k x 1 vector');
        end
        results.iprior = 1;
        results.bmean = c;
    elseif strcmp(fields{i},'bcov')
        TI = prior.bcov;
        if size(TI,2) ~= nvar
            error('sem_panel_FE_g: wrong size prior variance-covariance, must be k x k matrix');
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


bmean = results.beta;
sige = results.sige;
% step 4) find fixed effects and their t-values
if model==1
    intercept=mean(y)-mean(x)*bmean;
    results.con=intercept;
    results.sfe=meanny-meannx*bmean-kron(et,intercept);
    xhat=x*bmean+kron(en,results.sfe)+kron(ent,intercept);
    results.tsfe=results.sfe./sqrt(sige/T*ones(N,1)+diag(sige*meannx*(xwith'*xwith)*meannx'));
    results.tcon=results.con/sqrt(sige/nobs+sige*mean(x)*(xwith'*xwith)*mean(x)');
    tnvar=N;
elseif model==2
    intercept=mean(y)-mean(x)*bmean;
    results.con=intercept;
    results.tfe=meanty-meantx*bmean-kron(en,intercept); 
    xhat=x*bmean+kron(results.tfe,et)+kron(ent,intercept);
    results.ttfe=results.tfe./sqrt(sige/N*ones(T,1)+diag(sige*meantx*(xwith'*xwith)*meantx'));
    results.tcon=results.con/sqrt(sige/nobs+sige*mean(x)*(xwith'*xwith)*mean(x)');
    tnvar=T;
elseif model==3
    intercept=mean(y)-mean(x)*bmean; 
    results.con=intercept;
    results.sfe=meanny-meannx*bmean-kron(et,intercept);
    results.tfe=meanty-meantx*bmean-kron(en,intercept);
    results.tsfe=results.sfe./sqrt(sige/T*ones(N,1)+diag(sige*meannx*(xwith'*xwith)*meannx'));
    results.ttfe=results.tfe./sqrt(sige/N*ones(T,1)+diag(sige*meantx*(xwith'*xwith)*meantx'));
    results.tcon=results.con/sqrt(sige/nobs+sige*mean(x)*(xwith'*xwith)*mean(x)');
    xhat=x*bmean+kron(en,results.sfe)+kron(results.tfe,et)+kron(ent,intercept);
    tnvar=N+T;
else
    xhat=x*bmean;
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

function results = sem_g(y,x,W,ndraw,nomit,prior)
% PURPOSE: Bayesian estimates of the spatial error model
%          y = XB + u, u = rho*W + e
%          e = N(0,sige*V), V = diag(v1,v2,...vn) 
%          r/vi = ID chi(r)/r, r = Gamma(m,k)
%          B = N(c,T), 
%          1/sige = Gamma(nu,d0), 
%          rho = Uniform(rmin,rmax), or rho = beta(a1,a2); 
%-------------------------------------------------------------
% USAGE: results = sem_g(y,x,W,ndraw,nomit,prior)
% where: y = dependent variable vector (nobs x 1)
%        x = independent variables matrix (nobs x nvar)
%        W = spatial weight matrix (standardized, row-sums = 1)
%    ndraw = # of draws
%    nomit = # of initial draws omitted for burn-in            
%    prior = a structure variable with:
%            prior.beta  = prior means for beta,   c above (default 0)
%            priov.bcov  = prior beta covariance , T above (default 1e+12)
%            prior.novi  = 1 turns off sampling for vi, producing homoscedastic model            
%            prior.rval  = r prior hyperparameter, default=4
%            prior.m     = informative Gamma(m,k) prior on r
%            prior.k     = (default: not used)
%            prior.nu    = informative Gamma(nu,d0) prior on sige
%            prior.d0    = default: nu=0,d0=0 (diffuse prior)
%            prior.a1    = parameter for beta(a1,a2) prior on rho see: 'help beta_prior'
%            prior.a2    = (default = 1.0, a uniform prior on rmin,rmax) 
%            prior.rmin  = (optional) min rho used in sampling (default = -1)
%            prior.rmax  = (optional) max rho used in sampling (default = +1)  
%            prior.eigs  = 0 to compute rmin/rmax using eigenvalues, (1 = don't compute default)
%            prior.lflag = 0 for full lndet computation (default = 1, fastest)
%                        = 1 for MC approx (fast for large problems)
%                        = 2 for Spline approx (medium speed)
%            prior.dflag = 1 for Metropolis-Hastings sampling for rho (default)
%                        = 0 for griddy gibbs with univariate numerical integration
%            prior.order = order to use with prior.lflag = 1 option (default = 50)
%            prior.iter  = iters to use with prior.lflag = 1 option (default = 30)   
%            prior.lndet = a matrix returned by sar, sar_g, sarp_g, etc.
%                          containing log-determinant information to save time
%            prior.mlog  = 0 for no log-marginal likelihood, 
%                        = 1 for log-marginal likelihood, default = 1
%-------------------------------------------------------------
% RETURNS:  a structure:
%          results.meth   = 'sem_g'
%          results.beta   = posterior mean of bhat
%          results.rho    = posterior mean of rho
%          results.sige   = posterior mean of sige
%          results.bdraw  = bhat draws (ndraw-nomit x nvar)
%          results.pdraw  = rho  draws (ndraw-nomit x 1)
%          results.sdraw  = sige draws (ndraw-nomit x 1)
%          results.vmean  = mean of vi draws (nobs x 1) 
%          results.rdraw  = r draws (ndraw-nomit x 1) (if m,k input)
%          results.bmean  = b prior means, prior.beta from input
%          results.bstd   = b prior std deviations sqrt(diag(prior.bcov))
%          results.r      = value of hyperparameter r (if input)
%          results.nobs   = # of observations
%          results.nvar   = # of variables in x-matrix
%          results.ndraw  = # of draws
%          results.nomit  = # of initial draws omitted
%          results.y      = y-vector from input (nobs x 1)
%          results.yhat   = mean of posterior predicted (nobs x 1)
%          results.resid  = residuals, based on posterior means
%          results.rsqr   = r-squared based on posterior means
%          results.rbar   = adjusted r-squared
%          results.nu     = nu prior parameter
%          results.d0     = d0 prior parameter
%          results.a1     = a1 parameter for beta prior on rho from input, or default value
%          results.a2     = a2 parameter for beta prior on rho from input, or default value
%          results.time1  = time for eigenvalue calculation
%          results.time2  = time for log determinant calcluation
%          results.time3  = time for sampling
%          results.time   = total time taken  
%          results.rmax   = 1/max eigenvalue of W (or rmax if input)
%          results.rmin   = 1/min eigenvalue of W (or rmin if input)          
%          results.tflag  = 'plevel' (default) for printing p-levels
%                         = 'tstat' for printing bogus t-statistics 
%          results.lflag  = lflag from input
%          results.iter   = prior.iter option from input
%          results.order  = prior.order option from input
%          results.limit  = matrix of [rho lower95,logdet approx, upper95] 
%                           intervals for the case of lflag = 1
%          results.lndet = a matrix containing log-determinant information
%                          (for use in later function calls to save time)
%          results.mlike = log marginal likelihood for model comparisons,
%                          (a vector ranging over rho-values from rmin to rmax that can be
%                          integrated for model comparison)
%          results.acc   = acceptance rate for M-H sampling (ndraw x 1) vector
% --------------------------------------------------------------
% NOTES: - use either improper prior.rval 
%          or informative Gamma prior.m, prior.k, not both of them
% - for n < 1000 you should use lflag = 0 to get exact results  
% - use a1 = 1.0 and a2 = 1.0 for uniform prior on rho
% --------------------------------------------------------------
% SEE ALSO: (sem_gd, sem_gd2 demos) prt
% --------------------------------------------------------------
% REFERENCES: James P. LeSage, `Bayesian Estimation of Spatial Autoregressive
%             Models',  International Regional Science Review, 1997 
%             Volume 20, number 1\&2, pp. 113-129.
% For lndet information see: Ronald Barry and R. Kelley Pace, 
% "A Monte Carlo Estimator of the Log Determinant of Large Sparse Matrices", 
% Linear Algebra and its Applications", Volume 289, Number 1-3, 1999, pp. 41-54.
% and: R. Kelley Pace and Ronald P. Barry 
% "Simulating Mixed Regressive Spatially autoregressive Estimators", 
% Computational Statistics, 1998, Vol. 13, pp. 397-418.
%----------------------------------------------------------------

% written by:
% James P. LeSage, 12/2001, updated 7/2003
% Dept of Finance & Economics
% Texas State University
% 601 University Drive
% San Marcos, TX 78666
% jlesage@spatial-econometrics.com


% NOTE: some of the speed for large problems comes from:
% the use of methods pioneered by Pace and Barry.
% R. Kelley Pace was kind enough to provide functions
% lndetmc, and lndetint from his spatial statistics toolbox
% for which I'm very grateful.

timet = clock;

% error checking on inputs
[n junk] = size(y);
results.y = y;
[n1 k] = size(x);
[n3 n4] = size(W);
time1 = 0;
time2 = 0;
time3 = 0;

if n1 ~= n
error('sem_g: x-matrix contains wrong # of observations');
elseif n3 ~= n4
error('sem_g: W matrix is not square');
elseif n3~= n
error('sem_g: W matrix is not the same size at y,x');
end;

if nargin == 5
    prior.lflag = 1;
end;

[nu,d0,rval,mm,kk,rho,sige,rmin,rmax,detval,ldetflag,eflag,order,iter,novi_flag,c,T,cc,metflag,a1,a2,inform_flag,mlog] = sem_parse(prior,k);

if novi_flag == 1
    results.rval = 0;
end

results.order = order;
results.iter = iter;

% error checking on prior information inputs
[checkk,junk] = size(c);
if checkk ~= k
error('sem_g: prior means are wrong');
elseif junk ~= 1
error('sem_g: prior means are wrong');
end;

[checkk junk] = size(T);
if checkk ~= k
error('sem_g: prior bcov is wrong');
elseif junk ~= k
error('sem_g: prior bcov is wrong');
end;

V = ones(n,1); in = ones(n,1); % initial value for V   
ys = y.*sqrt(V);
vi = in;
          
bsave = zeros(ndraw-nomit,1);    % allocate storage for results
ssave = zeros(ndraw-nomit,1);
vmean = zeros(n,1);
yhat = zeros(n,1);

if mm~= 0                        % storage for draws on rvalue
rsave = zeros(ndraw-nomit,1);
end;

[rmin,rmax,time1] = sem_eigs(eflag,W,rmin,rmax,n);


results.rmin = rmin;
results.rmax = rmax;
results.lflag = ldetflag;

[detval,time2] = sem_lndet(ldetflag,W,rmin,rmax,detval,order,iter);


% storage for draws
          bsave = zeros(ndraw-nomit,k);
          if mm~= 0
          rsave = zeros(ndraw-nomit,1);
          end;
          psave = zeros(ndraw-nomit,1);
          ssave = zeros(ndraw-nomit,1);
          vmean = zeros(n,1);
          acc_rate = zeros(ndraw,1);

% ====== initializations
% compute this stuff once to save time
TI = inv(T);
TIc = TI*c;
iter = 1;
in = ones(n,1);
V = in;
Wy = sparse(W)*y;
Wx = sparse(W)*x;
vi = in;
V = vi;

switch (novi_flag) 
    
case{0} % we do heteroscedastic model
    
% hwait = waitbar(0,'sem\_g: MCMC sampling ...');
t0 = clock;                  
iter = 1;
acc = 0;
          while (iter <= ndraw); % start sampling;
                  
          % update beta   
          xs = matmul(sqrt(V),x);
          ys = sqrt(V).*y;;
          Wxs = W*xs;
          Wys = W*ys;
          xss = xs - rho*Wxs;
          AI = (xss'*xss + sige*TI)\eye(k);
          yss = ys - rho*Wys;
          b = xss'*yss + sige*TIc;
          b0 = AI*b;
          bhat = norm_rnd(sige*AI) + b0; 
            
          % update sige
          nu1 = n + 2*nu; 
          e = yss-xss*bhat;
          ed = e - rho*sparse(W)*e;
          d1 = 2*d0 + ed'*ed;
          chi = chis_rnd(1,nu1);
          sige = d1/chi;

          % update vi
          ev = ys - xs*bhat; 
          chiv = chis_rnd(n,rval+1);   
%           chiv = chis_rnd(m,rval+1);
          vi = ((ev.*ev/sige) + in*rval)./chiv;
          V = in./vi; 
              
          % update rval
          if mm ~= 0           
          rval = gamm_rnd(1,1,mm,kk);  
          end;

if metflag == 0
          % update rho using numerical integration          
          rho = draw_rho(detval,y,x,Wy,Wx,V,n,k,rmin,rmax,rho);
else
          % update rho using metropolis-hastings
          % numerical integration is too slow here
          xb = x*bhat;
          rhox = c_sem(rho,y,x,bhat,sige,W,detval,ones(n,1),a1,a2);
          accept = 0;
          rho2 = rho + cc*randn(1,1);
          while accept == 0
           if ((rho2 > rmin) & (rho2 < rmax)); 
           accept = 1;  
           else
           rho2 = rho + cc*randn(1,1);
           end; 
          end; 
           rhoy = c_sem(rho2,y,x,bhat,sige,W,detval,ones(n,1),a1,a2);
          ru = unif_rnd(1,0,1);
          if ((rhoy - rhox) > exp(1)),
          p = 1;
          else,          
          ratio = exp(rhoy-rhox);
          p = min(1,ratio);
          end;
              if (ru < p)
              rho = rho2;
              acc = acc + 1;
              end;
      acc_rate(iter,1) = acc/iter;
      % update cc based on std of rho draws
       if acc_rate(iter,1) < 0.4
       cc = cc/1.1;
       end;
       if acc_rate(iter,1) > 0.6
       cc = cc*1.1;
       end;
end; % end of if metflag
                                                         
    if iter > nomit % if we are past burn-in, save the draws
    bsave(iter-nomit,1:k) = bhat';
    ssave(iter-nomit,1) = sige;
    psave(iter-nomit,1) = rho;
    vmean = vmean + vi;

    if mm~= 0
        rsave(iter-nomit,1) = rval;
    end;         
    end;
                    

iter = iter + 1; 
% waitbar(iter/ndraw);         
end; % end of sampling loop
% close(hwait);

time3 = etime(clock,t0);

% compute posterior means and evaluate the log-marginal
vmean = vmean/(ndraw-nomit);
bmean = mean(bsave);
bmean = bmean';
rho = mean(psave);

V = in./vmean;
ys = y.*sqrt(V);
xs = matmul(x,sqrt(V));
Wys = sparse(W)*ys;
Wxs = sparse(W)*xs;
[nobs,nvar] = size(xs);
% if mlog == 1
%     % compute log marginal likelihood for model comparisions
%     if inform_flag == 0
%     mlike = sem_marginal(detval,ys,xs,Wys,Wxs,nobs,nvar,a1,a2);
%     else
%     mlike = sem_marginal2(detval,ys,xs,Wys,Wxs,nobs,nvar,a1,a2,c,TI,sige);
%     end;
% end;
% 
[n nvar] = size(x);
yhat = x*bmean;
y = results.y;
n = length(y);
e = y-yhat;
eD = e - rho*sparse(W)*e;
epe = eD'*eD;

sigu = epe;
sige = sigu/(n-nvar);
ym = y - mean(y);
rsqr1 = sigu;
rsqr2 = ym'*ym;
rsqr = 1.0 - rsqr1/rsqr2; % conventional r-squared
rsqr1 = rsqr1/(n-nvar);
rsqr2 = rsqr2/(n-1.0);

time = etime(clock,timet);

results.meth  = 'sem_g';
results.beta = bmean;
results.rho = rho;
results.sige = sige;
results.bdraw = bsave;
results.pdraw = psave;
results.sdraw = ssave;
results.vmean = vmean;
results.yhat  = yhat;
results.bmean = c;
results.bstd  = sqrt(diag(T));
results.rsqr  = rsqr;
results.rbar = 1 - (rsqr1/rsqr2); % rbar-squared
results.sige = sige;
results.nobs  = n;
results.nvar  = nvar;
results.ndraw = ndraw;
results.nomit = nomit;
results.time  = time;
results.time1 = time1;
results.time2 = time2;
results.time3 = time3;
results.acc = acc_rate;
results.dflag = metflag;
results.nu = nu;
results.d0 = d0;
results.a1 = a1;
results.a2 = a2;
% results.mlike = mlike;
results.tflag = 'plevel';
results.novi = novi_flag;
results.lndet = detval;
results.priorb = inform_flag;
% if mm~= 0
% results.rdraw = rsave;
% results.m     = mm;
% results.k     = kk;
% else
% results.r     = rval;
% results.rdraw = 0;
% end;

case{1} % we do homoscedastic model 
    

% hwait = waitbar(0,'sem\_g: MCMC sampling ...');
t0 = clock;                  
iter = 1;
acc = 0;
          while (iter <= ndraw); % start sampling;
                  
          % update beta   
          xs = x - rho*Wx;
          AI = inv(xs'*xs + sige*TI);
          ys = y - rho*Wy;
          b = xs'*ys + sige*TIc;
          b0 = AI*b;
          bhat = norm_rnd(sige*AI) + b0; 
            
          % update sige
          nu1 = n + 2*nu; 
          e = (y - x*bhat);
          es = e - rho*W*e;
          d1 = 2*d0 + es'*es;
          chi = chis_rnd(1,nu1);
          sige = d1/chi;
              
if metflag == 0
          % update rho using numerical integration          
          rho = draw_rho(detval,y,x,Wy,Wx,V,n,k,rmin,rmax,rho);

else 
          % update rho using metropolis-hastings
           xb = x*bhat;
           rhox = c_sem(rho,y,x,bhat,sige,W,detval,ones(n,1),a1,a2);
          accept = 0;
          rho2 = rho + cc*randn(1,1);
          while accept == 0
           if ((rho2 > rmin) & (rho2 < rmax)); 
           accept = 1;  
           else
           rho2 = rho + cc*randn(1,1);
           end; 
          end; 
           rhoy = c_sem(rho2,y,x,bhat,sige,W,detval,ones(n,1),a1,a2);
          ru = unif_rnd(1,0,1);
          if ((rhoy - rhox) > exp(1)),
          p = 1;
          else,          
          ratio = exp(rhoy-rhox);
          p = min(1,ratio);
          end;
              if (ru < p)
              rho = rho2;
              acc = acc + 1;
              end;
      acc_rate(iter,1) = acc/iter;
      % update cc based on std of rho draws
       if acc_rate(iter,1) < 0.4
       cc = cc/1.1;
       end;
       if acc_rate(iter,1) > 0.6
       cc = cc*1.1;
       end;       
end; % end of if metflag
                                                      
    if iter > nomit % if we are past burn-in, save the draws
    bsave(iter-nomit,1:k) = bhat';
    ssave(iter-nomit,1) = sige;
    psave(iter-nomit,1) = rho;
%     if mm~= 0
%         rsave(iter-nomit,1) = rval;
%     end;         
    end;
                    

iter = iter + 1; 
% waitbar(iter/ndraw);         
end; % end of sampling loop
% close(hwait);

time3 = etime(clock,t0);

% find posterior means and compute log-marginal
bmean = mean(bsave);
bmean = bmean';
rho = mean(psave);
sige = mean(ssave);

[nobs,nvar] = size(x);

% mlike = 0;
% if mlog == 1
%     % compute log marginal likelihood for model comparisions
%     if inform_flag == 0
%     mlike = sem_marginal(detval,y,x,Wy,Wx,nobs,nvar,a1,a2);
%     else
%     mlike = sem_marginal2(detval,y,x,Wy,Wx,nobs,nvar,a1,a2,c,TI,sige);
%     end;
% end;

yhat = x*bmean;
y = results.y;
n = length(y);
e = y-yhat;
eD = e - rho*sparse(W)*e;
epe = eD'*eD;
sigu = epe;
sige = sigu/(nobs-nvar);
ym = y - mean(y);
rsqr1 = sigu;
rsqr2 = ym'*ym;
rsqr = 1.0 - rsqr1/rsqr2; % conventional r-squared
rsqr1 = rsqr1/(nobs-nvar);
rsqr2 = rsqr2/(nobs-1.0);

time = etime(clock,timet);

results.meth  = 'sem_g';
results.beta = bmean;
results.rho = rho;
results.sige = sige;
results.bdraw = bsave;
results.pdraw = psave;
results.sdraw = ssave;
results.vmean = vmean;
results.yhat  = yhat;
results.bmean = c;
results.bstd  = sqrt(diag(T));
results.rsqr  = rsqr;
results.rbar = 1 - (rsqr1/rsqr2); % rbar-squared
results.sige = sige;
results.nobs  = n;
results.nvar  = nvar;
results.ndraw = ndraw;
results.nomit = nomit;
results.time  = time;
results.time1 = time1;
results.time2 = time2;
results.time3 = time3;
results.acc = acc_rate;
results.dflag = metflag;
results.nu = nu;
results.d0 = d0;
results.a1 = a1;
results.a2 = a2;
% results.mlike = mlike;
results.tflag = 'plevel';
results.novi = novi_flag;
results.lndet = detval;
results.priorb = inform_flag;
% if mm~= 0
% results.rdraw = rsave;
% results.m     = mm;
% results.k     = kk;
% else
% results.r     = rval;
% results.rdraw = 0;
% end;


otherwise
error('sem_g: unrecognized novi_flag value on input');
% we should never get here

end; % end of homoscedastic vs. heteroscedastic options

% =========================================================================
% support functions are below
% =========================================================================

function rho = draw_rho(detval,y,x,Wy,Wx,V,n,k,rmin,rmax,rho)
% update rho via univariate numerical integration
% for the heteroscedastic model case

nmk = (n-k)/2;
nrho = length(detval(:,1));
rgrid = rmin+0.01:0.01:rmax-0.01;
ng = length(rgrid);
iota = ones(nrho,1);
rvec = detval(:,1);
epet = zeros(ng,1);
detxt = zeros(ng,1);
for i=1:ng;
xs = x - rgrid(i)*Wx;
xs = matmul(xs,sqrt(V));
ys = y - rgrid(i)*Wy;
ys = ys.*sqrt(V);
bs = (xs'*xs)\(xs'*ys);
e = ys - xs*bs;
epet(i,1) = e'*e;
% detxt(i,1) = det(xs'*xs);
end;

% interpolate a finer grid
epe = interp1(rgrid',epet,detval(:,1),'spline');
% detx = interp1(rgrid',detxt,detval(:,1),'spline');

den = detval(:,2) - nmk*log(epe);
adj = max(den);
den = den - adj;
den = exp(den);

n = length(den);
y = detval(:,1);
x = den;

% trapezoid rule
isum = sum((y(2:n,1) + y(1:n-1,1)).*(x(2:n,1) - x(1:n-1,1))/2);

z = abs(x/isum);
den = cumsum(z);

rnd = unif_rnd(1,0,1)*sum(z);
ind = find(den <= rnd);
idraw = max(ind);
if (idraw > 0 & idraw < nrho)
rho = detval(idraw,1);
end;


function cout = c_sem(rho,y,x,b,sige,W,detval,vi,a1,a2);
% PURPOSE: evaluate the conditional distribution of rho given sige
%  spatial autoregressive model using sparse matrix algorithms
% ---------------------------------------------------
%  USAGE:cout = c_sar(rho,y,x,b,sige,W,detval,a1,a2)
%  where:  rho  = spatial autoregressive parameter
%          y    = dependent variable vector
%          W    = spatial weight matrix
%        detval = an (ngrid,2) matrix of values for det(I-rho*W) 
%                 over a grid of rho values 
%                 detval(:,1) = determinant values
%                 detval(:,2) = associated rho values
%          sige = sige value
%          a1    = (optional) prior parameter for rho
%          a2    = (optional) prior parameter for rho
% ---------------------------------------------------
%  RETURNS: a conditional used in Metropolis-Hastings sampling
%  NOTE: called only by sar_g
%  --------------------------------------------------
%  SEE ALSO: sar_g, c_far, c_sac, c_sem
% ---------------------------------------------------

gsize = detval(2,1) - detval(1,1);
% Note these are actually log detvalues
i1 = find(detval(:,1) <= rho + gsize);
i2 = find(detval(:,1) <= rho - gsize);
i1 = max(i1);
i2 = max(i2);
index = round((i1+i2)/2);
if isempty(index)
index = 1;
end;
detm = detval(index,2); 

[n,k] = size(x);
nmk = (n-k)/2;
z = speye(n) - rho*sparse(W);
xs = z*x;
ys = z*y;

% detx = 0.5*log(det(xs'*xs));

n = length(y);
e = ys - xs*b;
ev = e.*sqrt(vi);
epe = nmk*log(ev'*ev);

cout =   detm - epe;



function [nu,d0,rval,mm,kk,rho,sige,rmin,rmax,detval,ldetflag,eflag,order,iter,novi_flag,c,T,cc,metflag,a1,a2,inform_flag,mlog] = sem_parse(prior,k)
% PURPOSE: parses input arguments for far, far_g models
% ---------------------------------------------------
%  USAGE: [nu,d0,rval,mm,kk,rho,sige,rmin,rmax,detval, ...
%          ldetflag,eflag,order,iter,novi_flag,c,T,prior_beta,cc,metflag,a1,a2,inform_flag] = 
%                           sem_parse(prior,k)
% where info contains the structure variable with inputs 
% and the outputs are either user-inputs or default values
% ---------------------------------------------------

% set defaults

eflag = 1;     % default to not computing eigenvalues
ldetflag = 1;  % default to 1999 Pace and Barry MC determinant approx
order = 50;    % there are parameters used by the MC det approx
iter = 30;     % defaults based on Pace and Barry recommendation
rmin = -1;     % use -1,1 rho interval as default
rmax = 1;
detval = 0;    % just a flag
rho = 0.5;
sige = 1.0;
rval = 4;
mm = 0;
kk = 0;
nu = 0;
d0 = 0;
a1 = 1.0;
a2 = 1.0;
c = zeros(k,1);   % diffuse prior for beta
T = eye(k)*1e+12;
novi_flag = 0;    % default is do vi-estimates
metflag = 1;      % default to Metropolis-Hasting sampling
cc = 0.2;         % initial tuning parameter for M-H sampling
inform_flag = 0;  % flag for diffuse prior on beta
mlog = 1;

fields = fieldnames(prior);
nf = length(fields);
if nf > 0
 for i=1:nf
    if strcmp(fields{i},'nu')
        nu = prior.nu;
    elseif strcmp(fields{i},'d0')
        d0 = prior.d0;  
    elseif strcmp(fields{i},'rval')
        rval = prior.rval; 
    elseif strcmp(fields{i},'eigs')
        eflag = prior.eigs;
    elseif strcmp(fields{i},'dflag')
        metflag = prior.dflag;
    elseif strcmp(fields{i},'a1')
       a1 = prior.a1; 
    elseif strcmp(fields{i},'a2')
       a2 = prior.a2; 
    elseif strcmp(fields{i},'m')
        mm = prior.m;
        kk = prior.k;
        rval = gamm_rnd(1,1,mm,kk);    % initial value for rval   
    elseif strcmp(fields{i},'beta')
        c = prior.beta; inform_flag = 1; % flag for informative prior on beta
    elseif strcmp(fields{i},'bcov')
        T = prior.bcov; inform_flag = 1; % flag for informative prior on beta
    elseif strcmp(fields{i},'rmin')
        rmin = prior.rmin;  
    elseif strcmp(fields{i},'rmax')
        rmax = prior.rmax; 
    elseif strcmp(fields{i},'lndet')
    detval = prior.lndet;
    ldetflag = -1;
    rmin = detval(1,1);
    nr = length(detval);
    rmax = detval(nr,1);
    elseif strcmp(fields{i},'lflag')
        tst = prior.lflag;
        if tst == 0,
        ldetflag = 0; 
        elseif tst == 1,
        ldetflag = 1; 
        elseif tst == 2,
        ldetflag = 2; 
        else
        error('sem_g: unrecognizable lflag value on input');
        end;
    elseif strcmp(fields{i},'order')
        order = prior.order;  
    elseif strcmp(fields{i},'iter')
        iter = prior.iter; 
    elseif strcmp(fields{i},'novi')
        novi_flag = prior.novi;
    elseif strcmp(fields{i},'mlog')
        mlog = prior.mlog;
    end;
 end;
 
else, % the user has input a blank info structure
      % so we use the defaults
end; 


function [rmin,rmax,time2] = sem_eigs(eflag,W,rmin,rmax,n);
% PURPOSE: compute the eigenvalues for the weight matrix
% ---------------------------------------------------
%  USAGE: [rmin,rmax,time2] = far_eigs(eflag,W,rmin,rmax,W)
% where eflag is an input flag, W is the weight matrix
%       rmin,rmax may be used as default outputs
% and the outputs are either user-inputs or default values
% ---------------------------------------------------


if eflag == 0
t0 = clock;
opt.tol = 1e-3; opt.disp = 0;
lambda = eigs(sparse(W),speye(n),1,'SR',opt);  
rmin = 1/lambda;   
rmax = 1;
time2 = etime(clock,t0);
else
time2 = 0;
end;


function [detval,time1] = sem_lndet(ldetflag,W,rmin,rmax,detval,order,iter);
% PURPOSE: compute the log determinant |I_n - rho*W|
% using the user-selected (or default) method
% ---------------------------------------------------
%  USAGE: detval = far_lndet(lflag,W,rmin,rmax)
% where eflag,rmin,rmax,W contains input flags 
% and the outputs are either user-inputs or default values
% ---------------------------------------------------


% do lndet approximation calculations if needed
if ldetflag == 0 % no approximation
t0 = clock;    
out = lndetfull(W,rmin,rmax);
time1 = etime(clock,t0);
tt=rmin:.001:rmax; % interpolate a finer grid
outi = interp1(out.rho,out.lndet,tt','spline');
detval = [tt' outi];
    
elseif ldetflag == 1 % use Pace and Barry, 1999 MC approximation

t0 = clock;    
out = lndetmc(order,iter,W,rmin,rmax);
time1 = etime(clock,t0);
results.limit = [out.rho out.lo95 out.lndet out.up95];
tt=rmin:.001:rmax; % interpolate a finer grid
outi = interp1(out.rho,out.lndet,tt','spline');
detval = [tt' outi];

elseif ldetflag == 2 % use Pace and Barry, 1998 spline interpolation

t0 = clock;
out = lndetint(W,rmin,rmax);
time1 = etime(clock,t0);
tt=rmin:.001:rmax; % interpolate a finer grid
outi = interp1(out.rho,out.lndet,tt','spline');
detval = [tt' outi];

elseif ldetflag == -1 % the user fed down a detval matrix
    time1 = 0;
        % check to see if this is right
        if detval == 0
            error('sem_g: wrgon lndet input argument');
        end;
        [n1,n2] = size(detval);
        if n2 ~= 2
            error('sem_g: wrong sized lndet input argument');
        elseif n1 == 1
            error('sem_g: wrong sized lndet input argument');
        end;          
end;

function  out = sem_marginal(detval,y,x,Wy,Wx,nobs,nvar,a1,a2)
% PURPOSE: returns a vector of the log-marginal over a grid of rho-values
% -------------------------------------------------------------------------
% USAGE: out = sem_marginal(detval,y,x,Wy,Wx,nobs,nvar,a1,a2)
% where:       detval = an ngrid x 2 matrix with rho-values and lndet values
%                  y = y-vector
%                  x = x-matrix
%                 Wy = W*y-vector
%                 Wx = W*x-matrix
%               nobs = # of observations
%               nvar = # of explanatory variables
%                 a1 = parameter for beta prior on rho
%                 a2 = parameter for beta prior on rho
% -------------------------------------------------------------------------
% RETURNS: out = a structure variable
%        out.log = log marginal, a vector the length of detval
%        out.lik = concentrated log-likelihood vector the length of detval
% -------------------------------------------------------------------------
% NOTES: works only for homoscedastic SEM model
% we must feed in ys = sqrt(V)*y, xs = sqrt(V)*X 
% as well as logdetx = log(xs'*xs) for heteroscedastic model
% -------------------------------------------------------------------------

% written by:
% James P. LeSage, 7/2003
% Dept of Economics
% University of Toledo
% 2801 W. Bancroft St,
% Toledo, OH 43606
% jlesage@spatial-econometrics.com

nmk = (nobs-nvar)/2;

nrho = length(detval(:,1));
iota = ones(nrho,1);
rvec = detval(:,1);
epe = zeros(nrho,1);
rgrid = detval(1,1)+0.001:0.1:detval(end,1)-0.001;
rgrid = rgrid';

epetmp = zeros(length(rgrid),1);
detxtmp = zeros(length(rgrid),1);
for i=1:length(rgrid);
xs = x - rgrid(i,1)*Wx;
ys = y - rgrid(i,1)*Wy;
bs = (xs'*xs)\(xs'*ys);
e = ys - xs*bs;
epetmp(i,1) = e'*e;
detxtmp(i,1) = det(xs'*xs);
end;

% spline interpolate epetmp
tt=rvec; % interpolate a finer grid
epe = interp1(rgrid,epetmp,rvec,'spline'); 
detx = interp1(rgrid,detxtmp,rvec,'spline'); 


bprior = beta_prior(detval(:,1),a1,a2);
% C is a constant of integration that can vary with nvars, so for model
% comparisions involving different nvars we need to include this
C = log(bprior) + gammaln(nmk) - nmk*log(2*pi) ;
den = detval(:,2) - 0.5*log(detx) - nmk*log(epe);
den = real(den);
out = den + C;

function  out = sem_marginal2(detval,y,x,Wy,Wx,nobs,nvar,a1,a2,c,TI,sige)
% PURPOSE: returns a vector of the log-marginal over a grid of rho-values
%          for the case of an informative prior on beta
% -------------------------------------------------------------------------
% USAGE: out = sem_marginal2(detval,y,x,Wy,Wx,nobs,nvar,a1,a2,c,TI,sige)
% where:       detval = an ngrid x 2 matrix with rho-values and lndet values
%                  y = y-vector
%                  x = x-matrix
%                 Wy = W*y-vector
%                 Wx = W*x-matrix
%               nobs = # of observations
%               nvar = # of explanatory variables
%                 a1 = parameter for beta prior on rho
%                 a2 = parameter for beta prior on rho
% -------------------------------------------------------------------------
% RETURNS: out = a structure variable
%        out.log = log marginal, a vector the length of detval
%        out.lik = concentrated log-likelihood vector the length of detval
% -------------------------------------------------------------------------
% NOTES: works only for homoscedastic SEM model
% we must feed in ys = sqrt(V)*y, xs = sqrt(V)*X 
% as well as logdetx = log(xs'*xs) for heteroscedastic model
% -------------------------------------------------------------------------

% written by:
% James P. LeSage, 7/2003
% Dept of Economics
% University of Toledo
% 2801 W. Bancroft St,
% Toledo, OH 43606
% jlesage@spatial-econometrics.com

nmk = (nobs-nvar)/2;

nrho = length(detval(:,1));
iota = ones(nrho,1);
rvec = detval(:,1);
epe = zeros(nrho,1);

rgrid = detval(1,1)+0.001:0.1:detval(end,1)-0.001;
rgrid = rgrid';
epetmp = zeros(length(rgrid),1);
detxtmp = zeros(length(rgrid),1);
Q1 = zeros(length(rgrid),1);
Q2 = zeros(length(rgrid),1);
sTI = sige*TI;

for i=1:length(rgrid);
xs = x - rgrid(i,1)*Wx;
ys = y - rgrid(i,1)*Wy;
bs = (xs'*xs)\(xs'*ys);
beta = (xs'*xs + sTI)\(xs'*ys + sTI*c);
e = ys - xs*bs;
epetmp(i,1) = e'*e;
detxtmp(i,1) = det(xs'*xs);
Q1(i,1) = (c - beta)'*sTI*(c - beta);
Q2(i,1) = (bs - beta)'*(xs'*xs)*(bs - beta);
end;

% spline interpolate epetmp
tt=rvec; % interpolate a finer grid
epe = interp1(rgrid,epetmp,rvec,'spline'); 
detx = interp1(rgrid,detxtmp,rvec,'spline'); 
Q1 = interp1(rgrid,Q1,rvec,'spline'); 
Q2 = interp1(rgrid,Q2,rvec,'spline'); 
bprior = beta_prior(detval(:,1),a1,a2);
% C is a constant of integration that can vary with nvars, so for model
% comparisions involving different nvars we need to include this
C = log(bprior) + gammaln(nmk) - nmk*log(2*pi) ;
den = detval(:,2) - 0.5*log(detx) - nmk*log(epe + Q1 + Q2);
den = real(den);
out = den + C;


function rho = olddraw_rho(detval,y,x,Wy,Wx,V,n,k,rmin,rmax,rho)
% update rho via univariate numerical integration
% for the heteroscedastic model case

nmk = (n-k)/2;
nrho = length(detval(:,1));
iota = ones(nrho,1);
rvec = detval(:,1);
epe = zeros(nrho,1);

for i=1:nrho;
xs = x - rvec(i,1)*Wx;
xs = matmul(xs,sqrt(V));
ys = y - rvec(i,1)*Wy;
ys = ys.*sqrt(V);
bs = (xs'*xs)\(xs'*ys);
e = ys - xs*bs;
epe(i,1) = e'*e;
end;

den = detval(:,2) - nmk*log(epe);
adj = max(den);
den = den - adj;
den = exp(den);


n = length(den);
y = detval(:,1);
x = den;

% trapezoid rule
isum = sum((y(2:n,1) + y(1:n-1,1)).*(x(2:n,1) - x(1:n-1,1))/2);

z = abs(x/isum);
den = cumsum(z);

rnd = unif_rnd(1,0,1)*sum(z);
ind = find(den <= rnd);
idraw = max(ind);
if (idraw > 0 & idraw < nrho)
rho = detval(idraw,1);
end;
