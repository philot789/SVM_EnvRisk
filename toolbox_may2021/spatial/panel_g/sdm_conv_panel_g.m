function results = sdm_conv_panel_g(y,x,Wmatrices,N,T,ndraw,nomit,prior)
% PURPOSE: MCMC estimates for SDM convex combination of W model
%          using a convex combination of m different W-matrices
%          y = rho*Wc*y + X*beta + W1*X*c1 + W2*X*c2 ... + e, e = N(0,sige*I_n), 
%          Wc = g1*W1 + g2*W2 + ... + (1-g1-g2- ... -gm)*Wm
%          no prior for beta, c1,c2, ..., sigma, rho, gamma
%-------------------------------------------------------------
% USAGE: results = sdm_conv_panel_g(y,x,Wmatrices,N,T,ndraw,nomit,prior)
% where: y = dependent variable vector (N*T x 1)
%        x = independent variables matrix (N*T x nvar), 
% N.B. the function creates [x W1*x W2*x ...] matrix for you
%        Wmatrices = (N,m*N)
%        e.g., Wmatrices = [W1 W2 ... Wm]
%        where each W1, W2, ... Wm are (N x N) row-normalized weight matrices
%    Or, you can enter a set of m LARGE matrices that are N*T x N*T, which
%    allows for different W-matrices for each time period of the panel
%    ndraw = # of draws (use lots of draws, say 25,000 to 50,000
%    nomit = # of initial draws omitted for burn-in  (probably around 5,000
%    prior = a structure variable with:
%       prior.model = 0 pooled model without fixed effects (default, x may contain an intercept)
%                   = 1 spatial fixed effects (x may not contain an intercept)
%                   = 2 time period fixed effects (x may not contain an intercept)
%                   = 3 spatial and time period fixed effects (x may not contain an intercept)
%       prior.fe    = report fixed effects and their t-values in prt_panel() printout
%                     (default=0=not reported; info.fe=1=report) 
%       prior.thin  = a thinning parameter for use in analyzing
%                          posterior distributions, default = 1 (no thinning of draws)
%                          recommended value for ndraw > 20,000 is 10
%            NOTE: thin is NOT used to determine how many times to MH-MC sample the
%            log-posterior using Monte Carlo integration, which is sampled (ndraw-nomit) times
%       prior.plt_flag = 0,1 default = 0 no plotting of MCMC draws
%                      = 1 will plot the MCMC draws during estimation
%    prior.T2 = pre-calculated 2nd order Taylor series traces fed to the function
%    prior.T3 = pre-calculated 2nd order Taylor series traces fed to the function
%    prior.T4 = pre-calculated 2nd order Taylor series traces fed to the function
% these will be calculated by this function, but can be pre-calculated and
% fed to the function in the case of Monte Carlo studies to save some time
%-------------------------------------------------------------
% RETURNS:  a structure:
%    results.meth  = 'sdmp_conv_g'    if prior.model=0, no fixed effects
%                  = 'sdmsfe_conv_g'  if prior.model=1, spatial fixed effects
%                  = 'sdmtfe_conv_g'  if prior.model=2, time fixed effects
%                  = 'sdmstfe_conv_g' if prior.model=3, both space and time fixed effects
%          results.beta     = posterior mean of bhat based on draws, includes c1,c2, ... above
%          results.rho      = posterior mean of rho based on draws
%          results.sige     = posterior mean of sige based on draws
%          results.gamma    = m x 1 vector of posterior means for g1,g2, ... gm
%                             where m is the number of weight matrices used on input
%          results.nobs   = # of observations
%          results.nvar   = # of variables in x-matrix
%          results.p      = # of variables in x-matrix
%          results.beta_std = std deviation of beta draws, includes c1,c2, ... above
%          results.sige_std = std deviation of sige draws
%          results.rho_std  = std deviation of rho draws
%          results.g_std    = m x 1 vector of posterior std deviations
%          results.sigma    = posterior mean of sige based on (e'*e)/(n-k)
%          results.bdraw    = bhat draws (1:thin:ndraw-nomit x nvar)
%          results.pdraw    = rho  draws (1:thin:ndraw-nomit x 1)
%          results.sdraw    = sige draws (1:thin:ndraw-nomit x 1)
%          results.gdraw    = gamma draws (1:thin:ndraw-nomit x m)
%          results.direct   = nvar x 5 matrix with direct effect, t-stat, t-prob, lower05, upper95
%          results.indirect = nvar x 5 matrix with indirect effect, t-stat, t-prob, lower05, upper95
%          results.total    = nvar x 5 matrix with total effect, t-stat, t-prob, lower05, upper95
%          results.ndraw  = # of draws
%          results.nomit  = # of initial draws omitted
%          results.thin     = thinning value from input
%          results.total_draws    = a (1:thin:ndraw-nomit,p) total x-impacts
%          results.direct_draws   = a (1:thin:ndraw-nomit,p) direct x-impacts
%          results.indirect_draws = a (1:thin:ndraw-nomit,p) indirect x-impacts
%          results.lmarginal= a scalar log-marginal likelihood, from MH-MC
%                             (Metropolis-Hastings Monte Carlo) integration of the log-posterior
%          results.drawpost = a vector of log-posterior draws (ndraw-nomit)x1
%          results.betapost = a matrix of log-posterior draws for beta (ndraw-nomit x k)
%          results.sigpost  = a vector of log-posterior draws for sige (ndraw-nomit x 1)
%          results.rho_gamma = (ndraw-nomit) x m+1 contains
%                              [rhodraws gammdraws] used to evaluate the log-posterior
%          results.rho_mode   = modal value of rho from the log-posterior
%          results.gamma_mode = modal values of gamma vector from the log-posterior
%          results.beta_mode  = modal values of beta vector from the log-posterior
%          results.sig_mode   = modal values of beta vector from the log-posterior
%          results.logC_sar   = constants associated with log-marginal
%                               logC_sar = gammaln(dof) - dof*log(2*pi)  -0.5*lndetx_sar
%                               dof = (n - m)/2; lndetx_sar = log(det(xp'*x));
%          results.logm_profile = a profile of the log-margainal over [rho_gamma isum];
%                                 where: isum = exp(logp -adj) + adj; adj = max(logp)
%          results.y      = y-vector from input (nobs x 1)
%          results.yhat   = mean of posterior predicted (nobs x 1)
%          results.resid  = residuals, based on posterior means
%          results.rsqr   = r-squared based on posterior means
%          results.rbar   = adjusted r-squared
%          results.taylor = taylor from input;
%          results.sampling_time = time for MCMC sampling
%          results.effects_time  = time to calculate effects estimates
%          results.trace_time    = time to calculate traces for taylor/chebyshev approximations
%          results.rmax   =  1  
%          results.rmin   = -1       
%          results.cflag  = 0 for intercept term, 1 for no intercept term
%          results.T2 = 2nd order Taylor series traces
%          results.T3 = 3rd order Taylor series traces
%          results.T4 = 4th order Taylor series traces
%          results.sfe   = spatial fixed effects (if info.model=1 or 3)
%          results.tfe   = time period fixed effects (if info.model=2 or 3)
%          results.tsfe  = t-values spatial fixed effects (if info.model=1 or 3)
%          results.ttfe  = t-values time period fixed effects (if info.model=2 or 3)
%          results.con   = intercept 
%          results.tcon  = t-value intercept
% --------------------------------------------------------------
% NOTES: - the intercept term (if you have one)
%          must be in the first column of the matrix x
% --------------------------------------------------------------
% SEE ALSO: (sdm_conv_panel_gd.m demos) 
% --------------------------------------------------------------
% REFERENCES: Debarsy and LeSage (2018) 
% Flexible dependence modeling using convex combinations of different
% types of connectivity structures, Regional Science & Urban Economics,
% Volume 69, pp. 46-68.
% Debarsy and LeSage (2020) 
% Bayesian model averaging for spatial autoregressive models
% based on convex combinations of different types of connectivity matrices
% Journal of Businesss & Economic Statistics
%----------------------------------------------------------------

% written by:
% James P. LeSage, last updated 11/2020
% Dept of Finance & Economics
% Texas State University-San Marcos
% 601 University Drive
% San Marcos, TX 78666
% james.lesage@txstate.edu

% error checking on inputs
[n junk] = size(y);
[nc, k] = size(x);
if (nc ~= n)
       error('sdm_conv_panel_g: wrong sized x-matrix');
end

% check if the user handled the intercept term okay
    n = length(y);
    if sum(x(:,1)) ~= n
        tst = sum(x); % we may have no intercept term
        ind = find(tst == n); % we do have an intercept term
        if length(ind) > 0
            error('sar_panel_FE_g: intercept term must be in first column of the x-matrix');
        elseif length(ind) == 0 % case of no intercept term
            cflag = 0;
            pp = size(x,2);
        end
    elseif sum(x(:,1)) == n % we have an intercept in the right place
        cflag = 1;
        pp = size(x,2)-1;
    end
    results.cflag = cflag;
    results.p = pp;
      


if ndraw <= 1000
       error('sdm_conv_panel_g: ndraw<=1000, increase ndraw to at least 10000');
end

if nargin == 7
            error('sdm_conv_panel_g: you must add prior.model=0,1,2,3');

elseif nargin == 8

fields = fieldnames(prior);
nf = length(fields);
if nf > 0
    for i=1:nf
        if strcmp(fields{i},'model') model = prior.model;
        elseif strcmp(fields{i},'fe') fe = prior.fe;
         elseif strcmp(fields{i},'thin') thin = prior.thin;
           
        end
    end
end

    [n1,n2] = size(Wmatrices);
    m = n2/n1;
    
    if n1 == N*T
        Wlarge = 1;
    else
        Wlarge = 0;
    end
    
    if m < 2
        error('sdm_conv_panel_g: only one W-matrix');
    end
    

     [rho,sige,rmin,rmax,gamma,thin,ccmin,ccmax,ddmin,ddmax,T2,T3,T4,tr_flag,plt_flag,fe] = sar_parse(prior,m);


% check the thinning parameter
 eff = (ndraw -nomit)/thin;
if eff < 999
    warning('sdm_conv_panel_g: < 1000 draws after thining');
end

else
    error('sdm_conv_panel_g: wrong # of input arguments to sar_conv_panel_g');
end

if Wlarge == 0
    
    begi = 1;
    endi = N;
    Wout = [];
    
    for ii=1:m
        Wout = [Wout kron(eye(T),Wmatrices(:,begi:endi))];
        begi = begi + N;
        endi = endi + N;
    end
Wmatrices = Wout;
    
end


[ywith,xwith,meanny,meannx,meanty,meantx]=demean(y,x,N,T,model);

prior.plt_flag = plt_flag;
prior.thin = thin;
prior.fe = fe;
prior.rmin = rmin;
prior.rmax = rmax;

result0 = sdm_conv_g(ywith,xwith,Wmatrices,ndraw,nomit,prior);

results = result0;
results.m = m;

results.fe = fe;
thin = round(thin);
results.thin = thin;
results.ndraw = ndraw;
results.nomit = nomit;
results.rmin = rmin;
results.rmax = rmax;
results.ccmin = ccmin;
results.ccmax = ccmax;
results.ddmin = ddmin;
results.ddmax = ddmax;

results.model = model;
results.N = N;
results.T = T;
results.nobs  = n;
results.nvar  = k;
nvar = k;
results.y = y;


if model==0
    results.meth='psdm_conv_g';
elseif model==1
    results.meth='sdmsfe_conv_g';
elseif model==2
    results.meth='sdmtfe_conv_g';
elseif model==3
    results.meth='sdmstfe_conv_g';
else
    error('sdm_conv_panel_g: wrong input number of prior.model');
end

% =========================================================
% calculate SFE and TFE based on posterior mean estimates
[n,nvar] = size(xwith);

en=ones(T,1);
et=ones(N,1);
ent=ones(n,1);

Wx = x;
begi = 1;
endi = n;
for ii=1:m
    if results.cflag == 0
    for jj=1:nvar
        Wx = [Wx Wmatrices(:,begi:endi)*x(:,jj)];
    end
    elseif results.cflag == 1
     for jj=2:nvar
        Wx = [Wx Wmatrices(:,begi:endi)*x(:,jj)];
    end
    end       
    begi = begi + n;
    endi = endi + n;
end

nvarw = size(Wx,2);

gmean = mean(results.gdraw);

begi = 1;
endi = n;
Wch = sparse(n,n);
for ii=1:m
    Wch = Wch + gmean(1,ii)*sparse(Wmatrices(:,begi:endi));
    begi = begi + n;
    endi = endi + n;
end

Wy = sparse(Wch)*y;
        
if (model==1 | model==3);
meanny=zeros(N,1);
meannwy=zeros(N,1);
meannx=zeros(N,nvarw);
for i=1:N
    ym=zeros(T,1);
    wym=zeros(T,1);
    xm=zeros(T,nvarw);
    for t=1:T
        ym(t)=y(i+(t-1)*N,1);
        wym(t)=Wy(i+(t-1)*N,1);
        xm(t,:)=Wx(i+(t-1)*N,:);
    end
    meanny(i)=mean(ym);
    meannwy(i)=mean(wym);
    meannx(i,:)=mean(xm);
end
clear ym wym xm;
end % if statement

if ( model==2 | model==3)
meanty=zeros(T,1);
meantwy=zeros(T,1);
meantx=zeros(T,nvarw);
for i=1:T
    t1=1+(i-1)*N;t2=i*N;
    ym=y([t1:t2],1);
    wym=Wy([t1:t2],1);
    xm=Wx([t1:t2],:);
    meanty(i)=mean(ym);
    meantwy(i)=mean(wym);
    meantx(i,:)=mean(xm);
end
clear ym wym xm;
end % if statement

rho = mean(results.pdraw);
bhat = mean(results.bdraw)';
gmean = mean(results.gdraw)';
sige = mean(results.sdraw);
nobs = N*T;


              if model==1
                  intercept=mean(y)-mean(Wy)*rho-mean(Wx)*bhat;
                  con=intercept;
                  sfe=meanny-meannwy*rho-meannx*bhat-kron(et,intercept);
                  xhat=Wx*bhat+kron(en,sfe)+kron(ent,intercept);
%                   tmp = [1
%                          -rho*gmean];
%                   Wys = Wy2*tmp;
%                   e =Wys - xhat;
%                   sige = (e'*e)/N*T;
                  tsfe=sfe./sqrt(sige/T*ones(N,1)+diag(sige*meannx*(Wx'*Wx)*meannx'));
                  tcon=con/sqrt(sige/nobs+sige*mean(Wx)*(Wx'*Wx)*mean(Wx)');
                  results.con = con;
                  results.sfe = sfe;
                  results.tsfe = tsfe;
                  results.tcon = tcon;
                  results.sige = sige;

              elseif model==2
                  intercept=mean(y)-mean(Wy)*rho-mean(Wx)*bhat;
                  con=intercept;
                  tfe=meanty-meantwy*rho-meantx*bhat-kron(en,intercept);
                  xhat=Wx*bhat+kron(tfe,et)+kron(ent,intercept);
%                    tmp = [1
%                          -rho*gmean];
%                   Wys = Wy2*tmp;
%                   e =Wys - xhat;
%                   sige = (e'*e)/N*T;
                  ttfe=tfe./sqrt(sige/N*ones(T,1)+diag(sige*meantx*(Wx'*Wx)*meantx'));
                  tcon=con/sqrt(sige/nobs+sige*mean(Wx)*(Wx'*Wx)*mean(Wx)');
                  results.con = con;
                  results.tfe = tfe;
                  results.ttfe = ttfe;
                  results.tcon = tcon;
%                   results.sige = sige;

              elseif model==3
                  intercept=mean(y)-mean(Wy)*rho-mean(Wx)*bhat;
                  con=intercept;
                  sfe=meanny-meannwy*rho-meannx*bhat-kron(et,intercept);
                  tfe=meanty-meantwy*rho-meantx*bhat-kron(en,intercept);
                  xhat=Wx*bhat+kron(en,sfe)+kron(tfe,et)+kron(ent,intercept);
%                     tmp = [1
%                          -rho*gmean];
%                   Wys = Wy2*tmp;
%                   e =Wys - xhat;
%                   sige = (e'*e)/N*T;                
                  tsfe=sfe./sqrt(sige/T*ones(N,1)+diag(sige*meannx*(Wx'*Wx)*meannx'));
                  ttfe=tfe./sqrt(sige/N*ones(T,1)+diag(sige*meantx*(Wx'*Wx)*meantx'));
                  tcon=con/sqrt(sige/nobs+sige*mean(Wx)*(Wx'*Wx)*mean(Wx)');
                  results.con = con;
                  results.sfe = sfe;
                  results.tfe = tfe;
                  results.tsfe = tsfe;
                  results.ttfe = ttfe;
                  results.tcon = tcon;
%                   results.sige = sige;

              else
                  xhat=Wx*bhat;
%                   tmp = [1
%                          -rho*gmean];
%                   Wys = Wy2*tmp;
%                   e =Wys - xhat;
%                   results.sige = (e'*e)/N*T;                
              end
 
              
results.resid = y - rho*Wy - xhat; 
yme=y-mean(y);
rsqr2=yme'*yme;
rsqr1 = results.resid'*results.resid;
results.rsqr=1.0-rsqr1/rsqr2; %rsquared

 % calculate correlation-squared estimate
yhat=zeros(nobs,1);
ywithhat=zeros(nobs,1);
for t=1:T
    t1=1+(t-1)*N;t2=t*N;
    ywithhat(t1:t2,1)=(speye(N) - rho*Wch(t1:t2,t1:t2))\Wx(t1:t2,:)*bhat;
    yhat(t1:t2,1)=(speye(N) - rho*Wch(t1:t2,t1:t2))\xhat(t1:t2,1);
end
res1=ywith-mean(ywith);
res2=ywithhat-mean(ywith);
rsq1=res1'*res2;
rsq2=res1'*res1;
rsq3=res2'*res2;
results.corr2=rsq1^2/(rsq2*rsq3); %corr2
results.yhat=yhat;
           

% =========================================================================
% support functions below
% =========================================================================
% 

function [rho,sige,rmin,rmax,gamma,thin,ccmin,ccmax,ddmin,ddmax,T2,T3,T4,tr_flag,plt_flag,fe] = sar_parse(prior,m)
% PURPOSE: parses input arguments for sar_conv_g models
% ---------------------------------------------------
%  USAGE: [rval,rho,sige,rmin,rmax,novi_flag] =  sar_parse(prior,k)
% returns values set by user or default values 
% ---------------------------------------------------

% written by:
% James P. LeSage, last updated 1/2018
% Dept of Finance & Economics
% Texas State University-San Marcos
% 601 University Drive
% San Marcos, TX 78666
% james.lesage@txstate.edu


% set defaults
gamma = ones(m,1)*(1/m);
rmin = -1;     % use -1,1 rho interval as default
rmax = 1;
rho = 0.7; % starting values
sige = 1;
thin = 1; % default to no thinning
ccmin = 0.4;
ccmax = 0.6;
ddmin = 0.1;
ddmax = 0.4;
tr_flag = 0;
plt_flag = 0;
T2 = [];
T3 = [];
T4 = [];
fe = 0; % default to NOT printing fixed effects estimates

fields = fieldnames(prior);
nf = length(fields);
if nf > 0
    for i=1:nf
        if strcmp(fields{i},'thin')
            thin = prior.thin;
        elseif strcmp(fields{i},'fe')
            fe = prior.fe;
        elseif strcmp(fields{i},'ccmin')
            ccmin = prior.ccmin;
         elseif strcmp(fields{i},'ccmax')
            ccmax = prior.ccmax;
         elseif strcmp(fields{i},'ddmin')
            ddmin = prior.ddmin;
         elseif strcmp(fields{i},'ddmax')
            ddmax = prior.ddmax;          
        elseif strcmp(fields{i},'T2') && strcmp(fields{i},'T3') && strcmp(fields{i},'T4')
            T2 = prior.T2;    
            T3 = prior.T3;    
            T4 = prior.T4;    
            tr_flag = 1;
        elseif strcmp(fields{i},'plt_flag')
            plt_flag = prior.plt_flag;          
        end
    end
    
    
else % the user has input a blank info structure
    % so we use the defaults
end

function results = sdm_conv_g(y,x,Wmatrices,ndraw,nomit,prior)
% PURPOSE: Bayesian estimates of the spatial autoregressive model
%          using a convex combination of m different W-matrices
%          y = rho*Wc*y + X*beta + W1*X*c1 + W2*X*c2 + ... + e, e = N(0,sige*I_n), 
%          Wc = g1*W1 + g2*W2 + ... + (1-g1-g2- ... -gm)*Wm
%          no priors for beta 
%          no priors for sige
%          uniform (-1,1) prior for rho
%          uniform (0,1) prior for g1,g2, ... gm
%-------------------------------------------------------------
% USAGE: results = sar_conv_g(y,x,Wmatrices,ndraw,nomit,prior)
% where: y = dependent variable vector (nobs x 1)
%        x = independent variables matrix (nobs x nvar), 
% N.B. the function creates [x W1*x W2*x ...] matrix for you
%        Wmatrices = (nobs,m*nobs)
%        e.g., Wmatrices = [W1 W2 ... Wm]
%        where each W1, W2, ... Wm are (nobs x nobs) row-normalized weight matrices
%    ndraw = # of draws (use lots of draws, say 25,000 to 50,000
%    nomit = # of initial draws omitted for burn-in  (probably around 5,000
%    prior.plt = 1 for plotting of MCMC draws, 0 for no plots, default = 0
%    prior = a structure variable with:
%            prior.thin  = a thinning parameter for use in analyzing
%                          posterior distributions, default = 1 (no thinning of draws)
%                          recommended value for ndraw > 20,000 is 10
%                          default = 1
%            NOTE: thin is NOT used to determine how many times to MH-MC sample the
%            log-posterior using Monte Carlo integration, which is sampled
%            (ndraw-nomit) times
%    prior.T2 = pre-calculated 2nd order Taylor series traces fed to the function
%    prior.T3 = pre-calculated 2nd order Taylor series traces fed to the function
%    prior.T4 = pre-calculated 2nd order Taylor series traces fed to the function
% these will be calculated by this function, but can be pre-calculated and
% fed to the function in the case of Monte Carlo studies to save some time
%-------------------------------------------------------------
% RETURNS:  a structure:
%          results.meth     = 'sar_conv_g'
%          results.beta     = posterior mean of bhat based on draws
%          results.rho      = posterior mean of rho based on draws
%          results.sige     = posterior mean of sige based on draws
%          results.gamma    = m x 1 vector of posterior means for g1,g2, ... gm
%                             where m is the number of weight matrices used on input
%          results.nobs   = # of observations
%          results.nvar   = # of variables in x-matrix
%          results.p      = # of variables in x-matrix (excluding constant
%                           term if used)
%          results.beta_std = std deviation of beta draws
%          results.sige_std = std deviation of sige draws
%          results.rho_std  = std deviation of rho draws
%          results.g_std    = m x 1 vector of posterior std deviations
%          results.sigma    = posterior mean of sige based on (e'*e)/(n-k)
%          results.bdraw    = bhat draws (1:thin:ndraw-nomit x nvar)
%          results.pdraw    = rho  draws (1:thin:ndraw-nomit x 1)
%          results.sdraw    = sige draws (1:thin:ndraw-nomit x 1)
%          results.gdraw    = gamma draws (1:thin:ndraw-nomit x m)
%          results.thin     = thinning value from input
%          results.direct   = nvar x 5 matrix with direct effect, t-stat, t-prob, lower01, upper99
%          results.indirect = nvar x 5 matrix with indirect effect, t-stat, t-prob, lower01, upper99
%          results.total    = nvar x 5 matrix with total effect, t-stat, t-prob, lower01, upper99
%          results.total_draws    = a (1:thin:ndraw-nomit,p) total x-impacts
%          results.direct_draws   = a (1:thin:ndraw-nomit,p) direct x-impacts
%          results.indirect_draws = a (1:thin:ndraw-nomit,p) indirect x-impacts
%          results.lmarginal= a scalar log-marginal likelihood, from MH-MC
%                             (Metropolis-Hastings Monte Carlo) integration of the log-posterior
%          results.drawpost = a vector of log-posterior draws (ndraw-nomit)x1
%          results.betapost = a matrix of log-posterior draws for beta (ndraw-nomit x k)
%          results.sigpost  = a vector of log-posterior draws for sige (ndraw-nomit x 1)
%          results.rho_gamma = (ndraw-nomit) x m+1 contains
%                              [rhodraws gammdraws] used to evaluate the log-posterior
%          results.rho_mode   = modal value of rho from the log-posterior
%          results.gamma_mode = modal values of gamma vector from the log-posterior
%          results.beta_mode  = modal values of beta vector from the log-posterior
%          results.sig_mode   = modal values of beta vector from the log-posterior
%          results.logC_sar   = constants associated with log-marginal
%                               logC_sar = gammaln(dof) - dof*log(2*pi)  -0.5*lndetx_sar
%                               dof = (n - m)/2; lndetx_sar = log(det(xp'*x));
%          results.logm_profile = a profile of the log-margainal over [rho_gamma isum];
%                                 where: isum = exp(logp -adj) + adj; adj = max(logp)
%          results.ndraw  = # of draws
%          results.nomit  = # of initial draws omitted
%          results.y      = y-vector from input (nobs x 1)
%          results.yhat   = mean of posterior predicted (nobs x 1)
%          results.resid  = residuals, based on posterior means
%          results.rsqr   = r-squared based on posterior means
%          results.rbar   = adjusted r-squared
%          results.taylor = taylor from input;
%          results.sampling_time = time for MCMC sampling
%          results.effects_time  = time to calculate effects estimates
%          results.trace_time    = time to calculate traces for taylor/chebyshev approximations
%          results.rmax   =  1  
%          results.rmin   = -1       
%          results.cflag  = 0 for intercept term, 1 for no intercept term
%          results.T2 = 2nd order Taylor series traces
%          results.T3 = 3rd order Taylor series traces
%          results.T4 = 4th order Taylor series traces
% --------------------------------------------------------------
% NOTES: - the intercept term (if you have one)
%          must be in the first column of the matrix x
% --------------------------------------------------------------
% SEE ALSO: (sar_conv_g_demo.m demos) 
% --------------------------------------------------------------
% REFERENCES: Debarsy and LeSage (2017) 
% Flexible dependence modeling using convex combinations of different
% types of connectivity structures, Regional Science & Urban Economics,
% Volume 69, pp. 46-68.

% Debarsy and LeSage (2018) 
% Bayesian model averaging for spatial autoregressive models
% based on convex combinations of different types of connectivity matrices
% unpublished manuscript
%----------------------------------------------------------------

% written by:
% James P. LeSage, last updated 12/2020
% Dept of Finance & Economics
% Texas State University-San Marcos
% 601 University Drive
% San Marcos, TX 78666
% james.lesage@txstate.edu

% error checking on inputs
[n junk] = size(y);
[nc, k] = size(x);
nvar = k;
if (nc ~= n)
       error('sdm_conv_panel_g: wrong sized x-matrix');
end

if ndraw <= 1000
       error('sdm_conv_panel_g: ndraw<=1000, increase ndraw to at least 10000');
end


results.nobs  = n;
results.nvar  = k;
results.y = y; 

    [n1,n2] = size(Wmatrices);
    m = n2/n1;
    if m~= round(m)
        error('sdm_conv_panel_g: wrong sized W-matrices');
    elseif n1 ~= n
        error('sdm_conv_panel_g: wrong sized W-matrices');
    elseif m < 2
        error('sdm_conv_panel_g: only one W-matrix');
    end;
   results.nmat=m;
if nargin == 5
    % use default arguments
    thin = 1;
    rmin = -1;     % use -1,1 rho interval as default
    rmax = 1;
    rho = 0.7; % starting values
    sige = 1;
    gamma = ones(m,1)*(1/m);
    results.rmin = -1;
    results.rmax = 1;
    ccmin = 0.4;
    ccmax = 0.6;
    ddmin = 0.1;
    ddmax = 0.4;
    tr_flag=0;
    
     results.ccmin = ccmin;
     results.ccmax = ccmax;
     results.ddmin = ddmin;
     results.ddmax = ddmax;
    
elseif nargin == 6

     [rho,sige,rmin,rmax,gamma,thin,ccmin,ccmax,ddmin,ddmax,T2,T3,T4,tr_flag,plt_flag] = sar_parse2(prior,m);

% check the thinning parameter
 eff = (ndraw -nomit)/thin;
if eff < 999
    warning('sdm_conv_panel_g: < 1000 draws after thining');
end

     thin = round(thin);
     results.thin = thin;
     results.rmin = rmin;
     results.rmax = rmax;
     results.ccmin = ccmin;
     results.ccmax = ccmax;
     results.ddmin = ddmin;
     results.ddmax = ddmax;
     

else
    error('sdm_conv_panel_g: wrong # of input arguments to sar_conv_panel_g');
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
            results.cflag = cflag;
            pp = size(x,2);
            results.p = pp;
        end
    elseif sum(x(:,1)) == n % we have an intercept in the right place
        cflag = 1;
        results.cflag = cflag;
        pp = size(x,2)-1;
        results.p = pp;
    end
    results.cflag = cflag;
    results.p = pp;
      
results.rmin = rmin;
results.rmax = rmax;


   % ====== initializations
ccmin = 0.4;
ccmax = 0.6;
ddmin = 0.10;
ddmax = 0.40;
acc_rate = zeros(ndraw,1);
gacc_rate = zeros(ndraw,1);
cc_save = zeros(ndraw,1);
gflag = zeros(ndraw,1);
cc = 0.1;
acc = 0;
gacc = 0;
dd = 3.0;
dd_save = zeros(ndraw,1);

Wy = y;
Wx = x;
begi = 1;
endi = n;
for ii=1:m
    Wy = [Wy Wmatrices(:,begi:endi)*y];
    if cflag == 0
        Wx = [Wx Wmatrices(:,begi:endi)*x];
    elseif cflag == 1
        Wx = [Wx Wmatrices(:,begi:endi)*x(:,2:end)];
    end
    begi = begi + n;
    endi = endi + n;
end

nvarw = size(Wx,2);

if (tr_flag == 0)
[T2,T3,T4,ctime] = calc_taylor_approx4(Wmatrices);
results.taylor_time = ctime;
elseif (tr_flag == 1)
results.taylor_time = 0;
end


results.T2=T2;
results.T3=T3;
results.T4=T4;

% storage for draws
bsave = zeros(ndraw,nvarw);
psave = zeros(ndraw,1);
ssave = zeros(ndraw,1);
gsave = zeros(ndraw,m);
drawpost = zeros(ndraw-nomit,1);
rho_gamma = zeros(ndraw-nomit,m+1);
betapost = zeros(ndraw-nomit,nvarw);
sigpost = zeros(ndraw-nomit,1);

timet = clock; % start the timer


    noo = 1000;
    
    for iter=1:ndraw
        xpx = Wx'*Wx;
        
        tmp = [1
            -rho*gamma];
        bd = xpx\(Wx'*Wy*tmp);
        
        AI = (xpx)\eye(nvarw);
        
        Wys = Wy*tmp;
        
        b = Wx'*Wys;
        b0 = AI*b;
        bhat = norm_rnd(sige*AI) + b0;
        xb = Wx*bhat;
        
        % update sige
        V = sum((Wys -xb).^2)/2;
        sige=1/gamrand(n/2,V);
        
        
        % ======================
        % M-H sample rho
        % ======================
        
        % Anonymous function
        lp_rho = @(r) cond_rho4(r,gamma,Wy,bd,Wx,T2,T3,T4,n);   % evaluates rho-conditional on gamma
        % obtain random-walk (tuned) proposed rho values
        rho2 = rho + cc*randn(1,1); % proposal for rho
        
        accept = 0;
        while accept == 0
            if ((rho2 > rmin) && (rho2 < rmax))
                accept = 1;
            else
                rho2 = rho + cc*randn(1,1);
            end
        end
        
        alpMH =  lp_rho(rho2) - lp_rho(rho);
        
        ru = unif_rnd(1,0,1);
        
        if alpMH > 0
            p = 1;
        else
            ratio = exp(alpMH);
            p = min(1,ratio);
        end
        if (ru < p)
            rho = rho2;
            acc = acc + 1;
        end
        
        acc_rate(iter,1) = acc/iter;
        
        % update cc based on std of rho draws
        if acc_rate(iter,1) < ccmin
            cc = cc/1.1;
        end
        if acc_rate(iter,1) > ccmax
            cc = cc*1.1;
        end
        
        if cc < 0.001
            cc = 0.1;
        end
        if cc > 1.000
            cc = 0.1;
        end
        
        cc_save(iter,1) = cc;
        % ======================
        % M-H sample gamma
        % ======================
        % Anonymous function
        lp_gamma = @(g) cond_gamma4(g,rho,Wy,bd,Wx,T2,T3,T4,n); % evaluates gamma conditional on rho
        if iter < noo
            % obtain a block of gamma candidates using reversible jump
            gtst = zeros(m-1,1);
            % flip a coin
            coin = unif_rnd(m-1,0,1);
            for jj=1:m-1
                if coin(jj,1) <= (1/3)
                    gtst(jj,1) = unif_rnd(1,0,gamma(jj,1));
                elseif (coin(jj,1)  > (1/3) && coin(jj,1)  <= (2/3))
                    gtst(jj,1) = gamma(jj,1);
                elseif coin(jj,1)  > (2/3)
                    gtst(jj,1) = unif_rnd(1,gamma(jj,1),1);
                end
            end
            
            gnew = [gtst
                1-sum(gtst)];
            accept = 0;
            cntr = 1;
            while (accept == 0)
                if all(gnew >= 0)
                    accept = 1;
                else
                    gtst = zeros(m-1,1);
                    % flip a 3-headed coin
                    coin = unif_rnd(m-1,0,1);
                    for jj=1:m-1
                        if coin(jj,1) <= (1/3)
                            gtst(jj,1) = unif_rnd(1,0,gamma(jj,1));
                        elseif (coin(jj,1)  > (1/3) && coin(jj,1)  <= (2/3))
                            gtst(jj,1) = gamma(jj,1);
                        elseif coin(jj,1)  > (2/3)
                            gtst(jj,1) = unif_rnd(1,gamma(jj,1),1);
                        end
                    end
                    gnew = [gtst
                        1-sum(gtst)];
                    cntr = cntr+1;
                end
            end
            
        elseif iter == noo
            rmp = std(gsave(iter-1000+1:iter-1,:));
            gam_std = rmp';
            ind = find(gam_std < 0.001);
            if length(ind) > 0
                gam_std(ind,1) = 0.01;
            end
            
            % obtain a block of gamma candidates using reversible jump
            gtst = zeros(m-1,1);
            % flip a coin
            coin = unif_rnd(m-1,0,1);
            for jj=1:m-1
                if coin(jj,1) <= (1/3)
                    gtst(jj,1) = unif_rnd(1,gamma(jj,1) - dd*gam_std(jj,1),gamma(jj,1));
                elseif (coin(jj,1)  > (1/3) && coin(jj,1)  <= (2/3))
                    gtst(jj,1) = gamma(jj,1);
                elseif coin(jj,1)  > (2/3)
                    gtst(jj,1) = unif_rnd(1,gamma(jj,1),gamma(jj,1) + dd*gam_std(jj,1));
                end
            end
            
            gnew = [gtst
                1-sum(gtst)];
            accept = 0;
            cntr = 1;
            while (accept == 0)
                if all(gnew >= 0)
                    accept = 1;
                else
                    gtst = zeros(m-1,1);
                    % flip a 3-headed coin
                    coin = unif_rnd(m-1,0,1);
                    for jj=1:m-1
                        if coin(jj,1) <= (1/3)
                            gtst(jj,1) = unif_rnd(1,gamma(jj,1) - dd*gam_std(jj,1),gamma(jj,1));
                        elseif (coin(jj,1)  > (1/3) && coin(jj,1)  <= (2/3))
                            gtst(jj,1) = gamma(jj,1);
                        elseif coin(jj,1)  > (2/3)
                            gtst(jj,1) = unif_rnd(1,gamma(jj,1),gamma(jj,1) + dd*gam_std(jj,1));
                        end
                    end
                    gnew = [gtst
                        1-sum(gtst)];
                    cntr = cntr+1;
                end
            end
            
        else % use gamma std to produce draws
            
            % obtain a block of gamma candidates using reversible jump
            gtst = zeros(m-1,1);
            % flip a coin
            coin = unif_rnd(m-1,0,1);
            for jj=1:m-1
                if coin(jj,1) <= (1/3)
                    gtst(jj,1) = unif_rnd(1,gamma(jj,1) - dd*gam_std(jj,1),gamma(jj,1));
                elseif (coin(jj,1)  > (1/3) && coin(jj,1)  <= (2/3))
                    gtst(jj,1) = gamma(jj,1);
                elseif coin(jj,1)  > (2/3)
                    gtst(jj,1) = unif_rnd(1,gamma(jj,1),gamma(jj,1) + dd*gam_std(jj,1));
                end
            end
            
            gnew = [gtst
                1-sum(gtst)];
            accept = 0;
            cntr = 1;
            while (accept == 0)
                if all(gnew >= 0)
                    accept = 1;
                else
                    gtst = zeros(m-1,1);
                    % flip a 3-headed coin
                    coin = unif_rnd(m-1,0.0001,0.9999);
                    for jj=1:m-1
                        if coin(jj,1) <= (1/3)
                            gtst(jj,1) = unif_rnd(1,gamma(jj,1) - dd*gam_std(jj,1),gamma(jj,1));
                        elseif (coin(jj,1)  > (1/3) && coin(jj,1)  <= (2/3))
                            gtst(jj,1) = gamma(jj,1);
                        elseif coin(jj,1)  > (2/3)
                            gtst(jj,1) = unif_rnd(1,gamma(jj,1),gamma(jj,1) + dd*gam_std(jj,1));
                        end
                    end
                    gnew = [gtst
                        1-sum(gtst)];
                    cntr = cntr+1;
                end
            end
        end
        
        % evaluate conditional for the block of gamma proposals
        
        alpMH =  lp_gamma(gnew) - lp_gamma(gamma);
        
        ru = unif_rnd(1,0,1);
        
        if alpMH > 0
            p = 1;
        else
            ratio = exp(alpMH);
            p = min(1,ratio);
        end
        
        if (ru < p)
            gflag(iter,1) = 1;
            gacc = gacc + 1;
            gamma = gnew;
        end
        
        gacc_rate(iter,1) = gacc/iter;
        
        % update dd based on std of gamma draws
        if gacc_rate(iter,1) < ddmin
            dd = dd/1.1;
        end
        if gacc_rate(iter,1) > ddmax
            dd = dd*1.1;
        end
        
        if dd>3.0
            dd=3.0;
        elseif dd<1.0
            dd=1.0;
        end
        
        
        dd_save(iter,1) = dd;
        
        % save draws
        bsave(iter,:) = bhat';
        ssave(iter,1) = sige;
        psave(iter,1) = rho;
        gsave(iter,:) = gamma';
        
        % Anonymous function
        log_post = @(r,g) joint_post4(r,g,Wy,bd,Wx,T2,T3,T4,n); % evaluates log posterior for both rho and gamma
        
        
        if ( mod(iter, noo) ==0 )
            rmp = std(gsave(iter-noo+1:iter,:));
            gam_std = rmp';
            ind = find(gam_std < 0.001);
            if length(ind) > 0
                gam_std(ind,1) = 0.01;
            end
            
            if plt_flag == 1
                subplot(2,2,1),
                tt=iter-noo+1:iter-1;
                subplot(2,2,1),
                plot(tt,gsave(iter-noo+1:iter-1,:));
                xlabel('gammas draws');
                subplot(2,2,2),
                plot(tt,psave(iter-noo+1:iter-1,1));
                xlabel('rho draws');
                subplot(2,2,3),
                plot(tt,bsave(iter-noo+1:iter-1,:));
                xlabel('beta draws');
                subplot(2,2,4),
                plot(tt,ssave(iter-noo+1:iter-1,:));
                xlabel('sigma draws');
                drawnow;
            end
            
        end
        
        
        if iter > nomit
            
            logpost = log_post(rho,gamma);
            drawpost(iter-nomit,1) = logpost;
            rho_gamma(iter-nomit,:) = [rho gamma'];
            tmp = [1
                -rho*gamma];
            btmp = bd;
            betapost(iter-nomit,:) = btmp';
            sigpost(iter-nomit,1) = ((Wy*tmp - Wx*btmp)'*(Wy*tmp - Wx*btmp))/n;
        end
        
        %
        
        
    end % end of draws loop
    
time = etime(clock,timet);
results.sampling_time = time;

results.gflag = gflag;

results.gacc_rate = gacc_rate;
results.cc = cc_save;
results.dd = dd_save;
results.acc_rate = acc_rate;

results.thin = thin;
results.bdraw = bsave(nomit+1:thin:ndraw,:);
results.pdraw = psave(nomit+1:thin:ndraw,1);
results.sdraw = ssave(nomit+1:thin:ndraw,1);
results.gdraw = gsave(nomit+1:thin:ndraw,:);
% results.acc_rate = acc_rate(nomit+1:ndraw,1);
% results.gacc_rate = gacc_rate(nomit+1:ndraw,1);
results.drawpost = drawpost; % we don't want to thin these
results.rho_gamma = rho_gamma;
results.betapost = betapost;
results.sigpost = sigpost;

% calculate log-marginal likelihood (using Mh-MC integration)
logp = results.drawpost;
rho_gamma = results.rho_gamma;
[adj,mind] = max(logp); 
results.rho_mode = rho_gamma(mind,1);
results.beta_mode = betapost(mind,:);
results.sig_mode = sigpost(mind,1);
results.gamma_mode = rho_gamma(mind,2:end);
isum = exp(logp -adj) + adj;
lndetx_sar = log(det(xpx));
% constant terms

dof = (n - m)/2; % we must include the # of weight matrices
D = (1 - 1/rmin); % from uniform prior on rho
logC_sar = -log(D) + gammaln(dof) - dof*log(2*pi)  -0.5*lndetx_sar;

results.logmarginal = mean(logp) + logC_sar;
results.logC_sar = logC_sar; % return constants
results.logm_profile = [rho_gamma betapost sigpost isum];
% results.logm_profile = [rho_gamma isum];


% compute posterior means for return arguments
bmean = mean(results.bdraw);
rho_mean = mean(results.pdraw);
smean = mean(results.sdraw);
gmean = mean(results.gdraw);

results.sige = smean;
results.beta = bmean';
results.rho = rho_mean;
results.gamma = gmean';

bdraw = results.bdraw;
if results.cflag == 1
    bdraw = bdraw(:,2:end);
end

pdraw = results.pdraw;

pp = results.p;

[niter,nvars] = size(bdraw);

begi = 1;
endi = n;
Wch = sparse(n,n);
for ii=1:m
    Wch = Wch + gmean(1,ii)*sparse(Wmatrices(:,begi:endi));
    begi = begi + n;
    endi = endi + n;
end

Wch = sparse(Wch);

uiter=150;
maxorderu=100;
nobs = n;
    rv=randn(nobs,uiter);
    tracew=zeros(maxorderu,1);
    wjjju=rv;
    for jjj=1:maxorderu
        wjjju=Wch*wjjju;
        tracew(jjj)=mean(mean(rv.*wjjju));
    end
    
    traces=[tracew];
    
        trs=[1;traces];
        ntrs=length(trs);
        
         trbig=trs';
         trbig2 = [trbig(1,2:end) trbig(1,end)];

         trmat = trbig;
         for s = 1:m
             trmat = [trmat
                 trbig2];
         end

    total = zeros(niter,pp,101);
    direct = zeros(niter,pp,101);
    indirect = zeros(niter,pp,101);
    
    for iter=1:niter
        
        ree = 0:1:ntrs-1;
        rmat = pdraw(iter,1).^ree;
        
        %         for i=1:niter
        rmat = pdraw(iter,1).^ree;
        for j=1:pp
            beta = bdraw(iter,j:pp:end);
             total(iter,j,:) = sum(beta,2)*rmat;
             direct(iter,j,:) = (beta*trmat).*rmat;
             indirect(iter,j,:) = total(iter,j,:) - direct(iter,j,:);
        end
    end

time = etime(clock,timet);
results.effects_time = time;

results.time = results.effects_time + results.sampling_time + results.taylor_time;
% 
% % ====================================================================
% 
total_draws = zeros(niter,pp);
direct_draws = zeros(niter,pp);
indirect_draws = zeros(niter,pp);
for i=1:pp
tmp = squeeze(total(:,i,:)); % an ndraw by 1 by ntraces matrix
total_draws(:,i) = (sum(tmp'))'; % an ndraw by 1 vector
tmp = squeeze(indirect(:,i,:)); % an ndraw by 1 by ntraces matrix
indirect_draws(:,i) = (sum(tmp'))'; % an ndraw by 1 vector
tmp = squeeze(direct(:,i,:)); % an ndraw by 1 by ntraces matrix
direct_draws(:,i) = (sum(tmp'))'; % an ndraw by 1 vector
end

results.total_draws = total_draws;
results.direct_draws = direct_draws;
results.indirect_draws = indirect_draws;

% Compute means, std deviation and upper and lower 0.95 intervals
total_out = zeros(pp,5);
total_save = zeros(niter,pp);
for i=1:pp
    tmp = squeeze(total(:,i,:)); % an ndraw by 1 by ntraces matrix
    total_mean = mean(tmp);
    total_std = std(tmp);
    % Bayesian 0.95 credible intervals
    % for the cumulative total effects
    total_sum = (sum(tmp'))'; % an ndraw by 1 vector
    cum_mean = cumsum(mean(tmp));
    cum_std = cumsum(std(tmp));
    total_save(:,i) = total_sum;
    bounds = cr_interval(total_sum,0.95);
    cmean = mean(total_sum);
    smean = std(total_sum);
    ubounds = bounds(1,1);
    lbounds = bounds(1,2);
    total_out(i,:) = [cmean cmean./smean tdis_prb(cmean./smean,nobs) lbounds ubounds];
end

% now do indirect effects
indirect_out = zeros(pp,5);
indirect_save = zeros(niter,pp);
for i=1:pp
    tmp = squeeze(indirect(:,i,:)); % an ndraw by 1 by ntraces matrix
    indirect_mean = mean(tmp);
    indirect_std = std(tmp);
    % Bayesian 0.95 credible intervals
    % for the cumulative indirect effects
    indirect_sum = (sum(tmp'))'; % an ndraw by 1 vector
    cum_mean = cumsum(mean(tmp));
    cum_std = cumsum(std(tmp));
    indirect_save(:,i) = indirect_sum;
    bounds = cr_interval(indirect_sum,0.95);
    cmean = mean(indirect_sum);
    smean = std(indirect_sum);
    ubounds = bounds(1,1);
    lbounds = bounds(1,2);
    indirect_out(i,:) = [cmean cmean./smean tdis_prb(cmean./smean,nobs) lbounds ubounds];
end


% now do direct effects
direct_out = zeros(pp,5);
direct_save = zeros(niter,pp);
for i=1:pp
    tmp = squeeze(direct(:,i,:)); % an ndraw by 1 by ntraces matrix
    direct_mean = mean(tmp);
    direct_std = std(tmp);
    % Bayesian 0.95 credible intervals
    % for the cumulative direct effects
    direct_sum = (sum(tmp'))'; % an ndraw by 1 vector
    cum_mean = cumsum(mean(tmp));
    cum_std = cumsum(std(tmp));
    direct_save(:,i) = direct_sum;
    bounds = cr_interval(direct_sum,0.95);
    cmean = mean(direct_sum);
    smean = std(direct_sum);
    ubounds = bounds(1,1);
    lbounds = bounds(1,2);
    direct_out(i,:) = [cmean cmean./smean tdis_prb(cmean./smean,nobs) lbounds ubounds];
end


results.direct = direct_out;
results.indirect = indirect_out;
results.total = total_out;


% results.meth  = 'sar_conv_g';
results.ndraw = ndraw;
results.nomit = nomit;
results.tflag = 'plevel';


% =========================================================================
% support functions below
% =========================================================================

function [logp] = joint_post4(rho,gamma,Wy,bd,x,T2,T3,T4,n)
% PURPOSE: evaluate the  joint distribution of rho and gamma
%  spatial autoregressive model using sparse matrix algorithms
% ---------------------------------------------------
%  USAGE: cout = joint_post(rho,gamma,Wy,bd,x,T2,T3,T4,n)
%  where:  rho  = spatial autoregressive parameter
%        gamma  = convex combination parameters
%          Wy   = AR dependence vector [y W1*y W2*y .. WL*y]
%          bd  = xpx\(x'*Wy);
%         T2, T3, T4 = traces used to find det(I-rho*W) 
%          n    =  nobs
% ---------------------------------------------------
%  RETURNS: a joint distribution value used in Monte Carlo integration
%           to produce the log-marginal likelihood
%  --------------------------------------------------

tmp = [1
     -rho*gamma];

Wys = Wy*tmp;
xb = x*bd;
ew = (Wys - xb);

epe = ew'*ew;


    a2=kron(gamma,gamma);
    aTa = a2'*T2;
    
    a3=kron(a2,gamma);
    aTTa =a3'*T3;
    
    a4=kron(a3,gamma);
    aTTTa =a4'*T4;
    
    wAw4 = (rho*rho*aTa/2) + (rho*rho*rho*aTTa/3) + (rho*rho*rho*rho*aTTTa/4);
    
    lndet = -wAw4;
    

logp =  lndet - (n/2)*log(epe);



function rhoc = cond_rho4(rho,gamma,Wy,bd,x,T2,T3,T4,n)
% PURPOSE: evaluate the  conditional distribution of rho given gamma
%  spatial autoregressive model using sparse matrix algorithms
% ---------------------------------------------------
%  USAGE:cout = cond_rho(rho,gamma,Wy,bd,x,T2,T3,T4,n)
%  where:  rho  = spatial autoregressive parameter
%        gamma  = convex combination parameters
%          Wy   = AR dependence vector [y W1*y W2*y .. WL*y]
%          bd  = xpx\(x'*Wy);
%         T2, T3, T4 = traces used to find det(I-rho*W) 
%                 using Chebyshev or Taylor series approximation 
%          n    =  nobs
% ---------------------------------------------------
%  RETURNS: a conditional used in Metropolis-Hastings sampling
%  NOTE: called only by sar_conv_g
%  --------------------------------------------------

tmp = [1
     -rho*gamma];
 
Wys = Wy*tmp;
xb = x*bd;
ew = (Wys - xb);

epe = ew'*ew;


    a2=kron(gamma,gamma);    
    aTa = a2'*T2;
    
    a3=kron(a2,gamma);
    aTTa =a3'*T3;
    
    a4=kron(a3,gamma);
    aTTTa =a4'*T4;
    
    wAw4 = (rho*rho*aTa/2) + (rho*rho*rho*aTTa/3) + (rho*rho*rho*rho*aTTTa/4);
    
    lndet = -wAw4;
    


rhoc =  lndet -(n/2)*log(epe);



% ====================================================
function gamc = cond_gamma4(gamma,rho,Wy,bd,x,T2,T3,T4,n)
% PURPOSE: evaluate the  conditional distribution of rho given gamma
%  spatial autoregressive model using sparse matrix algorithms
% ---------------------------------------------------
%  USAGE:cout = cond_gamma(gamma,rho,Wy,,bd,x,T2,T3,T4,n)
%  where:  rho  = spatial autoregressive parameter
%        gamma  = convex combination parameters
%          Wy   = AR dependence vector [y W1*y W2*y .. WL*y]
%          bd  = xpx\(x'*Wy);
%         T2, T3, T4 = traces used to find det(I-rho*W) 
%          n    =  nobs
% ---------------------------------------------------
%  RETURNS: a conditional used in Metropolis-Hastings sampling
%  NOTE: called only by sar_conv_g
%  --------------------------------------------------

tmp = [1
     -rho*gamma];

Wys = Wy*tmp;
xb = x*bd;
ew = (Wys - xb);

epe = ew'*ew;

    a2=kron(gamma,gamma);
    aTa = a2'*T2;
    
    a3=kron(a2,gamma);
    aTTa =a3'*T3;
    
    a4=kron(a3,gamma);
    aTTTa =a4'*T4;
    
    wAw4 = (rho*rho*aTa/2) + (rho*rho*rho*aTTa/3) + (rho*rho*rho*rho*aTTTa/4);
    
    lndet = -wAw4;
    

gamc =  lndet -(n/2)*log(epe);

% ====================================================

% ===========================================================================


function [rho,sige,rmin,rmax,gamma,thin,ccmin,ccmax,ddmin,ddmax,T2,T3,T4,tr_flag,plt_flag] = sar_parse2(prior,m)
% PURPOSE: parses input arguments for sar_conv_g models
% ---------------------------------------------------
%  USAGE: [rval,rho,sige,rmin,rmax,novi_flag] =  sar_parse(prior,k)
% returns values set by user or default values 
% ---------------------------------------------------

% written by:
% James P. LeSage, last updated 1/2018
% Dept of Finance & Economics
% Texas State University-San Marcos
% 601 University Drive
% San Marcos, TX 78666
% james.lesage@txstate.edu


% set defaults
gamma = ones(m,1)*(1/m);
rmin = -1;     % use -1,1 rho interval as default
rmax = 1;
rho = 0.7; % starting values
sige = 1;
thin = 1; % default to no thinning
ccmin = 0.4;
ccmax = 0.6;
ddmin = 0.1;
ddmax = 0.4;
tr_flag = 0;
plt_flag = 0;
T2 = [];
T3 = [];
T4 = [];


fields = fieldnames(prior);
nf = length(fields);
if nf > 0
    for i=1:nf
        if strcmp(fields{i},'thin')
            thin = prior.thin;
        elseif strcmp(fields{i},'ccmin')
            ccmin = prior.ccmin;
         elseif strcmp(fields{i},'ccmax')
            ccmax = prior.ccmax;
         elseif strcmp(fields{i},'ddmin')
            ddmin = prior.ddmin;
         elseif strcmp(fields{i},'ddmax')
            ddmax = prior.ddmax;          
        elseif strcmp(fields{i},'T2') && strcmp(fields{i},'T3') && strcmp(fields{i},'T4')
            T2 = prior.T2;    
            T3 = prior.T3;    
            T4 = prior.T4;    
            tr_flag = 1;
        elseif strcmp(fields{i},'plt_flag')
            plt_flag = prior.plt_flag;          
        end
    end
    
    
else % the user has input a blank info structure
    % so we use the defaults
end

function [tmat2,tmat3,tmat4,ctime] = calc_taylor_approx4(Wmatrices)
% PURPOSE: calculate 4th other trace matrices for Taylor Series
% approximation to the log-determinant

[n,nm] = size(Wmatrices);
m = nm/n;

% ========================================
% 2nd order
tmat2 = zeros(m*m,1);
% ========================================
% 3rd order
tmat3 = zeros(m*m*m,1);
% ========================================
% 4th order
tmat4 = zeros(m*m*m*m,1);


    tic;

    begi = 1;
    endi = n;
    cnti = 1;
    cntj = 1;
    cntk = 1;
    for ii=1:m;
        begj = 1;
        endj = n;
        Wi = sparse(Wmatrices(:,begi:endi));
        for jj=1:m;
            begk = 1;
            endk = n;
            Wj = sparse(Wmatrices(:,begj:endj));
            ijsave = (Wi*Wj);
            if (cnti <= m*m)
                tmat2(cnti,1) = sum(sum(Wi.*Wj'));
                cnti = cnti + 1;
            end;
            
            for kk=1:m;
                begl = 1;
                endl = n;
                Wk = sparse(Wmatrices(:,begk:endk));
                ijksave = ijsave*Wk;
                if (cntj <= m*m*m)
                    tmat3(cntj,1) = sum(sum((ijsave).*Wk'));
                    cntj = cntj + 1;
                end;
                for ll=1:m;
                    Wl = sparse(Wmatrices(:,begl:endl));
                    tmat4(cntk,1) = sum(sum((ijksave).*Wl'));
                    cntk = cntk+1;
                    begl = begl+n;
                    endl = endl+n;
                end;
                begk = begk + n;
                endk = endk + n;
            end;
            begj = begj+n;
            endj = endj+n;
        end;
        begi = begi + n;
        endi = endi + n;
    end;
    ctime = toc;
    


function x=gamrand(alpha,lambda)
% Gamma(alpha,lambda) generator using Marsaglia and Tsang method
% Algorithm 4.33
if alpha>1
    d=alpha-1/3; c=1/sqrt(9*d); flag=1;
    while flag
        Z=randn;
        if Z>-1/c
            V=(1+c*Z)^3; U=rand;
            flag=log(U)>(0.5*Z^2+d-d*V+d*log(V));
        end
    end
    x=d*V/lambda;
else
    x=gamrand(alpha+1,lambda);
    x=x*rand^(1/alpha);
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



