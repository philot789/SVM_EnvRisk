function results = sem_conv_panel_g(y,x,Wmatrices,N,T,ndraw,nomit,prior)
% PURPOSE: MCMC estimates for SEM convex combination of W model
%          using a convex combination of m different W-matrices
%          y = X*beta + u, u = rho*Wc*u + e, e = N(0,sige*I_n), 
%          Wc = g1*W1 + g2*W2 + ... + (1-g1-g2- ... -gm)*Wm
%          no prior for beta, sigma, rho, gamma
%-------------------------------------------------------------
% USAGE: results = sem_conv_panel_g(y,x,Wmatrices,N,T,ndraw,nomit,prior)
% where: y = dependent variable vector (N*T x 1)
%        x = independent variables matrix (N*T x nvar), 
%        Wmatrices = (nobs,m*nobs), nobs = N*T
%        e.g., Wmatrices = [W1 W2 ... Wm]
%        where each W1, W2, ... Wm are (nobs x nobs) row-normalized weight matrices
%    ndraw = # of draws (use lots of draws, say 25,000 to 50,000
%    nomit = # of initial draws omitted for burn-in  (probably around 2,500 to 5,000           
%    prior = a structure variable with:
%       prior.model = 0 pooled model without fixed effects (x may NOT contain an intercept)
%                   = 1 spatial fixed effects (x may NOT contain an intercept)
%                   = 2 time period fixed effects (x may NOT contain an intercept)
%                   = 3 spatial and time period fixed effects (x may NOT contain an intercept)
%       prior.fe    = report fixed effects and their t-values in prt_panel 
%                     (default=0=not reported; prior.fe=1=report) 
%            prior.thin  = a thinning parameter for use in analyzing
%                          posterior distributions, default = 1 (no thinning of draws)
%                          recommended value for ndraw > 20,000 is 100
%            prior.plt_flag = 0 for no plotting of MCMC draws, 1 for plotting every 1000 draws (default = 0)
%-------------------------------------------------------------
% RETURNS:  a structure:
%    results.meth  = 'semp_conv_g'    if prior.model=0
%                  = 'semsfe_conv_g'  if prior.model=1
%                  = 'semtfe_conv_g'  if prior.model=2
%                  = 'semstfe_conv_g' if prior.model=3          
%          results.beta     = posterior mean of bhat based on draws
%          results.rho      = posterior mean of rho based on draws
%          results.sige     = posterior mean of sige based on draws
%          results.gamma    = m x 1 vector of posterior means for g1,g2, ... gm
%                             where m is the number of weight matrices used on input
%          results.sfe   = spatial fixed effects (if info.model=1 or 3)
%          results.tfe   = time period fixed effects (if info.model=2 or 3)
%          results.tsfe  = t-values spatial fixed effects (if info.model=1 or 3)
%          results.ttfe  = t-values time period fixed effects (if info.model=2 or 3)
%          results.con   = intercept 
%          results.tcon   = t-value intercept
%          results.nobs   = # of observations
%          results.nvar   = # of variables in x-matrix
%          results.nobs   = # of observations
%          results.nvar   = # of variables in x-matrix
%          results.p      = # of variables in x-matrix (excluding constant
%                           term if used)
%          results.beta_std = std deviation of beta draws
%          results.sige_std = std deviation of sige draws
%          results.rho_std  = std deviation of rho draws
%          results.g_std    = m x 1 vector of posterior std deviations
%          results.sigma    = posterior mean of sige based on (e'*e)/(n-k)
%          results.bdraw    = bhat draws (ndraw-nomit x nvar)
%          results.pdraw    = rho  draws (ndraw-nomit x 1)
%          results.sdraw    = sige draws (ndraw-nomit x 1)
%          results.gdraw    = gamma draws (ndraw-nomit x m)
%          results.thin     = thinning value from input
%          results.ndraw  = # of draws
%          results.nomit  = # of initial draws omitted
%          results.y      = y-vector from input (nobs x 1)
%          results.yhat   = mean of posterior predicted (nobs x 1)
%          results.resid  = residuals, based on posterior means
%          results.rsqr   = r-squared based on posterior means
%          results.corr2  = goodness-of-fit between actual and fitted values
%          results.sampling_time = time for MCMC sampling
%          results.taylor_time   = time to calculate Taylor series terms
%          results.logmarginal   = log-marginal likelihood
%          results.rmax   = 1  
%          results.rmin   = -1       
%          results.tflag  = 'plevel' (default) for printing p-levels
%                         = 'tstat' for printing bogus t-statistics 
%          results.cflag  = 1 for intercept term, 0 for no intercept term
%          results.tflag  = 'plevel' (default) for printing p-levels
%                         = 'tstat' for printing bogus t-statistics 
%          results.lflag  = 0,1 for fast, slow log-marginal from input (or default = 0)
% --------------------------------------------------------------
% NOTES: - the intercept term (if you have one)
%          must be in the first column of the matrix x
% --------------------------------------------------------------
% SEE ALSO: (house_sem_demo.m, house_sem_demo2.m demos) 
% --------------------------------------------------------------
% REFERENCES: Debarsy and LeSage (2017) 
% Flexible dependence modeling using convex combinations of different
% types of connectivity structures
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
if (nc ~= n)
       error('sem_conv_panel_g: wrong sized x-matrix');
end; 

results.nobs  = n;
results.nvar  = k;
nvar = k;
results.y = y; 
results.N = N;
results.T = T;

    [n1,n2] = size(Wmatrices);
    m = n2/n1;
    
    if n1 == N*T
        Wlarge = 1;
    else
        Wlarge = 0;
    end
    
    if m < 2
        error('sem_conv_panel_g: only one W-matrix');
    end
    
   results.nmat=m;
   
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

fe = 0;

fields = fieldnames(prior);
nf = length(fields);
if nf > 0
    for i=1:nf
        if strcmp(fields{i},'model') model = prior.model;
        elseif strcmp(fields{i},'fe') fe = prior.fe;
        end
    end
end


[ywith,xwith,meanny,meannx,meanty,meantx]=demean(y,x,N,T,model);

results = sem_conv_gnew(ywith,xwith,Wmatrices,ndraw,nomit,prior);

if model==0
    results.meth='semp_conv_g';
elseif model==1
    results.meth='semsfe_conv_g';
elseif model==2
    results.meth='semtfe_conv_g';
elseif model==3
    results.meth='semstfe_conv_g';
else
    error('sem_conv_panel_g: wrong input number of prior.model');
end
results.model = model;
results.fe = fe;
results.nmat = m;


%          results.beta     = posterior mean of bhat based on draws
%          results.rho      = posterior mean of rho based on draws
%          results.sige     = posterior mean of sige based on draws
%          results.gamma    = m x 1 vector of posterior means for g1,g2, ... gm
%                             where m is the number of weight matrices used on input

sige = results.sige;
results.beta = results.beta;

en=ones(T,1);
et=ones(N,1);
ent=ones(N*T,1);
nobs = N*T;

% step 4) find fixed effects and their t-values
if model==1
    intercept=mean(y)-mean(x)*results.beta;
    results.con=intercept;
    results.sfe=meanny-meannx*results.beta-kron(et,intercept);
    xhat=x*results.beta+kron(en,results.sfe)+kron(ent,intercept);
    results.tsfe=results.sfe./sqrt(sige/T*ones(N,1)+diag(sige*meannx*(xwith'*xwith)*meannx'));
    results.tcon=results.con/sqrt(sige/nobs+sige*mean(x)*(xwith'*xwith)*mean(x)');
    tnvar=N;
elseif model==2
    intercept=mean(y)-mean(x)*results.beta;
    results.con=intercept;
    results.tfe=meanty-meantx*results.beta-kron(en,intercept); 
    xhat=x*results.beta+kron(results.tfe,et)+kron(ent,intercept);
    results.ttfe=results.tfe./sqrt(sige/N*ones(T,1)+diag(sige*meantx*(xwith'*xwith)*meantx'));
    results.tcon=results.con/sqrt(sige/nobs+sige*mean(x)*(xwith'*xwith)*mean(x)');
    tnvar=T;
elseif model==3
    intercept=mean(y)-mean(x)*results.beta; 
    results.con=intercept;
    results.sfe=meanny-meannx*results.beta-kron(et,intercept);
    results.tfe=meanty-meantx*results.beta-kron(en,intercept);
    results.tsfe=results.sfe./sqrt(sige/T*ones(N,1)+diag(sige*meannx*(xwith'*xwith)*meannx'));
    results.ttfe=results.tfe./sqrt(sige/N*ones(T,1)+diag(sige*meantx*(xwith'*xwith)*meantx'));
    results.tcon=results.con/sqrt(sige/nobs+sige*mean(x)*(xwith'*xwith)*mean(x)');
    xhat=x*results.beta+kron(en,results.sfe)+kron(results.tfe,et)+kron(ent,intercept);
    tnvar=N+T;
else
    xhat=x*results.beta;
    tnvar=0;
end
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

% =========================================================================
% support functions below
% =========================================================================


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


function results = sem_conv_gnew(y,x,Wmatrices,ndraw,nomit,prior)
% PURPOSE: Bayesian estimates of the spatial autoregressive error model
%          using a convex combination of m different W-matrices
%          y = X*beta + u, u = rho*Wc*u + e, e = N(0,sige*I_n), 
%          Wc = g1*W1 + g2*W2 + ... + (1-g1-g2- ... -gm)*Wm
%          b = N(c,T), priors for beta 
%          c = k x 1 vector or prior means
%          T = k x k prior variance-covariance
%          1/sige = Gamma(a,b), diffuse priors for sige
%          no prior for rho
%          no prior for gamma
%-------------------------------------------------------------
% USAGE: results = sem_conv_g(y,x,Wmatrices,ndraw,nomit,prior)
% where: y = dependent variable vector (nobs x 1)
%        x = independent variables matrix (nobs x nvar), 
%            the intercept term (if present) must be in the first column of the matrix x
%        Wmatrices = (nobs,m*nobs)
%        e.g., Wmatrices = [W1 W2 ... Wm]
%        where each W1, W2, ... Wm are (nobs x nobs) row-normalized weight matrices
%    ndraw = # of draws (use lots of draws, say 25,000 to 50,000
%    nomit = # of initial draws omitted for burn-in  (probably around 5,000           
%    prior = a structure variable with:
%            prior.c = k x 1 vector of prior means for beta (default = 0)
%            prior.T = k x k prior variance-covariance matrix (default = eye(k)*1e+12)
%            prior.a = scalar prior value for IG(a,b) for sigma (default = 0)
%            prior.b = scalar prior value for IG(a,b) for sigma (default = 0)
%            prior.thin  = a thinning parameter for use in analyzing
%                          posterior distributions, default = 1 (no thinning of draws)
%                          recommended value for ndraw > 20,000 is 100
%            prior.lflag = 0 for fast log-marginal approximation, 1 for exact log-marginal
%                          (default = 0)
%-------------------------------------------------------------
% RETURNS:  a structure:
%          results.meth     = 'sem_conv_g'
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
%          results.bdraw    = bhat draws (ndraw-nomit x nvar)
%          results.pdraw    = rho  draws (ndraw-nomit x 1)
%          results.sdraw    = sige draws (ndraw-nomit x 1)
%          results.gdraw    = gamma draws (ndraw-nomit x m)
%          results.thin     = thinning value from input
%          results.ndraw  = # of draws
%          results.nomit  = # of initial draws omitted
%          results.y      = y-vector from input (nobs x 1)
%          results.yhat   = mean of posterior predicted (nobs x 1)
%          results.resid  = residuals, based on posterior means
%          results.rsqr   = r-squared based on posterior means
%          results.rbar   = adjusted r-squared
%          results.bprior_mean = c;
%          results.bprior_bcov = T;
%          results.sprior_a = a;
%          results.sprior_b = b;
%          results.sampling_time = time for MCMC sampling
%          results.taylor_time   = time to calculate Taylor series terms
%          results.logm_time     = time to calculate log-marginal likelihoods
%          results.logm_sem_W    = m+1 vector of log-marginal likelihoods
%                                logm_sem_W(1,1) = logm for Wc matrix
%                                logm_sem_W(2,1) = logm for W1 matrix
%                                logm_sem_W(.,1) = logm for W. matrix
%                                logm_sem_W(m+1,1) = logm for Wm matrix
%          results.rmax   = 1  
%          results.rmin   = -1       
%          results.tflag  = 'plevel' (default) for printing p-levels
%                         = 'tstat' for printing bogus t-statistics 
%          results.cflag  = 1 for intercept term, 0 for no intercept term
%          results.tflag  = 'plevel' (default) for printing p-levels
%                         = 'tstat' for printing bogus t-statistics 
%          results.lflag  = 0,1 for fast, slow log-marginal from input (or default = 0)
% --------------------------------------------------------------
% NOTES: - the intercept term (if you have one)
%          must be in the first column of the matrix x
% --------------------------------------------------------------
% SEE ALSO: (house_sem_demo.m, house_sem_demo2.m demos) 
% --------------------------------------------------------------
% REFERENCES: Debarsy and LeSage (2017) 
% Flexible dependence modeling using convex combinations of different
% types of connectivity structures
%----------------------------------------------------------------

% written by:
% James P. LeSage, last updated 6/2017
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
       error('sem_conv_g: wrong sized x-matrix');
end; 

results.nobs  = n;
results.nvar  = k;
results.y = y; 

    [n1,n2] = size(Wmatrices);
    m = n2/n1;
    if n2 ~= m*n1
        error('sem_conv_g: wrong sized W-matrices');
    elseif n1 ~= n
        error('sem_conv_g: wrong sized W-matrices');
    elseif m < 1
        error('sem_conv_g: wrong # of W-matrices');
    end;
    

if nargin == 5
    % use default arguments
    nskip = 1;
    rmin = -0.9999;     % use -1,1 rho interval as default
    rmax = 0.9999;
    lam = 0.7;
    sige = 1;
    gamma = ones(m,1)*(1/m);
    nu = 0;
    d0 = 0;
    c = zeros(k,1);   % diffuse prior for beta
    T = eye(k)*1e+12;
    results.priorb = 0;
    lflag = 0;
    results.lflag = 0;
    results.rmin = -1;
    results.rmax = 1;
    
    
elseif nargin == 6
     [nu,d0,lam,sige,rmin,rmax,gamma,c,T,thin,ccmin,ccmax,ggmin,ggmax,priorb,lflag,plt_flag] = sar_parse(prior,k,m);
%[nu,d0,lam,sige,rmin,rmax,gamma,c,T,thin,ccmin,ccmax,ggmin,ggmax,priorb,lflag,plt_flag] = sar_parse(prior,k,m);

     results.thin = thin;
     results.rmin = rmin;
     results.rmax = rmax;
     results.bprior_mean = c;
     results.bprior_bcov = T;
     results.sprior_a = nu;
     results.sprior_b = d0;
     results.ccmin = ccmin;
     results.ccmax = ccmax;
     results.ggmin = ggmin;
     results.ggmax = ggmax;
     results.priorb = priorb;
     results.lflag = lflag;

else
    error('sem_conv_g: wrong # of input arguments to sem_conv_g');
end;

% check if the user handled the intercept term okay
    n = length(y);
    if sum(x(:,1)) ~= n
        tst = sum(x); % we may have no intercept term
        ind = find(tst == n); % we do have an intercept term
        if length(ind) > 0
            error('sem_conv_g: intercept term must be in first column of the x-matrix');
        elseif length(ind) == 0 % case of no intercept term
            cflag = 0;
            pp = size(x,2);
        end;
    elseif sum(x(:,1)) == n % we have an intercept in the right place
        cflag = 1;
        pp = size(x,2)-1;
    end;
     
    results.cflag = cflag;
    results.p = pp;
    

results.rmin = rmin;
results.rmax = rmax;


% storage for draws
          bsave = zeros(ndraw,k);
          psave = zeros(ndraw,1);
          ssave = zeros(ndraw,1);
          gsave = zeros(ndraw,m);
          drawpost = zeros(ndraw,1);

          
% ====== initializations
% compute this stuff once to save time

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

% ===================================================================


sym = 0;
[T2,T3,T4,ctime] = calc_taylor_approx(Wmatrices,sym);

results.taylor_time = ctime;


timet = clock; % start the timer
    
ys = zeros(n,(m+1));

ys(:,1) = y;
    begi = 1;
    endi = n;

for ii=1:m
    ys(:,ii+1) =  -Wmatrices(:,begi:endi)*y;
    begi = begi+n;
    endi = endi+n;
end

xs = zeros(n,(m+1),k);

for j=1:k
    xj = x(:,j);
    xs(:,1,j) = xj;
    begi = 1;
    endi = n;

    for ii=1:m
        xs(:,ii+1,j) = -Wmatrices(:,begi:endi)*xj;
        begi = begi+n;
        endi = endi+n;
    end
end

    noo = 1000;

for iter=1:ndraw
    
    % update beta
    % ===================================================
    
    omega = [1
         lam*gamma];
     
xo = zeros(n,k);
for j=1:k
    xo(:,j) = squeeze(xs(:,:,j))*omega;
end

AI = (xo'*xo)\eye(k);
b = (xo'*ys*omega);
    b0 = AI*b;
    bhat = norm_rnd(sige*AI) + b0;


    % update sige
    % ===================================================
    nu1 = n;
    e = ys*omega - xo*bhat;

    d1 = e'*e;
    chi = chis_rnd(1,nu1);
    sige = d1/chi;

    
    % update lambda
    % ==============================================================
    % Metropolis-Hastings sampling for lambda
    
    % Anonymous function
    lp_lam = @(r) cond_lam4(r,gamma,ys,xs,bhat,T2,T3,T4,n,k);   % evaluates rho-conditional on gamma
% function rhoc = cond_lam4(rho,gamma,ys,xs,bhat,T2,T3,T4,n)

    % evaluate lambda conditional at current lambda and gamma values

    lam2 = lam + cc*randn(1,1); % proposal for rho
    
    accept = 0;
    while accept == 0
        if ((lam2 > rmin) && (lam2 < rmax))
            accept = 1;
        else
            lam2 = lam + cc*randn(1,1);
        end
    end
    
     
    alpMH =  lp_lam(lam2) - lp_lam(lam);

    ru = unif_rnd(1,0,1);
    
    if alpMH > 0
        p = 1;
    else
        ratio = exp(alpMH);
        p = min(1,ratio);
    end
    if (ru < p)
        lam = lam2;
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
    
% ===================================================
% sample gamma values
% Anonymous function
    lp_gamma = @(g) cond_gamma4(g,lam,ys,xs,bhat,T2,T3,T4,n,k); % evaluates gamma conditional on rho
% function gamc = cond_gamma4(gamma,lam,ys,xs,bhat,T2,T3,T4,n)

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
      
    bsave(iter,:) = bhat';
    ssave(iter,1) = sige;
    psave(iter,1) = lam;
    gsave(iter,:) = gamma';

    if iter > nomit
    % Anonymous function    
    log_post = @(r,g) joint_post4(r,g,ys,bhat,xs,T2,T3,T4,n,k); % evaluates log posterior for both rho and gamma 
     logpost = log_post(lam,gamma);
     drawpost(iter,1) = logpost;
    end
% 

         
    if ( mod(iter, noo) ==0 )
        
% Anonymous function    
    log_post = @(r,g) joint_post4(r,g,ys,bhat,xs,T2,T3,T4,n,k); % evaluates log posterior for both rho and gamma 
     logpost = log_post(lam,gamma);
     drawpost(iter,1) = logpost;
% 
        rmp = std(gsave(iter-noo+1:iter,:));
        gam_std = rmp';
        ind = find(gam_std < 0.001);
        if length(ind) > 0
            gam_std(ind,1) = 0.01;
        end
     
        if plt_flag == 1
        subplot(2,2,1),       
        tt=iter-noo+1:iter-1;
        plot(tt,gsave(iter-noo+1:iter-1,:));
        xlabel('gammas draws');
        subplot(2,2,2),
        plot(tt,psave(iter-noo+1:iter-1,1));
        xlabel('lambda draws');
        subplot(2,2,3),
        plot(tt,bsave(iter-noo+1:iter-1,:));
        xlabel('beta draws');
        subplot(2,2,4),
        plot(tt,ssave(iter-noo+1:iter-1,:));
        xlabel('sigma draws');       
        drawnow;
        end
        
    end
    

end
    
    
    

time = etime(clock,timet);
results.sampling_time = time;

results.time = results.sampling_time + results.taylor_time;

results.gacc_rate = gacc_rate;
results.cc = cc_save;
results.acc_rate = acc_rate;

results.thin = thin;

% results.logm = logm;

results.bdraw = bsave(nomit+1:thin:ndraw,:);
results.pdraw = psave(nomit+1:thin:ndraw,1);
results.sdraw = ssave(nomit+1:thin:ndraw,1);
results.gdraw = gsave(nomit+1:thin:ndraw,:);
results.drawpost = drawpost(nomit+1:ndraw,1);


% compute posterior means for return arguments
bmean = mean(results.bdraw);
rho_mean = mean(results.pdraw);
smean = mean(results.sdraw);
gmean = mean(results.gdraw);

results.sige = smean;
results.beta = bmean';
results.rho = rho_mean;
results.gamma = gmean';

results.sige_mean = smean;
results.beta_mean = bmean';
results.rho_mean = rho_mean;
results.gamma_mean = gmean';



% calculate log-marginal likelihood (using Mh-MC integration)
logp = results.drawpost;
% rho_gamma = [results.pdraw results.gdraw];
% [adj,mind] = max(logp); 
% results.rho_mode = rho_gamma(mind,1);
% results.beta_mode = betapost(mind,:);
% results.sig_mode = sigpost(mind,1);
% results.gamma_mode = rho_gamma(mind,2:end);
% isum = exp(logp -adj) + adj;
% constant terms
omega = [1
         rho_mean*gmean'];
     
xo = zeros(n,nvar);
for j=1:nvar
    xo(:,j) = squeeze(xs(:,:,j))*omega;
end;

xpx = (xo'*xo);
lndetx_sar = log(det(xpx));

dof = (n - m)/2; % we must include the # of weight matrices
D = (1 - 1/rmin); % from uniform prior on rho
logC_sar = -log(D) + gammaln(dof) - dof*log(2*pi)  -0.5*lndetx_sar;

results.logmarginal = mean(logp) + logC_sar;


% calculate fit statistics using posterior means for gamma


[nobs,nvar] = size(x);

        omega= [1
            rho_mean*gmean'];
        
     
xo = zeros(n,k);
for j=1:k
    xo(:,j) = squeeze(xs(:,:,j))*omega;
end;

    e = ys*omega - xo*bhat;

    epe = e'*e;

% 
ym = y - mean(y);
rsqr1 = epe;
rsqr2 = ym'*ym;
results.rsqr = 1- rsqr1/rsqr2; % r-squared
rsqr1 = rsqr1/(nobs-nvar);
rsqr2 = rsqr2/(nobs-1.0);
results.rbar = 1 - (rsqr1/rsqr2); % rbar-squared
% 


results.meth  = 'sem_conv_g';
% results.total = total_save;         
% results.direct = direct_save;       
% results.indirect = indirect_save;   
% results.beta = beta;
% results.rho = rho;
% results.gamma = gamma;
results.beta_std = std(results.bdraw)';
results.sige_std = std(results.sdraw);
results.rho_std = std(results.pdraw);
results.gamma_std = std(results.gdraw)';
% results.vmean = vmean;
% results.yhat  = yhat;
% results.resid = e;
results.ndraw = ndraw;
results.nomit = nomit;
results.tflag = 'plevel';
% results.rval = rval;


% =========================================================================
% support functions below
% =========================================================================

function lamx = cond_lam4(lam,gamma,ys,xs,bhat,T2,T3,T4,n,k)

omega = [1
    lam*gamma];


xo = zeros(n,k);
for j=1:k
    xo(:,j) = squeeze(xs(:,:,j))*omega;
end


nmk = (n-k)/2;
xpx = xo'*xo;
detx = 0.5*log(det(xpx));

e = ys*omega - xo*bhat;

epe = e'*e;

epe = nmk*log(epe);

a2=kron(gamma,gamma);
aTa = a2'*T2;

a3=kron(a2,gamma);
aTTa =a3'*T3;

a4=kron(a3,gamma);
aTTTa =a4'*T4;

wAw4 = (lam*lam*aTa/2) + (lam*lam*lam*aTTa/3) + (lam*lam*lam*lam*aTTTa/4);

detm = -wAw4;
% conditional at current lambda value
% lamx =   detm - detx - epe;
lamx = detm -detx - epe;

function gamx = cond_gamma4(gamma,lam,ys,xs,bhat,T2,T3,T4,n,k)

% ===================================================
% evaluate conditional for the current lambda and gamma values
omega = [1
    lam*gamma];

xo = zeros(n,k);
for j=1:k
    xo(:,j) = squeeze(xs(:,:,j))*omega;
end;

e = ys*omega - xo*bhat;

epe = e'*e;

nmk = (n-k)/2;

epe = nmk*log(epe);

xpx = xo'*xo;

detx = 0.5*log(det(xpx));

a2=kron(gamma,gamma);
aTa = a2'*T2;

a3=kron(a2,gamma);
aTTa =a3'*T3;

a4=kron(a3,gamma);
aTTTa =a4'*T4;

wAw4 = (lam*lam*aTa/2) + (lam*lam*lam*aTTa/3) + (lam*lam*lam*lam*aTTTa/4);

detm = -wAw4;
% conditional at current values
gamx = detm -detx - epe;


function [nu,d0,rho,sige,rmin,rmax,gamma,c,T,thin,ccmin,ccmax,ggmin,ggmax,priorb,lflag,plt_flag] = sar_parse(prior,k,m)
% PURPOSE: parses input arguments for sem_conv_g models
% ---------------------------------------------------
%  USAGE: [rval,rho,sige,rmin,rmax,novi_flag] =  sar_parse(prior,k)
% returns values set by user or default values 
% ---------------------------------------------------

% written by:
% James P. LeSage, last updated 12/2020
% Dept of Finance & Economics
% Texas State University-San Marcos
% 601 University Drive
% San Marcos, TX 78666
% james.lesage@txstate.edu


% set defaults
plt_flag = 0;
gamma = ones(m,1)*(1/m);
rmin = -0.9999;     % use -1,1 rho interval as default
rmax = 0.9999;
rho = 0.7;
sige = 1;
nu = 0;
d0 = 0;
c = zeros(k,1);   % diffuse prior for beta
T = eye(k)*1e+12;
thin = 1; % default to no thinning
ccmin = 0.4;
ccmax = 0.6;
ggmin = 0.01;
ggmax = 0.2;
priorb = 0;
lflag = 0; % default to fast log-marginal approximation

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
        elseif strcmp(fields{i},'ggmin')
            ggmin = prior.ggmin;
        elseif strcmp(fields{i},'ggmax')
            ggmax = prior.ggmax;
        elseif strcmp(fields{i},'beta')
            c = prior.beta;
            priorb = 1;
        elseif strcmp(fields{i},'bcov')
            T = prior.bcov;
            priorb = 1;
        elseif strcmp(fields{i},'a')
            nu = prior.nu;
        elseif strcmp(fields{i},'b')
            d0 = prior.d0;
        elseif strcmp(fields{i},'lflag')
            lflag = prior.lflag;
        elseif strcmp(fields{i},'plt_flag')
            plt_flag = prior.plt_flag;

        end;
    end;
    
    
else % the user has input a blank info structure
    % so we use the defaults
end;


function [tmat2,tmat3,tmat4,ctime] = calc_taylor_approx(Wmatrices,sym)
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

if (sym == 1)
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
                tmat2(cnti,1) = sum(sum(Wi.*Wj));
                cnti = cnti + 1;
            end;
            
            for kk=1:m;
                begl = 1;
                endl = n;
                Wk = sparse(Wmatrices(:,begk:endk));
                ijksave = ijsave*Wk;
                if (cntj <= m*m*m)
                    %                  tmat3(cntj,1) = sum(sum((Wi*Wj).*Wk));
                    tmat3(cntj,1) = sum(sum((ijsave).*Wk));
                    cntj = cntj + 1;
                end;
                for ll=1:m;
                    Wl = sparse(Wmatrices(:,begl:endl));
                    %                 tmat4(cntk,1) = sum(sum((Wi*Wj*Wk).*Wl));
                    tmat4(cntk,1) = sum(sum((ijksave).*Wl));
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
  
elseif sym == 0 % case of asymmetric matrices
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
            ijsave = (Wi*Wj');
            if (cnti <= m*m)
                tmat2(cnti,1) = sum(sum(Wi.*Wj'));
                cnti = cnti + 1;
            end;
            
            for kk=1:m;
                begl = 1;
                endl = n;
                Wk = sparse(Wmatrices(:,begk:endk));
                ijksave = ijsave*Wk';
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
    
else
    tic;
    error('calc_taylor_approx: sym must be 0 or 1');
ctime = 0;
end;


function out=lndetfull(rvec,W)
% PURPOSE: computes Pace and Barry's grid for log det(I-rho*W) using sparse matrices
% -----------------------------------------------------------------------
% USAGE: out = lndetfull(W,lmin,lmax)
% where:    
%             W     = symmetric spatial weight matrix (standardized)
%             lmin  = lower bound on rho
%             lmax  = upper bound on rho
% -----------------------------------------------------------------------
% RETURNS: out = a structure variable
%          out.lndet = a vector of log determinants for 0 < rho < 1
%          out.rho   = a vector of rho values associated with lndet values
% -----------------------------------------------------------------------
% NOTES: should use 1/lambda(max) to 1/lambda(min) for all possible rho values
% -----------------------------------------------------------------------
% References: % R. Kelley Pace and  Ronald Barry. 1997. ``Quick
% Computation of Spatial Autoregressive Estimators'', Geographical Analysis
% -----------------------------------------------------------------------
 
% written by:
% James P. LeSage, Dept of Economics
% University of Toledo
% 2801 W. Bancroft St,
% Toledo, OH 43606
% jpl@jpl.econ.utoledo.edu

spparms('tight'); 
[n junk] = size(W);
z = speye(n) - 0.1*sparse(W);
p = colamd(z);
niter = length(rvec);
dettmp = zeros(niter,2);
for i=1:niter;
    rho = rvec(i);
    z = speye(n) - rho*sparse(W);
    [l,u] = lu(z(:,p));
    dettmp(i,1) = rho;
    dettmp(i,2) = sum(log(abs(diag(u))));
end;

out = dettmp(:,2);


function out=lndetmc(order,iter,wsw,xx)
% PURPOSE: computes Barry and Pace MC approximation to log det(I-rho*W)
% -----------------------------------------------------------------------
% USAGE: out = lndetmc(order,iter,W,rmin,rmax)
% where:      order = # of moments u'(wsw^j)u/(u'u) to examine (default = 50)
%              iter = how many realizations are employed (default = 30)
%                 W = symmetric spatial weight matrix (standardized)  
%              grid = increment for lndet grid (default = 0.01)
% -----------------------------------------------------------------------
% RETURNS: out = a structure variable
%          out.lndet = a vector of log determinants for -1 < rho < 1
%          out.rho   = a vector of rho values associated with lndet values
%          out.up95  = an upper 95% confidence interval on the approximation
%          out.lo95  = a lower  95% confidence interval on the approximation
% -----------------------------------------------------------------------
% NOTES: only produces results for a grid of 0 < rho < 1 by default
%        where the grid ranges by 0.01 increments
% -----------------------------------------------------------------------
% References: Ronald Barry and R. Kelley Pace, "A Monte Carlo Estimator
% of the Log Determinant of Large Sparse Matrices", Linear Algebra and
% its Applications", Volume 289, Number 1-3, 1999, pp. 41-54.
% -----------------------------------------------------------------------
 
 
% Written by Kelley Pace, 6/23/97 
% (named fmcdetnormgen1.m in the spatial statistics toolbox )
% Documentation modified by J. LeSage 11/99

[n,n]=size(wsw);

% Exact moments from 1 to oexact
td=full([0;sum(sum(wsw.^2))/2]);
oexact=length(td);

o=order;
% Stochastic moments

mavmomi=zeros(o,iter);
for j=1:iter;
u=randn(n,1);
v=u;
utu=u'*u;
for i=1:o;
v=wsw*v;
mavmomi(i,j)=n*((u'*v)/(i*utu));
end;
end;

mavmomi(1:oexact,:)=td(:,ones(iter,1));

%averages across iterations
avmomi=mean(mavmomi')';

clear u,v;

%alpha matrix

alpha=xx;
valpha=vander(alpha);
valphaf=fliplr(valpha);
alomat=-valphaf(:,(2:(o+1)));

%Estimated ln|I-aD| using mixture of exact, stochastic moments
%exact from 1 to oexact, stochastic from (oexact+1) to o

lndetmat=alomat*avmomi;

out = lndetmat;



% =========================================================================
% support functions below
% =========================================================================


function [logp] = joint_post4(rho,gamma,ys,bhat,xs,T2,T3,T4,n,nvar)
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

        omega = [1
         rho*gamma];
     

xo = zeros(n,nvar);
for j=1:nvar
    xo(:,j) = squeeze(xs(:,:,j))*omega;
end

    e = ys*omega - xo*bhat;

    epe = e'*e;
    
    nmk = (n - nvar)/2;
    
    epe = nmk*log(epe);

    xpx = xo'*xo;

    detx = 0.5*log(det(xpx));

a2=kron(gamma,gamma);
aTa = a2'*T2;

a3=kron(a2,gamma);
aTTa =a3'*T3;

a4=kron(a3,gamma);
aTTTa =a4'*T4;


    % evaluate conditional at proposed gamma values
    wAw4 = (rho*rho*aTa/2) + (rho*rho*rho*aTTa/3) + (rho*rho*rho*rho*aTTTa/4);
    
    detm = -wAw4;
    % conditional at current lambda value
logp = detm - (n/2)*log(epe);
    

