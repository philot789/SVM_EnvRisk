function results = sar_panel_FE_g(y,x,W,T,ndraw,nomit,prior)
% PURPOSE: MCMC SAR model estimates for static spatial panels 
%          (N regions*T time periods) with spatial fixed effects (sfe) 
%          and/or time period fixed effects (tfe)
%          y = rho*W*y + X*b + sfe(optional) + tfe(optional) + e, 
%          e = N(0,sige*V), 
%          V = diag(v_1,v_2,...v_N*T), r/vi = ID chi(r)/r, r = 5 (default)
%          b   = N(c,C), default c = 0, C = eye(k)*1e+12
%          sige = gamma(nu,d0), default nu=0, d0=0   
%          no prior for rho
% Supply data sorted first by time and then by spatial units, so first region 1,
% region 2, et cetera, in the first year, then region 1, region 2, et
% cetera in the second year, and so on
% sar_panel_FE_g transforms y and x to deviation of the spatial and/or time means
% ---------------------------------------------------
%  USAGE: results = sar_panel_FE_g(y,x,W,T,ndraw,nomit,prior)
%  where:  y = N*T x 1 dependent variable vector
%          x = N*T x k independent variables matrix
%          W = spatial weights matrix (standardized)
%          N. B. W-matrix can be N*T x N*T or N x N
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
%       prior.fe    = report fixed effects and their t-values in prt_panel
%                     (default=0=not reported; prior.fe=1=report) 
%               prior.beta, prior means for beta,   b (default (k x 1) vector = 0)
%               priov.bcov, prior beta covariance, C above (default eye(k)*1e+12)
%               prior.rval, rval prior hyperparameter, default=4
%               prior.nu,   informative Gamma(nu,d0) prior on sige
%               prior.d0    informative Gamma(nu,d0) prior on sige
%                           default for above: nu=0,d0=0 (diffuse prior)
%       prior.rmin  = (optional) minimum value of rho to use in search  
%       prior.rmax  = (optional) maximum value of rho to use in search    
%       prior.lflag = 0 for full lndet computation (default = 1, fastest)
%                  = 1 for MC lndet approximation (fast for very large problems)
%                  = 2 for Spline lndet approximation (medium speed)
%       prior.order = order to use with info.lflag = 1 option (default = 50)
%       prior.iter  = iterations to use with info.lflag = 1 option (default = 30)  
%       prior.lndet = a matrix returned in results.lndet containing log-determinant information to save time
% ---------------------------------------------------
%  RETURNS: a structure
%         results.meth  = 'sarsfe_g' if prior.model=1
%                       = 'sartfe_g' if prior.model=2
%                       = 'sarstfe_g' if prior.model=3
%         results.beta  = bhat
%         results.rho   = rho 
%         results.bdraw = (ndraw-nomit)xk matrix of MCMC draws for beta
%         results.pdraw = (ndraw-nomit)x1 vector of MCMC draws for rho
%         results.vmean = N*T x 1 vector of v_{it} means
%         results.sdraw = (ndraw-nomit)x1 vector of MCMC draws for sige
%         results.bmean = b prior means (prior.beta from input)
%         results.bstd  = b prior std deviation, sqrt(diag(prior.bcov))
%         results.nu    = prior nu-value for sige prior (default = 0)
%         results.d0    = prior d0-value for sige prior (default = 0)
%         results.iprior = 1 for informative prior on beta, 
%                        = 0 for default no prior on beta
%         results.direct   = nvar x 5 matrix with direct effect, t-stat, t-prob, lower05, upper95
%         results.indirect = nvar x 5 matrix with indirect effect, t-stat, t-prob, lower05, upper95
%         results.total    = nvar x 5 matrix with total effect, t-stat, t-prob, lower05, upper95
%         results.direct_draws   = ndraw x nvar matrix of direct effect draws
%         results.indirect_draws = ndraw x nvar matrix of indirect effect draws
%         results.total_draws    = ndraw x nvar matrix of total effect draws
%         results.cov   = asymptotic variance-covariance matrix of the parameters b(eta) and rho
%         results.tstat = asymp t-stat (last entry is rho=spatial autoregressive coefficient)
%         results.yhat  = [inv(y-p*W)]*[x*b+fixed effects] (according to prediction formula)
%         results.resid = y-p*W*y-x*b
%         results.sige  = (y-p*W*y-x*b)'*(y-p*W*y-x*b)/n
%         results.rsqr  = rsquared
%         results.corr2 = goodness-of-fit between actual and fitted values
%         results.sfe   = spatial fixed effects (if prior.model=1 or 3)
%         results.tfe   = time period fixed effects (if prior.model=2 or 3)
%         results.tsfe  = t-values spatial fixed effects (if prior.model=1 or 3)
%         results.ttfe  = t-values time period fixed effects (if prior.model=2 or 3)
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
%         You should use lflag = 0 to get exact results, 
%         Fixed effects and their t-values are calculated as the deviation
%         from the mean intercept
% ---------------------------------------------------

% written by James P LeSage (Texas State University)
% last updated 11/2020

time1 = 0; 
time2 = 0;
time3 = 0;
time4 = 0;

timet = clock; % start the clock for overall timing

W=sparse(W);

fe=0;
model=0;

k = size(x,2);

results.nvar = k;

fields = fieldnames(prior);
nf = length(fields);
novi_flag = 0;
rval = 5; % rval = 5 is default
results.rval = rval;
nu = 0; d0 = 0; % default to a diffuse prior on sige
c = zeros(k,1); 
Tj = eye(k)*1e+12;
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
        if size(c,1) ~= k
            error('sar_panel_FE_g: wrong size prior means, must be k x 1 vector');
        end
        results.iprior = 1;
        results.bmean = c;
    elseif strcmp(fields{i},'bcov')
        TI = prior.bcov;
        if size(TI,2) ~= k
            error('sar_panel_FE_g: wrong size prior variance-covariance, must be k x k matrix');
        end
        results.iprior = 1;
        results.bstd = diag(sqrt(TI));
        
        Q = inv(TI);
        Qpc = Q*c;
        
        
    end
end


fields = fieldnames(prior);
nf = length(fields);
if nf > 0
    for i=1:nf
        if strcmp(fields{i},'model') model = prior.model;
        elseif strcmp(fields{i},'fe') fe = prior.fe;
        end
    end
end
if model==0
    results.meth='psar_g';
elseif model==1
    results.meth='sarsfe_g';
elseif model==2
    results.meth='sartfe_g';
elseif model==3
    results.meth='sarstfe_g';
else
    error('sar_panel_FE_g: wrong input number of info.model');
end

% check size of user inputs for comformability
[nobs, nvar] = size(x);
N = nobs/T;

    [n1,n2] = size(W);
    
    if n1 == N*T
        Wlarge = 1;
    else
        Wlarge = 0;
        W = kron(eye(T),W);
    end

[nchk junk] = size(y);
if nchk ~= nobs
error('sar_panel_FE_g: wrong size vector y or matrix x');
end;

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
      
if (fe==1 & model==0 ) error('info.fe=1, but cannot compute fixed effects if info.model is set to 0 or not specified'); end

% parse input options
[rmin,rmax,detval,ldetflag,eflag,order,miter,sflag,p,cflag] = sar_parse(prior);

if novi_flag == 1
    results.homo = 1;
    results.hetero = 0;
elseif novi_flag == 0
    results.hetero = 1;
    results.homo = 0;
end

results.ndraw = ndraw;
results.nomit = nomit;

% compute eigenvalues or limits
[rmin,rmax,time2] = sar_eigs(eflag,W,rmin,rmax,N); % function of LeSage

% do log-det calculations
[detval,time1] = sar_lndet(ldetflag,W,rmin,rmax,detval,order,miter); % function of LeSage


for t=1:T
    t1=1+(t-1)*N;t2=t*N;
    Wtmp = W(t1:t2,t1:t2);
    Wy(t1:t2,1)=Wtmp*y(t1:t2,1);
end

% demeaning of the y and x variables, depending on (info.)model

if (model==1 | model==3);
meanny=zeros(N,1);
meannwy=zeros(N,1);
meannx=zeros(N,nvar);
for i=1:N
    ym=zeros(T,1);
    wym=zeros(T,1);
    xm=zeros(T,nvar);
    for t=1:T
        ym(t)=y(i+(t-1)*N,1);
        wym(t)=Wy(i+(t-1)*N,1);
        xm(t,:)=x(i+(t-1)*N,:);
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
meantx=zeros(T,nvar);
for i=1:T
    t1=1+(i-1)*N;t2=i*N;
    ym=y([t1:t2],1);
    wym=Wy([t1:t2],1);
    xm=x([t1:t2],:);
    meanty(i)=mean(ym);
    meantwy(i)=mean(wym);
    meantx(i,:)=mean(xm);
end
clear ym wym xm;
end % if statement
    
en=ones(T,1);
et=ones(N,1);
ent=ones(nobs,1);

if model==1
    ywith=y-kron(en,meanny);
    wywith=Wy-kron(en,meannwy);
    xwith=x-kron(en,meannx);
elseif model==2
    ywith=y-kron(meanty,et);
    wywith=Wy-kron(meantwy,et);
    xwith=x-kron(meantx,et);
elseif model==3
    ywith=y-kron(en,meanny)-kron(meanty,et)+kron(ent,mean(y));
    wywith=Wy-kron(en,meannwy)-kron(meantwy,et)+kron(ent,mean(Wy));
    xwith=x-kron(en,meannx)-kron(meantx,et)+kron(ent,mean(x));
else
    ywith=y;
    wywith=Wy;
    xwith=x;
end % if statement


% ====== initializations
% compute this stuff once to save time

in = ones(N*T,1);
V = in;
vi = in;
sige = 1;
rho = 0.5;
vmean = vi;


bsave = zeros(ndraw-nomit,nvar);
ssave = zeros(ndraw-nomit,1);
psave = zeros(ndraw-nomit,1);

Wbig = W;

% step 1) do regressions
t0 = clock;
          
          % do MCMC sampling 
switch novi_flag
    
case{0} % we do heteroscedastic model    

%hwait = waitbar(0,'sar: MCMC sampling ...');

t0 = clock;                  
iter = 1;
          while (iter <= ndraw) % start sampling;
                  
          % update beta   
          xs = matmul(xwith,sqrt(V));
          ys = sqrt(V).*ywith;
          Wys = sqrt(V).*(Wbig*ywith);
          AI = (xs'*xs + sige*Q)\eye(nvar);         
          yss = ys - rho*Wys;          
          xpy = xs'*yss;
          b = xs'*yss + sige*Qpc;
          b0 = (xs'*xs + sige*Q)\b;
          bhat = norm_rnd(sige*AI) + b0;  
          xb = xs*bhat;
                    
          % update sige
          nu1 = nu + N*T; 
          e = (yss - xb);
          d1 = d0 + e'*e;
          chi = chis_rnd(1,nu1);
          sige = d1/chi;

          % update vi
          ev = ywith - rho*Wbig*ywith - xwith*bhat; 
          chiv = chis_rnd(N*T,rval+1);  
          %chiv = chi2rnd(rval+1,n,1); % Statistics Toolbox function
          vi = ((ev.*ev/sige) + in*rval)./chiv;
          V = in./vi; 
                        
          
      % we use griddy Gibbs to perform rho-draw
          b0 = (xs'*xs + sige*Q )\(xs'*ys + sige*Qpc);
          bd = (xs'*xs + sige*Q)\(xs'*Wys + sige*Qpc);
          e0 = ys - xs*b0;
          ed = Wys - xs*bd;
          epe0 = e0'*e0;
          eped = ed'*ed;
          epe0d = ed'*e0;
          
          rho = draw_rho(detval,epe0,eped,epe0d,N,T,nvar,rho);   
          
          if iter > nomit % if we are past burn-in, save the draws
              bsave(iter-nomit,1:nvar) = bhat';
              ssave(iter-nomit,1) = sige;
              psave(iter-nomit,1) = rho;
              vmean = vmean + vi;
              

          end
          
iter = iter + 1; 
%waitbar(iter/ndraw);         
          end % end of sampling loop
%close(hwait);

time3 = etime(clock,t0);
results.time3 = time3;

case{1} % we do homoscedastic model 
    
%hwait = waitbar(0,'sar: MCMC sampling ...');

t0 = clock;                  
iter = 1;

          while (iter <= ndraw) % start sampling;
                  
          % update beta   
          ys = ywith - rho*Wbig*ywith;         
          b = xwith'*ys + sige*Qpc;
          AI = (xwith'*xwith + sige*Q)\eye(nvar);         
          b0 = AI*b;
          bhat = norm_rnd(sige*AI) + b0;  
          xb = xwith*bhat;
                   
          % update sige
          nu1 = nu + N*T; 
          e = (ys - xb);
          d1 = d0 + e'*e;
          chi = chis_rnd(1,nu1);
          sige = d1/chi;
          
          % update rho using griddy Gibbs
          xs = xwith;
          ys = ywith;
          Wys = Wbig*ys;
          b0 = (xs'*xs + sige*Q )\(xs'*ys + sige*Qpc);
          bd = (xs'*xs + sige*Q)\(xs'*Wys + sige*Qpc);
          e0 = ys - xs*b0;
          ed = Wys - xs*bd;
          epe0 = e0'*e0;
          eped = ed'*ed;
          epe0d = ed'*e0;

          rho = draw_rho(detval,epe0,eped,epe0d,N,T,nvar,rho);

          if iter > nomit % if we are past burn-in, save the draws
              bsave(iter-nomit,1:nvar) = bhat';
              ssave(iter-nomit,1) = sige;
              psave(iter-nomit,1) = rho;
              vmean = vmean + vi;
          end
                    
iter = iter + 1; 
%waitbar(iter/ndraw);         
          end % end of sampling loop
%close(hwait);

time3 = etime(clock,t0);
results.time3 = time3;


otherwise
error('sar_panel_FE_g: unrecognized prior.novi_flag value on input');
% we should never get here

end; % end of homoscedastic vs. heteroscedastic vs. log-marginal options

% calculate effects estimates
        
t0 = clock; 

% pre-calculate traces for the x-impacts calculations
uiter=50;
maxorderu=100;
nobs = N*T;
rv=randn(nobs,uiter);
tracew=zeros(maxorderu,1);
wjjju=rv;
for jjj=1:maxorderu
    wjjju=W*wjjju;
    tracew(jjj)=mean(mean(rv.*wjjju));
    
end

[ndrawsg,junk] = size(bsave);

p = results.p;

traces=[tracew];
traces(1)=0;
traces(2)=sum(sum(W'.*W))/nobs;
trs=[1;traces];
ntrs=length(trs);
trbig=trs';
                 
        bdraws = bsave;
        pdraws = psave;

        ree = 0:1:ntrs-1;

        rmat = zeros(1,ntrs);
        total = zeros(ndrawsg,p,ntrs);
        direct = zeros(ndrawsg,p,ntrs);
        indirect = zeros(ndrawsg,p,ntrs);
        
for i=1:ndrawsg
    rmat = pdraws(i,1).^ree;
    for j=1:p
            beta = [bdraws(i,j)];
            total(i,j,:) = beta(1,1)*rmat;
    direct(i,j,:) = (beta*trbig).*rmat;
    indirect(i,j,:) = total(i,j,:) - direct(i,j,:);
    end

end

% Compute means, std deviation and upper and lower 0.95 intervals
% p = nvar;
ndraw = ndraw - nomit;
total_out = zeros(p,5);
total_save = zeros(ndraw,p);
for i=1:p
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
indirect_out = zeros(p,5);
indirect_save = zeros(ndraw,p);
for i=1:p;
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
end;


% now do direct effects
direct_out = zeros(p,5);
direct_save = zeros(ndraw,p);
for i=1:p;
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
end;


results.direct = direct_out;
results.indirect = indirect_out;
results.total = total_out;

total_draws = zeros(ndraw,p);
direct_draws = zeros(ndraw,p);
indirect_draws = zeros(ndraw,p);
for i=1:p;
tmp = squeeze(total(:,i,:)); % an ndraw by 1 by ntraces matrix
total_draws(:,i) = (sum(tmp'))'; % an ndraw by 1 vector
tmp = squeeze(indirect(:,i,:)); % an ndraw by 1 by ntraces matrix
indirect_draws(:,i) = (sum(tmp'))'; % an ndraw by 1 vector
tmp = squeeze(direct(:,i,:)); % an ndraw by 1 by ntraces matrix
direct_draws(:,i) = (sum(tmp'))'; % an ndraw by 1 vector
end;

results.total_draws = total_draws;
results.direct_draws = direct_draws;
results.indirect_draws = indirect_draws;

% step 3) find b,sige MCMC estimates

bmean = mean(bsave);
results.beta = bmean';
results.bdraw = bsave;
rmean = mean(psave);
results.rho = rmean; 
results.pdraw = psave;
results.vmean = vmean/(ndraw-nomit);
smean = mean(ssave);
results.sige = smean;
results.sdraw = ssave;
sige = results.sige;

% step 4) find fixed effects and their t-values
if model==1
    intercept=mean(y)-mean(Wy)*results.rho-mean(x)*results.beta;
    results.con=intercept;
    results.sfe=meanny-meannwy*results.rho-meannx*results.beta-kron(et,intercept);
    xhat=x*results.beta+kron(en,results.sfe)+kron(ent,intercept);
    results.tsfe=results.sfe./sqrt(sige/T*ones(N,1)+diag(sige*meannx*(xwith'*xwith)*meannx'));
    results.tcon=results.con/sqrt(sige/nobs+sige*mean(x)*(xwith'*xwith)*mean(x)');
    tnvar=N; 
elseif model==2
    intercept=mean(y)-mean(Wy)*results.rho-mean(x)*results.beta;
    results.con=intercept;
    results.tfe=meanty-meantwy*results.rho-meantx*results.beta-kron(en,intercept); 
    xhat=x*results.beta+kron(results.tfe,et)+kron(ent,intercept);
    results.ttfe=results.tfe./sqrt(sige/N*ones(T,1)+diag(sige*meantx*(xwith'*xwith)*meantx'));
    results.tcon=results.con/sqrt(sige/nobs+sige*mean(x)*(xwith'*xwith)*mean(x)');
    tnvar=T;
elseif model==3
    intercept=mean(y)-mean(Wy)*results.rho-mean(x)*results.beta; 
    results.con=intercept;
    results.sfe=meanny-meannwy*results.rho-meannx*results.beta-kron(et,intercept);
    results.tfe=meanty-meantwy*results.rho-meantx*results.beta-kron(en,intercept);
    results.tsfe=results.sfe./sqrt(sige/T*ones(N,1)+diag(sige*meannx*(xwith'*xwith)*meannx'));
    results.ttfe=results.tfe./sqrt(sige/N*ones(T,1)+diag(sige*meantx*(xwith'*xwith)*meantx'));
    results.tcon=results.con/sqrt(sige/nobs+sige*mean(x)*(xwith'*xwith)*mean(x)');
    xhat=x*results.beta+kron(en,results.sfe)+kron(results.tfe,et)+kron(ent,intercept);
    tnvar=N+T;
else
    xhat=x*results.beta;
    tnvar=0;
end    


% r-squared and corr-squared between actual and fitted values
results.tnvar=tnvar;
results.resid = y - rmean*Wy - xhat; 
yme=y-mean(y);
rsqr2=yme'*yme;
rsqr1 = results.resid'*results.resid;
results.rsqr=1.0-rsqr1/rsqr2; %rsquared

yhat=zeros(nobs,1);
ywithhat=zeros(nobs,1);
for t=1:T
    t1=1+(t-1)*N;t2=t*N;
    ywithhat(t1:t2,1)=(speye(N) - rmean*W(t1:t2,t1:t2))\xwith(t1:t2,:)*results.beta;
    yhat(t1:t2,1)=(speye(N) - rmean*W(t1:t2,t1:t2))\xhat(t1:t2,1);
end
res1=ywith-mean(ywith);
res2=ywithhat-mean(ywith);
rsq1=res1'*res2;
rsq2=res1'*res1;
rsq3=res2'*res2;
results.corr2=rsq1^2/(rsq2*rsq3); %corr2
results.yhat=yhat;

% calculate t-stats
bstd = std(bsave);
pstd = std(psave);

mean_parms = [bmean'
              rmean];
std_parms = [bstd'
             pstd];
         
tstat = mean_parms./std_parms;

results.tstat = tstat;

parms = [bsave psave ssave];

results.cov = cov(parms);


% return stuff
results.nobs  = nobs; 
results.nvar  = nvar;
results.rmax  = rmax;      
results.rmin  = rmin;
results.lflag = ldetflag;
results.order = order;
results.miter = miter;
results.fe    = fe;
results.time  = etime(clock,timet);
results.time1 = time1;
results.time2 = time2;
results.time3 = time3;
results.time4 = time4;
results.lndet = detval;
results.N     = N;
results.T     = T;
results.model = model;


function [rmin,rmax,detval,ldetflag,eflag,order,iter,sflag,p,cflag] = sar_parse(info)
% PURPOSE: parses input arguments for sar model
% ---------------------------------------------------
%  USAGE: [rmin,rmax,convg,maxit,detval,ldetflag,eflag,order,iter,options] = sar_parse(info)
% where info contains the structure variable with inputs 
% and the outputs are either user-inputs or default values
% ---------------------------------------------------

% set defaults

eflag = 0;     % default to not computing eigenvalues
ldetflag = 1;  % default to 1999 Pace and Barry MC determinant approx
order = 50;    % there are parameters used by the MC det approx
iter = 30;     % defaults based on Pace and Barry recommendation
rmin = -1;     % use -1,1 rho interval as default
rmax = 1;
detval = 0;    % just a flag
sflag = 0;
p = 0;
cflag = 0;

fields = fieldnames(info);
nf = length(fields);
if nf > 0
    
 for i=1:nf
    if strcmp(fields{i},'rmin')
        rmin = info.rmin;  eflag = 0;
    elseif strcmp(fields{i},'rmax')
        rmax = info.rmax; eflag = 0;
    elseif strcmp(fields{i},'p')
        p = info.p;
    elseif strcmp(fields{i},'cflag')
        cflag = info.cflag;
    elseif strcmp(fields{i},'lndet')
    detval = info.lndet;
    ldetflag = -1;
    eflag = 0;
    rmin = detval(1,1);
    nr = length(detval);
    rmax = detval(nr,1);
    elseif strcmp(fields{i},'lflag')
        tst = info.lflag;
        if tst == 0,
        ldetflag = 0; % compute full lndet, no approximation
        elseif tst == 1,
        ldetflag = 1; % use Pace-Barry approximation
        elseif tst == 2,
        ldetflag = 2; % use spline interpolation approximation
        else
        error('sar_panel_g: unrecognizable lflag value on input');
        end
    elseif strcmp(fields{i},'order')
        order = info.order;  
    elseif strcmp(fields{i},'iter')
        iter = info.iter;  
    elseif strcmp(fields{i},'eig')
        eflag = info.eig;  
     elseif strcmp(fields{i},'sflag')
        sflag = info.sflag; 
    end
 end
 
else, % the user has input a blank info structure
      % so we use the defaults
end 

function [rmin,rmax,time2] = sar_eigs(eflag,W,rmin,rmax,n);
% PURPOSE: compute the eigenvalues for the weight matrix
% ---------------------------------------------------
%  USAGE: [rmin,rmax,time2] = far_eigs(eflag,W,rmin,rmax,W)
% where eflag is an input flag, W is the weight matrix
%       rmin,rmax may be used as default outputs
% and the outputs are either user-inputs or default values
% ---------------------------------------------------


if eflag == 1 % do eigenvalue calculations
t0 = clock;
opt.tol = 1e-3; opt.disp = 0;
lambda = eigs(sparse(W),speye(n),1,'SR',opt);  
rmin = real(1/lambda);   
rmax = 1.0;
time2 = etime(clock,t0);
else % use rmin,rmax arguments from input or defaults -1,1
time2 = 0;
end;


function [detval,time1] = sar_lndet(ldetflag,W,rmin,rmax,detval,order,iter);
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
            error('sar: wrong lndet input argument');
        end;
        [n1,n2] = size(detval);
        if n2 ~= 2
            error('sar: wrong sized lndet input argument');
        elseif n1 == 1
            error('sar: wrong sized lndet input argument');
        end;          
end;

function rho = draw_rho(detval,epe0,eped,epe0d,n,t,k,rho)
% PURPOSE: draws rho-values using griddy Gibbs and inversion
% ---------------------------------------------------
%  USAGE: rho = draw_rho(detval,epe0,eped,epe0d,n,k,rho,a1,a2)
% where info contains the structure variable with inputs 
% and the outputs are either user-inputs or default values
% ---------------------------------------------------
% REFERENCES: LeSage and Pace (2009) Introduction to Spatial Econometrics
% Chapter 5, pp 136-141 on Bayesian spatial regression models.


% written by:
% James P. LeSage, last updated 3/2010
% Dept of Finance & Economics
% Texas State University-San Marcos
% 601 University Drive
% San Marcos, TX 78666
% jlesage@spatial-econometrics.com


nmk = (n*t-k)/2;
nrho = length(detval(:,1));
iota = ones(nrho,1);

z = epe0*iota - 2*detval(:,1)*epe0d + detval(:,1).*detval(:,1)*eped;
z = -nmk*log(z);
den =  detval(:,2) + z;

% bprior = beta_prior(detval(:,1),a1,a2);
% den = den + log(bprior);
n = length(den);
y = detval(:,1);
adj = max(den);
den = den - adj;
x = exp(den);
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

% To see how this works, uncomment the following lines
% plot(detval(:,1),den/1000,'-');
% line([detval(idraw,1) detval(idraw,1)],[0 den(idraw,1)/1000]);
% hold on;
% line([detval(idraw,1) 0],[den(idraw,1)/1000 den(idraw,1)/1000]);
% drawnow;
% pause;


function H = hessian(f,x,varargin)
% PURPOSE: Computes finite difference Hessian
% -------------------------------------------------------
% Usage:  H = hessian(func,x,varargin)
% Where: func = function name, fval = func(x,varargin)
%           x = vector of parameters (n x 1)
%    varargin = optional arguments passed to the function
% -------------------------------------------------------
% RETURNS:
%           H = finite differnce hessian
% -------------------------------------------------------

% Code from:
% COMPECON toolbox [www4.ncsu.edu/~pfackler]
% documentation modified to fit the format of the Ecoometrics Toolbox
% by James P. LeSage, Dept of Economics
% University of Toledo
% 2801 W. Bancroft St,
% Toledo, OH 43606
% jlesage@spatial-econometrics.com

eps = 1e-6;

n = size(x,1);
fx = feval(f,x,varargin{:});
 
% Compute the stepsize (h)
h = eps.^(1/3)*max(abs(x),1e-2);
xh = x+h;
h = xh-x;    
ee = sparse(1:n,1:n,h,n,n);
 
% Compute forward step 
g = zeros(n,1);
for i=1:n
  g(i) = feval(f,x+ee(:,i),varargin{:});
end
   
H=h*h';
% Compute "double" forward step 
for i=1:n
for j=i:n
  H(i,j) = (feval(f,x+ee(:,i)+ee(:,j),varargin{:})-g(i)-g(j)+fx)/H(i,j);
  H(j,i) = H(i,j);
end
end


