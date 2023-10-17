function results = ols_panel_FE_g(y,x,T,ndraw,nomit,prior)
% PURPOSE: computes MCMC regression model estimates for static panels 
%          (N regions*T time periods) with spatial fixed effects (sfe) 
%          and/or time period fixed effects (tfe)
%          y = X*b + sfe(optional) + tfe(optional) + e, 
%          e = N(0,sige*V), 
%          V = diag(v_1,v_2,...v_N*T), r/vi = ID chi(r)/r, r = 5 (default)
%          b   = N(c,C), default c = 0, C = eye(2*k)*1e+12
%          sige = gamma(nu,d0), default nu=0, d0=0   
% Supply data sorted first by time and then by spatial units, so first region 1,
% region 2, et cetera, in the first year, then region 1, region 2, et
% cetera in the second year, and so on
% ols_panel_FE_g transforms y and x to deviation of the spatial and/or time means
% ---------------------------------------------------
%  USAGE: results = ols_panel_FE_g(y,x,T,ndraw,nomit,prior)
%  where:  y = N*T x 1 dependent variable vector
%          x = N*T x k independent variables matrix
%          T = number of points in time
%      ndraw = # of MCMC draws
%      nomit = # of draws to exclude for start-up
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
%       prior.beta, prior means for beta,   b,g above (default 2*(k x 1) vector = 0)
%       priov.bcov, prior beta covariance, T above (default eye(2*k)*1e+12)
%       prior.rval, rval prior hyperparameter, default=4
%       prior.nu,   informative Gamma(nu,d0) prior on sige
%       prior.d0    informative Gamma(nu,d0) prior on sige
%                   default for above: nu=0,d0=0 (diffuse prior)
% ---------------------------------------------------
%  RETURNS: a structure
%         results.meth  = 'olssfe' if prior.model=1
%                       = 'olstfe' if prior.model=2
%                       = 'olsstfe' if prior.model=3
%         results.beta  = bhat
%         results.tstat = asymp t-stat
%         results.bdraw = (ndraw-nomit) x k matrix of MCMC draws
%         results.vmean = mean of vi draws (nobs x 1)
%         results.bmean = b prior means (prior.beta from input)
%         results.bstd  = b prior std deviation, sqrt(prior.bcov)
%         results.nu    = prior nu-value for sige prior (default = 0)
%         results.d0    = prior d0-value for sige prior (default = 0)
%         results.iprior = 1 for informative prior on beta, 
%                        = 0 for default no prior on beta
%         results.yhat  = x*b+fixed effects] (according to prediction formula)
%         results.resid = y-x*b
%         results.sige  = (y-x*b)'*(y-x*b)/n
%         results.rsqr  = rsquared
%         results.corr2 = goodness-of-fit between actual and fitted values
%         results.sfe   = spatial fixed effects (if prior.model=1 or 3)
%         results.tfe   = time period fixed effects (if prior.model=2 or 3)
%         results.tsfe  = t-values spatial fixed effects (if prior.model=1 or 3)
%         results.ttfe  = t-values time period fixed effects (if prior.model=2 or 3)
%         results.con   = intercept 
%         results.con   = t-value intercept
%         results.nobs  = # of observations
%         results.nvar  = # of explanatory variables in x 
%         results.tnvar = # fixed effects
%         results.fe    = fe from input
%         results.time1 = time for MCMC sampling
%         results.time  = total time taken      
% --------------------------------------------------
%  NOTES: 
%         Fixed effects and their t-values are calculated as the deviation
%         from the mean intercept
% ---------------------------------------------------

% written by James P LeSage (Texas State University)
% last updated 11/2020


timet = clock; % start the clock for overall timing

fe=0;
model=0;

results.ndraw = ndraw;
results.nomit = nomit;

k = size(x,2);

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
        novi_flag = 0;
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
            error('ols_panel_FE_g: wrong size prior means, must be k x 1 vector');
        end
        results.iprior = 1;
        results.bmean = c;
    elseif strcmp(fields{i},'bcov')
        TI = prior.bcov;
        if size(TI,2) ~= k
            error('ols_panel_FE_g: wrong size prior variance-covariance, must be k x k matrix');
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
%         elseif strcmp(fields{i},'Nhes') Nhes = info.Nhes;
        end
    end
end
if model==0
    results.meth='pols_g';
elseif model==1
    results.meth='olssfe_g';
elseif model==2
    results.meth='olstfe_g';
elseif model==3
    results.meth='olsstfe_g';
else
    error('ols_panel_FE_g: wrong input number of info.model');
end

% check size of user inputs for comformability
[nobs nvar] = size(x);
N = nobs/T;
% [N Ncol] = size(W);
% if N ~= Ncol
% error('ols_panel_FE_g: wrong size weight matrix W');
if N ~= nobs/T
error('ols_panel_FE_g: wrong size matrix x');
end;
[nchk junk] = size(y);
if nchk ~= nobs
error('ols_panel_FE_g: wrong size vector y or matrix x');
end;

if (fe==1 & model==0 ) error('info.fe=1, but cannot compute fixed effects if info.model is set to 0 or not specified'); end


if novi_flag == 1
    results.homo = 1;
    results.hetero = 0;
    results.rval = 0;
    rval = 0;
elseif novi_flag == 0
    results.hetero = 1;
    results.homo = 0;
end


% demeaning of the y and x variables, depending on (info.)model

if (model==1 | model==3);
meanny=zeros(N,1);
meannwy=zeros(N,1);
meannx=zeros(N,nvar);
for i=1:N
    ym=zeros(T,1);
%     wym=zeros(T,1);
    xm=zeros(T,nvar);
    for t=1:T
        ym(t)=y(i+(t-1)*N,1);
%         wym(t)=Wy(i+(t-1)*N,1);
        xm(t,:)=x(i+(t-1)*N,:);
    end
    meanny(i)=mean(ym);
%     meannwy(i)=mean(wym);
    meannx(i,:)=mean(xm);
end
clear ym wym xm;
end % if statement

if ( model==2 | model==3)
meanty=zeros(T,1);
% meantwy=zeros(T,1);
meantx=zeros(T,nvar);
for i=1:T
    t1=1+(i-1)*N;t2=i*N;
    ym=y([t1:t2],1);
%     wym=Wy([t1:t2],1);
    xm=x([t1:t2],:);
    meanty(i)=mean(ym);
%     meantwy(i)=mean(wym);
    meantx(i,:)=mean(xm);
end
% clear ym wym xm;
end % if statement
    
en=ones(T,1);
et=ones(N,1);
ent=ones(nobs,1);

if model==1
    ywith=y-kron(en,meanny);
%     wywith=Wy-kron(en,meannwy);
    xwith=x-kron(en,meannx);
elseif model==2
    ywith=y-kron(meanty,et);
%     wywith=Wy-kron(meantwy,et);
    xwith=x-kron(meantx,et);
elseif model==3
    ywith=y-kron(en,meanny)-kron(meanty,et)+kron(ent,mean(y));
%     wywith=Wy-kron(en,meannwy)-kron(meantwy,et)+kron(ent,mean(Wy));
    xwith=x-kron(en,meannx)-kron(meantx,et)+kron(ent,mean(x));
else
    ywith=y;
%     wywith=Wy;
    xwith=x;
end % if statement

% ====== initializations

in = ones(N*T,1);
V = in;
vi = in;
sige = 1;
vmean = vi;

bsave = zeros(ndraw,nvar);
ssave = zeros(ndraw,1);

% step 1) do regressions
t0 = clock;
          AI = xwith'*xwith;
          b0 = AI\(xwith'*ywith);
%           bd = AI\(xwith'*wywith);
          e0 = ywith - xwith*b0;
%           ed = wywith - xwith*bd;
          epe0 = e0'*e0;
%           eped = ed'*ed;
%           epe0d = ed'*e0;
          
          % do MCMC sampling 
switch novi_flag
    
case{0} % we do heteroscedastic model    
results.hetero = 1;
%hwait = waitbar(0,'sar: MCMC sampling ...');

t0 = clock;                  
iter = 1;
          while (iter <= ndraw); % start sampling;
                  
          % update beta   
          xs = matmul(xwith,sqrt(V));
          ys = sqrt(V).*ywith;
          AI = (xs'*xs + sige*Q)\eye(nvar);         
          b = xs'*ys + sige*Qpc;
          b0 = (xs'*xs + sige*Q)\b;
          bhat = norm_rnd(sige*AI) + b0;  
          xb = xs*bhat;
                    
          % update sige
          nu1 = N*T + nu; 
          e = (ys - xb);
          d1 = d0 + e'*e;
          chi = chis_rnd(1,nu1);
          t2 = chi/d1;
          sige = 1/t2;

          % update vi
          ev = ywith  - xwith*bhat; 
          chiv = chis_rnd(N*T,rval+1);  
          %chiv = chi2rnd(rval+1,n,1); % Statistics Toolbox function
          vi = ((ev.*ev/sige) + in*rval)./chiv;
          V = in./vi; 
                        
          
              bsave(iter,1:nvar) = bhat';
              ssave(iter,1) = sige;
              vmean = vmean + vi;
              

iter = iter + 1; 
%waitbar(iter/ndraw);         
end; % end of sampling loop
%close(hwait);

time3 = etime(clock,t0);
results.time3 = time3;

case{1} % we do homoscedastic model 
    
%hwait = waitbar(0,'sar: MCMC sampling ...');
results.homo = 1;

t0 = clock;                  
iter = 1;

          while (iter <= ndraw); % start sampling;
                  
          % update beta   
          ys = ywith;        
          b = xwith'*ys + sige*Qpc;
          AI = (xwith'*xwith + sige*Q)\eye(nvar);         
          b0 = AI*b;
          bhat = norm_rnd(sige*AI) + b0;  
          xb = xwith*bhat;
                    
          % update sige
          nu1 = N*T + nu; 
          e = (ys - xb);
          d1 = d0 + e'*e;
          chi = chis_rnd(1,nu1);
          t2 = chi/d1;
          sige = 1/t2;          
          
%           if iter > nomit % if we are past burn-in, save the draws
              bsave(iter,1:nvar) = bhat';
              ssave(iter,1) = sige;
              vmean = vmean + vi;

%           end
          
iter = iter + 1; 
%waitbar(iter/ndraw);         
end % end of sampling loop
%close(hwait);


time3 = etime(clock,t0);
results.time3 = time3;


otherwise
error('ols_panel_FE_g: unrecognized prior.novi_flag value on input');
% we should never get here

end % end of homoscedastic vs. heteroscedastic vs. log-marginal options


% step 3) find b,sige MCMC estimates
% calculate mean of MCMC effects estimates

bmean = mean(bsave(nomit+1:end,:));
results.beta = bmean';
results.bdraw = bsave(nomit+1:end,:);
results.vmean = vmean/ndraw;
smean = mean(ssave(nomit+1:end,1));
results.sige = smean;
results.sdraw = ssave(nomit+1:end,:);
sige = results.sige;

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


% r-squared and corr-squared between actual and fitted values

results.tnvar=tnvar;
results.resid = y - xhat; 
yme=y-mean(y);
rsqr2=yme'*yme;
rsqr1 = results.resid'*results.resid;
results.rsqr=1.0-rsqr1/rsqr2; %rsquared

yhat=zeros(nobs,1);
ywithhat=zeros(nobs,1);
for t=1:T
    t1=1+(t-1)*N;t2=t*N;
    ywithhat(t1:t2,1)=xwith(t1:t2,:)*results.beta;
    yhat(t1:t2,1)=xhat(t1:t2,1);
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

mean_parms = bmean';
std_parms = bstd';
         
tstat = mean_parms./std_parms;

results.tstat = tstat;

% tvar = abs(diag(hessi));
% tmp = [results.beta
%        results.rho];
% results.tstat = tmp./sqrt(tvar(1:end-1,1));


parm = [results.beta
        results.sige];
% 
results.lik = f2_olspanel(parm,ywith,xwith,T); %LeSage
% 
% % Determination variance-covariance matrix

results.cov = cov(bsave);


% return stuff
results.nobs  = nobs; 
results.nvar  = nvar;
results.fe    = fe;
results.time  = etime(clock,timet);
results.N     = N;
results.T     = T;
results.model = model;




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


function llike = f2_olspanel(parm,y,x,T)
% PURPOSE: evaluates log-likelihood -- given estimates
%  ols panel model using sparse matrix algorithms
% ---------------------------------------------------
%  USAGE:llike = f2_olspanel(parm,y,x,,T)
%  where: parm = vector of maximum likelihood parameters
%                parm(1:k-1,1) = b, parm(k,1) = sige
%         y    = dependent variable vector (N x 1)
%         x    = explanatory variables matrix (N x k)
%         T    = number of time points
% ---------------------------------------------------
%  RETURNS: a  scalar equal to minus the log-likelihood
%           function value at the ML parameters
% ---------------------------------------------------

% Updated by J.P. Elhorst summer 2008
% REFERENCES: 
% Elhorst JP (2003) Specification and Estimation of Spatial Panel Data Models,
% International Regional Science Review 26: 244-268.
% Elhorst JP (2009) Spatial Panel Data Models. In Fischer MM, Getis A (Eds.) 
% Handbook of Applied Spatial Analysis, Ch. C.2. Springer: Berlin Heidelberg New York.

N = length(y)/T; 
k = length(parm);
b = parm(1:k-1,1);
sige = parm(k,1);


e=y-x*b;
for t=1:T
    t1=1+(t-1)*N;t2=t*N;
    e(t1:t2,1)=e(t1:t2,1);
end

epe = e'*e;
tmp2 = 1/(2*sige);
llike = -(N*T/2)*log(2*pi*sige)  - tmp2*epe;