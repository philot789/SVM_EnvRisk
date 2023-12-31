function results = sdm_panel_FE(y,x,W,T,info)
% PURPOSE: ML SDM model estimates for spatial panels 
%          (N regions*T time periods) with spatial fixed effects (u) 
%          and/or time period fixed effects (v)
%          y = p*W*y + X*b + W*X*c + u (optional) + v(optional) + e, using sparse matrix algorithms
% Supply data sorted first by time and then by spatial units, so first region 1,
% region 2, et cetera, in the first year, then region 1, region 2, et
% cetera in the second year, and so on
% sar_panel_FE computes y and x in deviation of the spatial and/or time means
% ---------------------------------------------------
%  USAGE: results = sdm_panel_FE(y,x,W,T,info)
%  where:  y = dependent variable vector
%          x = independent variables matrix
%          (the function adds W*x for you)
%          W = spatial weights matrix (standardized)
%          T = number of points in time
%       info = an (optional) structure variable with input options:
%       info.model = 0 pooled model without fixed effects (default, x may contain an intercept)
%                  = 1 spatial fixed effects (x may not contain an intercept)
%                  = 2 time period fixed effects (x may not contain an intercept)
%                  = 3 spatial and time period fixed effects (x may not contain an intercept)
%       info.fe    = report fixed effects and their t-values in prt_sp (default=0=not reported; info.fe=1=report) 
%       info.Nhes  = N =< Nhes asymptotic variance matrix is computed using analytical formulas,
%                    N > Nhes asymptotic variance matrix is computed using numerical formulas
%                    (Default NHes=500)
%       info.rmin  = (optional) minimum value of rho to use in search  
%       info.rmax  = (optional) maximum value of rho to use in search    
%       info.convg = (optional) convergence criterion (default = 1e-8)
%       info.maxit = (optional) maximum # of iterations (default = 500)
%       info.lflag = 0 for full lndet computation (default = 1, fastest)
%                  = 1 for MC lndet approximation (fast for very large problems)
%                  = 2 for Spline lndet approximation (medium speed)
%       info.order = order to use with info.lflag = 1 option (default = 50)
%       info.iter  = iterations to use with info.lflag = 1 option (default = 30)  
%       info.lndet = a matrix returned by sar containing log-determinant information to save time
% ---------------------------------------------------
%  RETURNS: a structure
%         results.meth  = 'psdm'   if info.model=0
%                       = 'sdmsfe' if info.model=1
%                       = 'sdmtfe' if info.model=2
%                       = 'sdmstfe' if info.model=3
%         results.beta  = bhat (includes coefficents c on the W*x variables)
%         results.rho   = rho (p above)
%         results.direct   = nvar x 5 matrix with direct effect, t-stat, t-prob, lower01, upper99
%         results.indirect = nvar x 5 matrix with indirect effect, t-stat, t-prob, lower01, upper99
%         results.total    = nvar x 5 matrix with total effect, t-stat, t-prob, lower01, upper99
%         results.direct_draws   = ndraw x nvar matrix of direct effect draws
%         results.indirect_draws = ndraw x nvar matrix of indirect effect draws
%         results.total_draws    = ndraw x nvar matrix of total effect draws
%         results.cov   = asymptotic variance-covariance matrix of the parameters beta and rho
%         results.tstat = asymp t-stat (last entry is rho=spatial autoregressive coefficient)
%         results.yhat  = [inv(y-p*W)]*[x*b+fixed effects] (according to prediction formula)
%         results.resid = y-p*W*y-x*b-W*x*c
%         results.sige  = (y-p*W*y-x*b-W*x*c)'*(y-p*W*y-x*b-W*x*c)/n
%         results.rsqr  = rsquared
%         results.corr2 = goodness-of-fit between actual and fitted values
%         results.sfe   = spatial fixed effects (if info.model=1 or 3)
%         results.tfe   = time period fixed effects (if info.model=2 or 3)
%         results.tsfe  = t-values spatial fixed effects (if info.model=1 or 3)
%         results.ttfe  = t-values time period fixed effects (if info.model=2 or 3)
%         results.con   = intercept 
%         results.con   = t-value intercept
%         results.nobs  = # of observations
%         results.nvar  = # of explanatory variables in x 
%         results.nvarw = # of explanatory variables in [x Wx]
%         results.tnvar = # fixed effects
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
%         results.time3 = time for hessian or information matrix calculation
%         results.time4 = time for optimization
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
%
% Updated by: J.Paul Elhorst summer 2008
% University of Groningen
% Department of Economics
% 9700AV Groningen
% the Netherlands
% j.p.elhorst@rug.nl
%
% REFERENCES: 
% Elhorst JP (2003) Specification and Estimation of Spatial Panel Data Models,
% International Regional Science Review 26: 244-268.
% Elhorst JP (2009) Spatial Panel Data Models. In Fischer MM, Getis A (Eds.) 
% Handbook of Applied Spatial Analysis, Ch. C.2. Springer: Berlin Heidelberg New York.

% This function is partly based on James. P LeSage's function SAR

time1 = 0; 
time2 = 0;
time3 = 0;
time4 = 0;

timet = clock; % start the clock for overall timing

W=sparse(W);

% if we have no options, invoke defaults
if nargin == 4
    info.lflag = 1;
    info.model=0;
    info.Nhes=500;
    fprintf(1,'default: pooled model without fixed effects \n');
end;

fe=0;
model=0;
Nhes=500;

fields = fieldnames(info);
nf = length(fields);
if nf > 0
    for i=1:nf
        if strcmp(fields{i},'model') model = info.model;
        elseif strcmp(fields{i},'fe') fe = info.fe;
        elseif strcmp(fields{i},'Nhes') Nhes = info.Nhes;
        end
    end
end
if model==0
    results.meth='psdm';
elseif model==1
    results.meth='sdmsfe';
elseif model==2
    results.meth='sdmtfe';
elseif model==3
    results.meth='sdmstfe';
else
    error('sdm_panel: wrong input number of info.model');
end

% check size of user inputs for comformability
[nobs nvar] = size(x);

[N Ncol] = size(W);
if N ~= Ncol
error('sdm: wrong size weight matrix W');
elseif N ~= nobs/T
error('sdm: wrong size weight matrix W or matrix x');
end;
[nchk junk] = size(y);
if nchk ~= nobs
error('sdm: wrong size vector y or matrix x');
end;

% check if the user handled the intercept term okay
    n = length(y);
    if sum(x(:,1)) ~= n
        tst = sum(x); % we may have no intercept term
        ind = find(tst == n); % we do have an intercept term
        if length(ind) > 0
            error('sdm_panel_FE: intercept term must be in first column of the x-matrix');
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
%[rmin,rmax,convg,maxit,detval,ldetflag,eflag,order,miter,options] = sar_parse(info); % function of LeSage
[rmin,rmax,convg,maxit,detval,ldetflag,eflag,order,miter,options,ndraw,sflag,p,cflag] = sar_parse(info);

% compute eigenvalues or limits
[rmin,rmax,time2] = sar_eigs(eflag,W,rmin,rmax,N); % function of LeSage

% do log-det calculations
[detval,time1] = sar_lndet(ldetflag,W,rmin,rmax,detval,order,miter); % function of LeSage

for t=1:T
    t1=1+(t-1)*N;t2=t*N;
    Wy(t1:t2,1)=W*y(t1:t2,1);
end

Wbig = kron(eye(T),W);

x = [x Wbig*x];
nvarw = size(x,2);
results.nvarw = nvarw;

% demeaning of the y and x variables, depending on (info.)model
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
meantx=zeros(T,nvarw);
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


if ( model==2 | model==3)
meanty=zeros(T,1);
meantwy=zeros(T,1);
meantx=zeros(T,nvarw);
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

% step 1) do regressions
t0 = clock;
          AI = xwith'*xwith;
          b0 = AI\(xwith'*ywith);
          bd = AI\(xwith'*wywith);
          e0 = ywith - xwith*b0;
          ed = wywith - xwith*bd;
          epe0 = e0'*e0;
          eped = ed'*ed;
          epe0d = ed'*e0;

% step 2) maximize concentrated likelihood function;
	options = optimset('fminbnd');
    [p,liktmp,exitflag,output] = fminbnd('f_sarpanel',rmin,rmax,options,detval,epe0,eped,epe0d,N,T);
   
time4 = etime(clock,t0);

if exitflag == 0
fprintf(1,'sdm_panel: convergence concentrated likelihood function not obtained in %4d iterations \n',output.iterations);
results.iter = 1;
end
results.iter = output.iterations;

% step 3) find b,sige maximum likelihood estimates
results.beta = b0 - p*bd; 
results.rho = p; 
bhat = results.beta;
results.sige = (1/nobs)*(e0-p*ed)'*(e0-p*ed); 
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
results.resid = y - p*Wy - xhat; 
yme=y-mean(y);
rsqr2=yme'*yme;
rsqr1 = results.resid'*results.resid;
results.rsqr=1.0-rsqr1/rsqr2; %rsquared

yhat=zeros(nobs,1);
ywithhat=zeros(nobs,1);
for t=1:T
    t1=1+(t-1)*N;t2=t*N;
    ywithhat(t1:t2,1)=(speye(N) - p*W)\xwith(t1:t2,:)*results.beta;
    yhat(t1:t2,1)=(speye(N) - p*W)\xhat(t1:t2,1);
end
res1=ywith-mean(ywith);
res2=ywithhat-mean(ywith);
rsq1=res1'*res2;
rsq2=res1'*res1;
rsq3=res2'*res2;
results.corr2=rsq1^2/(rsq2*rsq3); %corr2
results.yhat=yhat;

parm = [results.beta
        results.rho
        results.sige];

results.lik = f2_sarpanel(parm,ywith,xwith,W,detval,T); %Elhorst

% Determination variance-covariance matrix
if N <= Nhes % Analytically
t0 = clock;
B = speye(N) - p*W; 
BI = inv(B); WB = W*BI;
pterm = trace(WB*WB + WB'*WB);
xpx = zeros(nvarw+2,nvarw+2);               
% bhat,bhat
xpx(1:nvarw,1:nvarw) = (1/sige)*(xwith'*xwith);     
% bhat,rho
ysum=zeros(nvarw,1);
for t=1:T
    t1=1+(t-1)*N;t2=t*N;
    ysum=ysum+(1/sige)*xwith(t1:t2,:)'*WB*xwith(t1:t2,:)*bhat;
end
xpx(1:nvarw,nvarw+1) = ysum;
xpx(nvarw+1,1:nvarw) = xpx(1:nvarw,nvarw+1)'; 
% rho,rho
ysom=0;
for t=1:T
    t1=1+(t-1)*N;t2=t*N;
    ysom=ysom+(1/sige)*bhat'*xwith(t1:t2,:)'*WB'*WB*xwith(t1:t2,:)*bhat + pterm;
end
xpx(nvarw+1,nvarw+1) = ysom;
% sige, sige
xpx(nvarw+2,nvarw+2) = nobs/(2*sige*sige);     
% rho,sige
xpx(nvarw+1,nvarw+2) = (T/sige)*trace(WB);  
xpx(nvarw+2,nvarw+1) = xpx(nvarw+1,nvarw+2);
xpxi = xpx\eye(size(xpx));
results.cov=xpxi(1:nvarw+1,1:nvarw+1);
tmp = diag(xpxi(1:nvarw+1,1:nvarw+1));
results.hessi = results.cov;
bvec = [results.beta
        results.rho];
tmp = bvec./(sqrt(tmp));
results.tstat = tmp;
time3 = etime(clock,t0);

else  % asymptotic t-stats using numerical hessian
t0 = clock;
dhessn = hessian('f2_sarpanel',parm,ywith,xwith,W,detval,T); %Elhorst
hessi = invpd(-dhessn);
results.cov=hessi(1:nvarw+1,1:nvarw+1);
tvar = abs(diag(hessi));
tmp = [results.beta
       results.rho];
results.tstat = tmp./sqrt(tvar(1:end-1,1));
time3 = etime(clock,t0);
results.hessi = hessi;

end; % end of t-stat calculations


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

result0 = panel_effects_sdm(results,W,nvar,ywith,xwith);

results.direct = result0.direct;
results.indirect = result0.indirect;
results.total = result0.total;

results.direct_draws = result0.direct_draws;
results.indirect_draws = result0.indirect_draws;
results.total_draws = result0.total_draws;


% ===============================
% support functions below
% ===============================


function result0 = panel_effects_sdm(results,W,nvar,ywith,xwith)
% PURPOSE: computes and prints direct, indirect and total effects estimates
%          for Elhorst SAR spatial panel models using the LeSage and Pace code
%---------------------------------------------------
% USAGE: panel_effects_sar(results,W)
% Where: results    = a structure returned by a spatial panel regression
%        vnames     = a structure of variable names
%        W          = spatial weights matrix used to estimate model
%---------------------------------------------------
%
% Effects estimates added by Donald J. Lacombe
% Donald J. Lacombe
% Research Associate Professor
% Regional Research Institute
% 886 Chestnut Ridge Road
% PO Box 6825
% Morgantown, WV 26506-6825
% donald.lacombe@mail.wvu.edu
% http://www.rri.wvu.edu/lacombe/~lacombe.htm
%
% REFERENCES:
% Elhorst JP (2010) Matlab Software for Spatial Panels. Under review.

% LeSage and Pace code for calcualting effects in a SAR model
% Adapted for the Elhorst spatial panel models
ndraw=1000;
uiter=50;
maxorderu=100;
nobs = results.N;
rv=randn(nobs,uiter);
tracew=zeros(maxorderu,1);
wjjju=rv;
for jjj=1:maxorderu
    wjjju=W*wjjju;
    tracew(jjj)=mean(mean(rv.*wjjju));
    
end

traces=[tracew];
traces(1)=0;
traces(2)=sum(sum(W'.*W))/nobs;
trs=[1;traces];
ntrs=length(trs);
trbig=trs';
trbig2 = [trbig(1,2:end) trbig(1,end)];
trmat = [trbig
         trbig2];



% cheat here to fix the numerical hessian
% Use MCMC to get good results
sigmat = results.hessi - diag(diag(results.hessi)) + diag(diag(abs(results.hessi)));
% sigmatt = sigmat(1:end-1,1:end-1);
[R,posdef] = chol(sigmat);

if posdef ~= 0 % even cheating did not work, so punt with a kludge
    tmp = [xwith W*ywith]'*[xwith W*ywith];
    sigmat = sige*(inv(tmp));
end;

tmp = [results.beta
    results.rho];

bdraws = matadd(norm_rndmat(sigmat,ndraw),tmp);
draws = bdraws';

psave = draws(:,end);
ind = find(psave > 1); % find bad rho draws
psave(ind,1) = 0.99;   % replace them with 0.99


bsave = draws(:,1:end-1);
% cflag = 0;
% if cflag == 1
%     bdraws = bsave(:,2:end);
%     nvar = nvar-1;
% elseif cflag == 0
    bdraws = bsave;
% end;
pdraws = psave;

ree = 0:1:ntrs-1;

rmat = zeros(1,ntrs);
total = zeros(ndraw,nvar,ntrs);
direct = zeros(ndraw,nvar,ntrs);
indirect = zeros(ndraw,nvar,ntrs);


for i=1:ndraw;
    rmat = pdraws(i,1).^ree;
    for j=1:nvar;
        beta = [bdraws(i,j) bdraws(i,j+nvar)];
        total(i,j,:) = (beta(1,1) + beta(1,2))*rmat;
        direct(i,j,:) = (beta*trmat).*rmat;
        indirect(i,j,:) = total(i,j,:) - direct(i,j,:);
    end;
    
end;

total_draws = zeros(ndraw,nvar);
direct_draws = zeros(ndraw,nvar);
indirect_draws = zeros(ndraw,nvar);
for i=1:nvar;
tmp = squeeze(total(:,i,:)); % an ndraw by 1 by ntraces matrix
total_draws(:,i) = (sum(tmp'))'; % an ndraw by 1 vector
tmp = squeeze(indirect(:,i,:)); % an ndraw by 1 by ntraces matrix
indirect_draws(:,i) = (sum(tmp'))'; % an ndraw by 1 vector
tmp = squeeze(direct(:,i,:)); % an ndraw by 1 by ntraces matrix
direct_draws(:,i) = (sum(tmp'))'; % an ndraw by 1 vector
end;

result0.total_draws = total_draws;
result0.direct_draws = direct_draws;
result0.indirect_draws = indirect_draws;


% Compute means, std deviation and upper and lower 0.95 intervals
p = nvar;
total_out = zeros(p,5);
total_save = zeros(ndraw,p);
for i=1:p;
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
end;

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

result0.direct = direct_out;
result0.indirect = indirect_out;
result0.total = total_out;




function [rmin,rmax,convg,maxit,detval,ldetflag,eflag,order,iter,options,ndraw,sflag,p,cflag] = sar_parse(info)
% PURPOSE: parses input arguments for sar model
% ---------------------------------------------------
%  USAGE: [rmin,rmax,convg,maxit,detval,ldetflag,eflag,order,iter,options] = sar_parse(info)
% where info contains the structure variable with inputs 
% and the outputs are either user-inputs or default values
% ---------------------------------------------------

% set defaults
options = zeros(1,18); % optimization options for fminbnd
options(1) = 0; 
options(2) = 1.e-6; 
options(14) = 500;

eflag = 0;     % default to not computing eigenvalues
ldetflag = 1;  % default to 1999 Pace and Barry MC determinant approx
order = 50;    % there are parameters used by the MC det approx
iter = 30;     % defaults based on Pace and Barry recommendation
rmin = -1;     % use -1,1 rho interval as default
rmax = 1;
detval = 0;    % just a flag
convg = 0.0001;
maxit = 500;
ndraw = 1000;
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
    elseif strcmp(fields{i},'convg')
        options(2) = info.convg;
    elseif strcmp(fields{i},'maxit')
        options(14) = info.maxit;  
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
        error('sar: unrecognizable lflag value on input');
        end;
    elseif strcmp(fields{i},'order')
        order = info.order;  
    elseif strcmp(fields{i},'eig')
        eflag = info.eig;  
    elseif strcmp(fields{i},'iter')
        iter = info.iter; 
     elseif strcmp(fields{i},'ndraw')
        ndraw = info.ndraw; 
     elseif strcmp(fields{i},'sflag')
        sflag = info.sflag; 
    end;
 end;
 
else, % the user has input a blank info structure
      % so we use the defaults
end; 

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


