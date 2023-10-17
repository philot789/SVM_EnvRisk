function prt_panel(results,vnames,nnames,tnames,fid)
% PURPOSE: Prints output using spatial regression results structures
%---------------------------------------------------
% USAGE: prt_panel(results,vnames,se_names,te_names,fid)
% Where: results = a structure returned by a spatial panel regression 
%        vnames  = an optional vector of variable names
%        fid     = optional file-id for printing results to a file
%                  (defaults to the MATLAB command window)
%--------------------------------------------------- 
%  NOTES: e.g. vnames = strvcat('y','const','x1','x2');
%         e.g. fid = fopen('ols.out','wr');
%  use prt_sp(results,[],fid) to print to a file with no vnames               
% --------------------------------------------------
%  RETURNS: nothing, just prints the spatial panel regression results
% --------------------------------------------------

% written by James P LeSage (Texas State University)
% last updated 11/2020

if ~isstruct(results)
 error('prt_panel requires structure argument');
elseif nargin == 1
 nflag = 0; fid = 1; fflag = 0; tflag = 0;
elseif nargin == 2
 fid = 1; nflag = 1;  fflag = 0; tflag = 0;
elseif nargin == 3
 nflag = 0;  fflag = 1; tflag = 0; fid = 1;
 [vsize junk] = size(vnames); % user may supply a blank argument
   if vsize > 0
   nflag = 1;          
   end
   region_names = nnames;   
   tperiods = [];
   for i=1:results.T
       tperiods = strvcat(tperiods,['time' num2str(i)]);
   end
   
elseif nargin == 4 % user supplied names for units and time
   fflag = 1; fid = 1; nflag = 1; tflag = 1;
   region_names = nnames;
   tperiods = tnames;
elseif nargin == 5 % user a file name for printing
   tflag = 1; fflag = 1; nflag = 1; 
   tperiods = tnames;
   region_names = nnames;
   fid = fopen(fid,'w+');
else
 error('Wrong # of arguments to prt_panel');
end

effects_flag = 0;

nvar = results.nvar;
nobs = results.nobs;

% handling of vnames
Vname = 'Variable';
 for i=1:nvar
    tmp = ['variable ',num2str(i)];
    Vname = strvcat(Vname,tmp);
 end;

if (nflag == 1) % the user supplied variable names
[tst_n nsize] = size(vnames);
 if tst_n ~= nvar+1
 fprintf(fid,'Wrong # of variable names in prt_panel -- check vnames argument \n');
 fprintf(fid,'will use generic variable names \n');
 nflag = 0;
 else,
Vname = 'Variable';
 for i=1:nvar
    Vname = strvcat(Vname,vnames(i+1,:));
 end;
 end; % end of if-else
end; % end of nflag issue

fprintf(fid,'\n');

switch results.meth
    
case {'psem_conv_g','semp_conv_g','semsfe_conv_g','semtfe_conv_g','semstfe_conv_g'} % <==convex combination fixed effects spatial lag models

if strcmp(results.meth,'psem_conv_g')
fprintf(fid,'MCMC SEM convex combination W model with no fixed effects \n');
end
if strcmp(results.meth,'semsfe_conv_g')
fprintf(fid,'MCMC SEM convex combination W model with region fixed effects \n');
end
if strcmp(results.meth,'semtfe_conv_g')
fprintf(fid,'MCMC SEM convex combination W model with time period fixed effects \n');
end
if strcmp(results.meth,'semstfe_conv_g')
fprintf(fid,'MCMC SEM convex combination W model with both region and time period fixed effects \n');
end

% add spatial rho parameter name
    Vname = strvcat(Vname,'rho');
fprintf(fid,'Homoscedastic model \n');
nmat = results.nmat;

m = results.nmat;
    for ii=1:m
            Vname = strvcat(Vname,['gamma' num2str(ii)]);
    end;

    nobs = results.nobs;
    nvar = results.nvar;
    % extract posterior means
    bout = [results.beta_mean
            results.rho_mean
            results.gamma_mean];

    sige = results.sige_mean;
        tmp1 = std(results.bdraw);
        tmp2 = std(results.pdraw);
        tmp3 = std(results.gdraw);

        bstd = [tmp1'
                tmp2
                tmp3']; 

    tout = bout./bstd;
results.tflag = 'tstat';
    if strcmp(results.tflag,'tstat')
        tstat = bout./bstd;
        % find t-stat marginal probabilities
        tout = tdis_prb(tstat,results.nobs);
        results.tstat = bout./bstd; % trick for printing below
    else % find plevels
        bout = plims(results.bdraw);
        pout = plims(results.pdraw);
        gout = plims(results.gdraw);
        bout = [bout
                pout
                gout];


    end

    rsqr = results.rsqr;


    fprintf(fid,'\n');
    fprintf(fid,'Bayesian spatial error convex W model \n');
    if (nflag == 1)
        fprintf(fid,'Dependent Variable  = %16s \n',vnames(1,:));
    end;
    fprintf(fid,'Log-marginal likeli = %9.4f \n',results.logmarginal);

     cstats2 = chainstats(results.drawpost);
     fprintf(fid,'Log-marginal MCerror= %9.6f\n',cstats2(1,2));
    fprintf(fid,'R-squared           = %9.4f \n',rsqr);
    fprintf(fid,'Rbar-squared        = %9.4f \n',results.rbar);
    fprintf(fid,'mean of sige draws  = %9.4f \n',mean(results.sdraw));
%     fprintf(fid,'posterior mode sige = %9.4f \n',results.sig_mode);
    fprintf(fid,'Nobs, Nvars         = %6d,%6d \n',results.nobs,results.nvar);
    fprintf(fid,'ndraws,nomit        = %6d,%6d \n',results.ndraw,results.nomit);
fprintf(fid,'total time in secs = %9.4f   \n',results.time);
% fprintf(fid,'time for effects   = %9.4f \n',results.effects_time);
fprintf(fid,'time for sampling  = %9.4f \n',results.sampling_time);
fprintf(fid,'time for Taylor    = %9.4f \n',results.taylor_time);
% fprintf(fid,'log-marginal likelihood  = %9.4f \n',results.logmarginal);
fprintf(fid,'min and max lambda = %9.4f,%9.4f \n',results.rmin,results.rmax);

    fprintf(fid,'***************************************************************\n');
    nd = (results.ndraw-results.nomit)/results.thin;
    fprintf(fid,'      MCMC diagnostics ndraws = %d \n',nd);
    cstats=chainstats([results.bdraw results.pdraw results.gdraw]);
    in.cnames = strvcat('mean','MC error','tau','Geweke');
    in.rnames = Vname;
    in.fmt = strvcat('%12.4f','12.8f','%10.6f','%10.6f');
    in.width = 10000;
    out = [cstats];
    mprint(out,in);

    fprintf(fid,'***************************************************************\n');
    fprintf(fid,'      Posterior Estimates \n');

    effects_flag = 0;
   
% =============================== end of sem convex
   
case {'psdem_conv_g','sdemp_conv_g','sdemsfe_conv_g','sdemtfe_conv_g','sdemstfe_conv_g'} % <==convex combination fixed effects spatial lag models

if strcmp(results.meth,'psdem_conv_g')
fprintf(fid,'MCMC SDEM convex combination W model with no fixed effects \n');
end
if strcmp(results.meth,'sdemsfe_conv_g')
fprintf(fid,'MCMC SDEM convex combination W model with region fixed effects \n');
end
if strcmp(results.meth,'sdemtfe_conv_g')
fprintf(fid,'MCMC SDEM convex combination W model with time period fixed effects \n');
end
if strcmp(results.meth,'sdemstfe_conv_g')
fprintf(fid,'MCMC SDEM convex combination W model with both region and time period fixed effects \n');
end


% add W-xnames
nvar = results.nvar;
m = results.nmat;
    
if (nflag == 1) % the user supplied variable names
    Vname = 'Variable';
    for i=1:nvar
        Vname = strvcat(Vname,vnames(i+1,:));
    end
    for j=1:m
        if results.cflag == 0
            for i=1:nvar
                Vname = strvcat(Vname,['W' num2str(j) '*' vnames(i+1,:)]);
            end
        elseif results.cflag == 1
            for i=2:nvar
                Vname = strvcat(Vname,['W' num2str(j) '*' vnames(i+1,:)]);
            end
            
        end
        % add spatial rho parameter name
        %         Vname = strvcat(Vname,'rho');
    end
    
end % end of nflag issue



% add spatial rho parameter name
    Vname = strvcat(Vname,'rho');
fprintf(fid,'Homoscedastic model \n');
nmat = results.nmat;

m = results.nmat;
    for ii=1:m
            Vname = strvcat(Vname,['gamma' num2str(ii)]);
    end;

    nobs = results.nobs;
    % extract posterior means
    bout = [results.beta_mean
            results.rho_mean
            results.gamma_mean];

    sige = results.sige_mean;
        tmp1 = std(results.bdraw);
        tmp2 = std(results.pdraw);
        tmp3 = std(results.gdraw);

        bstd = [tmp1'
                tmp2
                tmp3']; 

    tout = bout./bstd;
results.tflag = 'tstat';
    if strcmp(results.tflag,'tstat')
        tstat = bout./bstd;
        % find t-stat marginal probabilities
        tout = tdis_prb(tstat,results.nobs);
        results.tstat = bout./bstd; % trick for printing below
    else % find plevels
        bout = plims(results.bdraw);
        pout = plims(results.pdraw);
        gout = plims(results.gdraw);
        bout = [bout
                pout
                gout];


    end

    rsqr = results.rsqr;


    fprintf(fid,'\n');
    fprintf(fid,'Bayesian spatial Durbin error convex W model \n');
    if (nflag == 1)
        fprintf(fid,'Dependent Variable  = %16s \n',vnames(1,:));
    end;
    fprintf(fid,'Log-marginal      = %9.4f \n',results.logmarginal);
     cstats2 = chainstats(results.drawpost);
     fprintf(fid,'Log-marginal MCerror= %9.6f\n',cstats2(1,2));
    fprintf(fid,'R-squared           = %9.4f \n',rsqr);
    fprintf(fid,'Rbar-squared        = %9.4f \n',results.rbar);
    fprintf(fid,'mean of sige draws  = %9.4f \n',mean(results.sdraw));
%     fprintf(fid,'posterior mode sige = %9.4f \n',results.sig_mode);
    fprintf(fid,'Nobs, Nvars         = %6d,%6d \n',results.nobs,results.nvar);
    fprintf(fid,'ndraws,nomit        = %6d,%6d \n',results.ndraw,results.nomit);
fprintf(fid,'total time in secs = %9.4f   \n',results.time);
% fprintf(fid,'time for effects   = %9.4f \n',results.effects_time);
fprintf(fid,'time for sampling  = %9.4f \n',results.sampling_time);
fprintf(fid,'time for Taylor    = %9.4f \n',results.taylor_time);
% fprintf(fid,'log-marginal time  = %9.4f \n',results.logm_time);
fprintf(fid,'min and max lambda = %9.4f,%9.4f \n',results.rmin,results.rmax);

    fprintf(fid,'***************************************************************\n');
    nd = (results.ndraw-results.nomit)/results.thin;
    fprintf(fid,'      MCMC diagnostics ndraws = %d \n',nd);
    cstats=chainstats([results.bdraw results.pdraw results.gdraw]);
    in.cnames = strvcat('mean','MC error','tau','Geweke');
    in.rnames = Vname;
    in.fmt = strvcat('%12.4f','12.8f','%10.6f','%10.6f');
    in.width = 10000;
    out = [cstats];
    mprint(out,in);

    fprintf(fid,'***************************************************************\n');
    fprintf(fid,'      Posterior Estimates \n');

    
total_out    = results.total;
indirect_out = results.indirect;
direct_out   = results.direct;

    effects_flag = 1;
   % end of sdem convex
        
case {'psar_conv_g','sarsfe_conv_g','sartfe_conv_g','sarstfe_conv_g'} % <==convex combination fixed effects spatial lag models

if strcmp(results.meth,'psar_conv_g')
fprintf(fid,'MCMC SAR convex combination W model with no fixed effects \n');
end
if strcmp(results.meth,'sarsfe_conv_g')
fprintf(fid,'MCMC SAR convex combination W model with region fixed effects \n');
end
if strcmp(results.meth,'sartfe_conv_g')
fprintf(fid,'MCMC SAR convex combination W model with time period fixed effects \n');
end
if strcmp(results.meth,'sarstfe_conv_g')
fprintf(fid,'MCMC SAR convex combination W model with both region and time period fixed effects \n');
end

% add spatial rho parameter name
    Vname = strvcat(Vname,'rho');
fprintf(fid,'Homoscedastic model \n');
m = results.nmat;
    for ii=1:m
            Vname = strvcat(Vname,['gamma' num2str(ii)]);
    end;

    nobs = results.nobs;
    nvar = results.nvar;
    % extract posterior means
    bout = [results.beta
            results.rho
            results.gamma];

    sige = results.sige;
        tmp1 = std(results.bdraw);
        tmp2 = std(results.pdraw);
        tmp3 = std(results.gdraw);

        bstd = [tmp1'
                tmp2
                tmp3']; 

    tout = bout./bstd;
results.tflag = 'tstat';
    if strcmp(results.tflag,'tstat')
        tstat = bout./bstd;
        % find t-stat marginal probabilities
        tout = tdis_prb(tstat,results.nobs);
        results.tstat = bout./bstd; % trick for printing below
    else % find plevels
        bout = plims(results.bdraw);
        pout = plims(results.pdraw);
        gout = plims(results.gdraw);
        bout = [bout
                pout
                gout];


    end

    rsqr = results.rsqr;


    fprintf(fid,'\n');
    fprintf(fid,'Bayesian spatial autoregressive convex W model \n');
    if (nflag == 1)
    fprintf(fid,'Dependent Variable  = %16s \n',vnames(1,:));
    end;
    fprintf(fid,'Log-marginal        = %9.4f \n',results.logmarginal);
    cstats2 = chainstats(results.drawpost);
    fprintf(fid,'Log-marginal MCerror= %9.6f\n',cstats2(1,2));
    fprintf(fid,'R-squared           = %9.4f \n',rsqr);
    fprintf(fid,'corr-squared        = %9.4f   \n',results.corr2);
    fprintf(fid,'mean of sige draws  = %9.4f \n',mean(results.sdraw));
    fprintf(fid,'posterior mode sige = %9.4f \n',results.sig_mode);
    fprintf(fid,'Nobs, Nvars         = %6d,%6d \n',results.nobs,results.nvar);
    fprintf(fid,'ndraws,nomit        = %6d,%6d \n',results.ndraw,results.nomit);
    % fprintf(fid,'total time in secs  = %9.4f   \n',results.time);
    fprintf(fid,'time for effects    = %9.4f \n',results.effects_time);
    fprintf(fid,'time for sampling   = %9.4f \n',results.sampling_time);
    fprintf(fid,'time for Taylor     = %9.4f \n',results.taylor_time);
    fprintf(fid,'thinning for draws  = %9d   \n',results.thin);
    fprintf(fid,'min and max rho     = %9.4f,%9.4f \n',results.rmin,results.rmax);

    fprintf(fid,'***************************************************************\n');
    nd = (results.ndraw-results.nomit)/results.thin;
    fprintf(fid,'      MCMC diagnostics ndraws = %d \n',nd);
    cstats=chainstats([results.bdraw results.pdraw results.gdraw]);
    in.cnames = strvcat('mode','mean','MC error','tau','Geweke');
    in.rnames = Vname;
    in.fmt = strvcat('%12.4f','%12.4f','12.8f','%10.6f','%10.6f');
    in.width = 10000;
    pmode = [results.beta_mode results.rho_mode results.gamma_mode]';
    out = [pmode cstats];
    mprint(out,in);

    fprintf(fid,'***************************************************************\n');
    fprintf(fid,'      Posterior Estimates \n');
    
    
total_out    = results.total;
indirect_out = results.indirect;
direct_out   = results.direct;

effects_flag = 1;

% ============= end of convex combination SAR models    
    
case {'psdm_conv_g','sdmsfe_conv_g','sdmtfe_conv_g','sdmstfe_conv_g'} % <==convex combination fixed effects spatial lag models

if strcmp(results.meth,'psdm_conv_g')
fprintf(fid,'MCMC SDM convex combination W model with no fixed effects \n');
end
if strcmp(results.meth,'sdmsfe_conv_g')
fprintf(fid,'MCMC SDM convex combination W model with region fixed effects \n');
end
if strcmp(results.meth,'sdmtfe_conv_g')
fprintf(fid,'MCMC SDM convex combination W model with time period fixed effects \n');
end
if strcmp(results.meth,'sdmstfe_conv_g')
fprintf(fid,'MCMC SDM convex combination W model with both region and time period fixed effects \n');
end

% add W-xnames
nvar = results.nvar;
m = results.m;
    
if (nflag == 1) % the user supplied variable names
    Vname = 'Variable';
    for i=1:nvar
        Vname = strvcat(Vname,vnames(i+1,:));
    end
    for j=1:m
        if results.cflag == 0
            for i=1:nvar
                Vname = strvcat(Vname,['W' num2str(j) '*' vnames(i+1,:)]);
            end
        elseif results.cflag == 1
            for i=2:nvar
                Vname = strvcat(Vname,['W' num2str(j) '*' vnames(i+1,:)]);
            end
            
        end
        % add spatial rho parameter name
        %         Vname = strvcat(Vname,'rho');
    end
    
end % end of nflag issue

% add spatial rho parameter name
    Vname = strvcat(Vname,'rho');
fprintf(fid,'Homoscedastic model \n');
m = results.nmat;
    for ii=1:m
            Vname = strvcat(Vname,['gamma' num2str(ii)]);
    end;

    nobs = results.nobs;
    nvar = results.nvar;
    % extract posterior means
    bout = [results.beta
            results.rho
            results.gamma];

    sige = results.sige;
        tmp1 = std(results.bdraw);
        tmp2 = std(results.pdraw);
        tmp3 = std(results.gdraw);

        bstd = [tmp1'
                tmp2
                tmp3']; 

    tout = bout./bstd;
results.tflag = 'tstat';
    if strcmp(results.tflag,'tstat')
        tstat = bout./bstd;
        % find t-stat marginal probabilities
        tout = tdis_prb(tstat,results.nobs);
        results.tstat = bout./bstd; % trick for printing below
    else % find plevels
        bout = plims(results.bdraw);
        pout = plims(results.pdraw);
        gout = plims(results.gdraw);
        bout = [bout
                pout
                gout];


    end
    
    rsqr = results.rsqr;

    fprintf(fid,'\n');
    fprintf(fid,'Bayesian spatial Durbin convex W model \n');
    if (nflag == 1)
    fprintf(fid,'Dependent Variable  = %16s \n',vnames(1,:));
    end;
    fprintf(fid,'Log-marginal        = %9.4f \n',results.logmarginal);
    cstats2 = chainstats(results.drawpost);
    fprintf(fid,'Log-marginal MCerror= %9.6f\n',cstats2(1,2));
    fprintf(fid,'R-squared           = %9.4f \n',rsqr);
    fprintf(fid,'corr-squared        = %9.4f   \n',results.corr2);
    fprintf(fid,'mean of sige draws  = %9.4f \n',mean(results.sdraw));
    fprintf(fid,'posterior mode sige = %9.4f \n',results.sig_mode);
    fprintf(fid,'Nobs, Nvars         = %6d,%6d \n',results.nobs,results.nvar);
    fprintf(fid,'ndraws,nomit        = %6d,%6d \n',results.ndraw,results.nomit);
    % fprintf(fid,'total time in secs  = %9.4f   \n',results.time);
    fprintf(fid,'time for effects    = %9.4f \n',results.effects_time);
    fprintf(fid,'time for sampling   = %9.4f \n',results.sampling_time);
    fprintf(fid,'time for Taylor     = %9.4f \n',results.taylor_time);
    fprintf(fid,'thinning for draws  = %9d   \n',results.thin);
    fprintf(fid,'min and max rho     = %9.4f,%9.4f \n',results.rmin,results.rmax);

    fprintf(fid,'***************************************************************\n');
    nd = (results.ndraw-results.nomit)/results.thin;
    fprintf(fid,'      MCMC diagnostics ndraws = %d \n',nd);
    cstats=chainstats([results.bdraw results.pdraw results.gdraw]);
    in.cnames = strvcat('mode','mean','MC error','tau','Geweke');
    in.rnames = Vname;
    in.fmt = strvcat('%12.4f','%12.4f','12.8f','%10.6f','%10.6f');
    in.width = 10000;
    pmode = [results.beta_mode results.rho_mode results.gamma_mode]';
    out = [pmode cstats];
    mprint(out,in);

    fprintf(fid,'***************************************************************\n');
    fprintf(fid,'      Posterior Estimates \n');

total_out    = results.total;
indirect_out = results.indirect;
direct_out   = results.direct;

effects_flag = 1;

    
% ============= end of convex combination SDM models    
    
case {'psem','semsfe','semtfe','semstfe'} % <=================== fixed effects Max Like spatial error models

if strcmp(results.meth,'psem')
fprintf(fid,'MaxLike SEM model with no fixed effects \n');
end
if strcmp(results.meth,'semsfe')
fprintf(fid,'MaxLike SEM model with region fixed effects \n');
end
if strcmp(results.meth,'semtfe')
fprintf(fid,'MaxLike SEM model with time period fixed effects \n');
end
if strcmp(results.meth,'semstfe')
fprintf(fid,'MaxLike SEM model with both region and time period fixed effects \n');
end

% add spatial rho parameter name
    Vname = strvcat(Vname,'rho');
fprintf(fid,'Homoscedastic model \n');
    
if (nflag == 1)
fprintf(fid,'Dependent Variable = %16s \n',vnames(1,:));
end;
fprintf(fid,'R-squared           = %9.4f   \n',results.rsqr);
fprintf(fid,'corr-squared        = %9.4f   \n',results.corr2);
fprintf(fid,'sigma^2             = %9.4f   \n',results.sige);
% fprintf(fid,'log-likelihood  = %16.8g  \n',results.lik);
fprintf(fid,'Nobs,Nvar,#FE       = %6d,%6d,%6d  \n',results.nobs,results.nvar,results.tnvar);
% fprintf(fid,'# iterations        = %6d     \n',results.iter);
fprintf(fid,'min and max rho     = %9.4f,%9.4f \n',results.rmin,results.rmax);
% print timing information
fprintf(fid,'total time in secs  = %9.4f \n',results.time);
% fprintf(fid,'time for optimiz   = %9.4f \n',results.time4);
if results.time1 ~= 0
fprintf(fid,'time for lndet      = %9.4f \n',results.time1);
end;
if results.time2 ~= 0
fprintf(fid,'time for eigs       = %9.4f \n',results.time2);
end;
% if results.time3 ~= 0
% fprintf(fid,'time for MCMC draws = %9.4f \n',results.time3);
% end;

if results.lflag == 0
fprintf(fid,'No lndet approximation used \n');
end;
% put in information regarding Pace and Barry approximations
if results.lflag == 1
fprintf(fid,'Pace and Barry, 1999 MC lndet approximation used \n');
fprintf(fid,'order for MC appr  = %6d  \n',results.order);
fprintf(fid,'iter  for MC appr  = %6d  \n',results.miter);
end;
if results.lflag == 2
fprintf(fid,'Pace and Barry, 1998 spline lndet approximation used \n');
end;

fprintf(fid,'***************************************************************\n');

bout = [results.beta
        results.rho];
        
% <=================== end of sem case

case {'psem_g','semsfe_g','semtfe_g','semstfe_g'} % <=================== fixed effects Bayesian spatial error models

if results.homo == 1
fprintf(fid,'Homoscedastic model \n');
results.rval = 0;
end
if  results.hetero == 1
fprintf(fid,'Heteroscedastic model \n');
end

if strcmp(results.meth,'psem_g')
fprintf(fid,'MCMC SEM model with no fixed effects \n');
end

if strcmp(results.meth,'semsfe_g')
fprintf(fid,'MCMC SEM model with region fixed effects \n');
end
if strcmp(results.meth,'semtfe_g')
fprintf(fid,'MCMC SEM model with time period fixed effects \n');
end
if strcmp(results.meth,'semstfe_g')
fprintf(fid,'MCMC SEM model with both region and time period fixed effects \n');
end


% add spatial rho parameter name
    Vname = strvcat(Vname,'rho');
    
if (nflag == 1)
fprintf(fid,'Dependent Variable = %16s \n',vnames(1,:));
end;
fprintf(fid,'R-squared           = %9.4f   \n',results.rsqr);
fprintf(fid,'corr-squared        = %9.4f   \n',results.corr2);
fprintf(fid,'sigma^2             = %9.4f   \n',results.sige);
% fprintf(fid,'log-likelihood  = %16.8g  \n',results.lik);
fprintf(fid,'Nobs,Nvar,#FE       = %6d,%6d,%6d  \n',results.nobs,results.nvar,results.tnvar);
fprintf(fid,'prior rvalue       = %6d   \n',results.rval);
% fprintf(fid,'# iterations        = %6d     \n',results.iter);
fprintf(fid,'min and max rho     = %9.4f,%9.4f \n',results.rmin,results.rmax);
fprintf(fid,'ndraws,nomit        = %6d,%6d \n',results.ndraw,results.nomit);
% print timing information
fprintf(fid,'total time in secs  = %9.4f \n',results.time);
% % fprintf(fid,'time for optimiz   = %9.4f \n',results.time4);
% if results.time1 ~= 0
% fprintf(fid,'time for lndet      = %9.4f \n',results.time1);
% end;
% if results.time2 ~= 0
% fprintf(fid,'time for eigs       = %9.4f \n',results.time2);
% end;
if results.time3 ~= 0
fprintf(fid,'time for MCMC draws = %9.4f \n',results.time3);
end;

if results.lflag == 0
fprintf(fid,'No lndet approximation used \n');
end;
% put in information regarding Pace and Barry approximations
if results.lflag == 1
fprintf(fid,'Pace and Barry, 1999 MC lndet approximation used \n');
fprintf(fid,'order for MC appr  = %6d  \n',results.order);
fprintf(fid,'iter  for MC appr  = %6d  \n',results.iter);
end;
if results.lflag == 2
fprintf(fid,'Pace and Barry, 1998 spline lndet approximation used \n');
end;

if results.iprior == 1
fprintf(fid,'***************************************************************\n');

    vstring = 'Variable';
    bstring = 'Prior Mean';
    tstring = 'Std Deviation';
    
    tmp = [results.bmean results.bstd];
    cnames = strvcat(bstring,tstring);
    rnames = vstring;
    nvarw = results.nvar;
%     Vname
    for i=1:nvarw
        rnames = strvcat(rnames,Vname(i+1,:));
    end;
%     rnames
    pin.fmt = '%16.6f';
    pin.fid = fid;
    pin.cnames = cnames;
    pin.rnames = rnames;
    mprint(tmp,pin);
end
fprintf(fid,'***************************************************************\n');

effects_flag = 0;

bout = [results.beta
        results.rho];
        
% <=================== end of sem_g case

case {'pols_g','olssfe_g','olstfe_g','olsstfe_g'} % <=================== fixed effects Bayesian ols models

    if results.homo == 1
fprintf(fid,'Homoscedastic model \n');
end
if  results.hetero == 1
fprintf(fid,'Heteroscedastic model \n');
end

if strcmp(results.meth,'pols_g')
fprintf(fid,'MCMC OLS model with no fixed effects \n');
end

if strcmp(results.meth,'olssfe_g')
fprintf(fid,'MCMC OLS model with region fixed effects \n');
end
if strcmp(results.meth,'olstfe_g')
fprintf(fid,'MCMC OLS model with time period fixed effects \n');
end
if strcmp(results.meth,'olsstfe_g')
fprintf(fid,'MCMC OLS model with both region and time period fixed effects \n');
end

    
if (nflag == 1)
fprintf(fid,'Dependent Variable = %16s \n',vnames(1,:));
end;
fprintf(fid,'R-squared           = %9.4f   \n',results.rsqr);
fprintf(fid,'corr-squared        = %9.4f   \n',results.corr2);
fprintf(fid,'sigma^2             = %9.4f \n',results.sige);
fprintf(fid,'Nobs,Nvar,#FE       = %6d,%6d,%6d  \n',results.nobs,results.nvar,results.tnvar); % +1 due to spatially lagged dependent variable
fprintf(fid,'log-likelihood     = %16.8g \n',results.lik);
fprintf(fid,'prior rvalue       = %6d   \n',results.rval);
% print timing information
fprintf(fid,'total time in secs  = %9.4f \n',results.time);
fprintf(fid,'ndraws,nomit        = %6d,%6d \n',results.ndraw,results.nomit);
if results.time3 ~= 0
fprintf(fid,'time for MCMC draws = %9.4f \n',results.time3);
end;



if results.iprior == 1
fprintf(fid,'***************************************************************\n');

    vstring = 'Variable';
    bstring = 'Prior Mean';
    tstring = 'Std Deviation';
    
    tmp = [results.bmean results.bstd];
    cnames = strvcat(bstring,tstring);
    rnames = vstring;
    nvarw = results.nvar;
%     Vname
    for i=1:nvarw
        rnames = strvcat(rnames,Vname(i+1,:));
    end;
%     rnames
    pin.fmt = '%16.6f';
    pin.fid = fid;
    pin.cnames = cnames;
    pin.rnames = rnames;
    mprint(tmp,pin);
end
fprintf(fid,'***************************************************************\n');

effects_flag = 0;

bout = results.beta;
        
% <=================== end of ols case

case {'pslx_g','slxsfe_g','slxtfe_g','slxstfe_g'} % <=================== fixed effects Bayesian spatial lag of X models


if results.homo == 1    
fprintf(fid,'Homoscedastic model \n');
end
if results.hetero ==1
 fprintf(fid,'Heterocedastic model \n');
end
   
if strcmp(results.meth,'pslx_g')
fprintf(fid,'MCMC SLX model with no fixed effects \n');
end
if strcmp(results.meth,'slxsfe_g')
fprintf(fid,'MCMC SLX model with region fixed effects \n');
end
if strcmp(results.meth,'slxtfe_g')
fprintf(fid,'MCMC SLX model with time period fixed effects \n');
end
if strcmp(results.meth,'slxstfe_g')
fprintf(fid,'MCMC SLX model with both region and time period fixed effects \n');
end
    
nvar = results.nvar;
    
if (nflag == 1) % the user supplied variable names
    Vname = 'Variable';
    for i=1:nvar
        Vname = strvcat(Vname,vnames(i+1,:));
    end
    for i=1:nvar
        Vname = strvcat(Vname,['W-' vnames(i+1,:)]);
    end
elseif (nflag == 0)
    for i=1:nvar
        Vname = strvcat(Vname,['W-' Vname(i+1,:)]);
    end
    
end
    

if (nflag == 1)
fprintf(fid,'Dependent Variable = %16s \n',vnames(1,:));
end;
fprintf(fid,'R-squared           = %9.4f   \n',results.rsqr);
fprintf(fid,'corr-squared        = %9.4f   \n',results.corr2);
fprintf(fid,'sigma^2             = %9.4f \n',results.sige);
fprintf(fid,'Nobs,Nvar,#FE       = %6d,%6d,%6d  \n',results.nobs,results.nvar+1,results.tnvar); % +1 due to spatially lagged dependent variable
fprintf(fid,'log-likelihood     = %16.8g \n',results.lik);
fprintf(fid,'prior rvalue       = %6d   \n',results.rval);
% fprintf(fid,'min and max rho     = %9.4f,%9.4f \n',results.rmin,results.rmax);
% print timing information
fprintf(fid,'total time in secs  = %9.4f \n',results.time);
% fprintf(fid,'time for optimiz   = %9.4f \n',results.time4);
% if results.time1 ~= 0
% fprintf(fid,'time for lndet      = %9.4f \n',results.time1);
% end;
% if results.time2 ~= 0
% fprintf(fid,'time for eigs       = %9.4f \n',results.time2);
% end;
if results.time3 ~= 0
fprintf(fid,'time for MCMC draws = %9.4f \n',results.time3);
end;



if results.iprior == 1
fprintf(fid,'***************************************************************\n');

    vstring = 'Variable';
    bstring = 'Prior Mean';
    tstring = 'Std Deviation';
    
    tmp = [results.bmean results.bstd];
    cnames = strvcat(bstring,tstring);
    rnames = vstring;
    nvarw = results.nvarw;
%     Vname
    for i=1:nvarw
        rnames = strvcat(rnames,Vname(i+1,:));
    end;
%     rnames
    pin.fmt = '%16.6f';
    pin.fid = fid;
    pin.cnames = cnames;
    pin.rnames = rnames;
    mprint(tmp,pin);
end

fprintf(fid,'***************************************************************\n');
    
total_out    = results.total;
indirect_out = results.indirect;
direct_out   = results.direct;

effects_flag = 1;

bout = results.beta;

    

case {'psar','sarsfe','sartfe','sarstfe'} % <=================== fixed effects Max Like spatial lag models

    
fprintf(fid,'Homoscedastic model \n');
if strcmp(results.meth,'psar')
fprintf(fid,'MaxLike SAR model with no fixed effects \n');
end
if strcmp(results.meth,'sarsfe')
fprintf(fid,'MaxLike SAR model with region fixed effects \n');
end
if strcmp(results.meth,'sartfe')
fprintf(fid,'MaxLike SAR model with time period fixed effects \n');
end
if strcmp(results.meth,'sarstfe')
fprintf(fid,'MaxLike SAR model with both region and time period fixed effects \n');
end

% add spatial rho parameter name
    Vname = strvcat(Vname,'rho');
    
if (nflag == 1)
fprintf(fid,'Dependent Variable = %16s \n',vnames(1,:));
end;
fprintf(fid,'R-squared           = %9.4f   \n',results.rsqr);
fprintf(fid,'corr-squared        = %9.4f   \n',results.corr2);
fprintf(fid,'sigma^2             = %9.4f \n',results.sige);
fprintf(fid,'Nobs,Nvar,#FE       = %6d,%6d,%6d  \n',results.nobs,results.nvar+1,results.tnvar); % +1 due to spatially lagged dependent variable
fprintf(fid,'log-likelihood     = %16.8g \n',results.lik);
fprintf(fid,'# of iterations    = %6d   \n',results.iter);
fprintf(fid,'min and max rho     = %9.4f,%9.4f \n',results.rmin,results.rmax);
% print timing information
fprintf(fid,'total time in secs  = %9.4f \n',results.time);
% fprintf(fid,'time for optimiz   = %9.4f \n',results.time4);
if results.time1 ~= 0
fprintf(fid,'time for lndet      = %9.4f \n',results.time1);
end;
if results.time2 ~= 0
fprintf(fid,'time for eigs       = %9.4f \n',results.time2);
end;
if results.time3 ~= 0
fprintf(fid,'time for MCMC draws = %9.4f \n',results.time3);
end;

if results.lflag == 0
fprintf(fid,'No lndet approximation used \n');
end;
% put in information regarding Pace and Barry approximations
if results.lflag == 1
fprintf(fid,'Pace and Barry, 1999 MC lndet approximation used \n');
fprintf(fid,'order for MC appr  = %6d  \n',results.order);
fprintf(fid,'iter  for MC appr  = %6d  \n',results.miter);
end;
if results.lflag == 2
fprintf(fid,'Pace and Barry, 1998 spline lndet approximation used \n');
end;

fprintf(fid,'***************************************************************\n');



bout = [results.beta
        results.rho];

total_out    = results.total;
indirect_out = results.indirect;
direct_out   = results.direct;

effects_flag = 1;

case {'psdm','sdmsfe','sdmtfe','sdmstfe'} % <=================== fixed effects Max Like SDM models

    
fprintf(fid,'Homoscedastic model \n');
if strcmp(results.meth,'psdm')
fprintf(fid,'MaxLike SDM model with no fixed effects \n');
end
if strcmp(results.meth,'sdmsfe')
fprintf(fid,'MaxLike SDM model with region fixed effects \n');
end
if strcmp(results.meth,'sdmtfe')
fprintf(fid,'MaxLike SDM model with time period fixed effects \n');
end
if strcmp(results.meth,'sdmstfe')
fprintf(fid,'MaxLike SDM model with both region and time period fixed effects \n');
end

nvar = results.nvar;
    
if (nflag == 1) % the user supplied variable names
    Vname = 'Variable';
    for i=1:nvar
        Vname = strvcat(Vname,vnames(i+1,:));
    end
    for i=1:nvar
        Vname = strvcat(Vname,['W-' vnames(i+1,:)]);
    end
    %     add spatial rho parameter name
    Vname = strvcat(Vname,'rho');
elseif (nflag == 0)
    for i=1:nvar
        Vname = strvcat(Vname,['W-' Vname(i+1,:)]);
    end
    %     add spatial rho parameter name
    Vname = strvcat(Vname,'rho');
    
end

if (nflag == 1)
fprintf(fid,'Dependent Variable = %16s \n',vnames(1,:));
end;
fprintf(fid,'R-squared           = %9.4f   \n',results.rsqr);
fprintf(fid,'corr-squared        = %9.4f   \n',results.corr2);
fprintf(fid,'sigma^2             = %9.4f \n',results.sige);
fprintf(fid,'Nobs,Nvar,#FE       = %6d,%6d,%6d  \n',results.nobs,results.nvarw,results.tnvar); % +1 due to spatially lagged dependent variable
fprintf(fid,'log-likelihood     = %16.8g \n',results.lik);
fprintf(fid,'# of iterations    = %6d   \n',results.iter);
fprintf(fid,'min and max rho     = %9.4f,%9.4f \n',results.rmin,results.rmax);
% print timing information
fprintf(fid,'total time in secs  = %9.4f \n',results.time);
% fprintf(fid,'time for optimiz   = %9.4f \n',results.time4);
if results.time1 ~= 0
fprintf(fid,'time for lndet      = %9.4f \n',results.time1);
end;
if results.time2 ~= 0
fprintf(fid,'time for eigs       = %9.4f \n',results.time2);
end;
% if results.time3 ~= 0
% fprintf(fid,'time for MCMC draws = %9.4f \n',results.time3);
% end;

if results.lflag == 0
fprintf(fid,'No lndet approximation used \n');
end;
% put in information regarding Pace and Barry approximations
if results.lflag == 1
fprintf(fid,'Pace and Barry, 1999 MC lndet approximation used \n');
fprintf(fid,'order for MC appr  = %6d  \n',results.order);
fprintf(fid,'iter  for MC appr  = %6d  \n',results.miter);
end;
if results.lflag == 2
fprintf(fid,'Pace and Barry, 1998 spline lndet approximation used \n');
end;

 
total_out    = results.total;
indirect_out = results.indirect;
direct_out   = results.direct;

effects_flag = 1;

fprintf(fid,'***************************************************************\n');

bout = [results.beta
        results.rho];


case {'psar_g','sarsfe_g','sartfe_g','sarstfe_g'} % <=================== fixed effects Bayesian spatial lag models

if results.homo == 1    
fprintf(fid,'Homoscedastic model \n');
end
if results.hetero ==1
 fprintf(fid,'Heterocedastic model \n');
end
   
if strcmp(results.meth,'psar_g')
fprintf(fid,'MCMC SAR model with no fixed effects \n');
end
if strcmp(results.meth,'sarsfe_g')
fprintf(fid,'MCMC SAR model with region fixed effects \n');
end
if strcmp(results.meth,'sartfe_g')
fprintf(fid,'MCMC SAR model with time period fixed effects \n');
end
if strcmp(results.meth,'sarstfe_g')
fprintf(fid,'MCMC SAR model with both region and time period fixed effects \n');
end

% add spatial rho parameter name
    Vname = strvcat(Vname,'rho');
    
if (nflag == 1)
fprintf(fid,'Dependent Variable = %16s \n',vnames(1,:));
end;
fprintf(fid,'R-squared           = %9.4f   \n',results.rsqr);
fprintf(fid,'corr-squared        = %9.4f   \n',results.corr2);
fprintf(fid,'sigma^2             = %9.4f \n',results.sige);
fprintf(fid,'Nobs,Nvar,#FE       = %6d,%6d,%6d  \n',results.nobs,results.nvar,results.tnvar); 
% fprintf(fid,'log-likelihood     = %16.8g \n',results.lik);
fprintf(fid,'ndraw,nomit         = %6d,%6d   \n',results.ndraw,results.nomit);
fprintf(fid,'rvalue              = %6d       \n',results.rval);
fprintf(fid,'min and max rho     = %9.4f,%9.4f \n',results.rmin,results.rmax);
% print timing information
fprintf(fid,'total time in secs  = %9.4f \n',results.time);
% fprintf(fid,'time for optimiz   = %9.4f \n',results.time4);
if results.time1 ~= 0
fprintf(fid,'time for lndet      = %9.4f \n',results.time1);
end;
if results.time2 ~= 0
fprintf(fid,'time for eigs       = %9.4f \n',results.time2);
end;
if results.time3 ~= 0
fprintf(fid,'time for MCMC draws = %9.4f \n',results.time3);
end;

if results.lflag == 0
fprintf(fid,'No lndet approximation used \n');
end;
% put in information regarding Pace and Barry approximations
if results.lflag == 1
fprintf(fid,'Pace and Barry, 1999 MC lndet approximation used \n');
fprintf(fid,'order for MC appr  = %6d  \n',results.order);
fprintf(fid,'iter  for MC appr  = %6d  \n',results.miter);
end;
if results.lflag == 2
fprintf(fid,'Pace and Barry, 1998 spline lndet approximation used \n');
end;


if results.iprior == 1
    fprintf(fid,'***************************************************************\n');

    vstring = 'Variable';
    bstring = 'Prior Mean';
    tstring = 'Std Deviation';
    
    tmp = [results.bmean results.bstd];
    cnames = strvcat(bstring,tstring);
    rnames = vstring;
    nvarw = results.nvar;
%     Vname
    for i=1:nvarw
        rnames = strvcat(rnames,Vname(i+1,:));
    end
%     rnames
    pin.fmt = '%16.6f';
    pin.fid = fid;
    pin.cnames = cnames;
    pin.rnames = rnames;
    mprint(tmp,pin);
end

bout = [results.beta
        results.rho];

fprintf(fid,'***************************************************************\n');

total_out    = results.total;
indirect_out = results.indirect;
direct_out   = results.direct;

effects_flag = 1;

%     
% <=================== end of sar_g case
case {'psdm_g','sdmsfe_g','sdmtfe_g','sdmstfe_g'} % <=================== fixed effects Bayesian spatial Durbin models


if results.homo == 1    
fprintf(fid,'Homoscedastic model \n');
end
if results.hetero ==1
 fprintf(fid,'Heterocedastic model \n');
end
   
if strcmp(results.meth,'psdm_g')
fprintf(fid,'MCMC SDM model with no fixed effects \n');
end
if strcmp(results.meth,'sdmsfe_g')
fprintf(fid,'MCMC SDM model with region fixed effects \n');
end
if strcmp(results.meth,'sdmtfe_g')
fprintf(fid,'MCMC SDM model with time period fixed effects \n');
end
if strcmp(results.meth,'sdmstfe_g')
fprintf(fid,'MCMC SDM model with both region and time period fixed effects \n');
end
    
nvar = results.nvar;
    
if (nflag == 1) % the user supplied variable names
    Vname = 'Variable';
    for i=1:nvar
        Vname = strvcat(Vname,vnames(i+1,:));
    end
    if results.cflag == 0
        for i=1:nvar
            Vname = strvcat(Vname,['W*' vnames(i+1,:)]);
        end
    elseif results.cflag == 1
        for i=2:nvar
            Vname = strvcat(Vname,['W*' vnames(i+1,:)]);
        end
    end
    
    %     add spatial rho parameter name
    Vname = strvcat(Vname,'rho');
end
    if (nflag == 0)
        for i=1:nvar
            Vname = strvcat(Vname,['W-' Vname(i+1,:)]);
        end
        %     add spatial rho parameter name
        Vname = strvcat(Vname,'rho');
        
    end
    

if (nflag == 1)
fprintf(fid,'Dependent Variable = %16s \n',vnames(1,:));
end;
fprintf(fid,'R-squared           = %9.4f   \n',results.rsqr);
fprintf(fid,'corr-squared        = %9.4f   \n',results.corr2);
fprintf(fid,'sigma^2             = %9.4f \n',results.sige);
fprintf(fid,'Nobs,Nvar,#FE       = %6d,%6d,%6d  \n',results.nobs,results.nvar,results.tnvar); 
fprintf(fid,'ndraw,nomit         = %6d,%6d   \n',results.ndraw,results.nomit);
fprintf(fid,'rvalue              = %6d       \n',results.rval);
% fprintf(fid,'log-likelihood     = %16.8g \n',results.lik);
% fprintf(fid,'# of iterations    = %6d   \n',results.iter);
fprintf(fid,'min and max rho     = %9.4f,%9.4f \n',results.rmin,results.rmax);
% print timing information
fprintf(fid,'total time in secs  = %9.4f \n',results.time);
% fprintf(fid,'time for optimiz   = %9.4f \n',results.time4);
%if results.time1 ~= 0
%fprintf(fid,'time for lndet      = %9.4f \n',results.time1);
%end;
%if results.time2 ~= 0
%fprintf(fid,'time for eigs       = %9.4f \n',results.time2);
%end;
if results.time3 ~= 0
fprintf(fid,'time for MCMC draws = %9.4f \n',results.time3);
end;

if results.lflag == 0
fprintf(fid,'No lndet approximation used \n');
end;
% put in information regarding Pace and Barry approximations
if results.lflag == 1
fprintf(fid,'Pace and Barry, 1999 MC lndet approximation used \n');
fprintf(fid,'order for MC appr  = %6d  \n',results.order);
fprintf(fid,'iter  for MC appr  = %6d  \n',results.miter);
end;
if results.lflag == 2
fprintf(fid,'Pace and Barry, 1998 spline lndet approximation used \n');
end;


if results.iprior == 1
    fprintf(fid,'***************************************************************\n');

    vstring = 'Variable';
    bstring = 'Prior Mean';
    tstring = 'Std Deviation';
    
    tmp = [results.bmean results.bstd];
    cnames = strvcat(bstring,tstring);
    rnames = vstring;
    nvarw = results.nvarw;
%     Vname
    for i=1:nvarw
        rnames = strvcat(rnames,Vname(i+1,:));
    end
%     rnames
    pin.fmt = '%16.6f';
    pin.fid = fid;
    pin.cnames = cnames;
    pin.rnames = rnames;
    mprint(tmp,pin);
end

fprintf(fid,'***************************************************************\n');

    
total_out    = results.total;
indirect_out = results.indirect;
direct_out   = results.direct;

effects_flag = 1;


bout = [results.beta
        results.rho];

case {'psdem_g','sdemsfe_g','sdemtfe_g','sdemstfe_g'} % <=================== fixed effects Bayesian spatial Durbin error models


if results.homo == 1    
fprintf(fid,'Homoscedastic model \n');
results.rval = 0;
end
if results.hetero ==1
 fprintf(fid,'Heterocedastic model \n');
end
   
if strcmp(results.meth,'psdem_g')
fprintf(fid,'MCMC SDEM model with no fixed effects \n');
end
if strcmp(results.meth,'sdemsfe_g')
fprintf(fid,'MCMC SDEM model with region fixed effects \n');
end
if strcmp(results.meth,'sdemtfe_g')
fprintf(fid,'MCMC SDEM model with time period fixed effects \n');
end
if strcmp(results.meth,'sdemstfe_g')
fprintf(fid,'MCMC SDEM model with both region and time period fixed effects \n');
end
    
nvar = results.nvar;
    
if (nflag == 1) % the user supplied variable names
    Vname = 'Variable';
    for i=1:nvar
        Vname = strvcat(Vname,vnames(i+1,:));
    end
        if results.cflag == 0
        for i=1:nvar
            Vname = strvcat(Vname,['W*' vnames(i+1,:)]);
        end
    elseif results.cflag == 1
        for i=2:nvar
            Vname = strvcat(Vname,['W*' vnames(i+1,:)]);
        end
    end

    %     add spatial rho parameter name
    Vname = strvcat(Vname,'rho');
end
if (nflag == 0)
    for i=1:nvar
        Vname = strvcat(Vname,['W-' Vname(i+1,:)]);
    end
    %     add spatial rho parameter name
    Vname = strvcat(Vname,'rho');
    
end
    

if (nflag == 1)
fprintf(fid,'Dependent Variable = %16s \n',vnames(1,:));
end;
fprintf(fid,'R-squared           = %9.4f   \n',results.rsqr);
fprintf(fid,'corr-squared        = %9.4f   \n',results.corr2);
fprintf(fid,'sigma^2             = %9.4f \n',results.sige);
fprintf(fid,'Nobs,Nvar,#FE       = %6d,%6d,%6d  \n',results.nobs,results.nvar,results.tnvar); 
fprintf(fid,'ndraw,nomit         = %6d,%6d   \n',results.ndraw,results.nomit);
fprintf(fid,'rvalue              = %6d       \n',results.rval);
% fprintf(fid,'log-likelihood     = %16.8g \n',results.lik);
% fprintf(fid,'# of iterations    = %6d   \n',results.iter);
fprintf(fid,'min and max rho     = %9.4f,%9.4f \n',results.rmin,results.rmax);
% print timing information
fprintf(fid,'total time in secs  = %9.4f \n',results.time);
% fprintf(fid,'time for optimiz   = %9.4f \n',results.time4);
if results.time1 ~= 0
fprintf(fid,'time for lndet      = %9.4f \n',results.time1);
end;
if results.time2 ~= 0
fprintf(fid,'time for eigs       = %9.4f \n',results.time2);
end;
if results.time3 ~= 0
fprintf(fid,'time for MCMC draws = %9.4f \n',results.time3);
end;

if results.lflag == 0
fprintf(fid,'No lndet approximation used \n');
end;
% put in information regarding Pace and Barry approximations
if results.lflag == 1
fprintf(fid,'Pace and Barry, 1999 MC lndet approximation used \n');
fprintf(fid,'order for MC appr  = %6d  \n',results.order);
fprintf(fid,'iter  for MC appr  = %6d  \n',results.iter);
end;
if results.lflag == 2
fprintf(fid,'Pace and Barry, 1998 spline lndet approximation used \n');
end;


if results.iprior == 1
    fprintf(fid,'***************************************************************\n');

    vstring = 'Variable';
    bstring = 'Prior Mean';
    tstring = 'Std Deviation';
    
    tmp = [results.bmean results.bstd];
    cnames = strvcat(bstring,tstring);
    rnames = vstring;
    nvarw = results.nvarw;
%     Vname
    for i=1:nvarw
        rnames = strvcat(rnames,Vname(i+1,:));
    end
%     rnames
    pin.fmt = '%16.6f';
    pin.fid = fid;
    pin.cnames = cnames;
    pin.rnames = rnames;
    mprint(tmp,pin);
end

fprintf(fid,'***************************************************************\n');

    
total_out    = results.total;
indirect_out = results.indirect;
direct_out   = results.direct;

effects_flag = 1;

bout = [results.beta
        results.rho];
% <=================== end of sdem_g case
    
case {'psdem','sdemsfe','sdemtfe','sdemstfe'} % <=================== fixed effects Max Like spatial Durbin error models

if strcmp(results.meth,'psdem')
fprintf(fid,'MaxLike SDEM model with no fixed effects \n');
end
if strcmp(results.meth,'sdemsfe')
fprintf(fid,'MaxLike SDEM model with region fixed effects \n');
end
if strcmp(results.meth,'sdemtfe')
fprintf(fid,'MaxLike SDEM model with time period fixed effects \n');
end
if strcmp(results.meth,'sdemstfe')
fprintf(fid,'MaxLike SDEM model with both region and time period fixed effects \n');
end

nvar = results.nvar;
    
if (nflag == 1) % the user supplied variable names
    Vname = 'Variable';
    for i=1:nvar
        Vname = strvcat(Vname,vnames(i+1,:));
    end
    for i=1:nvar
        Vname = strvcat(Vname,['W-' vnames(i+1,:)]);
    end
elseif (nflag == 0)
    for i=1:nvar
        Vname = strvcat(Vname,['W-' Vname(i+1,:)]);
    end
    
end

% add spatial rho parameter name
    Vname = strvcat(Vname,'rho');
    
    
fprintf(fid,'Homoscedastic model \n');
    
if (nflag == 1)
fprintf(fid,'Dependent Variable = %16s \n',vnames(1,:));
end;
fprintf(fid,'R-squared           = %9.4f   \n',results.rsqr);
fprintf(fid,'corr-squared        = %9.4f   \n',results.corr2);
fprintf(fid,'sigma^2             = %9.4f   \n',results.sige);
% fprintf(fid,'log-likelihood  = %16.8g  \n',results.lik);
fprintf(fid,'Nobs,Nvar,#FE       = %6d,%6d,%6d  \n',results.nobs,results.nvar,results.tnvar);
% fprintf(fid,'# iterations        = %6d     \n',results.iter);
fprintf(fid,'min and max rho     = %9.4f,%9.4f \n',results.rmin,results.rmax);
% print timing information
fprintf(fid,'total time in secs  = %9.4f \n',results.time);
% fprintf(fid,'time for optimiz   = %9.4f \n',results.time4);
if results.time1 ~= 0
fprintf(fid,'time for lndet      = %9.4f \n',results.time1);
end;
if results.time2 ~= 0
fprintf(fid,'time for eigs       = %9.4f \n',results.time2);
end;
% if results.time3 ~= 0
% fprintf(fid,'time for MCMC draws = %9.4f \n',results.time3);
% end;

if results.lflag == 0
fprintf(fid,'No lndet approximation used \n');
end;
% put in information regarding Pace and Barry approximations
if results.lflag == 1
fprintf(fid,'Pace and Barry, 1999 MC lndet approximation used \n');
fprintf(fid,'order for MC appr  = %6d  \n',results.order);
fprintf(fid,'iter  for MC appr  = %6d  \n',results.miter);
end;
if results.lflag == 2
fprintf(fid,'Pace and Barry, 1998 spline lndet approximation used \n');
end;

fprintf(fid,'***************************************************************\n');

bout = [results.beta
        results.rho];
        
% <=================== end of sdem case

  
otherwise
error('results structure not known by prt_panel function');
end

% now print coefficient estimates, t-statistics and probabilities
tout = norm_prb(results.tstat); % find asymptotic z (normal) probabilities, function of LeSage
tmp = [bout results.tstat tout];  % matrix to be printed
% column labels for printing results
bstring = 'Coefficient'; tstring = 'Asymptot t-stat'; pstring = 'z-probability';
cnames = strvcat(bstring,tstring,pstring);
in.cnames = cnames;
in.rnames = Vname;
in.fmt = '%16.6f';
in.fid = fid;
mprint(tmp,in); %function of LeSage


if effects_flag == 1
% now print x-effects estimates

bstring = 'Coefficient'; 
tstring = 't-stat'; 
pstring = 't-prob';
lstring = 'lower 05';
ustring = 'upper 95';
cnames = strvcat(bstring,tstring,pstring,lstring,ustring);
ini.cnames = cnames;
ini.width = 2000;

nvar = results.nvar;
if results.cflag == 1
 vnameso = strvcat(Vname(3:nvar+1,:));    
elseif results.cflag == 0
 vnameso = strvcat(Vname(2:nvar+1,:));    
end

% print effects estimates
ini.rnames = strvcat('Direct  ',vnameso);
ini.fmt = '%16.6f';
ini.fid = fid;

% set up print out matrix
printout = direct_out;
mprint(printout,ini);

printout = indirect_out;
ini.rnames = strvcat('Indirect',vnameso);
mprint(printout,ini);

printout = total_out;
ini.rnames = strvcat('Total   ',vnameso);
mprint(printout,ini);

end

% print fixed effects, t-statistics and probabilities at request of the researcher
if (results.fe==1)
model=results.model;
N=results.N;
T=results.T;
if (model==1)
    fprintf(fid,'Mean intercept and region fixed effects \n');
    FEname=strvcat('Variable','intercept');
    if fflag == 0
    for i=1:N
        tmp = ['sfe ',num2str(i)];
        FEname = strvcat(FEname,tmp);
    end
    bout=[results.con;results.sfe];
    tstat=[results.tcon;results.tsfe];
    tout = norm_prb(tstat);
    tmp = [bout tstat tout];
    in.rnames = FEname;
    mprint(tmp,in);
    elseif fflag == 1
    for i=1:N
        FEname = strvcat(FEname,region_names(i,:));
    end
    bout=[results.con;results.sfe];
    tstat=[results.tcon;results.tsfe];
    tout = norm_prb(tstat);
    tmp = [bout tstat tout];
    in.rnames = FEname;
    mprint(tmp,in);
    end        
elseif (model==2)
    fprintf(fid,'Mean intercept and time period fixed effects \n');
    FEname=strvcat('Variable','intercept');
    if tflag == 0
    for i=1:T
        tmp = ['tfe ',num2str(i)];
        FEname = strvcat(FEname,tmp);
    end
    bout=[results.con;results.tfe];
    tstat=[results.tcon;results.ttfe];
    tout = norm_prb(tstat);
    tmp = [bout tstat tout];
    in.rnames = FEname;
    mprint(tmp,in);
    elseif tflag == 1
    for i=1:T
        FEname = strvcat(FEname,tperiods(i,:));
    end
    bout=[results.con;results.tfe];
    tstat=[results.tcon;results.ttfe];
    tout = norm_prb(tstat);
    tmp = [bout tstat tout];
    in.rnames = FEname;
    mprint(tmp,in);
    end        
elseif (model==3)
    fprintf(fid,'Mean intercept, region and time period fixed effects \n');
    FEname=strvcat('Variable','intercept');
    if fflag == 0
    for i=1:N
        tmp = ['sfe ',num2str(i)];
        FEname = strvcat(FEname,tmp);
    end
    elseif fflag == 1
    for i=1:N
        FEname = strvcat(FEname,region_names(i,:));
    end
    end     
    if tflag == 0
    for i=1:T
        tmp = ['tfe ',num2str(i)];
        FEname = strvcat(FEname,tmp);
    end
    elseif tflag == 1
      for i=1:T
        FEname = strvcat(FEname,tperiods(i,:));
    end
    end      
    bout=[results.con;results.sfe;results.tfe];
    tstat=[results.tcon;results.tsfe;results.ttfe];
    tout = norm_prb(tstat);
    tmp = [bout tstat tout];
    in.rnames = FEname;
    mprint(tmp,in);
end
end


function stats=chainstats(chain)
% MCMC convergence statistics 
% USAGE: stats = chainstats(ndraws x nparms, matrix of mcmc draws)
% RETURNS:
% stats = [mean, mcerr, tau, p], nparms x 4 matrix
% mean = vector of means nparms x 1
% mcerr = Monte Carlo estimate of the Monte Carlo std of the estimates, nparms x 1 
% tau = integrated autocorrelation time, Sokal's adaptive truncated periodogram, nparms x 1 
% p = Geweke's MCMC diagnostic Test for equality of the means of the first 10% and last 50% draws

% $Revision: 1.4 $  $Date: 2009/08/13 15:47:35 $

mcerr = bmstd(chain)./sqrt(size(chain,1));

[z,p]  = geweke(chain);
tau    = iact(chain);
stats  = [mean(chain)',mcerr',tau', p'];
% mcerr = Monte Carlo estimate of the Monte Carlo std of the estimates 
% tau = integrated autocorrelation time using Sokal's adaptive truncated periodogram estimator.
% p = Geweke's MCMC convergence diagnostic Test for equality of the means of the first 10% and last 50% of a Markov chain.


function [z,p]=geweke(chain,a,b)
%GEWEKE Geweke's MCMC convergence diagnostic
% [z,p] = geweke(chain,a,b)
% Test for equality of the means of the first a% (default 10%) and
% last b% (50%) of a Markov chain.
% See:
% Stephen P. Brooks and Gareth O. Roberts.
% Assessing convergence of Markov chain Monte Carlo algorithms.
% Statistics and Computing, 8:319--335, 1998.

% ML, 2002
% $Revision: 1.3 $  $Date: 2003/05/07 12:22:19 $

[nsimu,npar]=size(chain);

if nargin<3
  a = 0.1;
  b = 0.5;
end

na = floor(a*nsimu);
nb = nsimu-floor(b*nsimu)+1;

if (na+nb)/nsimu >= 1
  error('Error with na and nb');
end

m1 = mean(chain(1:na,:));
m2 = mean(chain(nb:end,:));

%%% Spectral estimates for variance
sa = spectrum0(chain(1:na,:));
sb = spectrum0(chain(nb:end,:));

z = (m1-m2)./(sqrt(sa/na+sb/(nsimu-nb+1)));
p = 2*(1-nordf(abs(z)));


function s=bmstd(x,b)
%BMSTD standard deviation calculated from batch means
% s = bmstd(x,b) - x matrix - b length of the batch
% bmstd(x) gives an estimate of the Monte Carlo std of the 
% estimates calculated from x

% Marko Laine <marko.laine@fmi.fi>
% $Revision: 1.2 $  $Date: 2012/09/27 11:47:34 $

[n,p] = size(x);

if nargin<2
  b = max(10,fix(n/20));
end

inds = 1:b:(n+1);
nb = length(inds)-1;
if nb < 2
  error('too few batches');
end

y = zeros(nb,p);

for i=1:nb
  y(i,:)=mean(x(inds(i):inds(i+1)-1,:));
end

% calculate the estimated std of MC estimate
s = sqrt( sum((y - repmat(mean(x),nb,1)).^2)/(nb-1)*b );

function s=spectrum0(x)
%SPECTRUM0 Spectral density at frequency zero
% spectrum0(x) spectral density at zero for columns of x

% ML, 2002
% $Revision: 1.3 $  $Date: 2003/05/07 12:22:19 $

[m,n]= size(x);
s = zeros(1,n);
for i=1:n
  spec = spectrum(x(:,i),m);
  s(i) = spec(1);
end


function [y,f]=spectrum(x,nfft,nw)
%SPECTRUM Power spectral density using Hanning window
%  [y,f]=spectrum(x,nfft,nw) 

% See also: psd.m in Signal Processing Toolbox 

% Marko Laine <marko.laine@fmi.fi>
% $Revision: 1.4 $  $Date: 2012/09/27 11:47:40 $

if nargin < 2 | isempty(nfft)
  nfft = min(length(x),256);
end
if nargin < 3 | isempty(nw)
  nw = fix(nfft/4);
end
noverlap = fix(nw/2);

% Hanning window
w = .5*(1 - cos(2*pi*(1:nw)'/(nw+1)));
% Daniel
%w = [0.5;ones(nw-2,1);0.5];
n = length(x);
if n < nw
    x(nw)=0;  n=nw;
end
x = x(:);

k = fix((n-noverlap)/(nw-noverlap)); % no of windows
index = 1:nw;
kmu = k*norm(w)^2; % Normalizing scale factor
y = zeros(nfft,1);
for i=1:k
% xw = w.*detrend(x(index),'linear');
  xw = w.*x(index);
  index = index + (nw - noverlap);
  Xx = abs(fft(xw,nfft)).^2;
  y = y + Xx;
end

y = y*(1/kmu); % normalize

n2 = floor(nfft/2);
y  = y(1:n2);
f  = 1./n*(0:(n2-1));


function y=nordf(x,mu,sigma2)
% NORDF the standard normal (Gaussian) cumulative distribution.
% NORPF(x,mu,sigma2) x quantile, mu mean, sigma2 variance

% Marko Laine <marko.laine@fmi.fi>
% $Revision: 1.6 $  $Date: 2012/09/27 11:47:38 $

if nargin < 2, mu     = 0; end
if nargin < 3, sigma2 = 1; end

%y = 0.5*erf(-Inf,sqrt(2)*0.5*x);
y = 0.5+0.5*erf((x-mu)/sqrt(sigma2)/sqrt(2));


function [tau,m] = iact(dati)
%IACT estimates the integrated autocorrelation time
%   using Sokal's adaptive truncated periodogram estimator.

% Originally contributed by Antonietta Mira by name sokal.m

% Marko Laine <marko.laine@fmi.fi>
% $Revision: 1.2 $  $Date: 2012/09/27 11:47:37 $

if length(dati) == prod(size(dati))
  dati = dati(:);
end

[mx,nx] = size(dati);
tau = zeros(1,nx);
m   = zeros(1,nx);

x  = fft(dati);
xr = real(x);
xi = imag(x);
xr = xr.^2+xi.^2; %%%controllare questo
xr(1,:)=0;
xr=real(fft(xr));
var=xr(1,:)./length(dati)/(length(dati)-1);

for j = 1:nx
  if var(j) == 0
    continue
  end
  xr(:,j)=xr(:,j)./xr(1,j);
  sum=-1/3;
  for i=1:length(dati)
    sum=sum+xr(i,j)-1/6;
    if sum<0
      tau(j)=2*(sum+(i-1)/6);
      m(j)=i;
      break
    end
  end
end

