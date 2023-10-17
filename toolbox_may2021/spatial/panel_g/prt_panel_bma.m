function prt_panel_bma(results,vnames,nnames,tnames,fid)
% PURPOSE: Prints BMA output using results structures
%---------------------------------------------------
% USAGE: prt_sar_conv(results,vnames,fid)
% Where: results = a structure returned by a sar_conv_g model
%        vnames  = an optional vector of variable names
%        fid     = optional file-id for printing results to a file
%                  (defaults to the MATLAB command window)
%        fmt     = format string, e.g., '%12.4f' (default)
%--------------------------------------------------- 
%  NOTES: e.g. vnames = strvcat('y','const','x1','x2');
%         e.g. fid = fopen('ols.out','wr');
%  use prt_spat(results,[],fid) to print to a file with no vnames               
% --------------------------------------------------
%  RETURNS: nothing, just prints the sar_conv_g estimation results
% --------------------------------------------------
% SEE ALSO: prt, plt
%---------------------------------------------------   

% written by:
% James P. LeSage, 4/2018
% Dept of Finance & Economics
% Texas State University-San Marcos
% 601 University Drive
% San Marcos, TX 78666
% jlesage@spatial-econometrics.com

if ~isstruct(results)
 error('prt_panel_bma requires structure argument');
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
 error('Wrong # of arguments to prt_panel_bma');
end


% 
nvar = results.nvar;
nmat = results.nmat; % Number of considered connectivy matrices.

if (nflag == 1) % the user supplied variable names
[tst_n nsize] = size(vnames);
 if tst_n ~= nvar+1
 fprintf(fid,'Wrong # of variable names in prt__panel_bma -- check vnames argument \n');
 fprintf(fid,'will use generic variable names \n');
 nflag = 0;
 end
end;

% handling of vnames
Vname = 'Variable';
if nflag == 0 % no user-supplied vnames or an incorrect vnames argument
    
    for i=1:nvar
        tmp = ['variable ',num2str(i)];
        Vname = strvcat(Vname,tmp);
    end;
    
    
% add spatial rho parameter name
    Vname = strvcat(Vname,'rho');
    
% add gamma parameter names
for ii=1:nmat
        Vname = strvcat(Vname,['gamma' num2str(ii)]);
end


elseif (nflag == 1) % the user supplied variable names
    Vname = 'Variable';
    for i=1:nvar
        Vname = strvcat(Vname,vnames(i+1,:));
    end
    % add spatial rho parameter name
    Vname = strvcat(Vname,'rho');
    % add gamma parameter names
    for ii=1:nmat
        Vname = strvcat(Vname,['gamma' num2str(ii)]);
    end
    
end % end of nflag issue

switch results.meth
    
    case{'sar_conv_panel_bma_g'}
        
        % extract posterior means
        bout = [results.beta
            results.rho
            results.gamma];
        
        % sige = results.sige;
        tmp1 = std(results.bdraw);
        tmp2 = std(results.pdraw);
        tmp3 = std(results.gdraw);
        
        bstd = [tmp1'
            tmp2
            tmp3'];
        
        tout = bout./bstd;
        
        % if strcmp(results.tflag,'tstat')
        %     tstat = bout./bstd;
        %     % find t-stat marginal probabilities
        %     tout = tdis_prb(tstat,results.nobs);
        %     results.tstat = bout./bstd; % trick for printing below
        % find plevels
        bout = plims(results.bdraw);
        pout = plims(results.pdraw);
        gout = plims(results.gdraw);
        bout = [bout
            pout
            gout];
        
        % rsqr = results.rsqr;
        
        % do effects estimates
        % =======================================================
        % a set of draws for the effects/impacts distribution
        direct_out =  plims(results.direct);
        indirect_out =  plims(results.indirect);
        total_out =  plims(results.total);
        
        
        fprintf(fid,'\n');
        fprintf(fid,'Bayesian Model Average of SAR convex panel W models \n');
        if (nflag == 1)
            fprintf(fid,'Dependent Variable  = %16s \n',vnames(1,:));
        end;
        fprintf(fid,'BMA Log-marginal    = %9.4f \n',results.lmarginal);
        % cstats2 = chainstats(results.drawpost);
        % fprintf(fid,'Log-marginal MCerror= %9.6f\n',cstats2(1,2));
        % fprintf(fid,'R-squared           = %9.4f \n',rsqr);
        % fprintf(fid,'Rbar-squared        = %9.4f \n',results.rbar);
        % fprintf(fid,'mean of sige draws  = %9.4f \n',mean(results.sdraw));
        % fprintf(fid,'posterior mode sige = %9.4f \n',results.sig_mode);
        fprintf(fid,'Nobs, T, Nvars      = %6d,%6d,%6d \n',results.nobs,results.ntime,results.nvar);
        fprintf(fid,'# weight matrices   = %6d \n', results.nmat);
        fprintf(fid,'ndraws,nomit        = %6d,%6d \n',results.ndraw,results.nomit);
        % fprintf(fid,'total time in secs  = %9.4f   \n',results.time);
        fprintf(fid,'total time          = %9.4f \n',results.time);
        % fprintf(fid,'time for sampling   = %9.4f \n',results.sampling_time);
        % fprintf(fid,'time for Taylor     = %9.4f \n',results.taylor_time);
        fprintf(fid,'thinning for draws  = %6d   \n',results.thin);
        fprintf(fid,'min and max rho     = %9.4f,%9.4f \n',results.rmin,results.rmax);
        
        fprintf(fid,'***************************************************************\n');
        nd = (results.ndraw-results.nomit)/results.thin;
        fprintf(fid,'      MCMC diagnostics ndraws = %d \n',nd);
        cstats=chainstats([results.bdraw results.pdraw results.gdraw]);
        in.cnames = strvcat('Mean','MC error','tau','Geweke');
        in.rnames = Vname;
        in.fmt = strvcat('%12.4f','12.8f','%10.6f','%10.6f');
        in.width = 10000;
        % pmode = [results.beta_mode results.rho_mode results.gamma_mode]';
        out = [cstats];
        mprint(out,in);
        
        fprintf(fid,'***************************************************************\n');
        fprintf(fid,'      Posterior Estimates \n');
        
        
        % % column labels for printing results
        cnames = strvcat('lower 0.01','lower 0.05','median','upper 0.95','upper 0.99');
        in.cnames = cnames;
        in.rnames = Vname;
        
        in.fmt ='%16.6f';
        in.fid = fid;
        mprint(bout,in);
        
        % print fixed effects, t-statistics and probabilities at request of the researcher
         if (results.fe==1)
             fprintf(fid,'fixed effects estimates are not available for BMA functions \n');
         end

        % now print x-effects estimates
                
        in.width = 2000;
        nvar = results.nvar;
        if results.cflag == 1
            vnameso = strvcat(Vname(3:nvar+1,:));
        elseif results.cflag == 0
            vnameso = strvcat(Vname(2:nvar+1,:));
        end
        

        % print effects estimates
        in.rnames = strvcat('Direct  ',vnameso);
        in.fmt = '%16.6f';
        in.fid = fid;
        
        cnames = strvcat('lower 0.01','lower 0.05','median','upper 0.95','upper 0.99');
        in.cnames = cnames;
        in.rnames = strvcat('Direct',vnameso);
        
        in.fmt ='%16.6f';
        in.fid = fid;
        mprint(direct_out,in);

        
        in.rnames = strvcat('Indirect',vnameso);
        mprint(indirect_out,in);
        
        in.rnames = strvcat('Total   ',vnameso);
        mprint(total_out,in);
        
        
    case{'sdm_conv_panel_bma_g'}
        
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
m = results.nmat;
    for ii=1:m
            Vname = strvcat(Vname,['gamma' num2str(ii)]);
    end;

        
        % extract posterior means
        bout = [results.beta
            results.rho
            results.gamma];
        
        % sige = results.sige;
        tmp1 = std(results.bdraw);
        tmp2 = std(results.pdraw);
        tmp3 = std(results.gdraw);
        
        bstd = [tmp1'
            tmp2
            tmp3'];
        
        tout = bout./bstd;
        
        % if strcmp(results.tflag,'tstat')
        %     tstat = bout./bstd;
        %     % find t-stat marginal probabilities
        %     tout = tdis_prb(tstat,results.nobs);
        %     results.tstat = bout./bstd; % trick for printing below
        % find plevels
        bout = plims(results.bdraw);
        pout = plims(results.pdraw);
        gout = plims(results.gdraw);
        bout = [bout
            pout
            gout];
        
        % rsqr = results.rsqr;
        
        % do effects estimates
        % =======================================================
        % a set of draws for the effects/impacts distribution
        direct_out =  plims(results.direct);
        indirect_out =  plims(results.indirect);
        total_out =  plims(results.total);
        
        
        fprintf(fid,'\n');
        fprintf(fid,'Bayesian Model Average of SDM convex panel W models \n');
        if (nflag == 1)
            fprintf(fid,'Dependent Variable  = %16s \n',vnames(1,:));
        end;
        fprintf(fid,'BMA Log-marginal    = %9.4f \n',results.lmarginal);
        % cstats2 = chainstats(results.drawpost);
        % fprintf(fid,'Log-marginal MCerror= %9.6f\n',cstats2(1,2));
        % fprintf(fid,'R-squared           = %9.4f \n',rsqr);
        % fprintf(fid,'Rbar-squared        = %9.4f \n',results.rbar);
        % fprintf(fid,'mean of sige draws  = %9.4f \n',mean(results.sdraw));
        % fprintf(fid,'posterior mode sige = %9.4f \n',results.sig_mode);
        fprintf(fid,'Nobs, T, Nvars      = %6d,%6d,%6d \n',results.nobs,results.ntime,results.nvar);
        fprintf(fid,'# weight matrices   = %6d \n', results.nmat);
        fprintf(fid,'ndraws,nomit        = %6d,%6d \n',results.ndraw,results.nomit);
        % fprintf(fid,'total time in secs  = %9.4f   \n',results.time);
        fprintf(fid,'total time          = %9.4f \n',results.time);
        % fprintf(fid,'time for sampling   = %9.4f \n',results.sampling_time);
        % fprintf(fid,'time for Taylor     = %9.4f \n',results.taylor_time);
        fprintf(fid,'thinning for draws  = %6d   \n',results.thin);
        fprintf(fid,'min and max rho     = %9.4f,%9.4f \n',results.rmin,results.rmax);
        
        fprintf(fid,'***************************************************************\n');
        nd = (results.ndraw-results.nomit)/results.thin;
        fprintf(fid,'      MCMC diagnostics ndraws = %d \n',nd);
        cstats=chainstats([results.bdraw results.pdraw results.gdraw]);
        in.cnames = strvcat('Mean','MC error','tau','Geweke');
        in.rnames = Vname;
        in.fmt = strvcat('%12.4f','12.8f','%10.6f','%10.6f');
        in.width = 10000;
        % pmode = [results.beta_mode results.rho_mode results.gamma_mode]';
        out = [cstats];
        mprint(out,in);
        
        fprintf(fid,'***************************************************************\n');
        fprintf(fid,'      Posterior Estimates \n');
        
        cnames = strvcat('lower 0.01','lower 0.05','median','upper 0.95','upper 0.99');
        in.cnames = cnames;
        in.rnames = Vname;
        
        in.fmt ='%16.6f';
        in.fid = fid;
        mprint(bout,in);
        

        % print fixed effects, t-statistics and probabilities at request of the researcher
         if (results.fe==1)
             fprintf(fid,'fixed effects not available for BMA models \n');
         end;

        in.width = 2000;
        nvar = results.nvar;
        
        % print effects estimates
        if results.cflag == 1
            vnameso = strvcat(Vname(3:nvar+1,:));
        elseif results.cflag == 0
            vnameso = strvcat(Vname(2:nvar+1,:));
        end
                
        in.rnames = strvcat('Direct  ',vnameso);
        
        in.fmt = '%16.6f';
        in.fid = fid;
        
        cnames = strvcat('lower 0.01','lower 0.05','median','upper 0.95','upper 0.99');
        in.cnames = cnames;
        in.rnames = strvcat('Direct',vnameso);
                
        in.fmt ='%16.6f';
        in.fid = fid;
        mprint(direct_out,in);
        
        in.rnames = strvcat('Indirect',vnameso);
        mprint(indirect_out,in);
        
        in.rnames = strvcat('Total   ',vnameso);
        mprint(total_out,in);
        

            case{'sdem_conv_panel_bma_g'}
        
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
m = results.nmat;
    for ii=1:m
            Vname = strvcat(Vname,['gamma' num2str(ii)]);
    end;

        
        % extract posterior means
        bout = [results.beta
            results.rho
            results.gamma];
        
        % sige = results.sige;
        tmp1 = std(results.bdraw);
        tmp2 = std(results.pdraw);
        tmp3 = std(results.gdraw);
        
        bstd = [tmp1'
            tmp2
            tmp3'];
        
        tout = bout./bstd;
        
        % if strcmp(results.tflag,'tstat')
        %     tstat = bout./bstd;
        %     % find t-stat marginal probabilities
        %     tout = tdis_prb(tstat,results.nobs);
        %     results.tstat = bout./bstd; % trick for printing below
        % find plevels
        bout = plims(results.bdraw);
        pout = plims(results.pdraw);
        gout = plims(results.gdraw);
        bout = [bout
            pout
            gout];
        
        % rsqr = results.rsqr;
        
        % do effects estimates
        % =======================================================
        % a set of draws for the effects/impacts distribution
        direct_out =  plims(results.direct);
        indirect_out =  plims(results.indirect);
        total_out =  plims(results.total);
        
        
        fprintf(fid,'\n');
        fprintf(fid,'Bayesian Model Average of SDEM convex panel W models \n');
        if (nflag == 1)
            fprintf(fid,'Dependent Variable  = %16s \n',vnames(1,:));
        end;
        fprintf(fid,'BMA Log-marginal    = %9.4f \n',results.lmarginal);
        % cstats2 = chainstats(results.drawpost);
        % fprintf(fid,'Log-marginal MCerror= %9.6f\n',cstats2(1,2));
        % fprintf(fid,'R-squared           = %9.4f \n',rsqr);
        % fprintf(fid,'Rbar-squared        = %9.4f \n',results.rbar);
        % fprintf(fid,'mean of sige draws  = %9.4f \n',mean(results.sdraw));
        % fprintf(fid,'posterior mode sige = %9.4f \n',results.sig_mode);
        fprintf(fid,'Nobs, T, Nvars      = %6d,%6d,%6d \n',results.nobs,results.ntime,results.nvar);
        fprintf(fid,'# weight matrices   = %6d \n', results.nmat);
        fprintf(fid,'ndraws,nomit        = %6d,%6d \n',results.ndraw,results.nomit);
        % fprintf(fid,'total time in secs  = %9.4f   \n',results.time);
        fprintf(fid,'total time          = %9.4f \n',results.time);
        % fprintf(fid,'time for sampling   = %9.4f \n',results.sampling_time);
        % fprintf(fid,'time for Taylor     = %9.4f \n',results.taylor_time);
        fprintf(fid,'thinning for draws  = %6d   \n',results.thin);
        fprintf(fid,'min and max rho     = %9.4f,%9.4f \n',results.rmin,results.rmax);
        
        fprintf(fid,'***************************************************************\n');
        nd = (results.ndraw-results.nomit)/results.thin;
        fprintf(fid,'      MCMC diagnostics ndraws = %d \n',nd);
        cstats=chainstats([results.bdraw results.pdraw results.gdraw]);
        in.cnames = strvcat('Mean','MC error','tau','Geweke');
        in.rnames = Vname;
        in.fmt = strvcat('%12.4f','12.8f','%10.6f','%10.6f');
        in.width = 10000;
        % pmode = [results.beta_mode results.rho_mode results.gamma_mode]';
        out = [cstats];
        mprint(out,in);
        
        fprintf(fid,'***************************************************************\n');
        fprintf(fid,'      Posterior Estimates \n');
        
        cnames = strvcat('lower 0.01','lower 0.05','median','upper 0.95','upper 0.99');
        in.cnames = cnames;
        in.rnames = Vname;
        
        in.fmt ='%16.6f';
        in.fid = fid;
        mprint(bout,in);
        

        % print fixed effects, t-statistics and probabilities at request of the researcher
         if (results.fe==1)
             fprintf(fid,'fixed effects not available for BMA models \n');
         end;

        in.width = 2000;
        nvar = results.nvar;
        % print effects estimates
        if results.cflag == 1
            vnameso = strvcat(Vname(3:nvar+1,:));
        elseif results.cflag == 0
            vnameso = strvcat(Vname(2:nvar+1,:));
        end

        in.rnames = strvcat('Direct  ',vnameso);
        in.fmt = '%16.6f';
        in.fid = fid;
        
        cnames = strvcat('lower 0.01','lower 0.05','median','upper 0.95','upper 0.99');
        in.cnames = cnames;
        in.rnames = strvcat('Direct',vnameso);
                
        in.fmt ='%16.6f';
        in.fid = fid;
        mprint(direct_out,in);
        
        in.rnames = strvcat('Indirect',vnameso);
        mprint(indirect_out,in);
        
        in.rnames = strvcat('Total   ',vnameso);
        mprint(total_out,in);
        

    otherwise
        
end


