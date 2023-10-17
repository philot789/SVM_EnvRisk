% PURPOSE: An example of using sar on a large data set   
%          Gibbs sampling spatial autoregressive model                         
%---------------------------------------------------
% USAGE: sar_d3 (see sar_d for a small data set)
%---------------------------------------------------

clear all;
n = 30000;

W = xy2cont(rand(n,1),rand(n,1));

beta = ones(2,1);
x = randn(n,2);
rho = 0.6;
evec = randn(n,1);

y = (speye(n) - rho*W)\(x*beta + evec);


info.lflag = 0;           % use full lndet calculation
result = sar(y,x,W,info); % maximum likelihood estimates
prt(result);

                      % use default MC approximation for lndet calculation
result2 = sar(y,x,W); % maximum likelihood estimates
prt(result2);
