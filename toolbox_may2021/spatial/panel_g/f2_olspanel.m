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