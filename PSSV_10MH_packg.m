function [Sdr, hdr, sigu2dr, muhdr, psidr, count_psi] = PSSV_10MH_packg(y, ...
    W, R1, R2, a0, b0, mh0, iVmh0, tau, h, psi, sigu2, mu_h)
%--------------------------------------------------------------------------
% input:
% R1: burn-in's
% R2: held-in's
%--------------------------------------------------------------------------
%       # hyperparameters #
% a0  : shape parameter for the IG
% b0  : scale parameter for the IG
% muh0: mean of the prior for mu_h
% iVh0: precision of the prior for mu_h
% tau : uniform(-1/tau, 1/tau) (spectral radius of W)
%--------------------------------------------------------------------------
%        # initial values #
% h      : initial nT \times 1 vector of latent log-volatilities
% psi    : initial value of spatial auto-regressive parameters
% sigu2  : initial variance of u_i's
% mu_h   : initial n \times 1 vector of mean for h's
%--------------------------------------------------------------------------

N = length(y);
n = size(W, 1);
T = N/n;
pk = length(psi);
cc = 0.6/2.38^2;
R = R1 + R2;
count_psi = 0;

Sdr = zeros(N, R2);
hdr = zeros(N, R2);
sigu2dr = zeros(R2, 1);
muhdr = zeros(n, R2);
psidr = zeros(R2, pk);

% AMH initialization for psi
CovPsi = zeros(R+1, pk);
CovPsi(1, :) = psi';


% -- normal mixture components

pj = [.00609 .04775 .13057 .20674 .22715 .18842 .12047 .05591 .01575 .00115];
mj = [1.92677 1.34744 .73504 .02266 -.85173 -1.97278 -3.46788 -5.55246 -8.68384 -14.65000];
sigj2 = [.11265 .17788 .26768 .40611 .62699 .98583 1.57469 2.54498 4.16591 7.33342];
sigj = sqrt(sigj2);

ystar = log(y.^2);

for isim = 1:R   
    
    % -- sample mixture component indicators 
    % (from a 10-component discrete distribution)
    tmprd = rand([N 1]);
    q = repmat(pj, [N 1]).*normpdf(repmat(ystar, [1 10]),...
                                   repmat(h, [1 10]) + repmat(mj, [N 1]),...
                                   repmat(sigj, [N 1]));
    q = q./repmat(sum(q, 2), [1 10]);
    S = 10 - sum(repmat(tmprd, [1 10]) < cumsum(q, 2), 2) + 1;
    
    % -- sample h
    d_s = mj(S)';
    iSig_s = sparse(1:N, 1:N, 1./sigj2(S));
    krho = K_rho(psi, W, n);
    ki_rho = krho\speye(n);
    [tmp_i, tmp_m, tmp_l, tmp_o] = OMi(psi, W, ki_rho, n);
    tmp_a = blockdiag(tmp_i, kron(speye(T-2), tmp_m), tmp_l);
    tmp_b = kron(sparse(1:(T-1), 2:T, 1, T, T), tmp_o);
    Om = tmp_a + tmp_b + tmp_b';
    Omi = Om./sigu2;       % Omega^\star
    Bh = (iSig_s + Omi);
    Bmh = kron(ones([T 1]), mu_h);
    bhat = Bh\(Omi*Bmh + iSig_s*(ystar - d_s));
    h = bhat + chol(Bh, 'lower')'\randn([N 1]);  

    % -- sample mu_h
    s_Oi = ((tmp_i + (T-2).*tmp_m + tmp_l) + (T-1).*tmp_o + (T-1).*tmp_o')./sigu2; 
    Vmhat = (iVmh0 + s_Oi);
    Muhhat = Vmhat\(iVmh0*mh0 + sum(reshape(Omi*h, [n T]), 2));
    mu_h = Muhhat + chol(Vmhat,'lower')'\randn([n 1]);
    
    % -- sample sigu2
    hm = h - kron(ones([T 1]), mu_h);
    J_rho = kron(speye(T), (speye(n) - psi(1).*W)) +...
        kron(sparse(2:T,1:(T-1),1,T,T), -1.*(psi(2).*speye(n) + psi(3).*W));
    Pi_rho = blockdiag(ki_rho, speye(n*(T-1)));
    sigu2 = 1/gamrnd(a0 + N/2, 1/(b0 + ((hm'*(J_rho'*Pi_rho*J_rho)*hm)/2)));

    % -- sample psi via AMH step (rho in the paper)
    lpsi_i = llike_psi(psi, h, mu_h, sigu2, n, T, W);
    accept = 0;

    while (accept == 0) 
            if (isim <= 2*pk)
                psic = psi + sqrt(0.1^2/pk).*randn([pk 1]);
            else
                V_psi = cov(CovPsi(1:isim,:));
                % check pos-def ness
                [Vps, Dps] = eig(V_psi);
                dd = diag(Dps);
                for ij = 1:length(dd)
                    if dd(ij) <= 1e-10
                        dd(ij) = 1e-10;
                    end
                end
                V_psi = Vps*(diag(dd))*Vps';
                %
                uf = double(rand() < .95);
                psic = uf*(psi + chol(cc*(2.38^2/pk).*V_psi,'lower')*randn([pk 1])) +...
                    (1-uf)*(psi + sqrt(0.1^2/pk).*randn([pk 1]));
            end       
        % check the stability
        if (sum(abs(psic)) < (1/tau)) 
            accept = 1;
        end
    end
    %

    %
    lpsi_c = llike_psi(psic, h, mu_h, sigu2, n, T, W);
    if (lpsi_c - lpsi_i) > exp(1)
        pp = 1;
    else
        ratio = exp(lpsi_c - lpsi_i);
        pp = min(1, ratio);
    end
    if (rand() < pp)
        psi = psic;
        count_psi = count_psi + 1;
    end

    CovPsi(isim + 1, :) = psi';
    if (count_psi/isim < 0.4), cc = cc/1.05; end
    if (count_psi/isim > 0.6), cc = cc*1.05; end

    % -- store keepers
    if isim > R1
        Sdr(:, isim-R1) = S;
        hdr(:, isim-R1) = h;
        sigu2dr(isim-R1) = sigu2;
        muhdr(:, isim-R1) = mu_h;
        psidr(isim-R1, :) = psi';
    end

end
end
%


function Krho = K_rho(psi, W, n)
    rho1 = psi(1);
    rho2 = psi(2);
    rho3 = psi(3);
    Sin = (speye(n) - rho1.*W)\speye(n);
    An = rho2.*speye(n) + rho3.*W;
    tmp1 = An*Sin;
    tmp2 = tmp1*tmp1';
    for j = 2:15
        tmp1 = tmp1*(An*Sin);
        tmp2 = tmp2 + tmp1*tmp1';
    end
    Krho = tmp2;
end


function [tmp_i, tmp_m, tmp_l, tmp_o] = OMi(psi, W, ki_rho, n)
    rho1 = psi(1);
    rho2 = psi(2);
    rho3 = psi(3);
    Sn = (speye(n) - rho1.*W);
    An = rho2.*speye(n) + rho3.*W;
    tmp_i = (Sn'*ki_rho*Sn + An'*An);
    tmp_m = (Sn'*Sn + An'*An);
    tmp_l = (Sn'*Sn);
    tmp_o = -(An'*Sn);
end


function llk = llike_psi(PSI, h, mu_h, sigu2, n, T, W)
    Sn = (speye(n) - PSI(1).*W);
    krho = K_rho(PSI, W, n);
    ki_rho = krho\speye(n);
    [tmp_i, tmp_m, tmp_l, tmp_o] = OMi(PSI, W, ki_rho, n);
    tmp_a = blockdiag(tmp_i, kron(speye(T-2), tmp_m), tmp_l);
    tmp_b = kron(sparse(1:(T-1), 2:T, 1, T, T), tmp_o);
    Om = tmp_a + tmp_b + tmp_b';
    Omi = Om./sigu2;
    Bmh = kron(ones([T 1]), mu_h);
    llk = -(n*T)*log(sigu2)/2 + T*log(det(Sn)) - log(det(krho))/2 - ((h - Bmh)'*Omi*(h - Bmh))/2;  
end
























    
    
    
    
