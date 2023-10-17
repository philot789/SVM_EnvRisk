
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
    

