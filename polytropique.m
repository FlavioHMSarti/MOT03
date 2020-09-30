function [n,i] = polytropique(vecP, vecV ,p,V, tol)

if (nargin < 5)
    tol = 1/100;
end

temp = vecP-p;
[minVal,indice] = min(temp(temp>0));
i = 0;
n = 0;
ratio = 10*tol;

    while (ratio > tol && i <= length(vecP))
        i = i+1;
        n_old = n;
        n = log10(vecP(i+indice)/vecP(indice))/log10(vecV(indice)/vecV(indice+i));
        ratio = abs(n-n_old)/abs(n_old);
    end

end