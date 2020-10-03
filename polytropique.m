function n = polytropique(vecP, vecV ,p)

% fonction pour trouver le coef n d'une transformation 
% polytropique du type p*V^n = cste

% INPUTS
% vecP = vecteur de pression (bar) du cycle (partie bas du cycle si
% compression ou partie haute si detente)
% vecV = vecteur volume (V/Vmin) du cycle
% p = pression de reference (p3p si compression ou p5 si detente)
% tol 

% OUTPUT
% n = valeur du coef

[~,indice] = min(abs(vecP-p));

i = 0;
n = zeros(1,length(vecP)-indice);

    while (i <= length(vecP)-indice-1) %ratio > tol && 
        i = i+1;
        pol_new = polyfit(log10(vecV(indice:indice+i)),log10(vecP(indice:indice+i)),1);
        n(i) = pol_new(1);
    end

n=sort(abs(n))';
n(n<1.1 | n > 1.5) = []; % pas reel
n = median(n);

end