function [p5,V5,indices] = trouver_isotherme(vecteur_p, vecteur_V, p4pp, V4pp, tol)

% Fonction pour trouver le point apres la combustion isotherme (point 5)

% Entree
% vecteur_p -> vecteur des points de la courbe du moteur
% vecteur_V -> vecteur des points de la courbe du moteur (V/V_min)
% p4pp et V4pp -> pression et volume du point 4''
% tol -> tolerance 

% Sortie
% p5 et V5 -> pression (bar) et volume (V/V_min) du point 5
% indices -> indices representant la isotherme dans le vecteur_V

if (nargin < 5 && nargin >= 4)
    tol = 1/100; % 1% valeur par defaut
end

C = p4pp * V4pp; % isotherme, alors p*V = cte = C

indices = abs(1-vecteur_p.*vecteur_V./C) < tol;
p5 = vecteur_p(indices); V5 = vecteur_V(indices);
p5 = ceil(p5(end));
V5 = V5(end)*(1+tol);


end