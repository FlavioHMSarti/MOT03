function valeur_interpolee = interp_data(M, val ,col1, col2)

% M -> matrice avec la table thermodynamiques d'un gaz
% val -> la valeur de temperature pour l'interpolation
% col1 -> l'indice de la colonne de temperature
% col2 -> l'indice de la colonne ou on veut la valeur interpolee

colMin = M(:,col1)-val;
pos = length(colMin(colMin<=0));

if(pos == 1 || pos == length(colMin))
    valeur_interpolee = M(pos,col2);
else
    if (abs(M(pos,col1)-M(pos+1,col1))>abs(M(pos,col1)-M(pos-1,col1)))
        pos2 = pos-1;
    else
        pos2 = pos+1;
    end
    valeur_interpolee = M(pos,col2)+(M(pos2,col2)-M(pos,col2))/(M(pos2,col1)-M(pos,col1))*(val-M(pos,col1));
end
end
