**Ce répertoire contient 4 dossiers et 1 fichier .m:**


## fichier "codes"

Ce dossier contient 5 sous-fonctions de la fonction principale. Sont-ils:
  * code_matlab_modelisation.m --> script principal qui lance toutes les équations contenues dans le rapport et appelle les sous-fonctions.
  * interp_data.m       --> fonction d'interpolation des tables thermodynamiques pour le calcul des propriétés du gaz et de la température adiabatique de flamme.
  * polytropique.m      --> fonction qui trouve le coefficient de transformation polytropique (3'4 ou 56) à partir des valeurs du cycle d'un moteur spécifique.
  * proprietes_gaz.m    --> fonction qui trouve les propriétés de l'air, du gaz frais et du gaz brûlés à partir des tables thermodynamiques.
  * trouver.isotherme.m --> fonction qui trouve automatiquement le point de départ et d'arrivée de la combustion isotherme à partir d'un cycle PV d'un moteur.

## fichier "info_moteurs"
**Ce dossier contient deux autres dossiers (4a et 4b), qui contiennent les informations nécessaires pour chaque moteur. Choisissez l'un des moteurs, copiez tous les 
fichiers .csv et collez-les dans le dossier "input_variables".**

## fichier "input_variables"
Dossier vide qui recevra les valeurs du dossier "info_moteurs".

## fichier "tables"
Contient les tables thermodynamiques attachées à l'extrémité du poly

## archive "fonction.principale.m"
Fonction principale à lancer dans matlab. Il se chargera d'appeler toutes les autres fonctions / scripts seul et créera le dossier output_variables avec
les résultats trouvés.
