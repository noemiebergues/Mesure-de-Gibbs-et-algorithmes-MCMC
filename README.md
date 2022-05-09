# Mesure-de-Gibbs-et-algorithmes-MCMC
Les fichiers disponibles sur ce github permettent d'implémenter: 
  - L'algorithme du rejet sur le modèle du composé chimique, 
  - Les algorithmes de Metropolis-Hastings et Gibbs sur les modèles du composé chimique, d'Ising et de Potts,
  - L'algorithme du recuit simulé sur un jeu du Sudoku,
  - L'algorithme de Propp-Wilson sur les modèles du composé chimique et d'Ising, 
  - Les algorithmes de descente de gradient et d'échange bruité sur le modèle d'Ising. 

Chaque dossier correspond à un modèle et comprend : 
 - Un fichier "NOMDUMODELE_modele.R" qui contient les fonctions qui caractérisent le modèle et qui sont utilisées par les algorithmes étudiés,
 - Des fichiers "NOMDUMODELE_NOMDELALGORITHME.R" qui contiennent les algorithmes implémentés pour ce modèle,
 - Un fichier "NOMDUMODELE_exemple.R" qui contient un exemple d'utilisation des algorithmes implémentés pour ce modèle. 

Pour faire fonctionner un des algorithmes sur un modèle avec des paramètres spécifiques, procédez de la manière suivante : 
  - Ouvrir le dossier correspondant au modèle choisi, 
  - Exécuter le code présent dans le fichier "NOMDUMODELE_modele.R",
  - Exécuter le code présent dans le fichier correspondant à l'algorithme souhaité, 
  - Appeler la fonction-algorithme avec les paramètres choisis. 
