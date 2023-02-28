# Mesure-de-Gibbs-et-algorithmes-MCMC

La mesure de Gibbs est une mesure de probabilité définie sur un ensemble fini de grand cardinal. Elle tient son nom de Josiah Willard Gibbs, père fondateur de la mécanique statistique dont l’objectif est d’expliquer les lois de la thermodynamique à l’aide de modèles statistiques pour des grands ensembles de particules. Cette mesure permet de modéliser des systèmes complexes incluant par exemple une structure de dépendance spatiale. Dans ce cas, la mesure permet de spécifier le comportement global du système à partir de caractéristiques locales. Il est ainsi naturel de retrouver cette mesure en physique statistique. Des exemples classiques sont les modèles d’Ising et de Potts qui décrivent l’évolution à l’échelle macroscopique d’un ensemble de particules à partir des interactions entre les paires de particules. Mais cette mesure se retrouve également dans d’autres domaines d’application comme l’écologie, le traitement d’image ou la chimie.

Bien qu’utiles d’un point de vue modélisation, les mesures de Gibbs présentent néanmoins une spécificité rendant difficile les approches probabilistes et statistiques. En effet, une mesure de Gibbs forme une famille exponentielle dont la fonction de partition est difficile si ce n’est impossible à calculer. Ainsi on ne dispose pas d’une forme analytique pour la densité ou la vraisemblance qui font partie intégrante de ces approches et il est nécessaire de mettre en place des solutions adaptées.

Dans ce mémoire de M1, nous nous intéressons à deux questions fondamentales. La première est celle de la simulation suivant une mesure de Gibbs. Tout d’abord, les méthodes de de Monte-Carlo par chaîne de Markov (MCMC) et le résultat fondateur énoncé par énoncé par Metropolis et al. (1953) apportent une réponse à cette problématique. L’idée générale est de construire une chaîne de Markov dont la loi de probabilité invariante est la mesure de Gibbs souhaitée. Ces méthodes ne permettent néanmoins pas de faire de la simulation exacte. Propp et Wilson (1996) ont proposé une solution basée sur des techniques de couplage permettant de lever cette limitation lorsque le noyau de transition de la chaîne vérifie une propriété de monotonie. La seconde question est celle de l’estimation du paramètre d’une mesure de Gibbs. Que ce soit dans un cadre fréquentiste ou bayésien, la difficulté intrinsèque vient de l’absence de forme explicite pour la vraisemblance du modèle. Cette difficulté peut néanmoins être contournée dès lors que l’on sait simuler suivant la mesure. L’estimateur du maximum de vraisemblance ne peut par exemple pas être déterminé directement mais le gradient peut néanmoins être estimé par des méthodes de Monte- Carlo et le problème d’optimisation se résout alors à l’aide de l’algorithme de Robbins et Monro (1951). Dans le paradigme bayésien, les estimateurs sont usuellement construits via la distribution empirique d’un échantillon produit par un algorithme MCMC. Chaque itération inclut en général une étape d’acceptation rejet nécessitant de calculer un rapport de vraisemblance. Cette étape étant impossible pour une mesure de Gibbs, Murray et al. (2012) puis Alquier et al. (2016) ont proposé de remplacer, là encore, ce rapport par un estimateur de Monte-Carlo.

Le mémoire s’organise en trois chapitres. Le Chapitre 1 introduit la mesure de Gibbs et les problématiques spécifiques nous intéressant. Dans le Chapitre 2, nous présentons les algorithmes MCMC et l’algorithme de couplage de Propp-Wilson permettant de simuler suivant les mesures considérées. Enfin, dans le Chapitre 3 nous présentons comment ces algorithmes de simulation permettent de contruire des estimateurs fréquentistes ou bayésiens du paramètre de la mesure.

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
