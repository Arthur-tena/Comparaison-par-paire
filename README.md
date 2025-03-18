# Comparaison par paire
**Auteur** : Arthur TENA

**Contact** : arthurtena3@gmail.com 

## Contexte d'étude 

Cette étude est mon sujet de stage intitulé **Comparaison par paire généralisée pour plusieurs critères de jugement classés par ordre de priorité : win ratio vs. Net Benefit of treatment** proposé par Maïlis Amico chercheuse à l'Institut Desbrest d'Epidémiologie et de Santé Publique (IDESP).

## Explication de la comparaison par paire
La **comparaison par paire** est une technique visant à comparer 2 tableaux de données selon des critères de jugements. Le but de cette comparaison est de savoir si un nouveau traitement est meilleur que le traitement de contrôle (placebo ou ancien traitement), s'il est moins bon ou s'il ne change pas les résultats.

Pour cela , 3 méthodes seront mise en place, la méthode **GPC** (*generalize pairwise comparison*), les **WR** (*Win ratio*) et les **WO** (*Win odds*). 
Ces 3 techniques ont un point commun, c'est l'affectation de score permettant de calculer la statistique qui décidera si le nouveau traitement est meilleur, moins bon ou neutre vis à vis du traitement de contrôle. Et la différence entre ces méthodes vient du calcul de cette statistique.

Un document détaillant ces techniques est fournit dans sa première version, Stage_IDESP.pdf.

## Simulation
Lors d'essai clinique, il est important d'effectuer des simulations afin de vérifier le déroulement des techniques car nous pouvons contrôler les données avec lesquelles nous travaillons. Il est important alors d'écarter des scénarios pouvant se produire avec des vraies données. J'ai choisi d'écarter 3 scénarios, *le nouveau traitement est uniformément meilleur que le traitement de contrôle*, *le nouveau traitement est uniformément moins bon que le traitement de contrôle* et *il n'y a pas de différence claire entre les 2 traitements*.

## Contenance du repository

Ce repo contient :
 - simulation.Rmd : un R markdown où j'effectue mes simulations et où sont mise mes fonctions
 - simulation.pdf : le document résultant de mon R markdown
 - Stage_IDESP.pdf : le rapport de mon stage 
