GOULANCOURT Rebecca M2BI
# ReadMe "Calcul de la surface accessible au solvant d'une protéine"

## Informations générales

Ce script est écrit en langage Python et calcule la surface accessible d'une protéine par le solvant, ici nous avons pris l'eau comme référence. 
Il est compilable via la version 3.9 de Python et nécessite l'installation de certaines librairies qui peuvent ne pas être initialement présentes dans votre environnement.
Les bibliothèques requises sont : **
- pandas 
- numpy
- scipy
- matplotlib
- sys
- mpl_toolkits
- collections**
***
Vous pouvez télécharger au préalable ces librairies via la commande "*conda install nom_du_module*" dans le shell, mais il faut avoir installé au minimum miniconda sur votre ordinateur.
Lien vers miniconda : https://docs.conda.io/en/latest/miniconda.html

```
$ conda install pandas
```
***
Certains modules ne peuvent pas être installés de cette manière comme mpl_toolkits car il requiert de connaître un channel précis.
Pour le télécharger il faut aller directement sur le site d'hébergement de modules anaconda-forge et taper la ligne de commande donnée par le lien : https://anaconda.org/conda-forge/mpld3. La commande utilisée dans notre cas est :
```
$ conda install -c conda-forge mpld3
```
***
Avant de lancer le code dans votre terminal il faut récupérer la fiche PDB de votre protéine d'intérêt et la stocker dans le même dossier que le script à lancer.
Lien vers la Protein Data Bank : https://www.rcsb.org/.

## Fonctionnement du code

Le code suit un cheminement particulier présenté ci-dessous sous forme simplifiée.
1. Parsing du fichier .pdb de la protéine d'intérêt.
2. Construction d'une matrice de distance atomique et détermination des atomes voisins.
3. On place sur chaque atome une sphère, si un point de la pshère est en contact avec un atome voisin alors on l'élimine. Si le point n'est pas en contact avec un voisin alors il est considéré comme exposé au solvant.
4. On regarde la distance séparant un point de la sphère exposé avec un atome voisin.
5. Si cette distance est inférieure à la taille d'une molécule d'eau alors celle-ci peut passer et on calcule une surface accessible au solvant.

## Lancement du script

Pour faire tourner le script il faut se placer dans le même dossier que le script et le fichier .pdb de la molécule à analyser. La commande à taper est : "*python3 script.py id_de_la_molécule.pdb*" tel que :
```
$ python3 script.py 1bfz.pdb
```
***
Le script se lancera puis vous demandera le nombre de points voulu sur la sphère. Plus le nombre de points sera élevé plus le calcul de la surface accessible au solvant sera précis.





