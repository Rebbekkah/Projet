"""Projet 'Calcul de la surface accessible au solvant d'une protéine'

------------------------------------------------------------------
Rebecca GOULANCOURT
M2 BIOLOGIE - INFORMATIQUE
Numéro étudiant : 71602003
------------------------------------------------------------------

Ce code n'est pas aboutit mais peut être retravaillé dans l'optique de calculer la surface accessible au solvant d'une protéine.

"""

# Import des bibliothèques nécessaires au fonctionnement du code
import pandas as pd
import numpy as np
import sys
import mpl_toolkits.mplot3d
import matplotlib.pyplot as pp
from collections import defaultdict
from scipy.spatial import distance_matrix


def parsing(filename) :
	""" Récupération du fichier .pdb de la molécule biologique et parsing pour avoir les coordonnées atomique
	
	Parameters
	----------
	filename : str
		Nom du fichier .pdb à analyser

	Returns
	-------
	coord (DataFrame)
		Dataframe contenant entre autre les positions {x, y, z} de chaque atome, son nom ainsi que le nom et la position des résidus

	"""

	with open(filename, "r") as pdb_file :
		atom_name = []
		res_name = []
		res_num = []
		x = []
		y = []
		z = []
		for line in pdb_file :
			if line.startswith("ATOM") :
				atom_name.append(line[76:78].strip()) # on isole le symbole de l'atome 
				res_name.append(line[17:20].strip()) # on isole les colonnes correspondant au nom des résidus du fichier pdb
				res_num.append(int(line[22:26])) # on isole le numéro du résidu
				x.append(float(line[30:38])) # coordonnées x
				y.append(float(line[38:46])) # coordonnées y
				z.append(float(line[46:54])) # coordonnées z
		print(atom_name, res_name, res_num)

		coord = pd.DataFrame({'atom_name' : atom_name, 'res_name' : res_name, 'res_num' : res_num, \
			'x' : x, 'y' : y, 'z' : z}) #on fabrique une DataFrame
		print(coord)
		coord.to_csv("coord_"+filename+".txt", sep = "\t")
		return coord



def distance(coord) :
	"""Calcul de la matrice de distances atomique

	Paremeters
	----------
	coord : DataFrame
		Matrice contenant la position (coordonnées) de chaque atome

	Returns
	-------
	mat_dist  : matrix
		Matrice contenant la distance entre chaques atomes pris deux à deux

	"""
	mat_dist = pd.DataFrame(distance_matrix(coord.iloc[:, 4:], coord.iloc[:, 4:])) # Calcul de la différence de position entre les atomes = distance
	print(mat_dist)
	return mat_dist



def neighbor(mat_dist) :
	"""Recherche des atomes voisins

	Parameters 
	----------
	mat_dist : matrix
		Matrice des distances atomiques

	Returns
	-------
	mat_dist : matrix
		Matrice des voisins (1 = voisins l'un de l'autre et 0 = trop éloignés)

	"""

	radius_WdW = {'H' : 1.20, 'C' : 1.70, 'N' : 1.55, 'O' : 1.52, \
		'F' :1.47, 'P' : 1.80, 'S' : 1.80, 'Cl' : 1.75, 'Cu' : 1.4} # Rayon de Van der Waals de chaque atome
	fold = radius_WdW["P"]*2 + 1.4
	print(fold)
	mat_dist[mat_dist <= fold] = 1 # Si la distance est inférieur à 2*radius_WdW + 1.4 alors les atomes sont voisins
	mat_dist[mat_dist > fold] = 0 # Sinon ils ne le sont pas
	print(mat_dist)

	return(mat_dist)



def neigh_coord(mat_dist, coord) :
	"""Récupération pour chaque atome de ses voisins et leurs coordonnées
	
	Parameters
	----------
	mat_dist : matrix
		Matrice de 0 et 1 représentant les atomes et leurs voisins

	coord : DataFrame
		Matrice contenant la position (coordonnées) de chaque atome
	
	Returns
	-------
	dic_neigh : dictionnaire
		dicitionnaire contenant pour chaque atome une liste de ses voisins

	list_coor : array
		Liste des voisins et leurs coordonnées

	"""
	index_list = mat_dist.index.values.tolist() # on transforme les lignes de la matrice en liste
	row_list = mat_dist.values.tolist() # on transforme les éléments de la matrice en liste
	# Initialisation des dictionnaires vides
	dic_neigh = {} 
	dic_coor_neigh = {}

	for index in range(len(index_list)) : # on parcourt les lignes de mat_dist
		list_neigh = []
		for row in range(len(row_list[index])-1) : # on parcourt les valeurs de la ligne de mat_dist
			if row_list[index][row] == 1 : # Si les atomes sont voisins alors on récupère le voisin
				list_neigh.append(row)
		dic_neigh[index] = list_neigh

	for atom in dic_neigh : # on récupère les atomes voisins puis y associer leurs coordonnées
		list_coor_neigh = []
		list_coor_neigh.append(atom)
		for value in dic_neigh[atom] :
			list_coor_neigh.append(value)
		list_coor = []
		for row in coord.iloc[list_coor_neigh[1:], ].iterrows() :
			list_coor.append([row[1][3], row[1][4], row[1][5]]) # récupération des coordonnées et stockage
	print(list_coor)

	return dic_neigh, list_coor



def Sphere(number_points, coord) :
	"""Définition de la sphère à placer autour de chaque atome, de ses points et leurs coordonnées

	Parameters
	----------
	number_points : int
		Nombre de points désiré à positionner sur la sphère

	coord : DataFrame
		Matrice contenant la position (coordonnées) de chaque atome

	"""

	radius_WdW = {'H' : 1.20, 'C' : 1.70, 'N' : 1.55, 'O' : 1.52, \
		'F' :1.47, 'P' : 1.80, 'S' : 1.80, 'Cl' : 1.75, 'Cu' : 1.4} # dictionnaire contenant les rayons de Van der Waals de chaque atome

	list_sphere = [] # liste des points de la sphère

	# initialisation des variables pour positionner une sphère de number_points équidistants
	goldenRatio = (1 + 5**0.5)/2
	index = np.arange(0, number_points) 
	phi = np.arccos(1 - 2*(index+0.5)/number_points)
	theta = np.pi * 2 * index / goldenRatio
	
	for i, atom in coord.loc[:, 'atom_name'].iteritems() :
		array_sphere = np.dot(np.ones(number_points, coord.iloc[:,3:]), np.dot(coord.iloc[:,3:], number_points)).shape() # ajustement de la dimension des données
		# Centrage des sphères sur les atomes
		x_point = (np.cos(theta) * np.sin(phi))*radius_WdW[atom] + coord.iloc[:,3]
		y_point = (np.sin(theta) * np.sin(phi))*radius_WdW[atom] + coord.iloc[:,4]
		z_point = np.cos(phi)*radius_WdW[atom] + coord.iloc[:,5]
		array_sphere = [x_point, y_point, z_point] # coordonnées de chaque points
		list_sphere.append(array_sphere) # contenant les points de la sphère et ses coordonnées

	pp.figure().add_subplot(111, projection = '3d').scatter(array_sphere)
	pp.show()




if __name__ == "__main__" :
	filename = sys.argv[1] # récupération du nom de fichier passé en argument

	# Appel des fonctions
	pars = parsing(filename)
	dist_matrix = distance(pars)
	neigh = neighbor(dist_matrix)
	dico_neigh, coor_neigh = neigh_coord(neigh, pars)
	number_points = int(input("Rentrez le nombre de points souhaité pour la sphère : "))
	Sphere(number_points, pars)


