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


# Définition des variables utilisées
radius_WdW = {'H' : 1.20, 'C' : 1.70, 'N' : 1.55, 'O' : 1.52, \
		'F' :1.47, 'P' : 1.80, 'S' : 1.80, 'Cl' : 1.75, 'Cu' : 1.4} # Rayon de Van der Waals de chaque atome
radius_solvent = 1.4


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

	fold = radius_WdW["P"]*2 + radius_solvent
	print(fold)
	mat_dist[mat_dist <= fold] = 1 # Si la distance est inférieur à 2*radius_WdW + radius du solvant alors les atomes sont voisins
	mat_dist[mat_dist > fold] = 0 # Sinon ils ne le sont pas
	print(mat_dist)
	#print(type(mat_dist))
	return(mat_dist)



#def neigh_coord(mat_dist, coord) :
def neigh_coord(mat_dist) :
	"""Récupération pour chaque atome de ses voisins et leurs coordonnées
	
	Parameters
	----------
	mat_dist : DataFrame
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
	
	dico_neigh = {}
	for atom in mat_dist :
		match = mat_dist[atom][mat_dist[atom] == 1]
		dico_neigh[atom] = match.index.tolist()

	#print(dico_neigh[0])

	return dico_neigh


def Fibonacci(number_points) :
	"""Visualisation de la sphère à placer autour de chaque atomes

	Parameters
	----------
	number_points : int
		Nombre de points désiré à positionner sur la sphère

	"""

	index = np.arange(0, number_points, dtype = float) + 0.5
	phi = np.arccos(1 - 2*index/number_points)
	theta = np.pi * (1 + 5**0.5) * index
	x, y, z = np.cos(theta) * np.sin(phi), np.sin(theta) * np.sin(phi), np.cos(phi)
	pp.figure().add_subplot(111, projection = '3d').scatter(x, y, z)
	pp.show()


def Sphere(number_points, coord, dico_neigh) :
	"""Définition de la sphère à placer autour de chaque atome, de ses points et leurs coordonnées

	Parameters
	----------
	number_points : int
		Nombre de points désiré à positionner sur la sphère

	coord : DataFrame
		Matrice contenant la position (coordonnées) de chaque atome

	"""

	list_sphere = []
	coord_sphere = {}
	coord_sphere['x'] = []
	coord_sphere['y'] = []
	coord_sphere['z'] = []

	goldenRatio = (1 + 5**0.5)/2
	index = np.arange(0, number_points) 
	phi = np.arccos(1 - 2*(index+0.5)/number_points)
	theta = np.pi * 2 * index / goldenRatio

	for atom in dico_neigh :
		x_point = (np.cos(theta) * np.sin(phi))*radius_WdW[coord.iloc[atom,0]] + coord.iloc[atom,3]
		y_point = (np.sin(theta) * np.sin(phi))*radius_WdW[coord.iloc[atom,0]] + coord.iloc[atom,4]
		z_point = np.cos(phi)*radius_WdW[coord.iloc[atom,0]] + coord.iloc[atom,5]
		coord_sphere['x'].append(x_point)
		coord_sphere['y'].append(y_point)
		coord_sphere['z'].append(z_point)

	print(coord_sphere)

	return coord_sphere



if __name__ == "__main__" :
	filename = sys.argv[1] # récupération du nom de fichier passé en argument
	number_points = int(sys.argv[2]) # nombre de points sur la sphère

	# Appel des fonctions
	pars = parsing(filename)
	dist_matrix = distance(pars)
	neigh = neighbor(dist_matrix)

	dico_neighbor = neigh_coord(neigh)
	Fibonacci(number_points)
	points_coord = Sphere(number_points, pars, dico_neighbor)


