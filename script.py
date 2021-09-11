import pandas as pd
import numpy as np
import sys
#from scipy.spatial.distance import pdist
from scipy.spatial import distance_matrix
from numpy import pi, cos, sin, arccos, arange
import mpl_toolkits.mplot3d
import matplotlib.pyplot as pp
from collections import defaultdict

def parsing(filename) : #récupération des coordonnées
	with open(filename, "r") as pdb_file :
		#res_count = 0 
		atom_name = []
		res_name = []
		res_num = []
		x = []
		y = []
		z = []
		for line in pdb_file :
			if line.startswith("ATOM") :
				#peut être pas nécessaire pour l'atome
				atom_name.append(line[76:78].strip()) # on isole le symbole de l'atome 
				res_name.append(line[17:20].strip()) # on isole les colonnes correspondant au nom des résidus du fichier pdb
				res_num.append(int(line[22:26])) # on isole le numéro du résidu
				#res_count .append(1
				x.append(float(line[30:38])) # coordonnées x
				y.append(float(line[38:46])) # coordonnées y
				z.append(float(line[46:54])) # coordonnées z
				#print(coord)
				'''
				with open("coord_"+filename+".txt", "w") as filout :
				   coord.to_csv('')
					'''
		print(atom_name, res_name, res_num)

		coord = pd.DataFrame({'atom_name' : atom_name, 'res_name' : res_name, 'res_num' : res_num, \
			'x' : x, 'y' : y, 'z' : z}) #on fabrique une DataFrame
		print(coord)
		# print(type(coord))
		coord.to_csv("coord_"+filename+".txt", sep = "\t")
		return coord


def distance(coord) :
	# coord_val = coord[['x', 'y', 'z']]
	mat_dist = pd.DataFrame(distance_matrix(coord.iloc[:, 4:], coord.iloc[:, 4:]))
	print(mat_dist)
	return mat_dist

def neighbor(mat_dist) :
	radius_WdW = {'H' : 1.20, 'C' : 1.70, 'N' : 1.55, 'O' : 1.52, \
		'F' :1.47, 'P' : 1.80, 'S' : 1.80, 'Cl' : 1.75, 'Cu' : 1.4}
	fold = radius_WdW["P"]*2 + 1.4
	print(fold)
	mat_dist[mat_dist <= fold] = 1
	mat_dist[mat_dist > fold] = 0
	# print(type(mat_dist))
	# matrix_distance[matrix_distance < fold] = 1
	# matrix_distance[matrix_distance > fold] = 0
	# matrix_neighbor = matrix_distance[matrix_distance == 1]
	print(mat_dist)

	# print(len(mat_dist))
	'''
	for i in range(len(mat_dist)) :
		index_list_+i = []
		for row in mat_dist.iterrows() :
			if (mat_dist[mat_dist == 1]) :
				index_list_+i.append(i)
	'''

	'''
	dic_list = defaultdict(list)
	for i in range(len(mat_dist)) :
		all_lists = "index_list_{}".format(i)
		#dic_list[all_lists] = {}
		for row in mat_dist.iterrows() :
			if (mat_dist[mat_dist == 1]) :
				dic_list[all_lists].append(i)
	'''

	# print(all_lists)
	return(mat_dist)
	# print(matrix_distance)
	# return matrix_neighbor

	# mat_dist = pdist(coord_val, 'euclidean')
	# mat_dist = pd.DataFrame(coord.pairwise(coord.to_numpy()))
	# mat_dist = dist.pairwise(coord)
	# print(mat_dist)
	# coord_val = coord[['x', 'y', 'z']]
	# coord_bis = coord_val.copy()
	# mat_dist = np.linalg.norm(coord_val - coord_bis)
	# print(mat_dist)
	# with open("mat_dist"+filename+".txt", "w") as filout :
	#     for value in mat_dist :
	#         filout.write(str(value))
	# print(coord_bis)
	
	# for column in coord :
	#     for value in column :
			
	
	# for column in coord :
	#     mat.coord = column[3:5]
	#     print(mat.coord)
	#     if column.startswith("x") || column.startswith("y") || column.startswith("z") :

'''
#distance(pars)
			if atom_name == "CA": # squelette carbone
				res_count += 1
				x = float(line[30:38]) # coordonnées x
				y = float(line[38:46]) # coordonnées y
				z = float(line[46:54]) # coordonnées z [13]: from scipy.spatial.distance import squareform, pdist
				print(res_name, res_num, x, y, z)
'''

def neigh_coord(neigh) :
	index_list = mat_dist.index.values.tolist()
	row_list = mat_dist.values.tolist()
	dic_neigh = {}
	for index in range(len(index_list)) :
		list_neigh = []
		for row in range(len(row_list[index])-1) :
			if row_list[index][row] == 1 :
				list.neigh.append(row)
		dic_neigh[index] = list_neigh
	print(dic_list, dic_list.items())

	return dic_neigh

'''
	dic_list = defaultdict(list)
	for i in range(len(mat_dist)) :
		all_lists = "index_list_{}".format(i)
		#dic_list[all_lists] = {}
		for row in mat_dist.iterrows() :
			if (mat_dist[mat_dist == 1]) :
				dic_list[all_lists].append(i)	
	print(dic_list, dic_list.items())
'''
def Sphere(number_points) :
	number_points = int(input("Rentrez le nombre de points souhaité pour la sphère"))
	index = arange(0, number_points, dtype = float) + 0.5
	phi = arccos(1 - 2*index/number_points)
	theta = pi * (1 + 5**0.5) * index
	x, y, z = cos(theta) * sin(phi), sin(theta) * sin(phi), cos(phi)
	pp.figure().add_subplot(111, projection = '3d').scatter(x, y, z)
	pp.show()



filename = sys.argv[1]
pars = parsing(filename)
#distance(pars)
dist_matrix = distance(pars)
neigh = neighbor(dist_matrix)
neigh_2 = neigh_coord(neigh)
Sphere(30)
