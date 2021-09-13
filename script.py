import pandas as pd
import numpy as np
import sys
#from scipy.spatial.distance import pdist
from scipy.spatial import distance_matrix
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
	# mat_dist.to_csv("mat_dist"+filename+".txt", sep = "\t")
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

def neigh_coord(mat_dist, coord) :
	index_list = mat_dist.index.values.tolist()
	row_list = mat_dist.values.tolist()
	dic_neigh = {}

	dic_coor_neigh = {}

	for index in range(len(index_list)) :
		list_neigh = []
		for row in range(len(row_list[index])-1) :
			if row_list[index][row] == 1 :
				list_neigh.append(row)
		dic_neigh[index] = list_neigh
	# print(dic_neigh, dic_neigh.items())
	# print(dic_neigh)
	'''
	with open("dic_neigh_"+filename+".txt", "w") as dic_file :
		for item in dic_neigh :
			dic_file.write(str(item))
	'''

	for atom in dic_neigh :
		list_coor_neigh = []
		list_coor_neigh.append(atom)
		for value in dic_neigh[atom] :
			list_coor_neigh.append(value)
		list_coor = []
		for row in coord.iloc[list_coor_neigh[1:], ].iterrows() :
			list_coor.append([row[1][3], row[1][4], row[1][5]])
	print(list_coor)

	'''
	for atom in dic_neigh.keys() :
		for index in dic_neigh.values() :
			list_neigh_coor = []
			for row in coord.iloc[:, 4:] :
				list_neigh_coor.append(row) 
			dic_coor_neigh[index] = list_neigh_coor
	print(dic_coor_neigh)
	'''

	return dic_neigh, list_coor

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
def Sphere(number_points, coord) :

	radius_WdW = {'H' : 1.20, 'C' : 1.70, 'N' : 1.55, 'O' : 1.52, \
		'F' :1.47, 'P' : 1.80, 'S' : 1.80, 'Cl' : 1.75, 'Cu' : 1.4}

	list_sphere = []
	goldenRatio = (1 + 5**0.5)/2
	#index = np.arange(0, number_points, dtype = float)
	index = np.arange(0, number_points) 
	phi = np.arccos(1 - 2*(index+0.5)/number_points)
	theta = np.pi * 2 * index / goldenRatio
	print(coord)
	for i, atom in coord.loc[:, 'atom_name'].iteritems() :
		array_sphere = np.zeros((number_points, 3))
		x_point = (np.cos(theta) * np.sin(phi))*radius_WdW[atom] + coord.loc[:, 'x']
		print(x_point)
		y_point = (np.sin(theta) * np.sin(phi))*radius_WdW[atom] + coord.loc[:, 'y']
		print(y_point)
		z_point = np.cos(phi)*radius_WdW[atom] + coord.loc[:, 'z']
		print(z_point)
		array_sphere = [x_point, y_point, z_point]
		print(array_sphere)
		list_sphere.append(array_sphere)
		print(list_sphere)
		#pp.figure().add_subplot(111, projection = '3d').scatter(x_point, y_point, z_point)
		#pp.show()



	'''
	for point in 
		sphere_point = []
		for point in range(len(number_points)) :
	
	x_point = []
	y_point = []
	z_point = []
	points = []

	radius_WdW = {'H' : 1.20, 'C' : 1.70, 'N' : 1.55, 'O' : 1.52, \
		'F' :1.47, 'P' : 1.80, 'S' : 1.80, 'Cl' : 1.75, 'Cu' : 1.4}

	for atom in coord.loc[:, 'atom_name'] :
		if atom == 'H' :
			radius = 1.20
		elif atom == 'C' :
			radius = 1.70
		elif atom == 'N' :
			radius = 1.70
		elif atom == 'O' :
			radius = 1.70
		elif atom == 'F' :
			radius = 1.70
		elif atom == 'P' :
			radius = 1.70
		elif atom == 'S' :
			radius = 1.70
		elif atom == 'Cl' :
			radius = 1.70
		elif atom == 'Cu' :
			radius = 1.70

		for point in range(number_points) :
			points.append(int(point))

			for coordinate in points :
				x_point.append((list(list_coor + radius))
				y_point.append((list(list_coor + radius))
				z_point.append((list(list_coor + radius))

				
				x_point.append(list_coor.values([0]) + radius)
				y_point.append(list_coor.values([1]) + radius)
				z_point.append(list_coor.values([2]) + radius)
				'''
'''		
				x_point.append(list_coor.values()[0] + radius)
				y_point.append(list_coor.values()[1] + radius)
				z_point.append(list_coor.values()[2] + radius)
				#x_points.append(point[coordinate]*)
				
		sphere = pd.DataFrame({'xSphere' : x_point, 'ySphere' : y_point, 'zSphere' : z_point})
	print(sphere)
	#print(points)
'''



if __name__ == "__main__" :
	filename = sys.argv[1]
	pars = parsing(filename)
	dist_matrix = distance(pars)
	neigh = neighbor(dist_matrix)
	dico_neigh, coor_neigh = neigh_coord(neigh, pars)
	number_points = int(input("Rentrez le nombre de points souhaité pour la sphère : "))
	Sphere(number_points, pars)

