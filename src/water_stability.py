# Author: Anthony Kai Kwang Ma
# Email: anthonyma27@gmail.com, akma327@stanford.edu
# Date: 07/13/17
# water_stability.py

import vmd, molecule, sys
from vmd import *
from atomsel import *
import datetime 
import glob 
import re
import numpy as np
import time 

USAGE_STR = """

# Purpose
# Aligns every frame of trajectory to topology. For each of these
# frames determine whether there are water molecules within a 
# specified distance from a particular crystallographic water position. 

# Usage 
# python water_stability.py <CRYS_REF> <TOP> <TRAJ> <OUTPUT> <DISTANCE_THRESH> <WATER_INDEX> <RESIDUE_INDICES>

# Arguments
# <CRYS_REF> Absolute path to topology of crystal structure pre-aligned to simulation topology
# <TOP> Absolute path to simulation topology 
# <TRAJ> Absolute path to trajectory file
# <OUTPUT> Output log 
# <DISTANCE_THRESH> Distance threshold to search for neighboring waters near a certain region. In nanometers
# <WATER_INDEX> The index of water molecule that serves as crystallographic reference 
# <RESIDUE_INDICES> Comma delimited indices of residues that surround the reference crystallographic water molecule 


# Example

CRYS_REF="/Users/anthony/Desktop/sherlock/MIF-waters/data/simulation/3DJH_wb/3djh_crys.pdb"
TOP="/Users/anthony/Desktop/sherlock/MIF-waters/data/simulation/3DJH_wb/3DJH_wb.pdb"
TRAJ="/Users/anthony/Desktop/sherlock/MIF-waters/data/simulation/3DJH_wb/3DJH_run.1.dcd"
OUTPUT="/Users/anthony/Desktop/sherlock/MIF-waters/data/water_stability/water_stability_1Ang.txt"
DISTANCE_THRESH=1.0
cd /Users/anthony/Desktop/dror/local-dev/github/StillWaters/src
python water_stability.py $CRYS_REF $TOP $TRAJ $OUTPUT $DISTANCE_THRESH

"""

K_MIN_ARG = 6

def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    return [ atoi(c) for c in re.split('(\d+)', text) ]

def euc_dist(pt1, pt2):
	a, b = np.array(pt1), np.array(pt2)
	return np.linalg.norm(a-b)



def map_water_to_coordinate(molid, selection, option="resid"):
	"""
		Takes in a molid and water selection. Returns a mapping from water resid 
		to coordinates. 
		option = "resid" to get resid key, "index" to get index as keys
	"""

	waters = evaltcl("set waters [atomselect " + str(molid) + " \" " + selection + " \" frame 0]")
	resid = map(str, (evaltcl("$waters get resid")).split(" "))
	indices = map(str, (evaltcl("$waters get index")).split(" "))
	xcoord = map(float, evaltcl("$waters get x").split(" "))
	ycoord = map(float, evaltcl("$waters get y").split(" "))
	zcoord = map(float, evaltcl("$waters get z").split(" "))


	### Create map from crystal water resid to coordinates 
	water_resid_to_coord = {}
	if(option == "resid"):
		for idx, resi in enumerate(resid):
			water_resid_to_coord[resi] = np.array([xcoord[idx], ycoord[idx], zcoord[idx]])
	elif(option == "index"):
		for idx, index in enumerate(indices):
			water_resid_to_coord[index] = np.array([xcoord[idx], ycoord[idx], zcoord[idx]])

	return water_resid_to_coord


def map_crystal_to_sim(CRYS_PDB, SIM_PDB):
	"""
		Map resid of crystallographic water to resid of corresponding simulation water 
	"""

	crys_molid = molecule.load('pdb', CRYS_PDB)
	sim_molid = molecule.load('pdb', SIM_PDB)

	# coordinate_selection = "noh and water and x > -15 and x < 25 and y > -25 and y < 25 and z > -17 and z < 17"
	coordinate_selection = "noh and water within 4 of protein"

	crys_resid_to_coord = map_water_to_coordinate(crys_molid, coordinate_selection, "resid")
	sim_resid_to_coord = map_water_to_coordinate(sim_molid, coordinate_selection, "index")

	### Map crystal water resid to simulation water resid
	crys_to_sim_index = {}

	for crys_water in crys_resid_to_coord:
		crys_coord = crys_resid_to_coord[crys_water]

		closest_sim_water, closest_sim_dist = "", float('Inf')
		### Find nearest simulation water 
		for sim_water in sim_resid_to_coord:
			sim_coord = sim_resid_to_coord[sim_water]
			dist = euc_dist(crys_coord, sim_coord)
			if(dist < closest_sim_dist):
				closest_sim_water = sim_water
				closest_sim_dist = dist
		crys_to_sim_index[crys_water] = closest_sim_water


	molecule.delete(crys_molid)
	molecule.delete(sim_molid)

	return crys_to_sim_index


def get_water_and_residue_indices(TOP, crys_to_sim_index):
	"""
		Input: Topology, mapping between crystal water resid and simulation water resid
		Output: List of tuples (water_index, residue_indices)
	"""

	structure_molid = molecule.load('pdb', TOP)
	water_index_list, neighboring_resids_list = [], []
	for crys_water_index in crys_to_sim_index:
		WATER_INDEX = crys_to_sim_index[crys_water_index]
		selection = "same residue as protein within 5 of (noh and water and index " + WATER_INDEX + ")"
		neighbor_residues_atoms = evaltcl("set neighbor_residues_atoms [atomselect " + str(structure_molid) + " \" " + selection + " \" frame 0]")
		resid = str(evaltcl("$neighbor_residues_atoms get resid")).split(" ")

		if(resid == ['']):
			print("No neighboring residues within 5 A of water " + str(WATER_INDEX))
			continue

		neighbor_residues = set()
		for resi in resid:
			neighbor_residues.add(resi)


		NEIGHBORING_RESIDS = " ".join(map(str, sorted(map(int, list(neighbor_residues)))))

		water_index_list.append(WATER_INDEX)
		neighboring_resids_list.append(NEIGHBORING_RESIDS)

	return water_index_list, neighboring_resids_list




def calc_water_occupancy(input_args):
	"""
		Worker method for calculating water stability of a set of positions

		BUG: Make this work generically. Currently printing out all 0's 
	"""

	traj_idx, TOP, TRAJ, DISTANCE_THRESH, water_index_list, neighboring_resids_list = input_args
	print("Processing fragment %s  Path: %s" % (traj_idx, TRAJ))


	### Initialize output variables
	num_waters = len(water_index_list)
	frames_with_water_density = np.array([0]*num_waters)

	traj_molid = molecule.load('pdb', TOP)
	molecule.read(traj_molid, 'netcdf', TRAJ, beg=0, end=-1, waitfor=-1)

	nFrames = molecule.numframes(traj_molid) - 1
	total_frames = np.array([nFrames]*num_waters)


	### Get reference crystallographic water position

	water_ref_coord = []
	for idx, WATER_INDEX in enumerate(water_index_list):
		# print(WATER_INDEX)

		water_ref = evaltcl("set water_ref [atomselect top \" noh and water and index " + WATER_INDEX + " and within 5 of protein\" frame 0]")
		resname = str(evaltcl("$water_ref get resname")).split(" ")
		name = evaltcl("$water_ref get name")
		index = evaltcl("$water_ref get index")
		xcoord, ycoord, zcoord = float(evaltcl("$water_ref get x")), float(evaltcl("$water_ref get y")), float(evaltcl("$water_ref get z"))
		print(xcoord, ycoord, zcoord)
		evaltcl("$water_ref delete")

		water_ref_coord.append([xcoord, ycoord, zcoord])

	water_ref_coord = np.array(water_ref_coord)

	### Align frames

	for idx, NEIGHBORING_RESIDS in enumerate(neighboring_resids_list):
		# if(idx != 2): continue
		refsel = atomsel("protein and resid " + NEIGHBORING_RESIDS, molid = traj_molid, frame = 0)
		for i in range(nFrames + 1):
			if(i == 0): continue

			print("Aligning frame " + str(i) + " ... for Water Index " + str(water_index_list[idx]) + " and resids " + NEIGHBORING_RESIDS)
			molecule.set_frame(traj_molid, i)
			b = atomsel("protein and resid " + NEIGHBORING_RESIDS, molid = traj_molid, frame = i)
			T = b.fit(refsel)
			atomsel("all", molid = traj_molid, frame = i).move(T)

			neighbor_waters = evaltcl("set sel [atomselect top \" noh and water within 5 of (protein and resid " + NEIGHBORING_RESIDS + ") \" frame " + str(i) + "]")
			resname = str(evaltcl("$sel get resname")).split(" ")


			### Skip if no waters to be found
			xstr, ystr, zstr = evaltcl("$sel get x"), evaltcl("$sel get y"), evaltcl("$sel get z")
			if(len(xstr) == 0 or len(ystr) == 0 or len(zstr) == 0): continue

			xcoord = map(float, xstr.split(" "))
			ycoord = map(float, ystr.split(" "))
			zcoord = map(float, zstr.split(" "))


			### Check distance of each candidate neighbor water 
			water_in_region = False
			for i in range(len(xcoord)):
				neighbor_water_coord = np.array([xcoord[i], ycoord[i], zcoord[i]])
				dist = euc_dist(water_ref_coord[idx], neighbor_water_coord)
				if(dist < DISTANCE_THRESH):
					# print(neighbor_water_coord, dist)
					water_in_region = True
					break

			if(water_in_region):
				frames_with_water_density[idx] += 1

			evaltcl("$sel delete")

	molecule.delete(traj_molid)

	print("\n\n\n%s/%s frames with water occupancy" % (frames_with_water_density, total_frames))
	return (traj_idx, total_frames, frames_with_water_density)





def compute_water_stability(CRYS_REF, TOP, TRAJ, OUTPUT, DISTANCE_THRESH):
	"""
		Compute the water stability of specified water index position
	"""

	### Find mapping between crystal structure water resid and simulation resid 
	crys_to_sim_index = map_crystal_to_sim(CRYS_REF, TOP)
	sim_to_crys_index = {crys_to_sim_index[k]: k for k in crys_to_sim_index}

	water_index_list, neighboring_resids_list = get_water_and_residue_indices(TOP, crys_to_sim_index)


	### Preparing input arguments 
	traj_fragments_list = [TRAJ]
	total_frames = 0
	frames_with_water_density = 0

	### Iterate through each trajectory fragment 
	output = []
	for traj_idx, TRAJ in enumerate(traj_fragments_list):
		# if(traj_idx > 3): break
		input_args = (traj_idx, TOP, TRAJ, DISTANCE_THRESH, water_index_list, neighboring_resids_list)
		output.append(calc_water_occupancy(input_args))


	### Piece together the trajectory statistics to get overall water stability for each position
	num_waters = len(water_index_list)
	print("num_waters", num_waters)
	tot, occupied = np.array([0]*num_waters), np.array([0]*num_waters)
	for traj_idx, total_frames, frames_with_water_density in output:
		tot += total_frames
		occupied += frames_with_water_density


	OUT_DIR = "/".join(OUTPUT.split("/")[:-1])
	if not os.path.exists(OUT_DIR):
		os.makedirs(OUT_DIR)
	f = open(OUTPUT, 'w')
	print("Water Stability\n\n\n")
	for idx, sim_water_index in enumerate(water_index_list):
		crys_water_resid = sim_to_crys_index[sim_water_index]
		water_stability = float(occupied[idx])/tot[idx]
		print(crys_water_resid + "\t" + str(water_stability))
		f.write(crys_water_resid + "\t" + str(water_stability) + "\n")




if __name__ == "__main__":
	if(len(sys.argv) < K_MIN_ARG):
		print(USAGE_STR)
		exit(1)

	(CRYS_REF, TOP, TRAJ, OUTPUT, DISTANCE_THRESH) = (sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], float(sys.argv[5]))
	tic = datetime.datetime.now()
	compute_water_stability(CRYS_REF, TOP, TRAJ, OUTPUT, DISTANCE_THRESH)
	toc = datetime.datetime.now()
	print("Serial Computation Time ... Start ", tic, " ... Finish ", toc, " Difference: ", (toc-tic))






