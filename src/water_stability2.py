# Author: Anthony Kai Kwang Ma
# Email: anthony.ma@yale.edu
# Date: 02/02/18
# water_stability2.py

from utils import *

USAGE_STR = """

# Purpose
# Aligns every frame of trajectory to topology. For each of these
# frames determine whether there are water molecules within a 
# specified distance from a particular crystallographic water position. 

# Usage 
# python water_stability.py <crys_struc> <top> <traj> <output> <distance_thresh>

# Arguments
# <crys_struc> Absolute path to topology of crystal structure pre-aligned to simulation topology
# <top> Absolute path to simulation topology 
# <traj> Absolute path to trajectory file
# <output> output file with stability values for each water in crystal structure 
# <distance_thresh> Distance threshold to search for neighboring waters near a certain region. In nanometers

# Example

crys_struc="/Users/anthony/Desktop/sherlock/MIF-waters/data/simulation/3DJH_wb/3djh_crys_aligned.pdb"
top="/Users/anthony/Desktop/sherlock/MIF-waters/data/simulation/3DJH_wb/3DJH_wb.pdb"
traj="/Users/anthony/Desktop/sherlock/MIF-waters/data/simulation/3DJH_wb/3DJH_run.1.dcd"
output="/Users/anthony/Desktop/sherlock/MIF-waters/data/water_stability/water_stability_1Ang.txt"
distance_thresh=1.0
cd /Users/anthony/Desktop/dror/local-dev/github/StillWaters/src
python water_stability2.py $crys_struc $top $traj $output $distance_thresh


# TODO: Bugs
- Accept any vmd input file 
- Check if the simulation pdb and crystal pdb must match for alignment, can work around it. 
- Take in long trajectory and perform internal chunking and parallelization 
- Verify output sanity check 

"""

K_MIN_ARG = 6

##############################################################################
# Global Variables
##############################################################################
TRAJ_FRAG_SIZE = 100

##############################################################################
# Functions
##############################################################################

def map_water_to_coordinate(structure, cutoff=4.0):
	"""
	Compute coordinates for every water molecule within cutoff distance from protein. 

	Parameters
	----------
	structure: string
		Path to structure to compute coordinates for waters near protein
	cutoff: float
		Cutoff distance for waters near protein to consider 

	Returns
	-------
	water_to_coord: dictionary 
		Mapping between water resid index to x,y,z coordinates 
	"""
	top_file_type = get_file_type(structure)
	crys_molid = molecule.load(top_file_type, structure)

	waters = evaltcl("set waters [atomselect %s \" noh and water within %s of protein \" frame 0]" % (crys_molid, cutoff))
	resid = map(int, (evaltcl("$waters get resid")).split(" "))
	indices = map(int, (evaltcl("$waters get index")).split(" "))
	xcoord = map(float, evaltcl("$waters get x").split(" "))
	ycoord = map(float, evaltcl("$waters get y").split(" "))
	zcoord = map(float, evaltcl("$waters get z").split(" "))

	water_to_coord = {}
	for idx, resi in enumerate(resid):
		water_to_coord[resi] = np.array([xcoord[idx], ycoord[idx], zcoord[idx]])

	molecule.delete(crys_molid)
	return water_to_coord



def map_waters_to_neighbor_resids(crys_struc, water_to_coord):
	"""
	Generate mapping of neighboring residues to every crystallographic water position 

	Parameters
	----------
	crys_struc: string
		Path to crystal structure to load into VMD
	water_to_coord: dictionary 
		Mapping between water resid index to x,y,z coordinates 

	Returns
	-------
	crys_water_neighboring_resids: dictionary 
		Mapping between water resid index to list of neighboring protein resid indices 
	"""

	top_file_type = get_file_type(crys_struc)
	crys_molid = molecule.load(top_file_type, crys_struc)

	### Generate mapping between water resid index to neighboring protein resid indices
	PROTEIN_BOUNDING_BOX = 5.0
	crys_water_neighboring_resids = {}
	for water_resid_idx in water_to_coord: 
		water_coord = water_to_coord[water_resid_idx]

		### Identify protein residues within 5 x 5 x 5 A bounding box of crystallographic water coordinates
		x_min = water_coord[0] - PROTEIN_BOUNDING_BOX/2
		x_max = water_coord[0] + PROTEIN_BOUNDING_BOX/2
		y_min = water_coord[1] - PROTEIN_BOUNDING_BOX/2
		y_max = water_coord[1] + PROTEIN_BOUNDING_BOX/2
		z_min = water_coord[2] - PROTEIN_BOUNDING_BOX/2
		z_max = water_coord[2] + PROTEIN_BOUNDING_BOX/2
		selection = "same residue as protein and (x > %s and x < %s and y > %s and y < %s and z > %s and z < %s)" % (x_min, x_max, y_min, y_max, z_min, z_max)
		neighbor_residue_atoms = evaltcl("set neighbor_residue_atoms [atomselect %s \" %s \" frame 0]" % (crys_molid, selection))
		neighbor_resid_list = str(evaltcl("$neighbor_residue_atoms get resid")).split(" ")

		if (neighbor_resid_list == ['']):
			crys_water_neighboring_resids[water_resid_idx] = []
		else:
			crys_water_neighboring_resids[water_resid_idx] = map(int, neighbor_resid_list)

	molecule.delete(crys_molid)
	return crys_water_neighboring_resids



def calc_water_occupancy(frag_idx, beg_frame, end_frame, top, traj, distance_thresh, water_to_coord, crys_water_neighboring_resids):
	"""
	Worker method for calculating water occupancy at the crystallographic water positions throughput simulation

	Parameters
	----------
	frag_idx: int 
		Fragment of simulation
	beg_frame: int
		Start frame of trajectory fragment
	end_frame: int
		End frame of trajectory fragment
	top: string
		Path to simulation topology
	traj: string
		Path to simulation trajectory
	distance_thresh: float
		Threshold distance to determine whether a water was within the crystallographic position during simulation 
	water_to_coord: dictionary 
		Mapping between water resid index to x,y,z coordinates
	crys_water_neighboring_resids: dictionary
		Mapping between water resid index to list of neighboring protein resid indices 

	Returns
	-------
	frag_idx: int
		Fragment of simulation
	frames_with_water_density: np.array
		For every water in the simulation, counter of how many frames a water was within cutoff distance to crystallographic position 
	total_frames: np.array
		Total number of frames for this simulation fragment 
	"""

	print (frag_idx, beg_frame, end_frame)

	### Read in simulation 
	traj_molid = load_traj(top, traj, beg_frame, end_frame, 1)
	num_frames = molecule.numframes(traj_molid) - 1

	### Initialize output variables
	num_waters = len(water_to_coord)
	frames_with_water_density = np.array([0]*num_waters)
	total_frames = np.array([num_frames]*num_waters)

	### Align frames
	water_index_list = sorted(water_to_coord.keys())
	for idx, water_resid_idx in enumerate(water_index_list):

		# if(idx != 2): continue # Testing

		neighboring_residues = " ".join(map(str, crys_water_neighboring_resids[water_resid_idx]))
		# If no neighboring residues ignore the water
		if(neighboring_residues == ''): continue
		refsel = atomsel("protein and resid " + neighboring_residues, molid = traj_molid, frame = 0)
		for i in range(num_frames + 1):
			if(i == 0): continue # Exclude first frame representing pdb 
			if(i > 2): continue # Testing

			print("Aligning frame %s ... for water index %s and protein resids %s" % (i, water_resid_idx, neighboring_residues))
			molecule.set_frame(traj_molid, i)
			b = atomsel("protein and resid " + neighboring_residues, molid = traj_molid, frame = i)
			T = b.fit(refsel)
			atomsel("all", molid = traj_molid, frame = i).move(T)

			neighbor_waters = evaltcl("set sel [atomselect top \" noh and water within 5 of (protein and resid %s ) \" frame %s]" % (neighboring_residues, i))
			resname = str(evaltcl("$sel get resname")).split(" ")

			### Skip if no waters to be found
			xstr = evaltcl("$sel get x")
			ystr = evaltcl("$sel get y")
			zstr = evaltcl("$sel get z")
			if(len(xstr) == 0 or len(ystr) == 0 or len(zstr) == 0): continue

			xcoord = map(float, xstr.split(" "))
			ycoord = map(float, ystr.split(" "))
			zcoord = map(float, zstr.split(" "))

			### Check distance of each candidate neighbor water 
			water_in_region = False
			for i in range(len(xcoord)):
				neighbor_water_coord = np.array([xcoord[i], ycoord[i], zcoord[i]])

				# Euclidean distance between crystallographic water position and neighboring waters
				dist = euc_dist(water_to_coord[water_resid_idx], neighbor_water_coord)
				if(dist < distance_thresh):
					water_in_region = True 
					break 
			if(water_in_region):
				frames_with_water_density[idx] += 1
			evaltcl("$sel delete")

	molecule.delete(traj_molid)

	print("\n\n\n%s/%s frames with water occupancy" % (frames_with_water_density, total_frames))
	return frag_idx, frames_with_water_density, total_frames

def calc_water_occupancy_helper(args):
	return calc_water_occupancy(*args)

def compute_water_stability(crys_struc, top, traj, output, distance_thresh = 1.0, cores=4):
	"""
	Compute water stability for all waters in the crystal structure throughout MD simulation

	Parameters
	----------
	crys_struc: string
		Path to crystal structure 
	top: string
		Path to simulation topology
	traj: string
		Path to simulation trajectory
	output: string
		Path to output file with water stability values for all crystallographic waters 
	distance_thresh: float
		Threshold distance to determine whether a water was within the crystallographic position during simulation 
	"""

	crys_water_coords = map_water_to_coordinate(crys_struc)
	crys_water_neighboring_resids = map_waters_to_neighbor_resids(crys_struc, crys_water_coords)

	estimated_sim_length = estimate_simulation_length(top, traj)
	input_args = []
	for frag_idx, beg_frame in enumerate(range(0, estimated_sim_length, TRAJ_FRAG_SIZE)):
		if(frag_idx > 1): break
		end_frame = beg_frame + TRAJ_FRAG_SIZE - 1
		print("Processing fragment %s beg_frame %s end_frame %s" % (frag_idx, beg_frame, end_frame))
		input_args.append((frag_idx, beg_frame, end_frame, top, traj, distance_thresh, crys_water_coords, crys_water_neighboring_resids))

	# Parallel computation
	pool = Pool(processes=cores)
	output_frags = pool.map(calc_water_occupancy_helper, input_args)
	pool.close()
	pool.join()

	### Piece together trajectory statistics to get overall water stability for each position
	output_frags = sorted(output_frags, key = lambda x: (x[0]))
	num_waters = len(crys_water_coords)
	occupied, total = np.array([0]*num_waters), np.array([0]*num_waters)
	for frag_idx, frames_with_water_density, total_frames in output_frags:
		occupied += frames_with_water_density
		total += total_frames

	f = open(output, 'w')
	water_index_list = sorted(crys_water_coords.keys())
	for idx, water_index in enumerate(water_index_list):
		water_stability = float(occupied[idx])/total[idx]
		f.write("%s\t%s\n" % (water_index, water_stability))


if __name__ == "__main__":
	if(len(sys.argv) < K_MIN_ARG):
		print (USAGE_STR)
		exit(1)
	(crys_struc, top, traj, output, distance_thresh) = (sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], float(sys.argv[5]))
	tic = datetime.datetime.now()
	compute_water_stability(crys_struc, top, traj, output, distance_thresh)
	toc = datetime.datetime.now()
	print("Computation Time: " + str((toc-tic).total_seconds()))



