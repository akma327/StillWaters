# Author: Anthony Kai Kwang Ma
# Email: anthonyma27@gmail.com, akma327@stanford.edu
# Date: 07/18/17
# color_water_stability.py

import pymol
from pymol import cmd, stored
import numpy as np
import math


USAGE_STR = """
# Purpose
# After computing water stability, map the occupancy of crystallographic water positions
# to the water molecules.

"""

def occupancy_to_color(occupancy):
	"""
		Return corresponding color from blue to red spectrum
	"""

	# color = [255*(1-isoval), 255*(1-isoval), 255] # white to blue
	# color = [255*isoval, 0, 255*(1-isoval)] # blue to red 
	color = "gray" + str(int(math.floor((1-occupancy)*100))).zfill(2)
	return color


def map_water_stability(TOP, WATER_STABILITY):
	top_name = TOP.split("/")[-1].split(".")[0]
	cmd.bg_color("white")
	cmd.load(TOP)
	cmd.hide()
	cmd.show("cartoon")
	cmd.cartoon("loop")
	cmd.do("set cartoon_transparency, 0.7")
	cmd.do("set sphere_scale, 0.75")
	cmd.color("white", top_name)

	f = open(WATER_STABILITY, 'r')
	for line in f:
		linfo = line.strip().split("\t")
		resi, occupancy = linfo[0], float(linfo[1])
		color = occupancy_to_color(occupancy)
		print(color, occupancy)
		cmd.do("show spheres, resn HOH and resi " + str(resi))
		cmd.color(color, "resn HOH and resi " + str(resi))



# DOR
# TOP="/Users/anthony/Desktop/dror/dynamic-networks/DynamicNetworks/data/water-density-with-neighbors/DOR-inactive-naltrindole-unpublished/condition-naltrindole-bound/4n6h_aligned_to_sim.pdb"
# WATER_STABILITY="/Users/anthony/Desktop/dror/dynamic-networks/DynamicNetworks/data/water-stability/dor_inactive_rep1_allwaters/occupancy-data/water_stability_1Ang.txt"

# DOR Inactive
# TOP="/Users/anthony/Desktop/dror/dynamic-networks/DynamicNetworks/data/water-stability/072717_a2a_mor/crystal_pdbs/4n6h_aligned.pdb"
# WATER_STABILITY="/Users/anthony/Desktop/dror/dynamic-networks/DynamicNetworks/data/water-stability/072717_a2a_mor/stability_values/dor_inactive_water_stability_1Ang.txt"


# # MOR Active
# TOP="/Users/anthony/Desktop/dror/dynamic-networks/DynamicNetworks/data/water-stability/072717_a2a_mor/crystal_pdbs/5c1m_aligned.pdb"
# WATER_STABILITY="/Users/anthony/Desktop/dror/dynamic-networks/DynamicNetworks/data/water-stability/072717_a2a_mor/stability_values/mor_active_water_stability_1Ang.txt"

# # A2A Inactive
TOP="/Users/anthony/Desktop/dror/dynamic-networks/DynamicNetworks/data/water-stability/072717_a2a_mor/crystal_pdbs/5iu4_aligned.pdb"
WATER_STABILITY="/Users/anthony/Desktop/dror/dynamic-networks/DynamicNetworks/data/water-stability/072717_a2a_mor/stability_values/a2a_inactive_water_stability_1Ang.txt"

map_water_stability(TOP, WATER_STABILITY)



