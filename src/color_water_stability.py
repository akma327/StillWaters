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
	i = 0
	for line in f:
		# if(i > 10): break
		# i += 1
		linfo = line.strip().split("\t")
		resi, occupancy = linfo[0], float(linfo[1])
		color = occupancy_to_color(occupancy)
		print(resi, color, occupancy)
		cmd.do("sele water_%s, resn HOH and resi %s" % (resi, resi))
		cmd.do("show spheres, water_%s" % (resi))
		# cmd.color(color, "water_%s" % (resi))
		# cmd.do("show spheres, resn HOH and resi " + str(resi))
		# cmd.color(color, "resn HOH and resi " + str(resi))



# MIF 
TOP="/Users/anthony/Desktop/sherlock/MIF-waters/data/simulation/3DJH_wb/3djh_crys_aligned.pdb"
WATER_STABILITY="/Users/anthony/Desktop/sherlock/MIF-waters/data/water_stability/water_stability_1Ang.txt"

map_water_stability(TOP, WATER_STABILITY)



