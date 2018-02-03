# Author: Anthony Kai Kwang Ma
# Email: anthonyma27@gmail.com, akma327@stanford.edu
# Date: 07/13/17
# utils.py

import os
import vmd, molecule, sys
from vmd import *
from atomsel import *
from multiprocessing import *
import datetime 
import glob 
import re
import numpy as np

def atoi(text):
	return int(text) if text.isdigit() else text

def natural_keys(text):
	return [ atoi(c) for c in re.split('(\d+)', text) ]

def euc_dist(pt1, pt2):
	a, b = np.array(pt1), np.array(pt2)
	return np.linalg.norm(a-b)

def get_file_type(file_name):
	"""
	Determine file type by extracting suffix of file_name
	"""
	if file_name is None:
		return None
	file_type = file_name.split(".")[-1].strip()
	if file_type == "nc":
		file_type = "netcdf"
	if file_type == "prmtop":
		file_type = "parm7"
	return file_type

def load_traj(TOP, TRAJ, beg_frame, end_frame, stride):
	"""
	Loads in topology and trajectory into VMD
	Parameters
	----------
	TOP: MD Topology
	TRAJ: MD Trajectory
	beg_frame: int
	end_frame: int
	stride: int
	Returns
	-------
	trajid: int
	simulation molid object
	"""
	top_file_type = get_file_type(TOP)
	traj_file_type = get_file_type(TRAJ)
	trajid = molecule.load(top_file_type, TOP)
	if TRAJ is not None:
		molecule.read(trajid, traj_file_type, TRAJ, beg=beg_frame, end=end_frame, skip=stride, waitfor=-1)
	else:
		molecule.read(trajid, top_file_type, TOP, beg=beg_frame, end=end_frame, skip=stride, waitfor=-1)
	return trajid


def estimate_simulation_length(TOP, TRAJ):
	"""
	Estimates an upper bound for simulation length for the purpose
	allocating chunks of frames for parallelization
	"""

	trajid = load_traj(TOP, TRAJ, 0, -1, 100)

	num_subsampled_frames = int(evaltcl("molinfo %s get numframes" % (trajid)))
	num_frames = (num_subsampled_frames - 1)*100
	return num_frames
