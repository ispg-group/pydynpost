#!/usr/bin/python
import os
import numpy as np


def getDirs(inp_dir):
    # This external function finds all subdirectories in a given 
    # input directory and returns them as a numpy array of chars.
    content_in_dir = os.listdir(inp_dir)
    tmp_mask = np.array([not("." in k) for k in content_in_dir])
    content_in_dir = np.array(content_in_dir)
    dirs_in_dir = content_in_dir[tmp_mask]
    return dirs_in_dir
