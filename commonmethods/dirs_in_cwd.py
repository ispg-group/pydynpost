#!/usr/bin/env python
import os
import numpy as np
def get_dirs(inp_dir):
    content_in_dir = os.listdir(inp_dir)
    tmp_mask = [not("." in k) for k in content_in_dir]
    content_in_dir = np.array(content_in_dir)
    dirs_in_dir = content_in_dir[tmp_mask]
    return dirs_in_dir
