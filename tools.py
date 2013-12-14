from __future__ import division
from iotbx import pdb
from cctbx.array_family import flex
import math, random
import sys, time
from scitbx import lbfgs
from mmtbx_nm_ext import *
from libtbx import adopt_init_args
from libtbx.test_utils import approx_equal
from cctbx import adptbx
from cctbx import xray
from libtbx.utils import user_plus_sys_time, Sorry
from libtbx.str_utils import line_breaker
from libtbx import group_args
import mmtbx.utils
import iotbx
import os

time_generate_evec = 0.0
time_convert_modes = 0.0

def show_time(out = None):
    if(out is None): out = sys.stdout
    total = time_generate_evec + time_convert_modes
    if(total > 0.01):
        print >> out, "NM refinement:"
        print >> out, " time_generate_evec = %-7.2f"% time_generate_evec
        print >> out, " time_convert_modes = %-7.2f"% time_convert_modes
#        print >> out, " time_init_nm_adp   = %-7.2f"% time_init_nm_adp
    return total

class nm_refinement(object):
  def __init__(self,
               fmodel,
               model,
               selections,
               selections_1d):
    pass

def read_nmval_file(evalin = "./eigenvalues.dat",
                    n_modes = 20,
                    zero_mode_input_flag = True,
                    zero_mode_flag = False):
    assert os.path.isfile(evalin) is True, "Cannot find: %s" % evalin 
    nmval = []
    if zero_mode_flag and not zero_mode_input_flag:
        nmode_start = 7
    else:
        nmode_start = 1
    print "reading eigenvalues..."
    for i in range(nmode_start - 1):
        nmval.append(0.0)
    i = 0
    with open(evalin, 'r') as file_eval:
        for line in file_eval:
            i += 1
            if i <= n_modes:
                nmval.append(float(line.split()[1]))
            else:
                break
    if zero_mode_flag:
        for i in range(6):
            if abs(nmval[i]) > 1.E-4:
                print "check modes: zero_modes may not be provided"
                nmval[i] = 0.0
    return nmval

def read_nmvec_file(pdb_hierarchy,
                    evecin = "eigenvectors.dat",
                    n_modes = 20,
                    zero_mode_input_flag = False,
                    zero_mode_flag = True):
    assert os.path.isfile(evecin) is True, "Cannot find: %s" % evecin 
    nmvec = []
    nmindex = []
    for model in pdb_hierarchy.models():
        for chain in model.chains():
            for residue_group in chain.residue_groups():
                for atom_group in residue_group.atom_groups():
                    for atom in atom_group.atoms():
                        nmvec.append([])
                        nmindex.append(0)
    if zero_mode_flag and not zero_mode_input_flag:
        nmode_start = 7
    else:
        nmode_start = 1
    print "reading eigenvectors..."
    atoms = pdb_hierarchy.atoms()
    with open("./eigenvectors.dat", "rb") as f:
        align_mark = 0
        for i in range(len(atoms)):
            nmvec_data = f.read(3*8)
            atom_index = struct.unpack('8s',nmvec_data[0:8])
            res_name = struct.unpack('3s',nmvec_data[8+5:2*8])
            atm_name = struct.unpack('3s',nmvec_data[2*8+5:3*8])
            print atom_index, res_name, atm_name
            for j in range(align_mark, len(atoms)):
                if res_name == atoms[j].parent().resname and atm_name == atoms[j].name:
                    nmindex[j] = i
                    align_mark = j+1
            if atom_index[0] == "END     ":
                break   

def generate_evec(selections,
                  xray_structure,
                  pdb_hierarchy,
                  filename,
                  n_modes,
                  zero_mode_input_flag = False,
                  zero_mode_flag = True):
    global time_generate_evec
    t1 = time.time()
    atoms = pdb_hierarchy.atoms()
    nm_init_manager = nm_init(filename = filename,
                            n_modes = n_modes, 
                            atoms = atoms, 
                            zero_mode_input_flag = zero_mode_input_flag, 
                            zero_mode_flag = zero_mode_flag)
    modes = []
    for i in range(n_modes):
        modes.append(nm_init_manager.return_modes(i))
    if zero_mode_flag == True and zero_mode_input_flag == False:
        count = 0
        for selection in selections:
            sites_cart_selected = xray_structure.sites_cart().select(selection)
            atomic_weights_selected = xray_structure.atomic_weights().select(selection)
            nm_init_manager.gen_zero_modes(sites_cart_selected, atomic_weights_selected)        
            padd = 6*count
            for i in range(6):
                selected_zero_modes = nm_init_manager.return_zero_modes(i+padd)
                modes[i].set_selected(selection, selected_zero_modes)
            count += 1
    t2 = time.time()
    time_generate_evec += (t2 - t1)
    return modes

def selected_modes_to_1D(modes,
                        n_modes,
                        selection):
    global time_convert_modes
    t1 = time.time()
    assert len(modes) == n_modes
    for i in range(n_modes):
        assert len(modes[i]) == len(modes[0])
    len_selected_modes = len(modes[0].select(selection))
#    modes1d = flex.vec3_double(len_selected_modes*n_modes, [0,0,0])
    modes1d = []
    for i in range(n_modes):
        modes_i_selected = modes[i].select(selection)
        for j in range(len_selected_modes):
            modes1d.append(modes_i_selected[j])
    modes1d = flex.vec3_double(modes1d)
    t2 = time.time()
    time_convert_modes += (t2 - t1)
    return modes1d

def init_nm_para(nmval, n_modes, zero_mode_flag = True):
    freq = []
    for i in range(n_modes):
        freq.append(0.0)
    if zero_mode_flag:
        nstart = 6
    else:
        nstart = 0
    for i in range(nstart, n_modes):
        if nmval[i] < 0.0:
            print "warning: non-positive modes"
            nmval[i] = abs(nmval[i])
        freq[i] = max(math.sqrt(nmval[i]), 1.e-2)
        
    ave_freq = 0.0
    for i in range(nstart, n_modes):
        ave_freq += freq[i]
    if nstart != n_modes:
        ave_freq = ave_freq/(n_modes - nstart)
    x = flex.double(n_modes*n_modes, 0.0)
    if zero_mode_flag:
        for i in range(6, n_modes):
            x[i+i*n_modes] = 1.0/freq[i]
        if n_modes < 6:
            mag = 1.0
        else:
            mag = .3/ave_freq
        for i in range(6):
            for j in range(i+1):
                x[j+i*n_modes] = mag*random.uniform(-1, 1)
        for i in range(7, n_modes):
            for j in range(6, i):
                x[j+i*n_modes] = mag*random.uniform(-1, 1)
    else:
        for i in range(n_modes):
            x[i+i*n_modes] = 1.0/freq[i]
        mag = .3/ave_freq
        for i in range(1, n_modes):
            for j in range(i):
                s[j+i*n_modes] = mag*random.uniform(-1, 1)
    return x

def convert_nm_adp(modes,
                    n_modes,
                    weights,
                    zero_mode_flag = True):
    global time_init_nm_adp
    t1 = time.time()
    adp_nma = init_nm_adp(modes = modes,
                        weights = weights,
                        n_mode = n_modes,
                        zero_mode_flag = zero_mode_flag)
    t2 = time.time()
    time_init_nm_adp  += (t2 - t1)
    return adp_nma
