from __future__ import division
#from mmtbx.nm import tools
from mmtbx_nm_ext import *
import mmtbx.f_model
import mmtbx.model
from mmtbx import monomer_library
import mmtbx.monomer_library.server
import mmtbx.monomer_library.pdb_interpretation
from cctbx.array_family import flex
from libtbx.test_utils import approx_equal
from libtbx.utils import format_cpu_times
from libtbx.utils import user_plus_sys_time
import time
import libtbx.load_env
import iotbx.pdb
import random
import sys, os
import struct

def run(args):
    #    if (len(args) == 0 ):
#        raise RuntimeError("Please specify one or more pdb file names.")
    mon_lib_srv = monomer_library.server.server()
    ener_lib = monomer_library.server.ener_lib()
#   pdb_file = libtbx.env.find_in_repositories(
#          relative_path="phenix_regression/pdb/1A32.pdb",
#          test=os.path.isfile)
    pdb_file = "./4KI8.pdb"
#   print pdb_file
    processed_pdb_file = monomer_library.pdb_interpretation.process(
           mon_lib_srv               = mon_lib_srv,
           ener_lib                  = ener_lib,
           file_name                 = pdb_file,
           raw_records               = None,
           force_symmetry            = True)
    xray_structure = processed_pdb_file.xray_structure()
    pdb_inp = iotbx.pdb.input(file_name=pdb_file)
    pdb_hierarchy = pdb_inp.construct_hierarchy()
    xray_structure.scattering_type_registry(table = "wk1995")
    xray_structure.convert_to_isotropic()
    u_iso_start = xray_structure.extract_u_iso_or_u_equiv()
    xray_structure.convert_to_anisotropic()
    atoms = pdb_hierarchy.atoms()
    atomsxyz = atoms.extract_xyz()
    sites_cart = xray_structure.sites_cart()
    atomic_weights = xray_structure.atomic_weights()
    t1 = time.time()
    nm_init_manager = nm_init(filename = "4ki8_evec.dat",
            n_modes = 50,
            atoms = atoms,
            zero_mode_input_flag = False,
            zero_mode_flag = True)
    eigenvector = nm_init_manager.return_modes(1)
    selections = []
    selection_strings = ["chain A", "chain B", "chain C"]
    for string in selection_strings:
        selections.append(processed_pdb_file.all_chain_proxies.selection(
                                                                string = string))
#    modes = tools.generate_evec(selections = selections,
#                                xray_structure = xray_structure,
#                                pdb_hierarchy = pdb_hierarchy,
#                                filename = "4ki8_evec.dat",
#                                n_modes = 50)
    eigenA = eigenvector.select(selections[0])
    for selection in selections:
        sites_cart_selected = xray_structure.sites_cart().select(selection)
        atomic_weights_selected = xray_structure.atomic_weights().select(selection)
        nm_init_manager.gen_zero_modes(sites_cart_selected, atomic_weights_selected)
#        zero_mode0 = nm_init_manager.return_zero_modes(0)
    t2 = time.time()
#    eigenvector.set_selected(selections[2], zero_mode0)
    t3 = time.time()
#    print t2 - t1
    print t3 - t1

if (__name__ == "__main__"):
    run(sys.argv[1:])
