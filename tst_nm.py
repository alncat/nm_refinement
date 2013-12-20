from __future__ import division
from mmtbx.nm import tools
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
from pprint import pprint
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
    u_cart = xray_structure.scatterers().extract_u_cart(xray_structure.unit_cell())
#    nm_init_manager = nm_init(filename = "4ki8_evec.dat",
#            n_modes = 50,
#            atoms = atoms,
#            zero_mode_input_flag = False,
#            zero_mode_flag = True)
#    eigenvector = nm_init_manager.return_modes(1)
    selections = []
    selection_strings = ["chain A", "chain B", "chain C", "chain D", "chain E", "chain F", "chain G"]
    for string in selection_strings:
        selections.append(processed_pdb_file.all_chain_proxies.selection(
                                                                string = string))
    modes = tools.generate_evec(selections = selections,
                                xray_structure = xray_structure,
                                pdb_hierarchy = pdb_hierarchy,
                                filename = "4ki8_evec.dat",
                                n_modes = 20)
    time_nm_from_uaniso = 0.0
    uanisos = flex.sym_mat3_double(xray_structure.sites_cart().size(), [0,0,0,0,0,0])
    nmval = tools.read_nmval_file(evalin = "4ki8_eval.dat",
                                  n_modes = 20,
                                  zero_mode_input_flag = False,
                                  zero_mode_flag = True)
    xs_initial = []
    for selection in selections:
        x = tools.init_nm_para(nmval = nmval,
                               n_modes = 20)
        modes1d = tools.selected_modes_to_1D(modes = modes, n_modes = 20, selection = selection)
        weights_selected = xray_structure.atomic_weights().select(selection)
        u_cart_selected = u_cart.select(selection)
        adp_nma = init_nm_adp(modes = modes1d,
                              weights = weights_selected,
                              n_modes = 20,
                              zero_mode_flag = True)
        uaniso_from_s_manager = uaniso_from_s(x = x,                            
                                              adp_nma = adp_nma, 
                                              weights = weights_selected,
                                              n_modes = 20,
                                              zero_mode_flag = True)
        u = uaniso_from_s_manager.u_cart()
        x_scaled = scale_x(x = x,
                           uanisos = u_cart_selected,
                           adp_all = u,
                           n_modes = 20,
                           zero_mode_flag = True)
        xs_initial.append(x_scaled)
        del adp_nma
    t1 = time.time()
    nm_from_uanisos = tools.nm_from_uanisos(xray_structure = xray_structure,
                                            selections = selections,
                                            modes      = modes,
                                            xs_initial = xs_initial,
                                            n_modes    = 20,
                                            number_of_macro_cycles = 5,
                                            verbose    = 1)
    t2 = time.time()
    print t2 - t1
    tools.show_time()

if (__name__ == "__main__"):
    run(sys.argv[1:])
