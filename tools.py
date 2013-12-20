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
time_nm_from_uanisos = 0.0

def show_time(out = None):
    if(out is None): out = sys.stdout
    total = time_generate_evec + time_convert_modes + time_nm_from_uanisos
    if(total > 0.01):
        print >> out, "NM refinement:"
        print >> out, " time_generate_evec     = %-7.2f"% time_generate_evec
        print >> out, " time_convert_modes     = %-7.2f"% time_convert_modes
        print >> out, " time_nm_from_uanisos   = %-7.2f"% time_nm_from_uanisos
    return total

class show_nm(object):
    def __init__(self, xs, n_modes, text="", zero_mode_flag=True, out=None):
        if(out is None): out = sys.stdout
        counter = 0
        print >> out, text
        for item in xs:
            counter += 1
            print >> out, "NM group number %d: " % counter
            nt = 0
            if zero_mode_flag == True:
                for i in range(6):
                    for j in range(i+1):
                        print >> out, '%4.4f ' % item[nt],
                        nt += 1
                    print >> out, ""
                for i in range(6, n_modes):
                    for j in range(i+1):
                        if j < 6:
                            print >> out, " "*10,
                        else:
                            print >> out, '%4.4f ' % item[nt],
                            nt += 1
                    print >> out, ""
            else:
                for i in range(n_modes):
                    for j in range(i+1):
                        print >> out, '%4.4f ' % item[nt],
                        nt += 1
                    print >> out, ""

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
    s = flex.double(n_modes*n_modes, 0.0)
    if zero_mode_flag:
        for i in range(6, n_modes):
            s[i+i*n_modes] = 1.0/freq[i]
        if n_modes < 6:
            mag = 1.0
        else:
            mag = .3/ave_freq
        for i in range(6):
            for j in range(i+1):
                s[j+i*n_modes] = mag*random.uniform(-1, 1)
        for i in range(7, n_modes):
            for j in range(6, i):
                s[j+i*n_modes] = mag*random.uniform(-1, 1)
    else:
        for i in range(n_modes):
            s[i+i*n_modes] = 1.0/freq[i]
        mag = .3/ave_freq
        for i in range(1, n_modes):
            for j in range(i):
                s[j+i*n_modes] = mag*random.uniform(-1, 1)
    if zero_mode_flag:
        n_nmpars = 21 + (n_modes - 5)*(n_modes - 6)/2
    else:
        n_nmpars = (n_modes + 1)*n_modes/2
    x = flex.double(int(n_nmpars), 0.0)
    s2x(s = s, x = x, n_modes = n_modes, zero_mode_flag = zero_mode_flag)
    return x

def s2x(s, x, n_modes, zero_mode_flag = True):
    nt = 0
    if zero_mode_flag:
        for i in range(6):
            for j in range(i+1):
                x[nt] = s[j + i*n_modes]
                nt += 1
        for i in range(6, n_modes):
            for j in range(6, i+1):
                x[nt] = s[j + i*n_modes]
                nt += 1
    else:
        for i in range(n_modes):
            for j in range(i+1):
                x[nt] = s[j + i*n_modes]
                nt += 1

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

class nm_from_uaniso_minimizer(object):
    def __init__(self,
                 uaniso,
                 x_initial,
                 adp_nma,
                 weights,
                 n_modes,
                 zero_mode_flag,
                 max_iterations):
        adopt_init_args(self, locals())
        assert self.uaniso.size() == self.weights.size()
        self.x     = self.x_initial
        self.x_min = self.x_initial
        self.n = self.x.size()
        t1 = time.time()
        self.minimizer = lbfgs.run(
                                  target_evaluator = self,
                                  termination_params = lbfgs.termination_parameters(
                                      max_iterations = max_iterations,
                                      max_calls = int(max_iterations*1.5)),
                                  exception_handling_params = 
                                  lbfgs.exception_handling_parameters(
                                      ignore_line_search_failed_step_at_lower_bound = True,
                                      ignore_line_search_failed_step_at_upper_bound = True,
                                      ignore_line_search_failed_maxfev              = True)
                                  )
        self.compute_functional_and_gradients()
        t2 = time.time()
        print t2 - t1

    def compute_functional_and_gradients(self):
        manager = nm_from_uaniso_target_and_grads(self.x,
                                                  self.weights,
                                                  self.adp_nma,
                                                  self.uaniso,
                                                  self.n_modes,
                                                  self.zero_mode_flag)
        self.f = manager.target()
        self.g = manager.grad_nm()
        return self.f, self.g

def nm_from_uanisos(xray_structure,
                     selections,
                     modes,
                     xs_initial,
                     n_modes,
                     number_of_macro_cycles = 3000,
                     max_iterations         = 1000,
                     zero_mode_flag         = True,
                     verbose                = -1,
                     out                    = None):
    global time_nm_from_uanisos
    t1 = time.time()
    if(out is None): out = sys.stdout
    if(verbose > 0):
        show_nm(xs = xs_initial,
                n_modes = n_modes,
                text = "NM from ADP: start NM values",
                zero_mode_flag = True,
                out = out)
    u_cart = xray_structure.scatterers().extract_u_cart(xray_structure.unit_cell())
    xs_min = []
    group_counter = 0
    for x_initial, selection in zip(xs_initial, selections):
        group_counter += 1
        stop_flag = 0
        target_stop = -1.0
        weights_selected = xray_structure.atomic_weights().select(selection)
        u_cart_selected = u_cart.select(selection)
        modes1d_selected = selected_modes_to_1D(modes = modes, n_modes = n_modes,
                                                selection = selection)
        assert len(modes1d_selected)/n_modes == len(weights_selected)
        adp_nma_selected = init_nm_adp(modes = modes1d_selected,
                              weights = weights_selected,
                              n_modes = n_modes,
                              zero_mode_flag = zero_mode_flag)
        for i in range(1, number_of_macro_cycles+1):
            target_start = target_stop
            minimized = nm_from_uaniso_minimizer(uaniso = u_cart_selected,
                                                 x_initial = x_initial,
                                                 adp_nma = adp_nma_selected,
                                                 weights = weights_selected,
                                                 n_modes = n_modes,
                                                 zero_mode_flag = zero_mode_flag,
                                                 max_iterations = max_iterations)
            x_initial = minimized.x_min
        if(verbose > 0):
            print >> out, "NM group %d: minimized target = " %(group_counter), minimized.f
            print >> out, "Macrocycle %d finished!" %i
        x_min_ = minimized.x_min
        xs_min.append(x_min_)
    if(verbose > 0):
        show_nm(xs = xs_min,
                n_modes = n_modes,
                text = "NM from ADP: final NM values", 
                zero_mode_flag = zero_mode_flag,
                out = out)
    t2 = time.time()
    time_nm_from_uanisos += (t2 - t1)
    return xs_min

class nm_xray_grads(object):

    def __init__(self, target_result, adp_nmas, selections, xs, n_modes, zero_mode_flag):
        self.grad_nm = []
        d_target_d_uaniso = target_result.gradients_wrt_atomic_parameters(
            u_aniso=True)
        for adp_nma, sel, x in zip(adp_nmas, selections, xs):
            d_target_d_nm_manager = d_target_d_nm(
                adp_nma = adp_nma,
                d_target_d_uaniso = d_target_d_uaniso.select(sel),
                x       = x,
                n_modes = n_modes,
                zero_mode_flag = zero_mode_flag)
            self.grad_nm.append(list(d_target_d_nm_manager.grad_nm()))

class nm_xray_target_minimizer(object):
    def __init__(self,
                 fmodel,
                 xs_initial,
                 modes,
                 selections,
                 selections_1d,
                 max_iterations,
                 n_modes,
                 run_finite_differences_test = False,
                 correct_adp = True,
                 zero_mode_flag = True):
        adopt_init_args(self, locals())
        fmodel.xray_structure.scatterers().flags_set_grads(state=False)
        xray.set_scatterer_grad_flags(scatterers = fmodel.xray_structure.scatterers(),
                                      u_aniso = True)
        if(self.run_finite_differences_test): self.correct_adp = False
        self.fmodel_copy = self.fmodel.deep_copy()
        self.target_functor = self.fmodel_copy.target_functor()
        self.run_finite_differences_test_counter = 0
        self.counter = 0
        self.n_groups = len(self.xs_initial)
        self.dim_x    = len(self.xs[0])
        self.xs_min    = self.xs_initial
        self.x = self.pack(self.xs_min)
        self.adp_nmas  = []
        for selection in selections:
            modes1d_selected = selected_modes_to_1D(modes = modes, n_modes = n_modes,
                                                selection = selection)
            assert len(modes1d_selected)/n_modes == len(weights_selected)
            adp_nma_selected = init_nm_adp(modes = modes1d_selected,
                              weights = weights_selected,
                              n_modes = n_modes,
                              zero_mode_flag = zero_mode_flag)
            self.adp_nmas.append[adp_nma_selected]
        self.minimizer = lbfgs.run(
                target_evaluator = self,
                core_params = lbfgs.core_parameters(maxfev = 10),
                termination_params = lbfgs.termination_parameters(
                    min_iterations = max_iterations,
                    max_calls = int(max_iterations*1.5)),
                exception_handling_params = lbfgs.exception_handling_parameters(
                    ignore_line_search_failed_step_at_lower_bound = True,
                    ignore_line_search_failed_step_at_upper_bound = True,
                    ignore_line_search_failed_maxfev              = True))
        self.compute_functional_and_gradients()
        self.xs_result = self.unpack()

        def pack(self, xs):
            v = []
            for x in xs:
                v += list(x)
            return flex.doule(tuple(v))

        def unpack(self):
            i = 0
            xs_min = []
            for j in xrange(self.n_groups):
                self.xs_min[j] = tuple(self.x)[i:i+self.dim_x]
                i += self.dim_x

        def compute_functional_and_gradients(self):
            self.counter += 1
            self.unpack()
            upate_xray_structure_with_nm(
                    xray_structure = self.fmodel_copy.xray_structure,
                    selections     = self.selections,
                    xs             = self.xs_min,
                    selections_1d  = self.selections_1d,
                    correct_adp    = correct_adp)
            self.fmodel_copy.update_xray_structure(update_f_calc=True)
            t_r = self.target_functor(compute_gradients=True)
            self.f = t_r.target_work()
            grad_manager = nm_xray_grads(
                    target_result = self.t_r,
                    adp_nmas      = self.adp_nams,
                    selections    = self.selections,
                    xs            = self.xs_min,
                    n_modes       = self.n_modes,
                    zero_mode_flag= self.zero_mode_flag)
            self.g = self.pack(grad_manager.grad_nm)
            if(self.run_finite_differences_test and
               self.run_finite_differences_test_counter < 2):
                pass
            return self.f, self.g

def upate_xray_structure_with_nm(xray_structure,
                                 selections,
                                 xs,
                                 selections_1d = None,
                                 correct_adp = True):
    global time_update_xray_structure_with_nm
    timer = user_plus_sys_time()
