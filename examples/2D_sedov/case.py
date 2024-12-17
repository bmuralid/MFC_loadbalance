#!/usr/bin/env python3

import json
import math

# Configuring case dictionary
print(json.dumps({
    # Logistics ================================================
    'run_time_info'                : 'T',
    # ==========================================================

    # Computational Domain Parameters ==========================
    'x_domain%beg'                 : -0.5,
    'x_domain%end'                 : 0.5,
    'y_domain%beg'                 : -0.5,
    'y_domain%end'                 : 0.5,
    'm'                            : 767,
    'n'                            : 767,
    'p'                            : 0,
    'dt'                           : 1e-07,
    't_step_start'                 : 0,
    't_step_stop'                  : 10000,
    't_step_save'                  : 250,
    # ==========================================================

    # Simulation Algorithm Parameters ==========================
    'num_patches'                  : 1,
    'model_eqns'                   : 2,
    'alt_soundspeed'               : 'F',
    'num_fluids'                   : 1,
    'mpp_lim'                      : 'F',
    'mixture_err'                  : 'T',
    'time_stepper'                 : 3,
    'mp_weno'                      : 'F',
    # 'recon_type'                   : 1,
    'weno_order'                   : 5,
    'weno_eps'                     : 1e-16,
    #'muscl_order'                  : 2,
    #'muscl_lim'                    : 1,
    'riemann_solver'               : 2,
    'wave_speeds'                  : 1,
    'avg_state'                    : 2,
    'bc_x%beg'                     : -3,
    'bc_x%end'                     : -3,
    'bc_y%beg'                     : -3,
    'bc_y%end'                     : -3,
    # ==========================================================

    # Formatted Database Files Structure Parameters ============
    'format'                       : 1,
    'precision'                    : 2,
    'prim_vars_wrt'                :'T',
    'parallel_io'                  :'T',
    # ==========================================================
                                                            
    # Patch 1: Base ============================================
    'patch_icpp(1)%geometry'       : 7,
    'patch_icpp(1)%hcid'           : 205,  
    'patch_icpp(1)%x_centroid'     : 0.0,
    'patch_icpp(1)%y_centroid'     : 0.0,
    'patch_icpp(1)%length_x'       : 1.0,
    'patch_icpp(1)%length_y'       : 1.0,
    'patch_icpp(1)%vel(1)'         : 0.0,
    'patch_icpp(1)%vel(2)'         : 0.0,
    'patch_icpp(1)%pres'           : 0.67e-4,
    'patch_icpp(1)%alpha_rho(1)'   : 1.0,
    'patch_icpp(1)%alpha(1)'       : 1,
    # ==========================================================
    # Fluids Physical Parameters ===============================
    'fluid_pp(1)%gamma'            : 1.E+00/(1.4E+00-1.E+00),
    'fluid_pp(1)%pi_inf'           : 0.E+00,
    # ==========================================================
}))
# ==============================================================================
