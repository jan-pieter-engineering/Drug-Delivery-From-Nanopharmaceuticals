# -----------------------------------------------------------------
# I N V E R S E   A N A L Y S I S
# -----------------------------------------------------------------
from model import run_model
from utility import calc_F_num_at_discrete_time_steps, check_cost_function
from plot_funcs import plot__F_exp_vs_F_num, print_opt_parameters, print_save__R_squared, print__R_squared_for_given_time_frame, print__residual_at_week_5
from scipy.optimize import least_squares
import numpy as np


def run_optimization(parameters_to_optimize, parameters_treated_as_constant,
                     t, t_break, t_reached_5, T__opt_end, T_explicit_Euler, dt, c_sat, c_s_t0, c_max, V_s, V_s_units, V_c_units, r_i, r_o, l_c, H,
                     A_m_t0, m_d_t0, const_distribution, parabolic_distribution, offset_parabola, U_0, U_1, beta, gamma, limiter_break, limiter_sat,
                     V, mesh, nx, nr, ds, degree, r, n_facet, num_steps, num_steps_explicit_Euler, F_num_komplex_n_Euler, runs_CN_high, fineness, tol_adaptiv_dt, safety, p, boundary_tol,
                     menu, distribution_menu, verbose, discretization_scheme, t_ref, c_ref, D_ref, l_ref, t_num, F_num,
                     t_exp, F_exp, opt_time_frame, opt_success, kind_of_1d_interpolation):
    """ Optimizes the unknown material parameters by fitting the numerical results (_t_num, _F_num) to the experimental data (t_exp, F_exp) """

    # 1. Define bounds
    if opt_time_frame == 'until_break':
        bnds = [[1e-2, 9.0, 3e-1], [1e4, 9e2, 3.0]]
    elif opt_time_frame == 'after_break':
        bnds = ([0.0, np.inf])


    # 2. Collect constants
    standard_constants = [t, t_break, t_reached_5, T__opt_end, T_explicit_Euler, dt, c_sat, c_s_t0, c_max, V_s, V_s_units, V_c_units, r_i, r_o, l_c, H, A_m_t0, m_d_t0, const_distribution, parabolic_distribution, offset_parabola, U_0, U_1, beta, gamma, limiter_break, limiter_sat,
                          V, mesh, nx, nr, ds, degree, r, n_facet, num_steps, num_steps_explicit_Euler, F_num_komplex_n_Euler, runs_CN_high, fineness, tol_adaptiv_dt, safety, p, boundary_tol,
                          menu, distribution_menu, verbose, discretization_scheme, t_ref, c_ref, D_ref, l_ref]


    # 3. Run least squares optimization
    optimized_parameters = least_squares(object_function, parameters_to_optimize, xtol=1e-08, bounds=bnds, method='trf', loss='linear', verbose=2,
                                         args=(parameters_treated_as_constant, standard_constants, t_exp, F_exp, opt_time_frame, opt_success, kind_of_1d_interpolation))


    return optimized_parameters



def object_function(parameters_to_optimize, parameters_treated_as_constant, standard_constants, t_exp, F_exp, opt_time_frame, opt_success, kind_of_1d_interpolation):
    """ Calcs weighted residuals by running the model and subtracting F_exp from the calculated _F_num and multiply that term by the weights """

    # 1. Extract needed constants
    T__opt_end        = standard_constants[3]
    A_m_t0            = standard_constants[16]
    l_ref             = standard_constants[-1]
    D_ref             = standard_constants[-2]
    c_ref             = standard_constants[-3]
    verbose           = standard_constants[-6]
    distribution_menu = standard_constants[-7]
    menu              = standard_constants[-8]


    # 2. Reset t_num & F_num (otherwise they contain the values from previous optimization steps)
    t_num = [0]
    F_num = [0]


    # 3. Calc _F_num by model
    if opt_time_frame == 'until_break':

        D            = parameters_to_optimize[0]
        alpha        = parameters_to_optimize[1]
        l_b          = parameters_to_optimize[2]
        A_m_degrad   = parameters_treated_as_constant[0]
        const_degrad = parameters_treated_as_constant[1]
        degrad_break = parameters_treated_as_constant[2]

    elif opt_time_frame == 'after_break':

        A_m_degrad   = parameters_to_optimize[0]
        D            = parameters_treated_as_constant[0]
        alpha        = parameters_treated_as_constant[1]
        l_b          = parameters_treated_as_constant[2]
        const_degrad = parameters_treated_as_constant[3]
        degrad_break = parameters_treated_as_constant[4]

    # Print parameters for the current optimization run
    print_opt_parameters(parameters_to_optimize, const_degrad, A_m_t0, c_ref, D_ref, l_ref, opt_time_frame, verbose)

    # Run model
    _t_num, _F_num = run_model(D, alpha, l_b, A_m_degrad, const_degrad, degrad_break, *standard_constants, t_num, F_num)


    # 4. Calc weighted residuals

    # Define weights for residuals
    t_opt_end = T__opt_end

    # Disable weights (_weight = 1.0)
    _weight = []
    for i in range(len([1 for i in t_exp if i <= t_opt_end])): # for T__opt_end=5/15 (no weights (1.0))
        _weight.append(1.0)

    # Choose weighting by removing #
    #_weight = [0.0001,0.0001,0.0001,0.0001,0.0001,1.0] # for T__opt_end=5
    #_weight = [1.0,0.01,0.01,0.5,0.5,1.0]              # for T__opt_end=5

    # Calc _F_num at discrete time steps (these are the time points, where measurements were taken during the experiment)
    _F_num_at_discrete_points = calc_F_num_at_discrete_time_steps(_t_num, _F_num, t_exp, T__opt_end, kind_of_1d_interpolation)

    # Calc Weighted residuals
    _weighted_residuals = [] # initialize
    for j in range(len(F_exp)):
        _weighted_residuals.append(_weight[j]*(F_exp[j] - _F_num_at_discrete_points[j]))

    # Convert list of weighted residuals to an array
    _weighted_residuals_arr = np.array(_weighted_residuals)


    # 5. Print/plot results
    #TODO: R(5,15) implementieren
    R_squared_0_T = print_save__R_squared(_F_num_at_discrete_points, F_exp, t_exp,0.0, T__opt_end, parameters_to_optimize, A_m_t0, D_ref, c_ref, l_ref, verbose, menu, distribution_menu, opt_time_frame)
    R_squared_0_1 = print__R_squared_for_given_time_frame(_F_num_at_discrete_points, F_exp, t_exp, 0.0, 1.0, verbose, opt_time_frame) # R² initial burst     (week:0-1)
    R_squared_3_5 = print__R_squared_for_given_time_frame(_F_num_at_discrete_points, F_exp, t_exp, 3.0, 5.0, verbose, opt_time_frame) # R² during stagnation (week:3-5)
    print__residual_at_week_5(_F_num_at_discrete_points, F_exp, verbose, opt_time_frame)
    plot__F_exp_vs_F_num(parameters_to_optimize, D_ref, c_ref, l_ref, t_exp, _t_num, F_exp, _F_num, T__opt_end, R_squared_0_T, kind_of_1d_interpolation, verbose, menu, opt_success, distribution_menu, opt_time_frame)
    check_cost_function(F_exp, _F_num_at_discrete_points, _weight, verbose)


    return _weighted_residuals_arr