# -----------------------------------------------------------------
# M E N U
# -----------------------------------------------------------------
from mesh_domains_boundaries import collect_mesh_parameters
from model import run_model
from convergence_analysis import run_convergence_analysis
from optimization import run_optimization
from utility import check_CFL_conditions, add_fitted_points_to_experimenta_data, adjust_t_exp_and_F_exp_to_optimization_time_frame
from testing import run_testing
import time
import sys


def compute_results(D, alpha, l_b, A_m_degrad, const_degrad, degrad_break,
                    t, t_break, t_reached_5, T, T_explicit_Euler, dt, dt_CN, c_sat, c_s_t0, c_max, V_s, V_s_units, V_c_units, r_i, r_o, l_c, H, A_m_t0, m_d_t0, const_distribution, parabolic_distribution, offset_parabola, U_0, U_1, beta, gamma, const_degrad_deactivated, limiter_break, limiter_sat,
                    nx, nr, h, degree, element_family, num_steps, num_steps_explicit_Euler, nt_CN, steps_convergence_analysis, L2_error_avg_list, F_num_komplex_n_Euler, t_num, F_num, runs_CN_high, fineness, tol_adaptiv_dt, safety, p, boundary_tol,
                    menu, distribution_menu, verbose, discretization_scheme, opt_success, t_ref, c_ref, D_ref, l_ref,
                    t_exp, F_exp, kind_of_1d_interpolation):

    if menu == 1: # S I M U L A T I O N

        # Mesh parameters
        mesh, V, ds, r, n_facet = collect_mesh_parameters(r_i, r_o, nx[0], degree[1], element_family, boundary_tol)

        # Check CFL-conditions
        dt[0] = check_CFL_conditions(dt[0], D, mesh, verbose, discretization_scheme, menu) # Set dt, that dt meets CFL-condition: dt<=dx^2/(2*D)

        # Run model
        _tnum, _Fnum = run_model(D, alpha, l_b, A_m_degrad, const_degrad, degrad_break,
                                 t, t_break, t_reached_5, T, T_explicit_Euler, dt[0], c_sat, c_s_t0, c_max, V_s, V_s_units, V_c_units, r_i, r_o, l_c, H, A_m_t0, m_d_t0, const_distribution, parabolic_distribution, offset_parabola, U_0, U_1, beta, gamma, limiter_break, limiter_sat,
                                 V, mesh, nx[0], nr[0], ds, degree[1], r, n_facet, num_steps, num_steps_explicit_Euler[0], F_num_komplex_n_Euler, runs_CN_high, fineness, tol_adaptiv_dt, safety, p, boundary_tol,
                                 menu, distribution_menu, verbose, discretization_scheme, t_ref, c_ref, D_ref, l_ref, t_num, F_num)


        return [_tnum, _Fnum, distribution_menu]



    elif menu == 2: # I N V E R S E   A N A L Y S I S

        # User chooses duration of optimization and the option to add fitted data points to experimental data
        parameters_to_optimize, parameters_treated_as_constant, T__opt_end, t_exp__opt, F_exp__opt, opt_time_frame = choose_optimization_time_frame(D, alpha, l_b, A_m_degrad, const_degrad, const_degrad_deactivated, degrad_break, kind_of_1d_interpolation, T, t_exp, F_exp)

        # Mesh parameters
        mesh, V, ds, r, n_facet = collect_mesh_parameters(r_i, r_o, nx[0], degree[1], element_family, boundary_tol)

        # Check CFL-conditions
        dt[0] = check_CFL_conditions(dt[0], D, mesh, verbose, discretization_scheme, menu) # Set dt, that dt meets CFL-condition: dt<=dx^2/(2*D)

        # Run inverse analysis (matrix break deactivated)
        _res_opt = run_optimization(parameters_to_optimize, parameters_treated_as_constant,
                                    t, t_break, t_reached_5, T__opt_end, T_explicit_Euler, dt[0], c_sat, c_s_t0, c_max,V_s, V_s_units, V_c_units, r_i, r_o, l_c, H, A_m_t0, m_d_t0, const_distribution, parabolic_distribution, offset_parabola, U_0, U_1, beta, gamma, limiter_break, limiter_sat,
                                    V, mesh, nx[0], nr[0], ds, degree[1], r, n_facet, num_steps, num_steps_explicit_Euler[0], F_num_komplex_n_Euler, runs_CN_high, fineness, tol_adaptiv_dt, safety, p, boundary_tol,
                                    menu, distribution_menu, verbose, discretization_scheme, t_ref, c_ref, D_ref, l_ref, t_num, F_num,
                                    t_exp__opt, F_exp__opt, opt_time_frame, opt_success, kind_of_1d_interpolation)

        if _res_opt.success == True:

            # Extract parameters
            if parameters_to_optimize == [D, alpha, l_b]: # Optimization for week 0-5

                D_opt     = _res_opt.x[0]
                alpha_opt = _res_opt.x[1]
                l_b_opt   = _res_opt.x[2]
                A_m_degrad_opt = A_m_degrad # A_m_degrad wasn't optimized

            elif parameters_to_optimize == [A_m_degrad]: # Optimization for week 5-15

                A_m_degrad_opt = _res_opt.x[0]
                D_opt = D         # D     wasn't optimized
                alpha_opt = alpha # alpha wasn't optimized
                l_b_opt = l_b     # l_b   wasn't optimized

        # Run model with optimized parameters
        _tnum, _Fnum = run_model(D_opt, alpha_opt, l_b_opt, A_m_degrad_opt, const_degrad, degrad_break,
                                 t, t_break, t_reached_5, T, T_explicit_Euler, dt[0], c_sat, c_s_t0, c_max, V_s, V_s_units, V_c_units, r_i, r_o, l_c, H, A_m_t0, m_d_t0, const_distribution, parabolic_distribution, offset_parabola, U_0, U_1, beta, gamma, limiter_break, limiter_sat,
                                 V, mesh, nx[0], nr[0], ds, degree[1], r, n_facet, num_steps, num_steps_explicit_Euler[0], F_num_komplex_n_Euler, runs_CN_high, fineness, tol_adaptiv_dt, safety, p, boundary_tol,
                                 menu, distribution_menu, verbose, discretization_scheme, t_ref, c_ref, D_ref, l_ref, t_num, F_num)



        return [_res_opt, _tnum, _Fnum, opt_time_frame, distribution_menu]



    elif menu == 4: # C O N V E R G E N C E   A N A L Y S I S

        _nx, _num_steps_explicit_Euler, \
        _L2_error_avg, _rate, _diff_F_num_total_avg, \
        _space_or_time, _discretization_scheme = run_convergence_analysis(D, alpha, l_b, A_m_degrad, const_degrad, degrad_break,
                                                                          t, t_break, t_reached_5, T, T_explicit_Euler, dt, dt_CN, c_sat, c_s_t0, c_max, V_s, V_s_units, V_c_units, r_i, r_o, l_c, H, A_m_t0, m_d_t0, const_distribution, parabolic_distribution, offset_parabola, U_0, U_1, beta, gamma, limiter_break, limiter_sat,
                                                                          degree, element_family, num_steps, nx, nr, h, steps_convergence_analysis, num_steps_explicit_Euler, nt_CN, L2_error_avg_list, F_num_komplex_n_Euler, runs_CN_high, fineness, tol_adaptiv_dt, safety, p, boundary_tol,
                                                                          menu, distribution_menu, verbose, t_ref, c_ref, D_ref, l_ref, t_num, F_num)


        return [_nx, _num_steps_explicit_Euler, _L2_error_avg, _rate, _diff_F_num_total_avg, _space_or_time, _discretization_scheme]



    elif menu == 5:  # T E S T I N G

        # Chose testing method
        menu_testing = choose_testing_method()

        # Run test
        _res_testing = run_testing(t, T, r_i, r_o, H, m_d_t0, l_ref, c_ref, menu, const_distribution, parabolic_distribution, U_0, U_1, offset_parabola, element_family, degree, nx, verbose, discretization_scheme, menu_testing, distribution_menu)

        return _res_testing



def main_menu():

    _menu = 0

    while _menu < 1 or _menu > 6:

        print('')
        print('  --------------------------------------------------------------')
        print('1. Simulation')
        print('2. Optimization')
        print('3' + '\u0336' + '.' + '\u0336' + ' ' + 'M' + '\u0336' + 'M' + '\u0336' + 'S' + '\u0336')
        print('4. Convergence analysis')
        print('5. Testing')
        print('')
        print('6. EXIT')
        print('  --------------------------------------------------------------')
        print('')

        _menu = int(input())

        if _menu < 1 or _menu > 6:
            print('You did not choose between "1-6"!')
            print('')


    if _menu == 6:
            sys.exit()


    return _menu



def u_t0_distribution_menu(D_constant, D_parabolic,
                           alpha_constant, alpha_parabolic,
                           l_b_constant, l_b_parabolic,
                           A_m_degrad_constant, A_m_degrad_parabolic,
                           const_degrad_constant, const_degrad_parabolic,
                           degrad_break_constant, degrad_break_parabolic,
                           menu):

    # Initialize distribution menu with dummy variable,
    # which is also used for convergence analysis, because there u(t=0) is given by the manufactured solution
    _distribution_menu = 2

    while _distribution_menu != 0 and _distribution_menu != 1:

        print('')
        print('  --------------------------------------------------------------')
        print('  G E N E R A L   O P T I O N S')
        print('  --------------------------------------------------------------')
        print('0. Constant  drug distribution [u(t=0)]')
        print('1. Parabolic drug distribution [u(t=0)]')
        print('  --------------------------------------------------------------')
        print('')

        _distribution_menu = int(input())

        if _distribution_menu != 0 and _distribution_menu != 1:
            print('You did not choose "1" or "0"!')
            print('')

    if _distribution_menu == 0:
        _D = D_constant
        _alpha = alpha_constant
        _l_b = l_b_constant
        _A_m_degrad = A_m_degrad_constant
        _const_degrad = const_degrad_constant
        _degrad_break = degrad_break_constant

    elif _distribution_menu == 1:
        _D = D_parabolic
        _alpha = alpha_parabolic
        _l_b = l_b_parabolic
        _A_m_degrad = A_m_degrad_parabolic
        _const_degrad = const_degrad_parabolic
        _degrad_break = degrad_break_parabolic


    return _D, _alpha, _l_b, _A_m_degrad, _const_degrad, _degrad_break, _distribution_menu



def verbose_mode():

    _verbose = 3 # initialize menu with dummy variables

    while _verbose != 0 and _verbose!= 1 and _verbose != 2:

        print('')
        print('  --------------------------------------------------------------')
        print('0. Verbose mode: SILENT')
        print('1. Verbose mode: BALANCED')
        print('2. Verbose mode: FULL')
        print('  --------------------------------------------------------------')
        print('')

        _verbose = int(input())

        if _verbose != 0 and _verbose!= 1 and _verbose != 2:
            print('You did not choose "0", "1" or "2"!')
            print('')


    return _verbose



def timer_menu():

    # Run menu
    _timer_menu = 2

    while _timer_menu != 0 and _timer_menu != 1:

        print('')
        print('  --------------------------------------------------------------')
        print('  Do you want to time the chosen calculation?')
        print('  --------------------------------------------------------------')
        print('0. No')
        print('1. Yes')
        print('  --------------------------------------------------------------')
        print('')


        _timer_menu = int(input())

        if _timer_menu != 0 and _timer_menu != 1:
                print('You did not choose "0" or "1"!')
                print('')

    if _timer_menu == 0:

        _timer_menu = 'timer_off'
        tic = 'dummy'

    elif _timer_menu == 1:

        _timer_menu = 'timer_on'
        tic = time.perf_counter()

    return _timer_menu, tic



def choose_optimization_time_frame(D, alpha, l_b, A_m_degrad, const_degrad, const_degrad_deactivated, degrad_break, kind_of_1d_interpolation, T, t_exp, F_exp):

    _menu_opt_time_frame = 0

    while _menu_opt_time_frame != 1 and _menu_opt_time_frame != 2:

        print('')
        print('    --------------------------------------------------------------')
        print('  1. Optimization from start (t=0) until matrix break (t=5)')
        print('  2. Optimization from matrix break (t=5) until end (t=15)')
        print('    --------------------------------------------------------------')
        print('')

        _menu_opt_time_frame = int(input())

        if _menu_opt_time_frame != 1 and _menu_opt_time_frame != 2:
            print('You did not choose between "1-2"!')
            print('')


    if _menu_opt_time_frame == 1:

        _opt_time_frame = 'until_break'
        _T__opt_end = 5.0
        _parameters_to_optimize = [D, alpha, l_b]
        _parameters_treated_as_constant = [A_m_degrad, const_degrad_deactivated, degrad_break]
        _t_exp__opt, _F_exp__opt = menu__add_20_fitted_points_to_experimenta_data(t_exp, F_exp, _T__opt_end, T, kind_of_1d_interpolation)

    elif _menu_opt_time_frame == 2:

        _opt_time_frame = 'after_break'
        _T__opt_end = T
        _parameters_to_optimize = [A_m_degrad]
        _parameters_treated_as_constant = [D, alpha, l_b, const_degrad, degrad_break]
        _t_exp__opt = t_exp
        _F_exp__opt = F_exp


    return _parameters_to_optimize, _parameters_treated_as_constant,_T__opt_end, _t_exp__opt, _F_exp__opt, _opt_time_frame



def menu__add_20_fitted_points_to_experimenta_data(t_exp, F_exp, T__opt_end, T, kind_of_1d_interpolation):


    # Adjust t_exp & F_exp to optimization time frame (t = 0-5)
    _t_exp__opt, _F_exp__opt = adjust_t_exp_and_F_exp_to_optimization_time_frame(t_exp, F_exp, T__opt_end)

    # Run menu
    _add_20_points_menu = 2

    while _add_20_points_menu != 0 and _add_20_points_menu != 1:

        print('')
        print('  --------------------------------------------------------------')
        print('  Do you want to add fitted points to the experimental data?')
        print('  --------------------------------------------------------------')
        print('0. No')
        print('1. Yes')
        print('  --------------------------------------------------------------')
        print('')

        _add_20_points_menu = int(input())

        if _add_20_points_menu != 0 and _add_20_points_menu!= 1:
                print('You did not choose "0" or "1"!')
                print('')

    if _add_20_points_menu == 0:

        return _t_exp__opt, _F_exp__opt



    elif _add_20_points_menu == 1:

        _add_20_points_menu__when = 0

        while _add_20_points_menu__when < 1 or _add_20_points_menu__when > 5:

            print('')
            print('    --------------------------------------------------------------')
            print('      When and how many?')
            print('    --------------------------------------------------------------')
            print('      1.  20 points @ t = 3.0 - 5.0')
            print('      2.  20 points @ t = 4.9 - 5.0')
            print('      3.  60 points @ t = 4.9 - 5.0')
            print('      4. 100 points @ t = 4.9 - 5.0')
            print('      5. 100 points @ t = 0.0 - 5.0')
            print('    --------------------------------------------------------------')
            print('')

            _add_20_points_menu__when = int(input())

            if _add_20_points_menu__when < 1 or _add_20_points_menu__when > 5:

                print('You did not choose between "1-5"!')
                print('')


        if _add_20_points_menu__when == 1:

            t_add_start = 3.0
            t_add_stop  = 5.0
            n_add = 20

        elif _add_20_points_menu__when == 2:

            t_add_start = 4.9
            t_add_stop  = 5.0
            n_add = 20

        elif _add_20_points_menu__when == 3:

            t_add_start = 4.9
            t_add_stop  = 5.0
            n_add = 60

        elif _add_20_points_menu__when == 4:

            t_add_start = 4.9
            t_add_stop  = 5.0
            n_add = 100

        elif _add_20_points_menu__when == 5:

            t_add_start = 0.0
            t_add_stop  = 5.0
            n_add = 100

        _t_exp_artificial, _F_exp_artificial = add_fitted_points_to_experimenta_data(t_add_start, t_add_stop, n_add, _t_exp__opt, _F_exp__opt, T__opt_end, T, kind_of_1d_interpolation)



        return _t_exp_artificial, _F_exp_artificial



def choose_testing_method():

    _menu_testing = 0

    while _menu_testing < 1 or _menu_testing > 2:

        print('')
        print('  --------------------------------------------------------------')
        print('1. Convergence analysis (TO BE IMPLEMENTED SOON)')
        print('2. Relative error of drug mass m_d depending on spatial discretization')
        print('  --------------------------------------------------------------')
        print('')

        _menu_testing = int(input())

        if _menu_testing < 1 or _menu_testing > 2:
            print('You did not choose between "1-2"!')
            print('')


    return _menu_testing
