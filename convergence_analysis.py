# -----------------------------------------------------------------
# C O N V E R G E N C E   A N A L Y S I S
# -----------------------------------------------------------------
from mesh_domains_boundaries import collect_mesh_parameters
from model import run_model
from utility import check_CFL_conditions, compute__convergence_rate
from save import save_error_to_csv, save_rate_to_csv


def run_convergence_analysis(D, alpha, l_b, A_m_degrad, const_degrad, degrad_break,
                             t, t_break, t_reached_5, T, T_explicit_Euler, dt, dt_CN, c_sat, c_s_t0, c_max, V_s, V_s_units, V_c_units, r_i, r_o, l_c, H, A_m_t0, m_d_t0, const_distribution, parabolic_distribution, offset_parabola, U_0, U_1, beta, gamma, limiter_break, limiter_sat,
                             degree, element_family, num_steps, nx, nr, h, steps_convergence_analysis, num_steps_explicit_Euler, nt_CN, L2_error_avg_list, F_num_komplex_n_Euler, runs_CN_high, fineness, tol_adaptiv_dt, safety, p, boundary_tol,
                             menu, distribution_menu, verboose, t_ref, c_ref, D_ref, l_ref, t_num, F_num):


    # User-menu
    space_or_time = menu_convergence_analysis__space_or_time()
    discretization_scheme = menu_convergence_analysis__discretization_scheme()

    # DUMMY parameter (used for collecting results of the 'discretization_scheme == 6/7'-case: )
    # TODO: ev. entfernen, falls case=6/7 nicht notwendig
    diff_F_num_total_avg = 0.0


    if discretization_scheme == 1: # explicit EULER

        nt = num_steps_explicit_Euler
        dt = dt
        T = T_explicit_Euler

    elif discretization_scheme == 2 or discretization_scheme == 3: # CRANK-NICOLSON (low and high precision)

        nt = nt_CN
        dt = dt_CN
        T_explicit_Euler = 0 # Choose T_explicit_Euler=0, that EULER loop will be skipped in model.py


    if space_or_time == 1: # spatial discretization

        # Set time step length
        nt = nt[0]
        dt = dt[0]

        # Run spatial convergence analysis
        for n in range(len(degree)):

            # Iterate through nx
            for i in range(len(nx)):

                print('')
                print(element_family, degree[n], ': nx =', nx[i])

                # Mesh parameters
                mesh, V, ds, r, n_facet = collect_mesh_parameters(r_i, r_o, nx[i], degree[n], element_family, boundary_tol)

                # Check dt for CFL-conditions (only for EULER's method)
                dt = check_CFL_conditions(dt, D, mesh, verboose, discretization_scheme, menu)

                # Compute average L2-error by running model
                L2_error_avg = run_model(D, alpha, l_b, A_m_degrad, const_degrad, degrad_break,
                                         t, t_break, t_reached_5, T, T_explicit_Euler, dt, c_sat, c_s_t0, c_max, V_s, V_s_units, V_c_units, r_i, r_o, l_c, H, A_m_t0, m_d_t0, const_distribution, parabolic_distribution, offset_parabola, U_0, U_1, beta, gamma, limiter_break, limiter_sat,
                                         V, mesh, nx[i], nr[i], ds, degree[n], r, n_facet, num_steps, nt, F_num_komplex_n_Euler, runs_CN_high, fineness, tol_adaptiv_dt, safety, p, boundary_tol,
                                         menu, distribution_menu, verboose, discretization_scheme, t_ref, c_ref, D_ref, l_ref, t_num, F_num)

                # Add L2_error to list
                L2_error_avg_list[n] += [L2_error_avg]


        # Compute and save convergence rate
        rate = compute__convergence_rate(L2_error_avg_list, degree, nx, nt, h, space_or_time)
        save_rate_to_csv(nx, nt, dt, rate, D, D_ref, T, space_or_time, discretization_scheme)

        # Save average L2-error for every spatial discretization
        save_error_to_csv(nx, nt, dt, L2_error_avg_list, D, D_ref, T, T_explicit_Euler, space_or_time, discretization_scheme)



    elif space_or_time == 2: # temporal discretization

        # Set amount of elements
        nx = nx[-1]
        nr = nr[-1] # only used for plots

        # Run temporal convergence analysis
        for n in range(len(degree)):

            # Iterate through dt
            for i in range(len(nt)):
                print('')
                print(element_family, degree[n], ': nt =', nt[i])

                # Mesh parameters
                mesh, V, ds, r, n_facet = collect_mesh_parameters(r_i, r_o, nx, degree[n], element_family, boundary_tol)

                # Check dt for CFL-conditions (only for EULER's method)
                dt[i] = check_CFL_conditions(dt[i], D, mesh, verboose, discretization_scheme, menu)

                # Compute average L2-error by running model
                L2_error_avg = run_model(D, alpha, l_b, A_m_degrad, const_degrad, degrad_break,
                                         t, t_break, t_reached_5, T, T_explicit_Euler, dt[i], c_sat, c_s_t0, c_max, V_s, V_s_units, V_c_units, r_i, r_o, l_c, H, A_m_t0, m_d_t0, const_distribution, parabolic_distribution, offset_parabola, U_0, U_1, beta, gamma, limiter_break, limiter_sat,
                                         V, mesh, nx, nr, ds, degree[n], r, n_facet, num_steps, nt[i], F_num_komplex_n_Euler, runs_CN_high, fineness, tol_adaptiv_dt, safety, p, boundary_tol,
                                         menu, distribution_menu, verboose, discretization_scheme, t_ref, c_ref, D_ref, l_ref, t_num, F_num)

                # Add L2_error to list
                L2_error_avg_list[n] += [L2_error_avg]


        # Compute and save convergence rate
        rate = compute__convergence_rate(L2_error_avg_list, degree, nx, nt, h, space_or_time)
        save_rate_to_csv(nx, nt, dt, rate, D, D_ref, T, space_or_time, discretization_scheme)

        # Save average L2-error for every temporal discretization
        save_error_to_csv(nx, nt, dt, L2_error_avg_list, D, D_ref, T, T_explicit_Euler, space_or_time, discretization_scheme)


    return nx, nt, L2_error_avg_list, rate, diff_F_num_total_avg, space_or_time, discretization_scheme



def menu_convergence_analysis__space_or_time():

    _space_or_time = 0 # initialize menu with dummy value

    while _space_or_time  != 1 and _space_or_time  != 2:

        print('')
        print('     ---------------------------------------------------------------------')
        print('      C O N V E R G E N C E   A N A L Y S I S')
        print('')
        print('   1. in SPACE')
        print('   2. in TIME')
        print('     ---------------------------------------------------------------------')
        print('')

        _space_or_time = int(input())

        if _space_or_time  != 1 and _space_or_time  != 2:
            print('   You did not choose "1" or "2"!')
            print('')


    return _space_or_time



def menu_convergence_analysis__discretization_scheme():

    _discretization_scheme = 0 # initialize menu with dummy value

    while _discretization_scheme < 1 or _discretization_scheme > 3:

        print('')
        print('        ---------------------------------------------------------------------')
        print('         D I S C R E T I Z A T I O N   S C H E M E   to test')
        print('')
        print('      1. MMS - EULER')
        print('      2. MMS - CRANK-NICOLSON (low  precision loop)')
        print('      3. MMS - CRANK-NICOLSON (high precision loop)')
        print('        ---------------------------------------------------------------------')
        print('')

        _discretization_scheme = int(input())

        if _discretization_scheme < 1 or _discretization_scheme > 3:

            print('      You did not choose between "1" - "3"!')
            print('')


    return _discretization_scheme