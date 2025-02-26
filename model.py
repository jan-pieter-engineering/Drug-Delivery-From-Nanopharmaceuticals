"""
Diffusion equation with Robin boundary conditions.
Test problem is chosen to give an exact solution at all nodes of the mesh.

  u'= D*Laplace(u)    on Interval mesh
"""
from fenics import *
from u_0 import define_u_t0, define_u_e
from boundary_conditions import h_modified
from degradation import degradation_model
from source_term import define_source_term
from surrounding_solution import calc_c_s, update_c_s_n, transfer_c_s_n_from_EULER_to_CN
from solver_Euler import solver_Euler
from solver_CN_low import solver_CN_low
from solver_CN_high import solver_CN_high
from utility import check__u_t0, choose_F_num, calc_F, adjust_dt_to_dt_fine, update_t_to_current_time
from save import collect_results, save_t_num_F_num, calc_plot_print_L2_error, create__save_files, save_solution
from plot_funcs import plot_results, show_plot, plot__h_mod, print__t_dt_num_steps
from plot_settings_inits import plot_settings

#set_log_active(False)         # mute console output
set_log_level(30)              # adjust level of console output (20 := solving linear var problem)


def run_model(D, alpha, l_b, A_m_degrad, const_degrad, degrad_break,
              t, t_break, t_reached_5, T, T_explicit_Euler, dt, c_sat, c_s_t0, c_max, V_s, V_s_units, V_c_units, r_i, r_o, l_c, H, A_m_t0, m_d_t0, const_distribution, parabolic_distribution, offset_parabola, U_0, U_1, beta, gamma, limiter_break, limiter_sat,
              V, mesh, nx, nr, ds, degree, r, n_facet, num_steps, num_steps_explicit_Euler, F_num_komplex_n_Euler, runs_CN_high, fineness, tol_adaptiv_dt, safety, p, boundary_tol,
              menu, distribution_menu, verbose, discretization_scheme, t_ref, c_ref, D_ref, l_ref, t_num, F_num):

    # Misc
    vtk_file, xdmf_file = create__save_files(menu, distribution_menu)                                         # create save files according to menu/verbose-mode user-choice
    ax1, ax2, ax3, ax4, ax5, ax6 = plot_settings(c_max, c_sat, m_d_t0, menu, discretization_scheme, verbose)  # set plot settings according to menu/verbose-mode user-choice
    L2_error = []                                                                                             # initialize list, that will save L2 errornorm for every time step, if convergence analysis is performed
    F_num_choice = choose_F_num(menu) # TODO: besser implementieren

    # Define t^0 values
    u_0  = define_u_t0(t, T, r_i, r_o, distribution_menu, const_distribution, parabolic_distribution, U_0, U_1, offset_parabola, degree, discretization_scheme, menu) # define u(t=0)
    u_t0 = interpolate(u_0, V)                                                     # interpolate u_0
    check__u_t0(u_t0, m_d_t0, r, H, l_ref, c_ref, verbose, menu) # check, if distribution-value is chosen correctly, that m_d_units = 7.5g

    # Define t^n values
    u_n        = interpolate(u_0, V)                      # u^n := u(t=0)

    # Define t^(n-1) values [DUMMY variables]
    u_n_1      = interpolate(u_0, V)                      # DUMMY variable
    c_s_n_1 = c_s_t0                                      # DUMMY variable
    h_mod_n_1  = D                                        # DUMMY variable
    F_num_komplex_n_1_Euler = 0                           # DUMMY variable

    # Define manufactured solution
    u_e      = define_u_e(t, T, U_0, U_1, degree+3+1)     # manufactured solution/exact solution; degree=element_degree+degree_rise+1 with degree_rise=3 (used for errornorm)


    # --------------------------------------------------------------
    #   Explicit   E U L E R
    # --------------------------------------------------------------
    # Initialize t^n Euler-values
    c_s_n_Euler = calc_c_s(t, u_n_1, c_s_n_1, c_s_t0, h_mod_n_1, T, dt, ds, H, gamma, U_0, U_1, r, r_i, r_o, degree, boundary_tol, discretization_scheme, menu)              # Konzentration der Umgebungslösung
    degrad_max_n_Euler = degradation_model(t, u_t0, u_n_1, const_degrad, r_i, l_c, nx, nr, verbose, ax2[0], menu)                                                            # Matrixdegradation
    h_mod_n_Euler = h_modified(t, D, alpha, l_b, c_s_n_Euler, c_sat, c_max, r, A_m_t0, A_m_degrad, degrad_max_n_Euler, degrad_break, U_1, T, discretization_scheme, degree, menu) # mod. Stoffübergangsparameter h'
    f_n_Euler = define_source_term(t, T, D, U_0, U_1, beta, degree, discretization_scheme, menu)                                                                             # source term
    # TODO: t=0 case implementieren
    F_num_simpel_n_Euler, F_num_komplex_n_Euler = calc_F(dt, D, u_n_1, c_s_n_Euler, F_num_komplex_n_1_Euler, c_max, V_s_units, V_c_units, H, beta, gamma, r, n_facet, ds, discretization_scheme, menu)

    # EULER time stepping
    while t<T_explicit_Euler-dt+(dt/100): # tol_dt = dt/100

        # Update num step counter
        num_steps += 1


        # Compute solution
        u_Euler = solver_Euler(t, dt, D, u_n, c_s_n_Euler, h_mod_n_Euler, f_n_Euler, beta, V, r, ds,
                               degree, linear_solver='direct')


        # Update t to current time
        t += dt

        # Update analytical solution
        u_e.t = t

        # Calc & Plot/Print L2-error (only during convergence analysis)
        L2_error = calc_plot_print_L2_error(t, T_explicit_Euler, u_e, u_Euler, L2_error, num_steps, num_steps_explicit_Euler, ax5, verbose, discretization_scheme)

        # Update t^n values for next timestep
        # TODO: muss hier nicht die update_c_s_n-function genutzt werden?
        c_s_n_Euler = calc_c_s(t, u_n, c_s_n_Euler, c_s_t0, h_mod_n_Euler, T, dt, ds, H, gamma, U_0, U_1, r, r_i, r_o, degree, boundary_tol, discretization_scheme, menu)
        degrad_max_n_Euler = degradation_model(t, u_t0, u_n, const_degrad, r_i, l_c, nx, nr, verbose, ax2[0], menu)
        h_mod_n_Euler = h_modified(t, D, alpha, l_b, c_s_n_Euler, c_sat, c_max, r, A_m_t0, A_m_degrad, degrad_max_n_Euler, degrad_break, U_1, T, discretization_scheme, degree, menu)
        F_num_simpel_n_Euler, F_num_komplex_n_Euler = calc_F(dt, D, u_n, c_s_n_Euler, F_num_komplex_n_Euler, c_max, V_s_units, V_c_units, H, beta, gamma, r, n_facet, ds, discretization_scheme, menu)
        print__t_dt_num_steps(t, dt, num_steps, T_explicit_Euler, menu, verbose)
        f_n_Euler.t=t
        u_n.assign(u_Euler)


        # Save
        save_solution(t, u_n, vtk_file, xdmf_file, menu)  # save to vtk- and XDMF-file

        # Plot
        plot_results(t, T_explicit_Euler, m_d_t0, V_s, H, r, ax1, ax2, ax3, ax6, verbose, discretization_scheme, menu,
                     u_n, c_s_n_Euler, degrad_max_n_Euler, F_num_simpel_n_Euler, F_num_komplex_n_Euler,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0)



    # --------------------------------------------------------------
    #   Adaptive C R A N K - N I C O L S O N - time-stepping
    # --------------------------------------------------------------
    # Initialize t^n CN-values by passing Euler values
    c_s_n_CN_low  = transfer_c_s_n_from_EULER_to_CN(c_s_n_Euler, degree, menu)
    c_s_n_CN_high = transfer_c_s_n_from_EULER_to_CN(c_s_n_Euler, degree, menu)
    h_mod_n_CN_low  = h_mod_n_Euler
    h_mod_n_CN_high = h_mod_n_Euler
    f_n_CN_low  = define_source_term(t, T, D, U_0, U_1, beta, degree, discretization_scheme, menu)
    f_n_CN_high = define_source_term(t, T, D, U_0, U_1, beta, degree, discretization_scheme, menu)
    F_num_komplex_n_CN_low  = F_num_komplex_n_Euler
    F_num_komplex_n_CN_high = F_num_komplex_n_Euler


    # CN time stepping
    num_steps_real = num_steps  # only time steps, that conform to Richardson interpolation
    while t>=T_explicit_Euler-(dt/100) and t<T: # tol_dt = dt/100

        # Update total step counter
        num_steps += 1

        # Adjust dt for special cases (t=15(end) and t=5(matrix break))
        if t+dt>T: # adjust last time-step, that t_num[last]=T
            dt=T-t
        elif t + dt > t_break and t_reached_5 == False and menu != 4: # adjust dt, that one exactly hits t=5 (not relevant for convergence analysis)
            dt = t_break - t
            t_reached_5 = True



        # Compute L O W precision solution
        u_CN_low, c_s_CN_low, degrad_max_CN_low, h_mod_CN_low, F_num_simpel_CN_low, F_num_komplex_CN_low = solver_CN_low(t, dt, u_n, c_s_n_CN_low, h_mod_n_CN_low, f_n_CN_low, F_num_komplex_n_CN_low,
                                                                                                                         u_t0, c_sat, c_max, c_s_t0, V_s_units, V_c_units, r_i, r_o, H, l_c, A_m_t0, A_m_degrad, beta, gamma, U_0, U_1, T, const_degrad, degrad_break, degree, nx, nr, V, n_facet, r, ds, boundary_tol, verbose, ax2[0], discretization_scheme, menu,
                                                                                                                         D, alpha, l_b,
                                                                                                                         degree, linear_solver='direct')
        # Adjust dt to dt_fine
        dt = adjust_dt_to_dt_fine(dt, fineness, discretization_scheme) # dt = dt/fineness [exception: convergence analysis for CN, then dt=dt]

        # Compute H I G H precision solution
        u_CN_high, c_s_CN_high, degrad_max_CN_high, h_mod_CN_high, F_num_simpel_CN_high, F_num_komplex_CN_high = solver_CN_high(t, dt, u_n, c_s_n_CN_high, h_mod_n_CN_high, f_n_CN_high, F_num_komplex_n_CN_high,
                                                                                                                                u_t0, c_sat, c_max, c_s_t0, V_s_units, V_c_units, r_i, r_o, H, l_c, A_m_t0, A_m_degrad, beta, gamma, U_0, U_1, T, const_degrad, degrad_break, runs_CN_high, degree, nx, nr, V, n_facet, r, ds, boundary_tol, verbose, ax2[0], discretization_scheme, menu,
                                                                                                                                D, alpha, l_b,
                                                                                                                                degree, linear_solver='direct')

        # Update t to current time
        t = update_t_to_current_time(t, dt, fineness, discretization_scheme) # t += fineness * dt [exception: convergence analysis for CN, then t += dt]


        # Adaptive CN-time-stepping
        richardson = sqrt(assemble(((u_CN_high - u_CN_low)**2)*dx)) / (2**p -1) # compute Richardson extrapolation

        if richardson < tol_adaptiv_dt and (discretization_scheme != 2 and discretization_scheme != 3):

            # Update step counter, if Richardson extrapolation is passed
            num_steps_real += 1

            # Update dt according to Richardson extrapolation
            dt = ((safety * tol_adaptiv_dt / richardson) ** (1 / p)) * dt

            # dt limiter
            if dt > 0.125:  # if dt grows too fast, reduce dt to dt=0.25
                dt = 0.125

            # Print dt
            print__t_dt_num_steps(t, dt, num_steps_real, T_explicit_Euler, menu, verbose)

            # Update analytical solution
            u_e.t = t

            # Reset dt to dt_coarse
            dt = fineness * dt


            # Update t^n values for next timestep (CN_low values are updated with CN_high results, because dt conformed to Richardson interpolation)
            c_s_n_CN_low  = update_c_s_n(t, c_s_CN_high, degree, menu)
            c_s_n_CN_high = update_c_s_n(t, c_s_CN_high, degree, menu)
            h_mod_n_CN_low  = h_mod_CN_high
            h_mod_n_CN_high = h_mod_CN_high
            f_n_CN_low.t=t
            f_n_CN_high.t=t
            F_num_komplex_n_CN_low  = F_num_komplex_CN_high
            F_num_komplex_n_CN_high = F_num_komplex_CN_high
            u_n.assign(u_CN_high)


            # Save
            save_solution(t, u_n, vtk_file, xdmf_file, menu)  # save to vtk- and XDMF-file
            t_num, F_num = save_t_num_F_num(t, t_num, F_num, F_num_simpel_CN_high, F_num_komplex_CN_high, F_num_choice)

            # Plot
            limiter_break, limiter_sat = plot__h_mod(t, h_mod_CN_high, degrad_max_CN_high, degrad_break, alpha, c_sat, c_s_CN_high, limiter_break, limiter_sat, ax4, verbose, menu, discretization_scheme)
            plot_results(t, T_explicit_Euler, m_d_t0, V_s, H, r, ax1, ax2, ax3, ax6, verbose, discretization_scheme, menu,
                        0, 0, 0, 0, 0,
                         u_CN_low, c_s_CN_low, degrad_max_CN_low, F_num_simpel_CN_low, F_num_komplex_CN_low,
                         u_CN_high, c_s_CN_high, degrad_max_CN_high, F_num_simpel_CN_high, F_num_komplex_CN_high)



        elif richardson > tol_adaptiv_dt and (discretization_scheme != 2 and discretization_scheme != 3):

            t -= fineness*dt # t wieder zurücksetzen, um wieder einen neuen Berechnungsschritt für selben Zeitraum durchzuführen, aber diesmal mit feinerem dt
            dt = dt/10       # adjust dt to higher precision/smaller timestep



        elif discretization_scheme == 2: # only active during convergence analysis (CN low precision)

            # Update
            num_steps_real += 1 # step counter
            u_e.t = t           # analytical solution

            # Calc & plot/print L2-errornorm (convergence analysis)
            L2_error = calc_plot_print_L2_error(t, T_explicit_Euler, u_e, u_CN_low, L2_error, num_steps_real, num_steps_explicit_Euler, ax5, verbose, discretization_scheme)

            # Update CN low variables for next timestep
            c_s_n_CN_low = update_c_s_n(t, c_s_CN_low, degree, menu)
            h_mod_n_CN_low  = h_mod_CN_low
            f_n_CN_low.t = t
            u_n.assign(u_CN_low)



        elif discretization_scheme == 3: # only active during convergence analysis (CN high precision)

            # Update
            num_steps_real += 1 # step counter
            u_e.t = t           # analytical solution

            # Calc & plot/print L2-errornorm (convergence analysis)
            L2_error = calc_plot_print_L2_error(t, T_explicit_Euler, u_e, u_CN_high, L2_error, num_steps_real, num_steps_explicit_Euler, ax5, verbose, discretization_scheme)

            # Update CN high variables for next timestep
            c_s_n_CN_high = update_c_s_n(t, c_s_CN_high, degree, menu)
            h_mod_n_CN_high = h_mod_CN_high
            f_n_CN_high.t = t
            u_n.assign(u_CN_high)

            # Reset dt to dt_coarse
            dt = fineness * dt



    
    show_plot(menu) # display final plot per optimization step for 5 sec
    results = collect_results(t_num, F_num, L2_error, menu, distribution_menu)  # collect results depending on menu choice

    return results