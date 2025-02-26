# -----------------------------------------------------------------
#  M A I N
# -----------------------------------------------------------------
from parameters import *
from menu import main_menu, u_t0_distribution_menu, timer_menu, verbose_mode, compute_results
from utility import stop_timer
from plot_funcs import print_plot_all_results
import matplotlib.pyplot as plt # version: 2.1.1


# -----------------------------------------------------------------
#   main M E N U
# -----------------------------------------------------------------
menu = main_menu()

# -----------------------------------------------------------------
#   D R U G   D I S T R I B U T I O N   u(t=0)
# -----------------------------------------------------------------
D, alpha, l_b, A_m_degrad, const_degrad, degrad_break, distribution_menu = u_t0_distribution_menu(D_constant, D_parabolic,
                                                                                                  alpha_constant, alpha_parabolic,
                                                                                                  l_b_constant, l_b_parabolic,
                                                                                                  A_m_degrad_constant, A_m_degrad_parabolic,
                                                                                                  const_degrad_constant, const_degrad_parabolic,
                                                                                                  degrad_break_constant, degrad_break_parabolic,
                                                                                                  menu)

# -----------------------------------------------------------------
#   toggle  V E R B O S E  mode
# -----------------------------------------------------------------
verbose = verbose_mode()

# -----------------------------------------------------------------
#   timer M E N U
# -----------------------------------------------------------------
timer, tic = timer_menu()

# -----------------------------------------------------------------
#   C O M P U T E   R E S U L T S   according to user choice
# -----------------------------------------------------------------
results = compute_results(D, alpha, l_b, A_m_degrad, const_degrad, degrad_break,
                          t, t_break, t_reached_5, T, T_explicit_Euler, dt, dt_CN, c_sat, c_s_t0, c_max, V_s, V_s_units, V_c_units, r_i, r_o, l_c, H, A_m_t0, m_d_t0, const_distribution, parabolic_distribution, offset_parabola, U_0, U_1, beta, gamma, const_degrad_deactivated, limiter_break, limiter_sat,
                          nx, nr, h, degree, element_family, num_steps, num_steps_explicit_Euler, nt_CN, steps_convergence_analysis, L2_error_avg_list, F_num_komplex_n_Euler, t_num, F_num, runs_CN_high, fineness, tol_adaptiv_dt, safety, p, boundary_tol,
                          menu, distribution_menu, verbose, discretization_scheme, opt_success, t_ref, c_ref, D_ref, l_ref,
                          t_exp, F_exp, kind_of_1d_interpolation)

stop_timer(tic, timer)

# -----------------------------------------------------------------
#   P R I N T   R E S U L T S
# -----------------------------------------------------------------
print_plot_all_results(results, menu, t_exp, F_exp, T, kind_of_1d_interpolation, element_family, degree, nx, h, D, alpha, l_b, const_degrad, degrad_break, A_m_t0, A_m_t0_units, A_m_degrad_constant_units, A_m_degrad_parabolic_units, opt_success, c_ref, D_ref, l_ref)
plt.show() # show plots
