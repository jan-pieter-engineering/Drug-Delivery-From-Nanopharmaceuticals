# -----------------------------------------------------------------
# P L O T
# -----------------------------------------------------------------
from fenics import *
from utility import calc_F_num_at_discrete_time_steps, adjust_t_exp_and_F_exp_to_optimization_time_frame
import matplotlib.pyplot as plt # version: 2.1.1
from scipy.interpolate import interp1d
import numpy as np
from sklearn.metrics import r2_score
import datetime
import os
import csv


def plot_results(t, T_explicit_Euler, m_d_t0, V_s, H, r, ax1, ax2, ax3, ax6, verbose, discretization_scheme, menu,
                 u_n, c_s_n_Euler, degrad_max_n_Euler, F_num_simpel_n_Euler, F_num_komplex_n_Euler,
                 u_CN_low, c_s_CN_low, degrad_max_CN_low, F_num_simpel_CN_low, F_num_komplex_CN_low,
                 u_CN_high, c_s_CN_high, degrad_max_CN_high, F_num_simpel_CN_high, F_num_komplex_CN_high):

    if (menu == 1 or menu == 2):

        if verbose == 1:

            if t < T_explicit_Euler:

                plot_EULER__F_num__simpel_vs_komplex(ax1, F_num_simpel_n_Euler, F_num_komplex_n_Euler, t)

            elif t > T_explicit_Euler:

                plot_CN__F_num__simpel_vs_komplex(t, F_num_simpel_CN_low, F_num_simpel_CN_high, F_num_komplex_CN_low, F_num_komplex_CN_high, ax1, verbose, menu, discretization_scheme)


        elif verbose == 2:

            if t < T_explicit_Euler:

                plot_EULER__F_num__simpel_vs_komplex(ax1, F_num_simpel_n_Euler, F_num_komplex_n_Euler, t)
                plot_EULER__degrad_max(ax2, t, degrad_max_n_Euler)
                plot_EULER__m_d(ax6, t, u_n, c_s_n_Euler, m_d_t0, V_s, H, r, dx)

            elif t > T_explicit_Euler:

                plot_CN__F_num__simpel_vs_komplex(t, F_num_simpel_CN_low, F_num_simpel_CN_high, F_num_komplex_CN_low, F_num_komplex_CN_high, ax1, verbose, menu, discretization_scheme)
                plot_CN__degrad_max(ax2, t, degrad_max_CN_low, degrad_max_CN_high)
                plot_CN__c_s(t, c_s_CN_high, ax3, verbose, menu, discretization_scheme)
                plot_CN__m_d(ax6, u_CN_low, u_CN_high, c_s_CN_low, c_s_CN_high, m_d_t0, V_s, H, r, dx, t, menu, discretization_scheme)



def plot_EULER__F_num__simpel_vs_komplex(ax1, F_num_simpel_Euler, F_num_komplex_Euler, t):

    # plot F_num
    ax1[0].scatter(t, F_num_simpel_Euler,                color='blue',   marker='o')
    ax1[1].scatter(t, F_num_komplex_Euler,               color='orange', marker='o')

    # irritierende Achsenskalierung fixen
    ax1[0].get_yaxis().get_major_formatter().set_useOffset(False)

    # Exponent der wissenschaftlichen Darstellung verschieben, um Überschneidung zu vermeiden
    ax1[0].get_yaxis().get_offset_text().set_x(-0.25)
    ax1[1].get_yaxis().get_offset_text().set_x(-0.25)



def plot_CN__F_num__simpel_vs_komplex(t, F_num_simpel_CN_low, F_num_simpel_CN_high, F_num_komplex_CN_low, F_num_komplex_CN_high, ax1, verbose, menu, discretization_scheme):

    if ((menu==1 and (verbose == 1 or verbose==2)) or (menu==2 and verbose==2)) and discretization_scheme != 1:

        ax1[0].scatter(t, F_num_simpel_CN_high,                color='purple',   marker='x')
        ax1[1].scatter(t, F_num_komplex_CN_high,               color='red', marker='x')

        if menu==1:
            plt.pause(0.005)



def plot_degrad(_degrad, nr, r_i, l_c, verbose, ax2, menu):
    """ Plots degradation of the matrix. """

    #if menu == 1 or menu == 2 and verbose == 2:

    r_plot = np.linspace(r_i, l_c+r_i, nr+1) # l_c = r_o-r_i
    ax2.plot(r_plot, _degrad)



def plot_EULER__m_d(ax6, t, u_Euler, c_s_Euler, m_d_t0, V_s, H, r, dx):

    # Calc
    m_d_s = c_s_Euler.val * V_s
    m_d_m = assemble(2.0 * np.pi * H * u_Euler * r * dx)
    m_d_total = m_d_s + m_d_m
    delta_m_d = (m_d_total - m_d_t0) / m_d_t0
    delta_m_d__prozent = delta_m_d * 100.0

    # plot
    ax6[0,0].plot(t, m_d_total, color='purple', marker='x')
    ax6[1,0].plot(t, delta_m_d__prozent, color='purple', marker='x')

    # irritierende Achsenskalierung fixen
    ax6[1,0].get_yaxis().get_major_formatter().set_useOffset(False)

    # Exponent der wissenschaftlichen Darstellung verschieben, um Überschneidung zu vermeiden
    ax6[0,0].get_yaxis().get_offset_text().set_x(0.0)
    ax6[0,0].get_xaxis().get_offset_text().set_x(0.9)



def plot_CN__m_d(ax6, u_CN_low, u_CN_high, c_s_CN_low, c_s_CN_high, m_d_t0, V_s, H, r, dx, t, menu, discretization_scheme):

    if discretization_scheme != 1:

        m_d_s_CN_low =  c_s_CN_low.val  * V_s
        m_d_s_CN_high = c_s_CN_high.val * V_s
        m_d_m_CN_low  = assemble(2.0 * np.pi * H * u_CN_low  * r * dx)
        m_d_m_CN_high = assemble(2.0 * np.pi * H * u_CN_high * r * dx)
        m_d_total_CN_low  = m_d_s_CN_low  + m_d_m_CN_low
        m_d_total_CN_high = m_d_s_CN_high + m_d_m_CN_high

        delta_m_d_CN_low  = (m_d_t0 - m_d_total_CN_low)  / m_d_t0
        delta_m_d_CN_high = (m_d_t0 - m_d_total_CN_high) / m_d_t0
        delta_m_d_CN_low__prozent  = delta_m_d_CN_low  * 100.0
        delta_m_d_CN_high__prozent = delta_m_d_CN_high * 100.0

        # plot m_d_total
        ax6[0,1].plot(t, m_d_total_CN_low,  color='purple', marker='x')
        ax6[0,2].plot(t, m_d_total_CN_high, color='purple', marker='x')
        ax6[1,1].plot(t, delta_m_d_CN_low__prozent, color='purple', marker='x')
        ax6[1,2].plot(t, delta_m_d_CN_high__prozent, color='purple', marker='x')

        # irritierende Achsenskalierung fixen
        ax6[1,1].get_yaxis().get_major_formatter().set_useOffset(False)
        ax6[1,2].get_yaxis().get_major_formatter().set_useOffset(False)

        if menu==1:
            plt.pause(0.005)



def plot_EULER__degrad_max(ax2, t, degrad_max_n_Euler):

    #plot
    ax2[1].plot(t, degrad_max_n_Euler, color='purple', marker='x')

    # irritierende Achsenskalierung fixen
    ax2[1].get_yaxis().get_major_formatter().set_useOffset(False)

    # Exponent der wissenschaftlichen Darstellung verschieben, um Überschneidung zu vermeiden
    ax2[1].get_yaxis().get_offset_text().set_x(0.0)
    ax2[1].get_xaxis().get_offset_text().set_x(0.9)



def plot_CN__degrad_max(ax2, t, degrad_max_CN_low, degrad_max_CN_high):

    #plot
    ax2[2].plot(t, degrad_max_CN_low, color='purple', marker='x')
    ax2[3].plot(t, degrad_max_CN_high, color='purple', marker='x')



def plot__F_exp_vs_F_num(parameters_to_optimize, D_ref, c_ref, l_ref, t_exp, tnum, F_exp, Fnum, T, R_squared_0_T, kind_of_1d_interpolation, verbose, menu, opt_success, distribution_menu, opt_time_frame):

    # Extract parameters
    if menu == 1:
        [D, alpha, l_b, A_m_degrad_units, A_m_t0_units, const_degrad, degrad_break] = parameters_to_optimize

    elif menu == 2:
        if opt_time_frame == 'until_break':
            [D_opt, alpha_opt, l_b_opt] = parameters_to_optimize

        elif opt_time_frame == 'after_break':
            A_m_degrad_opt = parameters_to_optimize

    # Initialize Plot
    fig, ax = plt.subplots()
    ax.set_xlabel(r"$\tilde{t}$", loc="right")
    ax.set_ylabel(r"$\tilde{F}(\tilde{t})$", loc="top", rotation="horizontal")
    ax.yaxis.set_label_coords(-0.05, 0.9)

    # Set title
    if menu == 1:
        plt.title(f"D = {D*D_ref:e} [mm²/s] , alpha = {alpha/c_ref:e} [mm³/g] , l_b = {l_b*l_ref:e} [mm], A_m_degrad = {A_m_degrad_units/A_m_t0_units}*A_m_t0, const_degrad = {const_degrad}, degrad_break = {degrad_break}", fontsize=6)

    elif menu == 2:
        if opt_time_frame == 'until_break':
            plt.title(f"D = {D_opt*D_ref:e} [mm²/s] , alpha = {alpha_opt/c_ref:e} [mm³/g] , l_b = {l_b_opt*l_ref:e} [mm]", fontsize=6)

        elif opt_time_frame == 'after_break':
            plt.title(f"A_m_degrad = {A_m_degrad_opt * l_ref**2} [mm²]", fontsize=6)

    if opt_success == True:
        plt.suptitle(f"$O P T I M I Z A T I O N   S U C C E S S F U L: \quad F\_num \; vs. \; F\_exp: \quad R²(0,{int(T)}) = {R_squared_0_T}$")
    else:
        plt.suptitle(f"$F\_num \; vs. \; F\_exp: \quad R²(0,{int(T)}) = {R_squared_0_T}$")

    # Set the range of axis
    plt.ylim(0, 1)

    # F_exp interpolieren
    _F_exp_interp1d = interp1d(t_exp, F_exp, kind_of_1d_interpolation)
    _t_plot = np.linspace(0, T, 100)
    _F_exp_interpolated = _F_exp_interp1d(_t_plot)

    # F_exp und F_exp interpoliert plotten
    ax.plot(_t_plot, _F_exp_interpolated, color='black', label=f'F_exp ({kind_of_1d_interpolation} interpol)')
    ax.scatter(t_exp, F_exp, color='black', marker='x', label='F_exp')

    # F_num interpolieren
    _F_num_interp1d = interp1d(tnum, Fnum, kind_of_1d_interpolation)
    _F_num_interpolated = _F_num_interp1d(_t_plot)

    # F_num plotten
    ax.plot(_t_plot, _F_num_interpolated, color='blue', label=f'F_num ({kind_of_1d_interpolation} interpol)')
    #ax.plot(tnum, Fnum, 'o')

    # Display legends
    ax.legend()

    # Save plot
    if menu == 2:

        # Set time stamp
        date = str(datetime.date.today())
        date_plus_time = str(datetime.datetime.now())

        # Check distribution
        if distribution_menu == 0:
            start_distribution = "constant_c0"
        elif distribution_menu == 1:
            start_distribution = "parabolic_c0"

        # Create path
        path_2 = str("results/02__optimization/" + date + "/" + start_distribution)

        # Check time frame, which is optimized
        if opt_time_frame == 'until_break':

            # Create folder, if necessary
            if not os.path.exists(path_2 + "/until_break"):
                os.makedirs(path_2 + "/until_break")

            # Save figure to .png
            plt.savefig(path_2 + "/until_break" + "/" + date_plus_time+".png", dpi=600)

        elif opt_time_frame == 'after_break':

            # Create folder, if necessary
            if not os.path.exists(path_2 + "/after_break"):
                os.makedirs(path_2 + "/after_break")

            # Save figure to .png
            plt.savefig(path_2 + "/after_break" + "/" + date_plus_time+".png", dpi=600)

    plt.axvline(x=5.0)

    if verbose != 0:
        plt.show()

    plt.close()



def plot_CN__c_s(t, c_s_CN_high, ax3, verbose, menu, discretization_scheme):

    if (verbose == 1 or verbose == 2) and discretization_scheme != 1:

        ax3.scatter(t, c_s_CN_high.val, color='purple', marker='x')

        if menu==1:
            plt.pause(0.005)



def plot__h_mod(t, h_mod_CN_high, degrad_max_CN_high, degrad_break, alpha, c_sat, c_s_CN_high, limiter_break, limiter_sat, ax4, verbose, menu, discretization_scheme):

    _limiter_break = limiter_break
    _limiter_sat = limiter_sat

    # TODO: später wieder auf orignale if Bedingung zurücksetzen
    if ((menu==1 or menu == 2) and (verbose == 1 or verbose == 2)) and discretization_scheme != 1:
    #if ((menu==1 or menu == 2) and verbose == 2) and discretization_scheme != 1:

        # Plot
        ax4.scatter(t, h_mod_CN_high, color='blue', marker='x')

        if degrad_max_CN_high > degrad_break and _limiter_break <= 1:

            _limiter_break+=1
            ax4.axvline(t, color='r', label='Matrixbruch')
            ax4.legend()

        if alpha * (c_sat-c_s_CN_high.val) < 1 and _limiter_sat <= 1:

            _limiter_sat+=1
            ax4.axvline(t, color='m', label='Saettigung')
            ax4.legend()

        if menu==1:
            plt.pause(0.005)


    return _limiter_break, _limiter_sat



def plot_print__L2_errornorm(ax5, t, T_explicit_Euler, num_steps, L2_errornorm):

    # Print
    if t < T_explicit_Euler:

        print(f'{num_steps}. time step (EULER): t = {t}: L2-Errornorm = {L2_errornorm}')

    else:

        print(f'{num_steps}. time step (CN):    t = {t}: L2-Errornorm = {L2_errornorm}')


    # Plot
    ax5.scatter(t, L2_errornorm, color='red')
    ax5.get_yaxis().get_major_formatter().set_useOffset(False)

    plt.pause(0.05)



def plot__convergence_analysis_results(num_steps_explicit_Euler, nx, L2_error_avg_list, degree, space_or_time):

    if space_or_time == 1:

        # Plot results
        for i in range(len(degree)):
            fig, ax = plt.subplots()
            fig.suptitle('P{}'.format(degree[i]), fontsize=15)
            plt.yscale("log")
            plt.xscale("log")
            ax.set_xlabel(r"$n_{\mathrm{s}}$", loc="right")
            ax.set_ylabel(r"$avg(L2-error)$", loc="center")
            ax.plot(nx, L2_error_avg_list[i])

    elif space_or_time == 2:

        # Plot results
        for i in range(len(degree)):
            fig, ax = plt.subplots()
            fig.suptitle(f'P{degree[i]}-Element (n_s = {nx})', fontsize=15)
            ax.set_xlabel(r"$n_{\mathrm{t}}$", loc="right", fontsize=15)
            ax.set_ylabel(r"$avg(L2-error)}$", loc="center", fontsize=15)
            ax.plot(num_steps_explicit_Euler, L2_error_avg_list[i])

    # Show plots
    plt.show()



def show_plot(menu):

    if menu == 1:
        plt.show()

    elif menu == 2:
        plt.show(block=False)
        ##plt.pause(5)
        #plt.close()
        #plt.close()



# ----------------------------------
# P R I N T
# ----------------------------------

def print__R_squared_for_given_time_frame(F_num_at_discrete_points, F_exp, t_exp, t_start, t_end, verbose, opt_time_frame):
    """ R² = 1-SRS/SST; SSR = sum of residuals squared, SST = total sum of squares"""
    # https://www.ncl.ac.uk/webtemplate/ask-assets/external/maths-resources/statistics/regression-and-correlation/coefficient-of-determination-r-squared.html

    if verbose != 0 and opt_time_frame == 'until_break':

        # Initialize internal parameters
        _t_exp = t_exp
        _F_exp = F_exp
        _F_num_at_discrete_points = list(F_num_at_discrete_points)

        # Consider time frame
        if t_start != 0:

            # Calc position of t_start in t_exp
            t_start_count = t_exp.index(t_start)

            # Remove elements from list
            _t_exp = _t_exp[t_start_count:]
            _F_exp = _F_exp[t_start_count:]
            _F_num_at_discrete_points = _F_num_at_discrete_points[t_start_count:]

        if t_end != 5:

            # Calc position of t_end in t_exp
            t_end_count = t_exp.index(t_end)

            # Remove elements from list
            _t_exp = _t_exp[:t_end_count+1]
            _F_exp = _F_exp[:t_end_count+1]
            _F_num_at_discrete_points = _F_num_at_discrete_points[:t_end_count+1]

        # Initialize and calc sum of squared residuals (SSR = sum of residuals squared
        _residuals_squared = []

        for i in range(len(_F_exp)):
            _residuals_squared += [(_F_exp[i]-_F_num_at_discrete_points[i])**2]

        _sum_of_residuals_squared = sum(_residuals_squared)

        # Initialize and calc sum of squared differences between F_exp and the average of F_exp (SST = total sum of squares)
        _avg_F_exp = sum(_F_exp)/len(_F_exp)
        _squared_differences_between_F_exp_and_avg_of_F_exp = []

        for j in range(len(_F_exp)):
            _squared_differences_between_F_exp_and_avg_of_F_exp += [(_F_exp[j]-_avg_F_exp)**2]

        _sum_of_squared_differences_between_F_exp_and_avg_of_F_exp = sum(_squared_differences_between_F_exp_and_avg_of_F_exp) # = SST

        # Calc R² = 1-SRS/SST
        _R_squared = 1.0 - (_sum_of_residuals_squared/_sum_of_squared_differences_between_F_exp_and_avg_of_F_exp)
        print(f' R²(week: {t_start}-{t_end})  = {_R_squared}')

        if verbose == 2:
            # Check R² by calc R² with sklearn.metrics r2_score package
            _r2 = r2_score(_F_exp, _F_num_at_discrete_points)
            print(f' R²(week: {t_start}-{t_end})  = {_r2} (calculated by sklearn.metrics.r2_score package)')



def print_save__R_squared(F_num_at_discrete_points, F_exp, t_exp, t_start, t_end, parameters_to_optimize, A_m_t0, D_ref, c_ref, l_ref, verbose, menu, distribution_menu, opt_time_frame):
    """ R² = 1-SRS/SST; SSR = sum of residuals squared, SST = total sum of squares"""
    # https://www.ncl.ac.uk/webtemplate/ask-assets/external/maths-resources/statistics/regression-and-correlation/coefficient-of-determination-r-squared.html

    # Initialize internal parameters
    _t_exp = t_exp
    _F_exp = F_exp
    _F_num_at_discrete_points = F_num_at_discrete_points

    """
    # Consider time frame
    if t_start != 0:

        # Calc position of t_start in t_exp
        t_start_count = t_exp.index(t_start)

        # Remove elements from list
        _t_exp = _t_exp[t_start_count:]
        _F_exp = _F_exp[t_start_count:]
        _F_num_at_discrete_points = _F_num_at_discrete_points[t_start_count:]
    """

    # Initialize and calc sum of squared residuals (SSR = sum of residuals squared)
    _residuals_squared = []

    for i in range(len(_F_exp)):
        _residuals_squared += [(_F_exp[i]-_F_num_at_discrete_points[i])**2]

    _sum_of_residuals_squared = sum(_residuals_squared)


    # Initialize and calc sum of squared differences between F_exp and the average of F_exp (SST = total sum of squares)
    _avg_F_exp = sum(_F_exp)/len(_F_exp)
    _squared_differences_between_F_exp_and_avg_of_F_exp = []

    for j in range(len(_F_exp)):
        _squared_differences_between_F_exp_and_avg_of_F_exp += [(_F_exp[j]-_avg_F_exp)**2]

    _sum_of_squared_differences_between_F_exp_and_avg_of_F_exp = sum(_squared_differences_between_F_exp_and_avg_of_F_exp) # = SST


    # Calc R² = 1-SRS/SST
    _R_squared = 1.0 - (_sum_of_residuals_squared/_sum_of_squared_differences_between_F_exp_and_avg_of_F_exp)

    # Save R² to CSV
    date = str(datetime.date.today())  # convert date to str

    if distribution_menu == 0:
        start_distribution = "constant_c0"
    elif distribution_menu == 1:
        start_distribution = "parabolic_c0"

    if menu == 1:

        [D, alpha, l_b, A_m_degrad, A_m_t0, const_degrad, degrad_break] = parameters_to_optimize

        # Create path
        path_1 = str("results/01__simulation/"+date+"/"+start_distribution+"/R_squared")  # create path to .csv-file

        # Create folder
        if not os.path.exists(path_1):
            os.makedirs(path_1)

        # Create filename
        filename_1 = "/"+date+"__"+start_distribution+"__"+"R_squared.csv"

        # Write column name to .csv-file
        if not os.path.isfile(path_1 + filename_1):
            with open(path_1 + filename_1, 'w') as f:
                writer = csv.writer(f, delimiter='\t')
                writer.writerow(['R**2 = 1-SRS/SST', 'D [mm**2/s]', 'alpha [mm**3/g]', 'l_b [mm]', 'A_m_degrad = x * A_m_t0 [-]', 'const_degrad [-]', 'degrad_break [-]' ])

        # Write results to .csv-file
        with open(path_1 + filename_1, 'a') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow([_R_squared, D * D_ref, alpha / c_ref, l_b * l_ref, A_m_degrad/A_m_t0, const_degrad, degrad_break])
            f.close()

    elif menu == 2:

        if t_start == 0 and opt_time_frame == 'until_break':

            [D, alpha, l_b] = parameters_to_optimize

            # Create path
            path_2a = str("results/02__optimization/" + date + "/" + start_distribution + "/until_break") # create path to .csv-file

            # Create folder
            if not os.path.exists(path_2a):
                os.makedirs(path_2a)

            # Create filename
            filename_2a = "/" + date + "__" + start_distribution + "__" + "R_squared_Optimization_until_break.csv"

            # Write column name to .csv-file
            if not os.path.isfile(path_2a + filename_2a):
                with open(path_2a + filename_2a, 'w') as f:
                    writer = csv.writer(f, delimiter='\t')
                    writer.writerow(['R² = 1-SRS/SST', 'D [mm²/s]', 'alpha [mm³/g]', 'l_b [mm]'])

            # Write results to .csv-file
            with open(path_2a + filename_2a, 'a') as f:
                writer = csv.writer(f, delimiter='\t')
                writer.writerow([_R_squared, D*D_ref, alpha/c_ref, l_b*l_ref])
                f.close()


        elif t_start == 0 and opt_time_frame == 'after_break':

            A_m_degrad = parameters_to_optimize[0]

            # Create path
            path_2b = str("results/02__optimization/" + date + "/" + start_distribution + "/after_break") # create path to .csv-file

            # Create folder
            if not os.path.exists(path_2b):
                os.makedirs(path_2b)

            # Create filename
            filename_2b = "/" + date + "__" + start_distribution + "__" + "R_squared_Optimization_after_break.csv"

            # Write column name to .csv-file
            if not os.path.isfile(path_2b + filename_2b):
                with open(path_2b + filename_2b, 'w') as f:
                    writer = csv.writer(f, delimiter='\t')
                    writer.writerow(['R^2 = 1-SRS/SST', 'A_m_degrad [mm^2]', 'A_m_degrad/A_m_t0 [-]'])

            # Write results to .csv-file
            with open(path_2b + filename_2b, 'a') as f:
                writer = csv.writer(f, delimiter='\t')
                writer.writerow([_R_squared, A_m_degrad*l_ref**2, A_m_degrad/A_m_t0])
                f.close()

    # Print
    if verbose != 0:

        # Print R²
        print('')
        print(f' R²(week: {t_start}-{t_end})  = {_R_squared}')

        if verbose == 2:

            # Check R² by calc R² with sklearn.metrics r2_score package
            _r2 = r2_score(_F_exp, _F_num_at_discrete_points)
            print(f' R²(week: {t_start}-{t_end})  = {_r2} (calculated by sklearn.metrics.r2_score package)')


    return _R_squared



def print__residual_at_week_5(Fnum_at_discrete_points, F_exp, verbose, opt_time_frame):

    if opt_time_frame == 'until_break' and verbose != 0:

        # Calc residuals at t=5
        print(' residual (week: 5) =', F_exp[5]-Fnum_at_discrete_points[5])
        print(" --------------------------------------------------------------")
        print('')



def print__convergence_rate(rate, nx, element_family, degree):

    for i in range(len(nx)-1):
        print('      nx=', nx[i], '--> nx=', nx[i+1])

        for n in range(len(degree)):

            print(element_family, degree[n], ': r =', rate[n][i])

        print('--------------------------')



def print_opt_parameters(parameters_to_optimize, const_degrad, A_m_t0, c_ref, D_ref, l_ref, opt_time_frame, verbose):

    if opt_time_frame == 'until_break':

        D_opt_units     = parameters_to_optimize[0] * D_ref
        alpha_opt_units = parameters_to_optimize[1] / c_ref
        l_b_opt_units   = parameters_to_optimize[2] * l_ref

        print('')
        print('#-------------------------------------------------------------')
        print('# D_opt = ', D_opt_units, '[mm^2/s] ,      alpha_opt = ', alpha_opt_units, '[mm^3/g] ,      l_b_opt = ', l_b_opt_units, '[mm^3/g]')
        print('#-------------------------------------------------------------')


    elif opt_time_frame == 'after_break':

        A_m_degrad_opt_units = parameters_to_optimize[0] * l_ref**2
        q_A_m = parameters_to_optimize[0]/A_m_t0

        print('')
        print('#-------------------------------------------------------------')
        print('# A_m_degrad_opt = ', A_m_degrad_opt_units, '[mm²]      mit const_degrad = ', const_degrad, "[-]")
        print('# A_m_degrad_opt = ', q_A_m, '* A_m_t0')
        print('#-------------------------------------------------------------')



def print_plot_all_results(results, menu, t_exp, F_exp, T, kind_of_1d_interpolation, element_family, degree, nx, h, D, alpha, l_b, const_degrad, degrad_break, A_m_t0, A_m_t0_units, A_m_degrad_constant_units, A_m_degrad_parabolic_units, opt_success, c_ref, D_ref, l_ref):

    if menu == 1: # simulation mode was choosen

        # Extract parameters
        _t_num = results[0]
        _F_num = results[1]
        _distribution_menu = results[2]

        # Define dummy variables to be able to use the same functions
        DUMMY_opt_time_frame = 'after_break'

        if _distribution_menu == 0:
            DUMMY_parameters_to_optimize = [D, alpha, l_b, A_m_degrad_constant_units, A_m_t0_units, const_degrad, degrad_break]
        elif _distribution_menu == 1:
            DUMMY_parameters_to_optimize = [D, alpha, l_b, A_m_degrad_parabolic_units, A_m_t0_units, const_degrad, degrad_break]


        # Post-processing
        # t = 0-5
        F_num_at_discrete_points_0_5  = calc_F_num_at_discrete_time_steps(_t_num, _F_num, t_exp, 5.0, kind_of_1d_interpolation)
        t_exp__opt_0_5, F_exp__opt_0_5 = adjust_t_exp_and_F_exp_to_optimization_time_frame(t_exp, F_exp, 5.0)
        R_squared_0_1 = print__R_squared_for_given_time_frame(F_num_at_discrete_points_0_5, F_exp__opt_0_5, t_exp__opt_0_5, 0.0, 1.0, 2,'until_break')  # R² initial burst     (week:0-1)
        R_squared_3_5 = print__R_squared_for_given_time_frame(F_num_at_discrete_points_0_5, F_exp__opt_0_5, t_exp__opt_0_5, 3.0, 5.0, 2,'until_break')  # R² during stagnation (week:3-5)
        # t = 0-15
        F_num_at_discrete_points_0_15 = calc_F_num_at_discrete_time_steps(_t_num, _F_num, t_exp, T, kind_of_1d_interpolation)
        R_squared_0_15 = print_save__R_squared(F_num_at_discrete_points_0_15, F_exp, t_exp, 0.0, T, DUMMY_parameters_to_optimize, A_m_t0, D_ref, c_ref, l_ref, 2, menu, _distribution_menu, DUMMY_opt_time_frame)
        plot__F_exp_vs_F_num(DUMMY_parameters_to_optimize, D_ref, c_ref, l_ref, t_exp, _t_num, F_exp, _F_num, T, R_squared_0_15, kind_of_1d_interpolation, 2, menu, opt_success, _distribution_menu, DUMMY_opt_time_frame)


    elif menu == 2: # inverse analysis

        # Extract optimized parameters
        _res_opt  = results[0]

        # Extract t_num & F_num for optimized parameters
        _t_num = results[1]
        _F_num = results[2]

        # Extract various
        _opt_time_frame    = results[3]
        _distribution_menu = results[4]

        print('')
        print('')
        print('')
        print('#-------------------------------------------------------------')
        print('#-------------------------------------------------------------')
        print('# S U C C E S S F U L L Y   O P T I M I Z E D')
        print('#-------------------------------------------------------------')
        print('#-------------------------------------------------------------')
        print('')

        # Post-processing
        print_opt_parameters(_res_opt.x, const_degrad, A_m_t0, c_ref, D_ref, l_ref, _opt_time_frame, 2)
        F_num_at_discrete_points =  calc_F_num_at_discrete_time_steps(_t_num, _F_num, t_exp, T, kind_of_1d_interpolation)
        R_squared_0_T = print_save__R_squared(F_num_at_discrete_points, F_exp, t_exp, 0, T, _res_opt.x, A_m_t0, D_ref, c_ref, l_ref, 2, 'dummy', _distribution_menu, _opt_time_frame) # menu = 'dummy' to avoid saving R²(0,15) to R²(0,5) resp. R²(5,15)
        plot__F_exp_vs_F_num(_res_opt.x, D_ref, c_ref, l_ref, t_exp, _t_num, F_exp, _F_num, T, R_squared_0_T, kind_of_1d_interpolation, 2, menu, opt_success, _distribution_menu, _opt_time_frame)


    elif menu == 4:  # convergence analysis

        # Extract parameters
        _nx                       = results[0]
        _nt                       = results[1]
        _L2_error_avg_list        = results[2]
        _rate                     = results[3]
        _diff_F_num_total_avg     = results[4]
        _space_or_time            = results[5]
        _discretization_scheme    = results[6]

        if 1 <= _discretization_scheme <= 5:

            plot__convergence_analysis_results(_nt, _nx, _L2_error_avg_list, degree, _space_or_time)

            if _space_or_time == 1:
                # Print convergence rate
                print__convergence_rate(_rate, _nx, element_family, degree)



        elif _discretization_scheme == 6 or _discretization_scheme == 7:

            plot__convergence_analysis_results(_nt, _nx, _diff_F_num_total_avg, _space_or_time)



    elif menu == 5:

        # Extract parameters
        rel_err_percent_list = results[0]
        distribution_menu    = results[1]

        # Plot
        plot__rel_err_m_d(nx, rel_err_percent_list, element_family, degree, distribution_menu)
        plt.show()



def print__t_dt_num_steps(t, dt, num_steps, T_explicit_Euler, menu, verbose):

    if menu == 1 or menu ==2:

        if verbose == 0 and menu == 2 and t > T_explicit_Euler:

            print(f'{num_steps}. CN-step: t = {t} (dt = {dt})')


        if verbose == 1:

            if t > T_explicit_Euler:

                print(f'{num_steps}. CN-step: t = {t} (dt = {dt})')


        elif verbose == 2:

            if t > T_explicit_Euler:
                print(f'{num_steps}. CN-step: t = {t} (dt = {dt})')

            else:
                print(f'{num_steps}. EULER-step: t = {t} (dt = {dt})')



def plot__rel_err_m_d(nx, rel_err_percent_list, element_family, degree, distribution_menu):

    # plots initialisieren
    fig, ax = plt.subplots(1, 2)

    if distribution_menu == 0:
        plt.suptitle(f"CONSTANT distribution for u(t=0)", fontsize=15)
    elif distribution_menu == 1:
        plt.suptitle(f"PARABOLIC distribution for u(t=0)", fontsize=15)

    # Achsenbezeichnung
    ax[0].set_xlabel(r"$n_\mathrm{s} [-]$", loc="right")
    ax[1].set_xlabel(r"$n_\mathrm{s} [-]$", loc="right")
    ax[0].set_ylabel(r"$err_{m_\mathrm{d}} [\%]$", loc="top", rotation="horizontal")
    ax[1].set_ylabel(r"$err_{m_\mathrm{d}} [\%]$", loc="top", rotation="horizontal")

    # Plot
    ax[0].plot(nx, rel_err_percent_list[0], color='orange', marker='x', label=f'{element_family}{degree[0]}-Element')
    ax[1].plot(nx, rel_err_percent_list[1], color='red'   , marker='x', label=f'{element_family}{degree[1]}-Element')

    # Show labels
    ax[0].legend()
    ax[1].legend()