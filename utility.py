# -----------------------------------------------------------------
#  U T I L I T Y   F U N C T I O N S
# -----------------------------------------------------------------
from fenics import *
from scipy.interpolate import interp1d
import numpy as np
import time
import sys


def check_CFL_conditions(dt, D, mesh, verbose, discretization_scheme, menu):
    """ Checks, if CFL=D * dt/dx < 1 is met and reduces dt by dt=dt/5 till condition is met. """

    _dt=dt
    element_size = mesh.hmin()
    CFL = D * _dt/element_size

    if menu == 1 or menu == 2:

        while CFL > 1.0:

            if CFL > 1.0:

                if verbose==1 or verbose==2:

                    print("#######")
                    print(" # ! #   CFL = D * dt/dx > 1    (CFL - Bedingungen nicht erfüllt)")
                    print("  # #")
                    print("   #")
                    print('#-----------------------------------------------------')

                _dt=_dt/5


            else:

                _dt=_dt

        if verbose == 1 or verbose == 2:
            print("CFL conditions satisfied:")
            print("CFL = ", CFL)
            print("dx = ", element_size)
        print('')

        return _dt


    elif menu == 4 and (discretization_scheme != 2 and discretization_scheme != 3):

        if CFL > 1.0:

            print("#######")
            print(" # ! #   CFL = D * dt/dx > 1    (CFL - Bedingungen nicht erfüllt)")
            print("  # #")
            print("   #")
            print('#-------------------------------------------------------------')

            sys.exit("Choose dt, that it will conform to CFL-conditions!")


        else:

            if verbose == 1 or verbose == 2:

                print("CFL conditions satisfied")

            return _dt


    elif menu == 4 and (discretization_scheme == 2 or discretization_scheme == 3):

        return _dt



def check__u_t0(u_t0, m_d_t0, r, H, l_ref, c_ref, verbose, menu):
    """ Checks, if distribution is chosen correctly, that m_d_units = 7.5g. """

    if (menu == 1 and (verbose == 1 or verbose == 2)) or (menu == 2 and verbose == 2) or menu == 5:

        # TODO: use m_d_t0_units from parameters_units.py
        m_d_t0_units             = 0.0075                                          # [g]      - Gesamtmasse GM, die in die Matrix eingebracht wurde (m_d = 7.5 mg)

        # Calc m_d_test_units [g]
        _m_d_test_units = 2.0 * np.pi * H * assemble(u_t0 * r * dx) * (c_ref * l_ref ** 3)

        # Calc rel error [%]
        _rel_err_percent = ((_m_d_test_units - m_d_t0_units) / m_d_t0_units) * 100


        if menu == 5:

            return _rel_err_percent



def choose_F_num(menu):

    if menu == 1:

        _F_num_option = 0 # initialize menu with dummy variables

        while _F_num_option != 1 and _F_num_option != 2:

            print('1. Using c_s to calculate F_num')
            print('2. Usinǵ F_num^complex')
            print('')

            _F_num_option = int(input())
            print('')

            if _F_num_option != 1 and _F_num_option != 2:
                print('You did not choose "1" or "2"!')
                print('')


        return _F_num_option


    elif menu == 2 or menu == 4:

        _F_num_option = 2

        # TODO: später wieder aktivieren, nur zur besseren Lesbarkeit bei der Optimierung deaktiviert
        """
        if menu == 2:
            print('')
            print('T O  -  D O: function "choose_F_num" muss für die Konvergenzanalyse noch sauber programmiert werden. Aktuell wird "_F_num_option = 1" automatisch gesetzt.')
            print('')
        """

        return _F_num_option



def calc_F(dt, D, u_n, c_s, F_num_komplex_n, c_max, V_s_units, V_c_units, H, beta, gamma, r, n_facet, ds, discretization_scheme, menu):
    """ Calc drug release value F_num. """

    _F_num_komplex = F_num_komplex_n - dt * ((2.0 * np.pi * H * gamma) / c_max) * assemble(D * inner(grad(u_n), n_facet) * r * ds)

    if menu == 1 or menu == 2:
        _F_num_simpel  = c_s.val / c_max

    elif menu == 4:
        _F_num_simpel = c_s / c_max


    return _F_num_simpel, _F_num_komplex



def calc_F_num_at_discrete_time_steps (_tnum, _Fnum, _t_exp, T, kind_of_1d_interpolation):
    """Calc _F_num at discrete time steps (these are the time points, where measurements were taken during the experiment)"""

    # Initialize
    _F_num_at_discrete_points = []

    # Interpolate _F_num
    _F_num_interpolated_1d = interp1d(_tnum, _Fnum, kind_of_1d_interpolation)

    # count user specified week i
    _t_i_count = _t_exp.index(T) + 1

    # Optimization until user specified week i
    for i in range(_t_i_count):
        _F_num_at_discrete_points += [_F_num_interpolated_1d(_t_exp[i])]


    return _F_num_at_discrete_points



def adjust_t_exp_and_F_exp_to_optimization_time_frame(t_exp, F_exp, T__opt_end):

    # Adjust t_exp & F_exp to optimization time frame (t = 0-5)
    t_i_count = t_exp.index(T__opt_end) + 1

    _t_exp__opt = []
    _F_exp__opt = []

    for i in range(t_i_count):
        _t_exp__opt += [t_exp[i]]
        _F_exp__opt += [F_exp[i]]


    return _t_exp__opt, _F_exp__opt



def adjust_dt_to_dt_fine(dt, fineness, discretization_scheme):
    """ Changes time step from coarse to fine. """

    # Disabled for convergence analysis
    if discretization_scheme == 2:
        return dt

    # Reduce time step
    else:
        return dt / fineness



def update_t_to_current_time(t, dt, fineness, discretization_scheme):

    _t=t # initialize internal variable

    if discretization_scheme == 2:
        _t += dt

    else:
        _t += fineness * dt


    return _t



def compute__convergence_rate(L2_error_avg_list, degree, nx, nt, h, space_or_time):
    """ Compute convergence rate. """

    if space_or_time == 1:
        for i in range(len(nx)):
            h += [1 / nx[i]]

    elif space_or_time == 2:
        for i in range(len(nt)):
            h += [1 / nt[i]]

    # rate = [ [[],[],[],[]], [[],[],[],[]], [[],[],[],[]]]
    rate = [[[] for _ in range(len(h) - 1)] for _ in range(len(degree))]

    for l in range(len(degree)):

        for m in range(len(h) - 1):
            rate[l][m] += [np.log(L2_error_avg_list[l][m + 1] / L2_error_avg_list[l][m]) / np.log(h[m + 1] / h[m])]


    return rate



def check_cost_function(F_exp, _F_num_at_discrete_points, _weight, verbose):
    """ Checks cost function during parameter optimization. """

    if verbose == 2:

        # Initialize weighted squared residuals
        _weighted_squared_residuals = []

        # Calc weighted squared residuals
        for k in range(len(F_exp)):
            _weighted_squared_residuals.append((_weight[k] * (F_exp[k] - _F_num_at_discrete_points[k]))**2)

        # Calc cost function
        _cost_function = 0.5 * sum(_weighted_squared_residuals)

        # Print cost function
        print('cost function (Check) = %.4e'% _cost_function)



def add_fitted_points_to_experimenta_data(t_add_start, t_add_end, n_add, t_exp, F_exp, T__opt_end, T, kind_of_1d_interpolation):

    # F_exp interpolieren
    _F_exp_interp1d = interp1d(t_exp, F_exp, kind_of_1d_interpolation)
    if t_add_start == 3.0:
        #n_add = 20
        n_plot = int(((T/(t_add_end - t_add_start))*n_add) + 1)
        #n_plot = 151             # 5-3 = 2    ; 15/2 = 7.5  ; 7.5*20 +1 = 151
    elif t_add_start == 4.9:
        #n_add = 20
        n_plot = int(((T/(t_add_end - t_add_start))*n_add) + 1)
        #n_plot = 3001            # 5-4.9 = 0.1; 15/0.1 = 150; 150 * 20 +1 = 3001
    elif t_add_start == 0.0:
        #n_add = 100
        n_plot = int(((T / (t_add_end - t_add_start)) * n_add) + 1)
        #n_plot = 301                # 5-0   = 5.0; 15/5.0 = 3  ; 3*100 +1    = 301

    _t_plot = np.linspace(0, T__opt_end, n_plot)
    _F_exp_interpolated = _F_exp_interp1d(_t_plot)

    # Create t_exp_artifical + F_exp_artificial
    _t_exp_artificial = []
    _F_exp_artificial = []

    for i in range(len(t_exp)):
        if t_exp[i] <= t_add_start:
            _t_exp_artificial.append(t_exp[i])
            _F_exp_artificial.append(F_exp[i])

    for i in range(len(_t_plot)):
        if _t_plot[i] > t_add_start and _t_plot[i] <= t_add_end:
            _t_exp_artificial.append(_t_plot[i])
            _F_exp_artificial.append(_F_exp_interpolated[i])

    for i in range(len(t_exp)):
        if t_exp[i] > t_add_end:
            _t_exp_artificial.append(t_exp[i])
            _F_exp_artificial.append(F_exp[i])


    return _t_exp_artificial, _F_exp_artificial



def stop_timer(tic, menu_timer):

    if menu_timer == 'timer_on':

        toc = time.perf_counter()

        print(f' Time needed to run calculation: {toc - tic:0.4f} seconds')
        print('')