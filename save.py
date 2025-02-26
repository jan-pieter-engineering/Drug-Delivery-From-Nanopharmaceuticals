# -----------------------------------------------------------------
#  S A V E
# -----------------------------------------------------------------
import csv
import math
import datetime
import os
from itertools import zip_longest
from dolfin import errornorm, File, XDMFFile
from plot_funcs import plot_print__L2_errornorm


def collect_results(t_num, F_num, L2_error, menu, distribution_menu):

    if menu == 1:

        save_tnum_Fnum_to_csv(t_num, F_num, distribution_menu)

        return [t_num, F_num]


    elif menu == 2:

        return [t_num, F_num]


    elif menu == 4:
        """
        For floating point numbers the numerical precision of sum (and np.add.reduce) is in general limited by directly adding each number individually
        to the result causing rounding errors in every step. However, often numpy will use a numerically better approach (partial pairwise summation)
        leading to improved precision in many use-cases. In contrast to NumPy, Python’s math.fsum function uses a slower but more precise approach to summation.
        
        https://numpy.org/doc/stable/reference/generated/numpy.sum.html
        """

        #Calc
        #_L2_error_avg = 1.0/len(L2_error) * np.sum(L2_error)
        _L2_error_avg = 1.0/len(L2_error) * math.fsum(L2_error)
        #_L2_error_avg = 1.0/len(L2_error) * sum(L2_error)

        # Print
        print('len(L2_error) =', len(L2_error))
        print('avg(L2_error) =', _L2_error_avg)

        return _L2_error_avg



# --------------------------------------------------------
#  save t_num, F_num
# --------------------------------------------------------
def save_t_num_F_num(t, t_num, F_num, F_num_simpel_CN_high, F_num_komplex_CN_high, F_num_choice):

    # Define internal variables
    _t = t
    _t_num = t_num
    _F_num = F_num
    _F_num_simpel_CN_high  = F_num_simpel_CN_high
    _F_num_komplex_CN_high = F_num_komplex_CN_high

    # Save t_num+F_num
    if _t - _t_num[len(_t_num) - 1] > 0.2 and _t < 1.0:  # exclude values during initial time-step finding

        _t_num += [_t]        # save t_num

        if F_num_choice == 1: # save F_num
            _F_num += [_F_num_simpel_CN_high]

        elif F_num_choice == 2:
            _F_num += [_F_num_komplex_CN_high]

    elif _t >= 1.0:  # includes values of time-step finding after matrix break

        _t_num += [_t]         # save t_num

        if F_num_choice == 1:  # save F_num
            _F_num += [_F_num_simpel_CN_high]

        elif F_num_choice == 2:
            _F_num += [_F_num_komplex_CN_high]


    return _t_num, _F_num


# --------------------------------------------------------
#  save L2-error
# --------------------------------------------------------
def calc_plot_print_L2_error(t, T_explicit_Euler, u_e, u, _L2_error, num_steps, nt, ax5, verbose, discretization_scheme):

    if discretization_scheme != 0:

        if num_steps <= nt:

            # Calc & collect L2-errornorm
            _L2_error += [errornorm(u_e, u, 'L2', degree_rise=3)]

            if verbose == 2:

                # Plot/print L2-errornorm
                plot_print__L2_errornorm(ax5, t, T_explicit_Euler, num_steps, _L2_error[-1])

        elif num_steps > nt:

            print('W A R N I N G: L2_error[', num_steps, '] was not added to L2_error')
            print('')


    return _L2_error



# --------------------------------------------------------
#  save pvd, XDMF
# --------------------------------------------------------
def create__save_files(menu, distribution_menu):

    # Read date
    date = str(datetime.date.today())  # convert date to str

    # Check distribution
    if distribution_menu == 0:
        start_distribution = "/constant_c0"
    elif distribution_menu == 1:
        start_distribution = "/parabolic_c0"

    # Create path
    path = str("results/01__simulation/" + date + start_distribution + "/solution")

    # Create folder
    if not os.path.exists(path):
        os.makedirs(path)

    # Save files
    if menu == 1:

        _vtk_file  =     File(path+"/pvd/results__simulation.pvd")  # create vtk-file
        _xdmf_file = XDMFFile(path+"/xdmf/results__simulation.xdmf") # create XDMF-file

    elif menu == 2 or menu == 4:

        _vtk_file  = 0 # dummy
        _xdmf_file = 0 # dummy


    return _vtk_file, _xdmf_file



def create__u_t0_save_files(distribution_menu, element_family, degree):

    # Create path
    if distribution_menu == 0:
        start_distribution = "constant_c0/"
    elif distribution_menu == 1:
        start_distribution = "parabolic_c0/"

    date = str(datetime.date.today())  # convert date to str
    path = str("results/05__testing/rel_err_m_d/" + start_distribution + date)  # create path to .csv-file

    # Create folder, if it doesn't exist
    if not os.path.exists(path+"/vtk"):
        os.makedirs(path+"/vtk")
    elif not os.path.exists(path+"/xdmf"):
        os.makedirs(path+"/xdmf")

    # Create file
    _vtk_file_u_t0  =     File(path+"/vtk/nx_ut0__"+str(element_family)+str(degree)+"_.pvd")   # create vtk-file
    _xdmf_file_u_t0 = XDMFFile(path+"/xdmf/nx_ut0__"+str(element_family)+str(degree)+"_.xdmf") # create XDMF-file


    return _vtk_file_u_t0, _xdmf_file_u_t0



def save_solution(t, u_n, vtk_file, xdmf_file, menu):

    # Initialize internal variables
    _vtk_file  = vtk_file
    _xdmf_file = xdmf_file

    if menu == 1:

        # Save u (numerical)
        _vtk_file << (u_n, t)    # (t, u_high) in vtk -file hinterlegen
        _xdmf_file.write(u_n, t) # (t, u_high) in XDMF-file hinterlegen


def save_u_t0(nx, u_t0, vtk_file, xdmf_file):

    # Initialize internal variables
    _vtk_file  = vtk_file
    _xdmf_file = xdmf_file

    # Save u_t0 (numerical)
    _vtk_file << (u_t0, nx)    # (nx, u_t0) in vtk -file hinterlegen
    _xdmf_file.write(u_t0, nx) # (nx, u_t0) in XDMF-file hinterlegen




# --------------------------------------------------------
#  save CSV
# --------------------------------------------------------
def save_error_to_csv(nx, nt, dt, L2_error_avg_list, D, D_ref, T, T_explicit_Euler, space_or_time, discretization_scheme):

    if discretization_scheme == 1 and space_or_time == 1:
        with open('results/04__convergence_analysis/011__space__EULER/results__convergence_analysis__L2_error__space__EULER.csv','w') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(['D', 'D_ref', 'T'])
            writer.writerow([D, D_ref, T_explicit_Euler])
            writer.writerow([''])
            writer.writerow(['n_s', 'err_n_s'])
            writer.writerow([''])
            writer.writerow(nx)
            writer.writerows(L2_error_avg_list)

    elif discretization_scheme == 1 and space_or_time == 2:
        with open('results/04__convergence_analysis/021__time__EULER/results__convergence_analysis__L2_error__time__EULER.csv','w') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(['D', 'D_ref', 'T', 'n_s'])
            writer.writerow([D, D_ref, T_explicit_Euler, nx])
            writer.writerow([''])
            writer.writerow(['n_t', 'err_n_t'])
            writer.writerow([''])
            writer.writerow(nt)
            writer.writerows(L2_error_avg_list)

    elif discretization_scheme == 2 and space_or_time == 1:
        with open('results/04__convergence_analysis/012__space__CRANK_NICOLSON__low_precision/results__convergence_analysis__L2_error__space__CRANK_NICOLSON__low_precision.csv','w') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(['D', 'D_ref', 'T'])
            writer.writerow([D, D_ref, T])
            writer.writerow([''])
            writer.writerow(['n_s', 'err_n_s'])
            writer.writerow([''])
            writer.writerow(nx)
            writer.writerows(L2_error_avg_list)

    elif discretization_scheme == 2 and space_or_time == 2:
        with open('results/04__convergence_analysis/022__time__CRANK_NICOLSON__low_precision/results__convergence_analysis__L2_error__time__CRANK_NICOLSON__low_precision.csv','w') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(['D', 'D_ref', 'T'])
            writer.writerow([D, D_ref, T])
            writer.writerow([''])
            writer.writerow(['n_t', 'err_n_t'])
            writer.writerow([''])
            writer.writerow(nt)
            writer.writerows(L2_error_avg_list)

    elif discretization_scheme == 3 and space_or_time == 1:
        with open('results/04__convergence_analysis/013__space__CRANK_NICOLSON__high_precision/results__convergence_analysis__L2_error__space__CRANK_NICOLSON__high_precision.csv','w') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(['D', 'D_ref', 'T'])
            writer.writerow([D, D_ref, T])
            writer.writerow([''])
            writer.writerow(['n_s', 'err_n_s'])
            writer.writerow([''])
            writer.writerow(nx)
            writer.writerows(L2_error_avg_list)

    elif discretization_scheme == 3 and space_or_time == 2:
        with open('results/04__convergence_analysis/023__time__CRANK_NICOLSON__high_precision/results__convergence_analysis__L2_error__time__CRANK_NICOLSON__high_precision.csv','w') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(['D', 'D_ref', 'T'])
            writer.writerow([D, D_ref, T])
            writer.writerow([''])
            writer.writerow(['n_t', 'dt', 'err_n_t'])
            writer.writerow([''])
            writer.writerow(nt)
            writer.writerows(L2_error_avg_list)



def save_rate_to_csv(nx, nt, dt, rate, D, D_ref, T, space_or_time, discretization_scheme):

        if discretization_scheme == 1 and space_or_time == 1:
            with open(
                    'results/04__convergence_analysis/011__space__EULER/results__convergence_analysis__rate__space__EULER.csv',
                    'w') as f:
                write = csv.writer(f, delimiter='\t')
                write.writerow(['D', 'D_ref', 'T'])
                write.writerow([D, D_ref, T])
                write.writerow([''])
                write.writerow(['n_s', 'Konvergenzrate'])
                write.writerow(nx)
                write.writerows(rate)


        elif discretization_scheme == 1 and space_or_time == 2:
            with open(
                    'results/04__convergence_analysis/021__time__EULER/results__convergence_analysis__rate__time__EULER.csv',
                    'w') as f:
                write = csv.writer(f, delimiter='\t')
                write.writerow(['D', 'D_ref', 'T'])
                write.writerow([D, D_ref, T])
                write.writerow([''])
                write.writerow(['n_t', 'Konvergenzrate'])
                write.writerow(nt)
                write.writerows(rate)


        elif discretization_scheme == 2 and space_or_time == 1:
            with open(
                    'results/04__convergence_analysis/012__space__CRANK_NICOLSON__low_precision/results__convergence_analysis__rate__space__CRANK_NICOLSON__low_precision.csv',
                    'w') as f:
                write = csv.writer(f, delimiter='\t')
                write.writerow(['D', 'D_ref', 'T'])
                write.writerow([D, D_ref, T])
                write.writerow([''])
                write.writerow(['n_s', 'Konvergenzrate'])
                write.writerow(nx)
                write.writerows(rate)


        elif discretization_scheme == 2 and space_or_time == 2:
            with open(
                    'results/04__convergence_analysis/022__time__CRANK_NICOLSON__low_precision/results__convergence_analysis__rate__time__CRANK_NICOLSON__low_precision.csv',
                    'w') as f:
                write = csv.writer(f, delimiter='\t')
                write.writerow(['D', 'D_ref', 'T'])
                write.writerow([D, D_ref, T])
                write.writerow([''])
                write.writerow(['n_t', 'Konvergenzrate'])
                write.writerow(nt)
                write.writerows(rate)


        elif discretization_scheme == 3 and space_or_time == 1:
            with open(
                    'results/04__convergence_analysis/013__space__CRANK_NICOLSON__high_precision/results__convergence_analysis__rate__space__CRANK_NICOLSON__high_precision.csv',
                    'w') as f:
                write = csv.writer(f, delimiter='\t')
                write.writerow(['D', 'D_ref', 'T'])
                write.writerow([D, D_ref, T])
                write.writerow([''])
                write.writerow(['n_s', 'Konvergenzrate'])
                write.writerow(nx)
                write.writerows(rate)


        elif discretization_scheme == 3 and space_or_time == 2:
            with open(
                    'results/04__convergence_analysis/023__time__CRANK_NICOLSON__high_precision/results__convergence_analysis__rate__time__CRANK_NICOLSON__high_precision.csv',
                    'w') as f:
                write = csv.writer(f, delimiter='\t')
                write.writerow(['D', 'D_ref', 'T'])
                write.writerow([D, D_ref, T])
                write.writerow([''])
                write.writerow(['n_t', 'Konvergenzrate'])
                write.writerow(nt)
                write.writerows(rate)



def save_tnum_Fnum_to_csv(t_num, F_num, distribution_menu):

    # Read date
    date = str(datetime.date.today())  # convert date to str

    # Check distribution
    if distribution_menu == 0:
        start_distribution = "constant_c0"
    elif distribution_menu == 1:
        start_distribution = "parabolic_c0"

    # Create path
    path = str("results/01__simulation/" + date + "/" + start_distribution + "/F_num")

    # Create folder
    if not os.path.exists(path):
        os.makedirs(path)

    # Save (t_num, F_num) in .csv file
    with open(path+'/'+date+'__'+start_distribution+'__'+'tnum_Fnum.csv', 'w') as f:
        write = csv.writer(f, delimiter='\t')
        write.writerow(['t_num', 'F_num'])
        write.writerow(t_num)
        write.writerow(F_num)


def save_rel_err_m_d_to_csv(nx, rel_err_m_d, degree, distribution_menu):

    # Create path
    if distribution_menu == 0:
        start_distribution = "constant_c0/"
    elif distribution_menu == 1:
        start_distribution = "parabolic_c0/"

    date = str(datetime.date.today()) # convert date to str
    path = str("results/05__testing/rel_err_m_d/"+start_distribution+date) # create path to .csv-file

    # Create folder, if it doesn't exist
    if not os.path.exists("results/05__testing/rel_err_m_d/"+start_distribution+date):
        os.makedirs("results/05__testing/rel_err_m_d/"+start_distribution+date)

    # Create csv-file, if it doesn't exist
    if not os.path.isfile("results/05__testing/rel_err_m_d/"+start_distribution+date+"/rel_err_m_d.csv"):
        with open(path + '/rel_err_m_d.csv', 'w') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(['n_s [-]', 'err_m_d [%] (P1)', 'err_m_d [%] (P2)', '','c(t=0)', start_distribution])

    # Prepare lists for processing
    rel_err_m_d_P1 = rel_err_m_d[0]
    rel_err_m_d_P2 = rel_err_m_d[1]
    combined_list = [nx, rel_err_m_d_P1, rel_err_m_d_P2]
    export_data = zip_longest(*combined_list, fillvalue='')

    # Save to csv
    with open(path+'/rel_err_m_d.csv', 'a') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerows(export_data)
        f.close()