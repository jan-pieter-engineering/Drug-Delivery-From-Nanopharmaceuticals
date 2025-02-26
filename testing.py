# -----------------------------------------------------------------
# T E S T I N G
# -----------------------------------------------------------------
from fenics import interpolate, plot
from mesh_domains_boundaries import define_mesh, define_r
from u_0 import define_u_t0
from utility import check__u_t0
from save import save_rel_err_m_d_to_csv, create__u_t0_save_files, save_u_t0
import matplotlib.pyplot as plt


def run_testing(t, T, r_i, r_o, H, m_d_t0, l_ref, c_ref, menu, const_distribution, parabolic_distribution, U_0, U_1, offset_parabola, element_family, degree, nx, verbose, discretization_scheme, menu_testing, distribution_menu):

    if menu_testing==1:
        print(f'CONVERGENCE ANALYSIS will be moved here')

    elif menu_testing==2: # Relative error of drug mass m_d depending on spatial discretization

        # Calc relative errors of m_d in relation to degree and amount of elements
        rel_err_m_d = calc_rel_err_m_d(t, T, r_i, r_o, H, m_d_t0, l_ref, c_ref, menu, const_distribution, parabolic_distribution, U_0, U_1, offset_parabola, element_family, degree, nx, verbose, discretization_scheme, distribution_menu)

        # Save relative errors
        save_rel_err_m_d_to_csv(nx, rel_err_m_d, degree, distribution_menu)


        return [rel_err_m_d, distribution_menu]



def calc_rel_err_m_d(t, T, r_i, r_o, H, m_d_t0, l_ref, c_ref, menu, const_distribution, parabolic_distribution, U_0, U_1, offset_parabola, element_family, degree, nx, verbose, discretization_scheme, distribution_menu):

    # Initialize list
    rel_err_percent_list = [[], []]

    # Create vtk/xdmf filea depending on distribution and element type
    vtk_file_u_t0__const_P1, xdmf_file_u_t0__const_P1         = create__u_t0_save_files(0, element_family, 1)
    vtk_file_u_t0__const_P2, xdmf_file_u_t0__const_P2         = create__u_t0_save_files(0, element_family, 2)
    vtk_file_u_t0__parabolic_P1, xdmf_file_u_t0__parabolic_P1 = create__u_t0_save_files(1, element_family, 1)
    vtk_file_u_t0__parabolic_P2, xdmf_file_u_t0__parabolic_P2 = create__u_t0_save_files(1, element_family, 2)

    # Iterate through element degrees
    for i in range(len(degree)):

        for j in range(len(nx)):  # Iterate through amount of elements

            # Calc relative error
            u_0 = 0 # reset u_0
            u_0 = define_u_t0(t, T, r_i, r_o, distribution_menu, const_distribution, parabolic_distribution, U_0, U_1, offset_parabola, degree[i], discretization_scheme, menu)
            mesh, _space_dim, V = define_mesh(nx[j], r_i, r_o, degree[i], element_family)
            r = define_r(degree[i])
            u_t0 = interpolate(u_0, V)  # interpolate u_0
            rel_err_percent = check__u_t0(u_t0, m_d_t0, r, H, l_ref, c_ref, verbose, menu)  # check, if distribution-value is chosen correctly, that m_d_units = 7.5g
            rel_err_percent_list[i] += [rel_err_percent]

            # Save u(t=0) to vtk/xdmf file
            if distribution_menu == 0:
                if degree[i] == 1:
                    save_u_t0(nx[j], u_t0, vtk_file_u_t0__const_P1, xdmf_file_u_t0__const_P1)
                elif degree[i] == 2:
                    save_u_t0(nx[j], u_t0, vtk_file_u_t0__const_P2, xdmf_file_u_t0__const_P2)

            elif distribution_menu == 1:
                if degree[i] == 1:
                    save_u_t0(nx[j], u_t0, vtk_file_u_t0__parabolic_P1, xdmf_file_u_t0__parabolic_P1)
                elif degree[i] == 2:
                    save_u_t0(nx[j], u_t0, vtk_file_u_t0__parabolic_P2, xdmf_file_u_t0__parabolic_P2)

            # Plot solution u
            if verbose == 2:
                plot(u_t0, title=f"u(t={t}) - {element_family}{degree[i]}-Element (nx = {nx[j]})")  # show plot
                plt.pause(10)  # hold plot for 0.5 s
                plt.clf()  # close plot


    return rel_err_percent_list