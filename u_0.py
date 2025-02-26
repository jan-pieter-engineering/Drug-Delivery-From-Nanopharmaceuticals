# -----------------------------------------------------------------
# u (t=0)
# -----------------------------------------------------------------
from fenics import Expression


def define_u_t0(t, T, r_i, r_o, distribution_menu, const_distribution, parabolic_distribution, U_0, U_1, offset_parabola, degree, discretization_scheme, menu):
    """ Define initial drug distribution of the matrix. """

    # Define initial drug distribution with real parameters
    if menu == 1 or menu == 2 or ((menu == 3 or menu == 4) and (discretization_scheme == 6 or discretization_scheme == 7)) or menu == 5:

        # Constant, initial drug distribution
        if distribution_menu == 0:                                                                                   # constant distribution of u(t=0)
            _u_0 = Expression('const_distribution', const_distribution=const_distribution, degree=degree)

        # Parabolic, initial drug distribution
        elif distribution_menu == 1:
            _u_0 = Expression("parabolic_distribution * (x[0] - ((r_i + r_o) / 2.0)) * (x[0] - ((r_i + r_o) / 2.0))" # parabolic distribution of u(t=0)
                              " + offset_parabola",
                              parabolic_distribution=parabolic_distribution, r_i=r_i, r_o=r_o, offset_parabola=offset_parabola, degree=degree)

        return _u_0

    # Define initial drug distribution for the convergence analysis based on the manufactured solution (MMS)
    elif ((menu == 3 or menu == 4) and (discretization_scheme != 6 or discretization_scheme != 7)):

        _u_0 = Expression('U_0 - U_1*t*sin(x[0])/T', U_0=U_0, U_1=U_1, t=t, T=T, degree=degree)

        return _u_0



def define_u_e(t, T, U_0, U_1, degree_u_e):
    """ Define the manufactured solution (MMS), To be used for the convergence analysis. """

    _u_e = Expression('U_0 - U_1*t*sin(x[0])/T', U_0=U_0, U_1=U_1, t=t, T=T, degree=degree_u_e)

    return _u_e