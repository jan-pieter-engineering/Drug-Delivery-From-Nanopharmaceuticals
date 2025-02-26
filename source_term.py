# -----------------------------------------------------------------
# S O U R C E   T E R M
# -----------------------------------------------------------------
from fenics import *

def define_source_term(t, T, D, U_0, U_1, beta, degree, discretization_scheme, menu):
    """ Source term, that is only !=0, if the convergence analysis is acitve. Then the source term is used as a forcing function. """

    # Source term is set to zero, because diffusion is decribed by a homogeneous diffusion equation
    if menu == 1 or menu == 2 or ((menu == 3 or menu == 4) and (discretization_scheme == 6 or discretization_scheme == 7)):

        return Expression('0', degree=degree)

    # Source term is used as a forcing function for the convergence analysis
    elif ((menu == 3 or menu == 4) and (discretization_scheme != 6 or discretization_scheme != 7)):

        return Expression('-U_1 * sin(x[0]) / T + beta*D*U_1*t*((cos(x[0])/x[0])-sin(x[0]))/T', t=t, T=T, U_0=U_0, U_1=U_1, beta=beta, D=D, degree=2)