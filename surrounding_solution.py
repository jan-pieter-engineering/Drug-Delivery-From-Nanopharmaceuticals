# -----------------------------------------------------------------
# S U R R O U N D I N G   S O L U T I O N
# -----------------------------------------------------------------
from fenics import *
import numpy as np


def calc_c_s(t, u_n, c_s_n, c_s_t0, h_mod_n, T_W, dt, ds, H, gamma, U_0, U_1, r, r_i, r_o, degree, boundary_tol, discretization_scheme, menu):
    """ Calc the surrounding solution concentration. """
    if menu == 1 or menu == 2:

        _c_s_n = c_s_SIM(t, u_n, c_s_n, c_s_t0, h_mod_n, dt, ds, H, gamma, r, degree)

    elif menu == 4:

        _c_s_n = C_s_MMS(t, T_W, U_0, U_1, r_i, r_o, boundary_tol)


    return _c_s_n



def c_s_SIM(t, u_n, c_s_n, c_s_t0, h_mod_n, dt, ds, H, gamma, r, degree):
    """ Calc the surrounding solution concentration for the simulation. Real parameters are used. """

    if t > 0:

        _c_s_val = c_s_n.val + dt * 2.0 * np.pi * H * gamma * h_mod_n * assemble((u_n - c_s_n.val) * r * ds)

        return Expression("val", degree=degree, val=_c_s_val)

    else:
        return Expression("val", degree=degree, val=c_s_t0)



class C_s_MMS(UserExpression):
    """ Calc the surrounding solution concentration for the convergence analysis. The manufactured solution (MMS) is used. """

    def __init__(self, t, T_W, U_0, U_1, r_i, r_o, boundary_tol, **kwargs):
        super().__init__(**kwargs)
        self.t = t
        self.T_W = T_W
        self.r_i = r_i
        self.r_o = r_o
        self.U_0 = U_0
        self.U_1 = U_1
        self.boundary_tol = boundary_tol


    #TODO: c_s.update(t,...) oder c_s.t=t, c_s.u=u, ... verwenden?
    def update(self, t):
        self.t = t


    def eval(self, values, x):

        if near(x[0], self.r_i, self.boundary_tol):

            values[0] = self.U_0 - self.U_1*self.t*sin(self.r_i)/self.T_W + 1.0


        elif near(x[0], self.r_o, self.boundary_tol):

            values[0] = self.U_0 - self.U_1*self.t*sin(self.r_o)/self.T_W - 1.0

        else:
            values[0] = 0


    def value_shape(self):
        return ()



def update_c_s_n(t, c_s, degree, menu):

    if menu == 1 or menu == 2:

        _c_s_n = Expression("val", degree=degree, val=c_s.val)  # Konzentration der Umgebungslösung

    elif menu == 4:

        c_s.t = t    # update c_s-Expression
        _c_s_n = c_s # pass updated c_s-Expression to internal _c_s_n variable, which will be used for the next time-step


    return _c_s_n



def transfer_c_s_n_from_EULER_to_CN(c_s_n_EULER, degree, menu):

    if menu == 1 or menu == 2:

        _c_s_n_CN = Expression("val", degree=degree, val=c_s_n_EULER.val)  # Konzentration der Umgebungslösung

    elif menu == 4:

        _c_s_n_CN = c_s_n_EULER


    return _c_s_n_CN