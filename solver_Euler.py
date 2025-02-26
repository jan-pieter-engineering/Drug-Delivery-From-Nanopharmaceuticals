# -----------------------------------------------------------------
# solver EULER
# -----------------------------------------------------------------
from fenics import *


def solver_Euler(t, dt, D, u_n, c_s_n_Euler, h_mod_n_Euler, f_n_Euler,
                 beta, V, r, ds,
                 element_degree=1,
                 linear_solver='Krylov',
                 abs_tol=1E-5,
                 rel_tol=1E-3,
                 max_iter=1000):

    """ Solve u'= D*Laplace(u) on the line interval [r_i, r_o] with Lagrange elements of specified degree and Robin boundary conditions by using the explicit EULER-method. """

    # Set linear solver parameters
    prm = LinearVariationalSolver.default_parameters()
    if linear_solver == 'Krylov':
        prm["linear_solver"] = 'gmres'
        prm["preconditioner"] = 'ilu'
        prm["krylov_solver"]["absolute_tolerance"] = abs_tol
        prm["krylov_solver"]["relative_tolerance"] = rel_tol
        prm["krylov_solver"]["maximum_iterations"] = max_iter
    else:
        prm["linear_solver"] = 'lu'


    # Initialize u_n
    _u_n_Euler = Function(V)
    _u_n_Euler.assign(u_n)


    # Define trial & test function
    _u_Euler = TrialFunction(V)
    v_Euler = TestFunction(V)

    # Weak Form
    F = (_u_Euler - _u_n_Euler) * v_Euler * r * dx \
      + dt * beta * h_mod_n_Euler * (_u_n_Euler - c_s_n_Euler) * v_Euler * r * ds \
      + dt * beta * D * Dx(_u_n_Euler, 0) * Dx(v_Euler, 0) * r * dx \
      - dt * f_n_Euler * v_Euler * r * dx

    a, L = lhs(F), rhs(F)

    # Compute solution
    _u_Euler = Function(V)
    solve(a == L, _u_Euler, solver_parameters=prm)


    return _u_Euler