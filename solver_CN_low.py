# -----------------------------------------------------------------
# solver CN (low precision/big timestep)
# -----------------------------------------------------------------
from fenics import *
from boundary_conditions import h_modified
from degradation import degradation_model
from source_term import define_source_term
from surrounding_solution import calc_c_s
from utility import calc_F


def solver_CN_low(t, dt, u_n, c_s_n_CN_low, h_mod_n_CN_low, f_n_CN_low, F_num_komplex_n_CN_low,
                  u_t0, c_sat, c_max, c_s_t0, V_s_units, V_c_units, r_i, r_o, H, l_c, A_m_t0, A_m_degrad, beta, gamma, U_0, U_1, T, const_degrad, degrad_break, degree, nx, nr, V, n_facet, r, ds, boundary_tol, verbose, ax2, discretization_scheme, menu,
                  D, alpha, l_b,
                  element_degree=1,
                  linear_solver='Krylov',
                  abs_tol=1E-5,
                  rel_tol=1E-3,
                  max_iter=10000):

    """ Solve u'= D*Laplace(u) on the line interval [r_i, r_o] with Lagrange elements of specified degree and Robin boundary conditions by using a low precision CRANK-NICOLSON-method. """

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
    _u_n_CN_low= Function(V)
    _u_n_CN_low.assign(u_n)

    # Update time
    t += dt

    # Initialize & update t^(n+1) variables
    c_s_CN_low = calc_c_s(t, _u_n_CN_low, c_s_n_CN_low, c_s_t0, h_mod_n_CN_low, T, dt, ds, H, gamma, U_0, U_1, r, r_i, r_o, degree, boundary_tol, discretization_scheme, menu) # Konzentration der Umgebungslösung
    degrad_max_CN_low = degradation_model(t, u_t0, _u_n_CN_low, const_degrad, r_i, l_c, nx, nr, verbose, ax2, menu)
    h_mod_CN_low = h_modified(t, D, alpha, l_b, c_s_CN_low, c_sat, c_max, r, A_m_t0, A_m_degrad, degrad_max_CN_low, degrad_break, U_1, T, discretization_scheme, degree, menu)
    f_CN_low = define_source_term(t, T, D, U_0, U_1, beta, degree, discretization_scheme, menu)


    # Define trial & test function
    _u_CN_low = TrialFunction(V)
    v_CN_low = TestFunction(V)

    # Weak form
    F = r * (_u_CN_low - _u_n_CN_low) * v_CN_low * dx \
      + 0.5 * beta * dt * D * r * (Dx(_u_CN_low, 0) + Dx(_u_n_CN_low, 0)) * Dx(v_CN_low, 0) * dx \
      + 0.5 * beta * dt * (h_mod_CN_low * (_u_CN_low - c_s_CN_low) + h_mod_n_CN_low * (_u_n_CN_low - c_s_n_CN_low)) * v_CN_low * r * ds \
      - 0.5 * dt * (f_CN_low + f_n_CN_low) * v_CN_low * r * dx

    a, L = lhs(F), rhs(F)

    # Compute solution
    _u_CN_low = Function(V)
    solve(a == L, _u_CN_low, solver_parameters=prm)


    # Calc F_num
    # TODO: 2.) ueberpruefen, ob _F_num_komplex_CN_low mit u oder u_n berechnet?
    F_num_simpel_CN_low, F_num_komplex_CN_low = calc_F(dt, D, _u_n_CN_low, c_s_CN_low, F_num_komplex_n_CN_low, c_max, V_s_units, V_c_units, H, beta, gamma, r, n_facet, ds, discretization_scheme, menu)


    return _u_CN_low, c_s_CN_low, degrad_max_CN_low, h_mod_CN_low, F_num_simpel_CN_low, F_num_komplex_CN_low