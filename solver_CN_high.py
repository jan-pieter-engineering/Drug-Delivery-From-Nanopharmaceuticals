# -----------------------------------------------------------------
# solver CN (high precision/small timestep)
# -----------------------------------------------------------------
from fenics import *
from boundary_conditions import h_modified
from degradation import degradation_model
from source_term import define_source_term
from surrounding_solution import calc_c_s, update_c_s_n
from utility import calc_F


def solver_CN_high(t, dt, u_n, c_s_n_CN_high, h_mod_n_CN_high, f_n_CN_high, F_num_komplex_n_CN_high,
                   u_t0, c_sat, c_max, c_s_t0, V_s_units, V_c_units, r_i, r_o, H, l_c, A_m_t0, A_m_degrad, beta, gamma, U_0, U_1, T, const_degrad, degrad_break, runs_CN_high, degree, nx, nr, V, n_facet, r, ds, boundary_tol, verbose, ax2, discretization_scheme, menu,
                   D, alpha, l_b,
                   element_degree=1,
                   linear_solver='Krylov',
                   abs_tol=1E-5,
                   rel_tol=1E-3,
                   max_iter=10000):

    """ Solve u'= D*Laplace(u) on the line interval [r_i, r_o] with Lagrange elements of specified degree and Robin boundary conditions by using a high precision CRANK-NICOLSON-method. """

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
    _u_n_CN_high = Function(V)
    _u_n_CN_high.assign(u_n)

    # Time stepping
    for n in range(runs_CN_high):

        # Update time
        t += dt

        # Initialize & update t^(n+1) values
        c_s_CN_high = calc_c_s(t, _u_n_CN_high, c_s_n_CN_high, c_s_t0, h_mod_n_CN_high, T, dt, ds, H, gamma, U_0, U_1, r, r_i, r_o, degree, boundary_tol, discretization_scheme, menu)
        degrad_max_CN_high = degradation_model(t, u_t0, _u_n_CN_high, const_degrad, r_i, l_c, nx, nr, verbose, ax2, menu)
        h_mod_CN_high = h_modified(t, D, alpha, l_b, c_s_CN_high, c_sat, c_max, r, A_m_t0, A_m_degrad, degrad_max_CN_high, degrad_break, U_1, T, discretization_scheme, degree, menu)
        f_CN_high = define_source_term(t, T, D, U_0, U_1, beta, degree, discretization_scheme, menu)
        F_num_simpel_CN_high, F_num_komplex_CN_high = calc_F(dt, D, _u_n_CN_high, c_s_CN_high, F_num_komplex_n_CN_high, c_max, V_s_units, V_c_units, H, beta, gamma, r, n_facet, ds, discretization_scheme, menu)

        # Define trial & test function
        _u_CN_high = TrialFunction(V)
        v_CN_high = TestFunction(V)

        # Weak form
        F = (_u_CN_high - _u_n_CN_high) * v_CN_high * r * dx \
          + 0.5 * dt * beta * (h_mod_CN_high * (_u_CN_high - c_s_CN_high) + h_mod_n_CN_high * (_u_n_CN_high - c_s_n_CN_high)) * v_CN_high * r * ds \
          + 0.5 * dt * beta * D * (Dx(_u_CN_high, 0) + Dx(_u_n_CN_high, 0)) * Dx(v_CN_high, 0) * r * dx \
          - 0.5 * dt * (f_CN_high + f_n_CN_high) * v_CN_high * r * dx

        a, L = lhs(F), rhs(F)

        # Compute solution
        _u_CN_high = Function(V)
        solve(a == L, _u_CN_high, solver_parameters=prm)


        # Update t^n values for next timestep
        c_s_n_CN_high = update_c_s_n(t, c_s_CN_high, degree, menu)
        h_mod_n_CN_high = h_mod_CN_high
        f_n_CN_high.t = t
        F_num_komplex_n_CN_high = F_num_komplex_CN_high
        _u_n_CN_high.assign(_u_CN_high)


    return _u_CN_high, c_s_CN_high, degrad_max_CN_high, h_mod_CN_high, F_num_simpel_CN_high, F_num_komplex_CN_high