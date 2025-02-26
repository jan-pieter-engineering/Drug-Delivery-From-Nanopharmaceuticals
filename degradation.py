# -----------------------------------------------------------------
# D E G R A D A T I O N   M O D E L
# -----------------------------------------------------------------
import numpy as np
from plot_funcs import plot_degrad


def degradation_model(t, u_t0, u_n, const_degrad, r_i, l_c, nx, nr, verbose, ax2, menu):
    """
    If the parabolic distribution of u(t=0) is choosen and the drug diffuses to the middle of the matrix,
    so u_t0<u_n applies, _degrad[i] is set to zero at these positions i.
    """

    # Reset _degrad
    _degrad = [0] * (nr + 1)

    # Calc degradation
    for i in range(nr + 1):

        r_coord = r_i + i * (l_c / nr)

        if u_n(r_coord) > u_t0(r_coord):
            _degrad[i] = 0.0

        else:
            _degrad[i] = 1.0-np.exp(-const_degrad*(((u_t0(r_coord)-u_n(r_coord))/u_t0(r_coord))*((u_t0(r_coord)-u_n(r_coord))/u_t0(r_coord))))

    # Calc minimum of degradation
    _degrad_min_no_zeros = np.min(_degrad[_degrad != 0])
    _degrad_min = np.min(_degrad)

    # Plot degradation
    if (menu == 1 or menu == 2) and verbose == 2:
        plot_degrad(_degrad, nr, r_i, l_c, verbose, ax2, menu)  # activate for debugging


    return _degrad_min_no_zeros
