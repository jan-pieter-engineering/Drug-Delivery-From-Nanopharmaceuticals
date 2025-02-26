# -----------------------------------------------------------------
# P L O T  ( S E T T I N G S  &  I N I T I A L I Z A T I O N )
# -----------------------------------------------------------------
from fenics import *
import matplotlib.pyplot as plt # version: 2.1.1


# ----------------------------------
# S E T T I N G S
# ----------------------------------
def plot_settings(c_max, c_sat, m_d_t0, menu, discretization_scheme, verbose):
    """ Defines which plots are shown depending on user menu & verbose choice. """
    if menu == 1 or menu == 2:

        if verbose == 0: # verbose mode OFF
            _fig1 = 0
            _fig2 = 0
            _ax1  = 0
            _ax2 = [0,0,0,0]
            _ax3 = 0
            _ax4 = 0
            _ax5 = 0
            _ax6 = 0

        elif verbose == 1: # verbose mode OFF
            _fig1, _ax1 = initialize_F_num_plots__simpel_vs_complex()  # initialize F_num plots
            _fig2 = 0
            _ax2 = [0,0,0,0]
            _ax3 = 0
            _fig4, _ax4 = intialize_h_mod_plot()
            _ax5 = 0
            _ax6 = 0

        elif verbose == 2: # verbose mode ON
            _vtk_file = File("results/pvd/u.pvd") # create vtk-file
            _xdmf_file = XDMFFile("results/xdmf/t_u.xdmf") # create XDMF-file
            _fig1, _ax1 = initialize_F_num_plots__simpel_vs_complex() # initialize F_num plots
            _fig2, _ax2 = initialize_plot__degrad() # initialize degradation plot
            _fig3, _ax3 = initialize_c_s_plot(c_max, c_sat)
            _fig4, _ax4 = intialize_h_mod_plot()
            _ax5 = 0
            _fig6, _ax6 = initialize_m_d_plot(m_d_t0)


    elif menu == 4:

        if verbose == 0 or verbose == 1: # verbose mode OFF
            _fig1 = 0
            _ax1 = 0
            _fig2 = 0
            _ax2  = [0,0,0,0]
            _ax3 = 0
            _ax4 = 0
            _ax5 = 0
            _ax6 = 0

        elif verbose == 2: # verbose mode ON
            _fig1 = 0
            _ax1 = 0
            _fig2 = 0
            _ax2  = [0,0,0,0]
            _ax3 = 0
            _ax4 = 0
            _ax5 = initialize_plot__L2_errornorm(discretization_scheme)
            _ax6 = 0


    return _ax1, _ax2, _ax3, _ax4, _ax5, _ax6



# ----------------------------------
# I N I T I A L I Z A T I O N S
# ----------------------------------
def initialize_F_num_plots__simpel_vs_complex():
    """ Initialization of (t_num, F_num_simpel), (t_num, F_num_complex) and (t_exp, F_exp). """

    # plots initialisieren
    fig1, ax1 = plt.subplots(1,2)

    # Achsenbezeichnung
    # x-Achse
    ax1[0].set_xlabel(r"$\tilde{t}$", loc="center")
    ax1[1].set_xlabel(r"$\tilde{t}$", loc="center")

    # y-Achse
    ax1[0].set_ylabel(r"$\tilde{F}_{num}$", loc="top", rotation="horizontal")
    ax1[1].set_ylabel(r"$\tilde{F}_{num}$", loc="top", rotation="horizontal")

    # Titel
    ax1[0].set_title(r"$\tilde{F}_\mathrm{num} = \frac{\tilde{c}_\mathrm{s}}{\tilde{c}_\mathrm{max}}$")
    ax1[1].set_title(r"$\tilde{F}_\mathrm{num}-Normalization$")

    # F_exp interpolieren
    import numpy as np
    from scipy.interpolate import interp1d

    t_exp = [0, 1, 2, 3, 4, 5, 7, 8, 9, 11, 15]
    F_exp = [0., 0.399, 0.517, 0.607, 0.623, 0.639, 0.792, 0.805, 0.814, 0.823, 0.844]
    T = 15.0

    kind_of_1d_interpolation = 'quadratic'  # [-]      - Wie soll F_exp interpoliert werden?
    _F_exp_interp1d = interp1d(t_exp, F_exp, kind_of_1d_interpolation)

    _t_plot = np.linspace(0, T, 100)
    _F_exp_interpolated = _F_exp_interp1d(_t_plot)

    # F_exp und F_exp interpoliert plotten
    ax1[1].plot(_t_plot, _F_exp_interpolated, color='black', label='F_exp (quad interpol)')
    ax1[1].scatter(t_exp, F_exp, color='black', marker='s', label='F_exp')


    return fig1, ax1



def initialize_plot__degrad():
    """ Initialization of degradation plot. """

    # Creating grid for subplots
    _fig2 = plt.figure()
    _fig2.set_figheight(6)
    _fig2.set_figwidth(6)

    # Initialize _ax2 with dummy values
    _ax2 = [0, 0, 0, 0]

    _ax2[0] = plt.subplot2grid(shape=(3, 3), loc=(0, 0), colspan=3, rowspan=2) # degrad
    _ax2[1] = plt.subplot2grid(shape=(3, 3), loc=(2, 0))                       # m_d Euler
    _ax2[2] = plt.subplot2grid(shape=(3, 3), loc=(2, 1))                       # m_d CN-low
    _ax2[3] = plt.subplot2grid(shape=(3, 3), loc=(2, 2))                       # m_d CN-high

    _ax2[0].set_xlabel(r"$\tilde{r}$", loc="right")
    _ax2[1].set_xlabel(r"$\tilde{t}$", loc="right")
    _ax2[2].set_xlabel(r"$\tilde{t}$", loc="right")
    _ax2[3].set_xlabel(r"$\tilde{t}$", loc="right")

    _ax2[0].set_ylabel(r"$degradation$", loc="top", rotation="horizontal")
    _ax2[1].set_ylabel(r"$degradation_\mathrm{max}$", loc="top", rotation="horizontal")
    _ax2[2].set_ylabel(r"$degradation_\mathrm{max}$", loc="top", rotation="horizontal")
    _ax2[3].set_ylabel(r"$degradation_\mathrm{max}$", loc="top", rotation="horizontal")


    return _fig2, _ax2



def initialize_c_s_plot(c_max, c_sat):
    """ Initialization of surrounding solution concentration plot. """

    # Plots initialisieren
    fig3, ax3 = plt.subplots()

    # Achsenbezeichnung
    ax3.set_xlabel(r"$\tilde{t}$", loc="right")
    ax3.set_ylabel(r"$\tilde{c}_{s}$", loc="top", rotation="horizontal")

    # Plot c_max for reference
    plt.axhline(c_max, label=r"$\tilde{c}_{max}$", color='r')

    # Plot c_sat
    plt.axhline(c_sat, label=r"$\tilde{c}_{sat}$", color='y')

    # Show labels
    ax3.legend()


    return fig3, ax3



def intialize_h_mod_plot():
    """ Initialization of mass transfer coefficient plot. """

    # Plots initialisieren
    fig4, ax4 = plt.subplots()

    # Achsenbezeichnung
    ax4.set_xlabel(r"$\tilde{t}$", loc="right")
    ax4.set_ylabel(r"$\tilde{h}_{mod}$", loc="top", rotation="horizontal")


    return fig4, ax4



def initialize_plot__L2_errornorm(discretization_scheme):
    """ Initialization of L2-norm errors between numerical and experimental release curve value F_num/F_exp. """

    # Initialize plots
    fig5, ax5 = plt.subplots()

    # Achsenbezeichnung
    ax5.set_xlabel(r"$\tilde{t}$", loc="center")
    ax5.set_ylabel("L2-Errornorm [-]", loc="top", rotation="horizontal")

    # Titel
    if discretization_scheme == 1:
        ax5.set_title('Explicit Euler')

    elif discretization_scheme == 2:
        ax5.set_title('CN (low precision)')

    elif discretization_scheme == 3:
        ax5.set_title('CN (high precision)')


    return ax5



def initialize_m_d_plot(m_d_t0):
    """ Initialization of total drug mass plot. """

    # Initialize plots
    _fig6, _ax6 = plt.subplots(2,3)

    # Title
    _fig6.suptitle(r"$\tilde{m}_{d, total} = \tilde{m}_{d, m} + \tilde{m}_{d, s}$", fontsize=15)

    # Achsenbezeichnung
    # x-Achse
    _ax6[0,0].set_xlabel(r"$\tilde{t}$", loc="right")
    _ax6[0,1].set_xlabel(r"$\tilde{t}$", loc="right")
    _ax6[0,2].set_xlabel(r"$\tilde{t}$", loc="right")

    _ax6[1,0].set_xlabel(r"$\tilde{t}$", loc="right")
    _ax6[1,1].set_xlabel(r"$\tilde{t}$", loc="right")
    _ax6[1,2].set_xlabel(r"$\tilde{t}$", loc="right")

    # y-Achse
    _ax6[0,0].set_ylabel(r"$\tilde{m}_{d, total}$", loc="top", rotation="horizontal")
    _ax6[0,1].set_ylabel(r"$\tilde{m}_{d, total}$", loc="top", rotation="horizontal")
    _ax6[0,2].set_ylabel(r"$\tilde{m}_{d, total}$", loc="top", rotation="horizontal")

    _ax6[1,0].set_ylabel(r"$\Delta \tilde{m}_{d, total} \; [\%]$", loc="top", rotation="horizontal")
    _ax6[1,1].set_ylabel(r"$\Delta \tilde{m}_{d, total} \; [\%]$", loc="top", rotation="horizontal")
    _ax6[1,2].set_ylabel(r"$\Delta \tilde{m}_{d, total} \; [\%]$", loc="top", rotation="horizontal")

    # Titel
    _ax6[0,0].set_title('Explicit Euler')
    _ax6[0,1].set_title('CN (low precision)')
    _ax6[0,2].set_title('CN (high precision)')

    # Plot m_d line for reference
    _ax6[0,0].axhline(m_d_t0, label=r"$\tilde{m}_{d}(\tilde{t}=0)$", color='y')
    _ax6[0,1].axhline(m_d_t0, label=r"$\tilde{m}_{d}(\tilde{t}=0)$", color='y')
    _ax6[0,2].axhline(m_d_t0, label=r"$\tilde{m}_{d}(\tilde{t}=0)$", color='y')

    # Show legend
    _ax6[0,0].legend()
    _ax6[0,1].legend()
    _ax6[0,2].legend()


    return _fig6, _ax6