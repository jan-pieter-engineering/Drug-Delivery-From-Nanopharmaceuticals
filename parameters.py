# -----------------------------------------------------------------
# P A R A M E T E R S  (unitless)
# -----------------------------------------------------------------
from parameters_units import D_expect_constant_units, alpha_expect_constant_units, l_b_expect_constant_units, D_expect_parabolic_units, alpha_expect_parabolic_units, l_b_expect_parabolic_units, r_i_units, r_o_units, H_units, A_m_t0_units, A_m_degrad_constant_units, A_m_degrad_parabolic_units, V_c_units, m_d_t0_units, c_max_units, V_s_units, D_ref, c_ref, l_ref, t_ref
import numpy as np

#-----------------------------------------------------------------
# U N B E K A N N T E   P A R A M E T E R (hier bereits optimiert)
#-----------------------------------------------------------------
D_constant               = D_expect_constant_units / D_ref       # [-]      - Diffusionskoeffizient (D = D_GM_exp/D_ref)
D_parabolic              = D_expect_parabolic_units / D_ref      # [-]      - Diffusionskoeffizient (D = D_GM_exp/D_ref)
alpha_constant           = alpha_expect_constant_units * c_ref   # [-]      - another transfer coefficient
alpha_parabolic          = alpha_expect_parabolic_units * c_ref  # [-]      - another transfer coefficient
l_b_constant             = l_b_expect_constant_units/l_ref       # [-]      - Diffusionsübergangsbereich
l_b_parabolic            = l_b_expect_parabolic_units/l_ref      # [-]      - Diffusionsübergangsbereich
A_m_degrad_constant      = A_m_degrad_constant_units/l_ref**2    # [-]      - Mantelfläche der Hohlzylindermatrix nachdem sie aufgebrochen ist
A_m_degrad_parabolic     = A_m_degrad_parabolic_units/l_ref**2   # [-]      - Mantelfläche der Hohlzylindermatrix nachdem sie aufgebrochen ist
const_degrad_constant    = 1.0                                   # [-]      - degrad = (1-exp(-const_degrad*(1-u_n)²))*100
const_degrad_parabolic   = 1.0                                   # [-]      - degrad = (1-exp(-const_degrad*(1-u_n)²))*100
degrad_break_constant    = 0.356                                 # [-]      - matrix breaks, if value [0,1] is reached
degrad_break_parabolic   = 0.503                                 # [-]      - matrix breaks, if value [0,1] is reached

# Zeitparameter
t                        = 0.0                                   # [-]      - Startzeitpunkt
T                        = 15.0                                  # [-]      - Gesamtversuchsdauer = 15 Wochen (T = 15 W[s]/t_ref = 15 W[s]/1 W[s] = 15)
num_steps                = 0                                     # [-]      - Startwert des Berechnungsschrittzählers

# Gentamicin
m_d_t0                   = m_d_t0_units/(c_ref*l_ref**3)         # [-]      - Gesamtmasse GM, die in die Matrix eingebracht wurde (m_d_units = 7.5 mg)

# Geometrie Hohlzylindermatrix
r_i                      = r_i_units/l_ref                       # [-]      - Innenradius r_i der Hohlzylindermatrix
r_o                      = r_o_units/l_ref                       # [-]      - Außenradius r_o der Hohlzylindermatrix
l_c                      = r_o-r_i                               # [-]      - Wandstaerke des Hohlzylinders
H                        = H_units/l_ref                         # [-]      - Höhe H der Hohlzylindermatrix
V_c                      = V_c_units/l_ref**3                    # [-]      - Volumen V_c der Hohlzylindermatrix
A_m_t0                   = A_m_t0_units/l_ref**2                 # [-]      - Mantelfläche der Hohlzylindermatrix

# Konzentration Hohlzylindermatrix
const_distribution       = m_d_t0/(np.pi*H*(r_o**2 - r_i**2))    # [-]      - u(t=0) = constant = m_d_t0 / (pi*H*r_o² - r_i²)
offset_parabola          = 0.0001                                # [-]      - u(t=0) = parabel * (x[0] - ((r_i + r_o) / 2.0))² + offset
parabolic_distribution   = ((-12.0*(m_d_t0-np.pi*H*offset_parabola*(r_o**2 - r_i**2)))/ # [-]      - Parameter, um parabelfoermige Medikamentenverteilung zu spezifizieren
                           (np.pi*H*(r_i-r_o)**3 * (r_i+r_o)))

# Parameter Umgebungslösung c_s
c_s_t0                   = 0.0                                   # [-]      - Medikamentenkonzentration in der Umgebungsloesung zum Startzeitpunkt t=0
c_max                    = c_max_units/c_ref                     # [-]      - maximal mögliche Konzentration innerhalb der Umgebungslösung
F_15                     = 0.844                                 # [-]      - F(t=15 W) = c_s/c_max = 0.844
c_sat                    = F_15*c_max                            # [-]      - Saettigungsloeslichkeit der Umgebungsloeslichkeit fuer Gentamicin
V_s                      = V_s_units/l_ref**3                    # [-]      - Volumen der Umgebungslösung (V_s_units = 15 ml)

# Weitere Degradationsparameter
t_break                  = 5.0                                   # [-]      - Matrix bricht in Woche 5 (t=5)
t_reached_5              = False                                 # [-]      - wird auf 'True' gesetztm sobald t=5 erreicht wird bzw. die Matrix gebrochen ist
const_degrad_deactivated = 0                                     # [-]      - matrix break deactivated (for optimization before break): degrad = (1-exp(-const_degrad*(1-u_n)²))*100 = 0
limiter_break            = 1                                     # [-]      - if limiter_break != 1: matrix is broken
limiter_sat              = 1                                     # [-]      - if limiter_sat   != 1: "saturation" started

# Hilfsparameter
beta                     = (D_ref*t_ref)/l_ref**2                # [-]      - Hilfskonstante zur besseren Lesbarkeit
gamma                    = D_ref*t_ref*l_ref/V_s_units           # [-]      - Hilfskonstante zur besseren Lesbarkeit (wird für c_s verwendet)

# Elementparameter
element_family           = 'P'
degree                   = [1,2]                                      # [-]      - Grad des gewaehlten Elementtyps
nx = nr                  = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]  # [-]      - Anzahl der Elemente der Intervalldomain
h                        = []                                         # [-]      - h = 1/nx bzw. h = 1/nt

# mesh Parameter
boundary_tol             = 1E-14                                 # [-]

# Expliziter EULER
T_explicit_Euler         = 70 * 10 ** -9
num_steps_explicit_Euler = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]

# adaptiver CN
fineness                 = 2                                     # [-]      - dt_coarse = fineness * dt_fine
runs_CN_high             = fineness                              # [-]      - solver Durchläufe pro Zeitschritt
p                        = 2                                     # [-]      - CN-Wert für Richardson-Extrapolation
tol_adaptiv_dt           = 10**-3                                # [-]      - besserer Wert? --> Bezug zur Messgenauigkeit?
safety                   = 0.9                                   # [-]      - 0 < safety < 1 (Quelle: Fenics hands-on)


# Inverse Analyse
F_num_komplex_n_Euler    = 0.0                                   # [-]      - Startwert: F_num(t=0)=0
kind_of_1d_interpolation = 'linear'                              # [-]      - Wie soll F_exp interpoliert werden?
t_exp                    = [0.,1.,   2.,   3.,   4.,   5.,   7.,   8.,   9.,   11.,  15.]
F_exp                    = [0.,0.399,0.517,0.607,0.623,0.639,0.792,0.805,0.814,0.823,0.844]
t_num                    = [0]                                   # [-]      - erster Wert: t=0
F_num                    = [0]                                   # [-]      - erster Wert: F(t=0) = 0
opt_success              = False                                 # [-]      - Status der Optimierung

# MMS
U_0                      = 2.0                                   # [-]      - Konzentration zur Berechnung der angenommenen Loesung
U_1                      = 1.0                                   # [-]      - Konzentration zur Berechnung der angenommenen Loesung

# Konvergenzanalyse
discretization_scheme      = 0                                   # [-]      - DUMMY-variable [only =! 0, if convergence analysis is active]
steps_convergence_analysis = len(num_steps_explicit_Euler)
L2_error_avg_list          = [[] for _ in range(len(degree))]

# Zeitschrittparameter
dt                         = len(num_steps_explicit_Euler)* [0]
for i in range(len(num_steps_explicit_Euler)):
    dt[i]= T_explicit_Euler/num_steps_explicit_Euler[i]

nt_CN                    = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]  # [-]  - Anzahl der Unterteilungen des untersuchten Zeitraums
dt_CN = len(nt_CN)* [0]                                               # [-]  - Wird nur für die Konvergenzanalyse verwendet
for j in range(len(nt_CN)):
    dt_CN[j] = T/nt_CN[j]