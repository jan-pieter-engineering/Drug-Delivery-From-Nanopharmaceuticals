# -----------------------------------------------------------------
# P A R A M E T E R S  (units)
# -----------------------------------------------------------------
import numpy as np

# Geometrieparameter Hohlzylindermatrix
r_i_units                    = 3.2                                             # [mm]     - Innenradius r_i der Hohlzylindermatrix
r_o_units                    = 3.4                                             # [mm]     - Außenradius r_o der Hohlzylindermatrix
H_units                      = 20.0                                            # [mm]     - Höhe H der Hohlzylindermatrix
V_c_units                    = np.pi * (r_o_units**2 - r_i_units**2) * H_units # [mm³]    - Volumen V_c der Hohlzylindermatrix
A_m_t0_units                 = 2.0*np.pi*H_units*(r_i_units+r_o_units)         # [mm²]    - Mantelfläche der Hohlzylindermatrix

# Gentamicin-Parameter
D_GM_exp_units               = 2.6 * 10**-7                                    # [mm²/s]  - erwartbare Groessenordung des Diffusionskoeffizientens
m_d_t0_units                 = 0.0075                                          # [g]      - Gesamtmasse GM, die in die Matrix eingebracht wurde (m_d = 7.5 mg)

# Umgebungslösungs-Parameter
V_s_units                    = 15000.0                                         # [mm³]    - Volumen der Umgebungslösung (V_s = 15 ml)
c_max_units                  = m_d_t0_units/V_s_units                          # [g/mm³]  - maximal mögliche Konzentration innerhalb der Umgebungslösung

#-----------------------------------------------------------------
# U N B E K A N N T E   P A R A M E T E R (hier bereits optimiert)
#-----------------------------------------------------------------
D_expect_constant_units      = 2.6032998653581 * 10**-7                        # [mm²/s]  - Diffusion coefficient
D_expect_parabolic_units     = 2.599321773496793 * 10**-7                        # [mm²/s]  - Diffusion coefficient
alpha_expect_constant_units  = 2.61919530791703 * 10**6                        # [mm³/g]  - transfer coefficient
alpha_expect_parabolic_units = 3.020448481655452 * 10**6                        # [mm³/g]  - transfer coefficient
l_b_expect_constant_units    = 3.19623827676102                                # [mm]     - Diffusionsübergangsbereich
l_b_expect_parabolic_units   = 3.200565154791571                                # [mm]     - Diffusionsübergangsbereich
A_m_degrad_constant_units    = 2.85545640764786 * A_m_t0_units                 # [mm²]    - Mantelfläche der Hohlzylindermatrix nachdem sie aufgebrochen ist
A_m_degrad_parabolic_units   = 4.5053678672349 * A_m_t0_units                  # [mm²]    - Mantelfläche der Hohlzylindermatrix nachdem sie aufgebrochen ist

#-----------------------------------------------------------------
# R E F E R E N Z K O N S T A N T E N
#-----------------------------------------------------------------
t_ref                        = 604800.0                                        # [s]      - = 1 Woche
c_ref                        = m_d_t0_units/V_c_units                          # [g/mm³]  - = c(t=0)
l_ref                        = r_i_units                                       # [mm]     - Innenradius der zylinderförmigen Matrix
D_ref                        = D_GM_exp_units                                  # [mm²/s]  - erwartbare Groessenordung des Diffusionskoeffizientens