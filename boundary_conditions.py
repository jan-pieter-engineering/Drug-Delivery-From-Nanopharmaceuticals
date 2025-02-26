# -----------------------------------------------------------------
# B O U N D A R Y   C O N D I T I O N S
# -----------------------------------------------------------------
from fenics import Expression


def h_modified(t, D, alpha, l_b, c_s, c_sat, c_max, r, A_m_t0, A_m_degrad, degrad_max, degrad_break, U_1, T, discretization_scheme, degree, menu):

    if menu == 1 or menu == 2:

        if t > 0:
            if alpha * (c_sat - c_s.val) >= 1: # "unsaturated"

                if degrad_max > degrad_break: # matrix breaks
                    print('M A T R I X B R U C H')
                    return (A_m_degrad/A_m_t0)*(D/l_b)

                else:
                    return D/l_b


            elif alpha * (c_sat - c_s.val) < 1: # "saturated"

                if degrad_max > degrad_break: # matrix breaks
                    print('S Ä T T I G U N G   +   M A T R I X B R U C H')
                    return (A_m_degrad/A_m_t0) * (D/l_b) * alpha * (c_sat - c_s.val)

                else:
                    return (D/l_b) * alpha * (c_sat - c_s.val)
                    print('S Ä T T I G U N G')


        else:
            return D/l_b


    elif menu == 4:

        return Expression('D*U_1*t*cos(r)/T', t=t, T=T, D=D, U_1=U_1, r=r, degree=degree)