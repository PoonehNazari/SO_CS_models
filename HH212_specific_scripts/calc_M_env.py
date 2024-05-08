import numpy as np
import matplotlib.pylab as plt

#Constants
d = 400. #pc
m_H = 1.67e-24
au = 1.496e13 #cm
###################################
###################################
def calc_M_env():
    r_inn = np.linspace(0.11,0.8,100) * d * au
    r_out = np.linspace(0.8,5,100) * d* au
    H_out = 0.8* d* au * (r_out/(0.8 * d* au))**0.221
    H_inn = 0.1* d* au * (r_inn/(0.3 * d* au))**2.124

    n_inn = 5e6 * (r_inn/(400* au))**-1.5
    n_out = 5e6 * (r_out/(400* au))**-1.5

    M_E = 0.
    for i_r in range(len(r_inn)-1):
        M_E_i = 2 * np.pi * r_inn[i_r] * 2 * H_inn[i_r] * n_inn[i_r] * (r_inn[i_r+1] - r_inn[i_r])
        M_E = M_E+M_E_i

    for i_r in range(len(r_out)-1):
        M_E_i = 2 * np.pi * r_out[i_r] * 2 * H_out[i_r] * n_out[i_r] * (r_out[i_r+1] - r_out[i_r])
        M_E = M_E+M_E_i

    M_E = M_E * 1.4 * m_H/1e33
    print(M_E)
###################################
###################################
calc_M_env()



# plt.figure()
# plt.plot(r_inn,H_inn)
# plt.plot(r_out, H_out)
# plt.xscale('log')
# #plt.xlim(0.01,6)
# plt.show()
