import numpy as np
import matplotlib.pyplot as plt

plotting = True
save = True
Uinf = 10
R = 50
NB = 3 
#Omega = Uinf*localTSR/R

#print(r'$C_T = $',np.sum(Faxial*chord * NB/ (0.5*1.225*Uinf**2*np.pi*R*R)))

# Read BEM data 
# keys [a, aline, r_R, fnorm, ftan, gamma, phi, AoA, Prandtl, Prandtltip, Prandtlroot, cl_chord, cd_chord, L_chord, D_chord, sigma]
BEM = np.load('BEM_results.npy') #[:, key, tsr]

# Define the expected column names
#header='Gammas,alpha,inflow,Cl,chord, points, FazimND, FaxialND'
column_names = ['GammasND', 'alpha', 'inflow','Cl', 'chord', 'points','FazimND', 'FaxialND']

# Read LL data with specified column names
LL_6 = np.genfromtxt('LL_6.csv', delimiter=',', skip_header=1, names=column_names, dtype=None, encoding=None)
LL_8 = np.genfromtxt('LL_8.csv', delimiter=',', skip_header=1, names=column_names, dtype=None, encoding=None)
LL_10 = np.genfromtxt('LL_10.csv', delimiter=',', skip_header=1, names=column_names, dtype=None, encoding=None)
LL_U_CS_8 = np.genfromtxt('LL_U_CS_8.csv', delimiter=',', skip_header=1, names=['AW', 'Ct', 'Cp'], dtype=None, encoding=None)

# Plot alphas from LL data
def plot_alpha(plotting, save):
    plt.figure()
    plt.plot(LL_6['points'], LL_6['alpha'],marker='o', color='r', label='LL:TSR = 6')
    plt.plot(LL_8['points'], LL_8['alpha'],marker='o', color='g', label='LL:TSR = 8')
    plt.plot(LL_10['points'], LL_10['alpha'],marker='o', color='b',  label='LL:TSR = 10')

    plt.plot(BEM[:,2,0], BEM[:,7,0], '--',marker='v', color='r', label='BEM:TSR = 6')
    plt.plot(BEM[:,2,1], BEM[:,7,1], '--',marker='v', color='g', label='BEM:TSR = 8')
    plt.plot(BEM[:,2,2], BEM[:,7,2], '--',marker='v', color='b', label='BEM:TSR = 10')

    # Adding labels and legend
    plt.grid()
    plt.xlabel('R_r')
    plt.ylabel(r'$\alpha$')
    plt.legend()
    if save:
        plt.savefig('./figures/alpha.pdf')
    if plotting:
        plt.show()
    plt.close()
# Plot inflow angles
def plot_inflow(plotting, save):
    plt.figure()
    plt.plot(LL_6['points'], LL_6['inflow'],marker='o', color='r', label='LL:TSR = 6')
    plt.plot(LL_8['points'], LL_8['inflow'],marker='o', color='g', label='LL:TSR = 8')
    plt.plot(LL_10['points'], LL_10['inflow'],marker='o', color='b',  label='LL:TSR = 10')

    plt.plot(BEM[:,2,0], BEM[:,6,0], '--',marker='v', color='r', label='BEM:TSR = 6')
    plt.plot(BEM[:,2,1], BEM[:,6,1], '--',marker='v', color='g', label='BEM:TSR = 8')
    plt.plot(BEM[:,2,2], BEM[:,6,2], '--',marker='v', color='b', label='BEM:TSR = 10')

    # Adding labels and legend
    plt.grid()
    plt.xlabel('R_r')
    plt.ylabel(r'$\phi$')
    plt.legend()
    if save:
        plt.savefig('./figures/phi.pdf')
    if plotting:
        plt.show()
    plt.close()
# Plot normalized gammas 
def plot_gammas(plotting, save):
    plt.figure()
    plt.plot(LL_6['points'], LL_6['GammasND'],marker='o', color='r', label='LL:TSR = 6')
    plt.plot(LL_8['points'], LL_8['GammasND'],marker='o', color='g', label='LL:TSR = 8')
    plt.plot(LL_10['points'], LL_10['GammasND'],marker='o', color='b',  label='LL:TSR = 10')

    plt.plot(BEM[:,2,0], BEM[:,5,0]/(Uinf*Uinf*np.pi/(NB*(Uinf*6/R))), '--',marker='v', color='r', label='BEM:TSR = 6')
    plt.plot(BEM[:,2,1], BEM[:,5,1]/(Uinf*Uinf*np.pi/(NB*(Uinf*8/R))), '--',marker='v', color='g', label='BEM:TSR = 8')
    plt.plot(BEM[:,2,2], BEM[:,5,2]/(Uinf*Uinf*np.pi/(NB*(Uinf*10/R))), '--',marker='v', color='b', label='BEM:TSR = 10')

    # Adding labels and legend
    plt.grid()
    plt.xlabel('R_r')
    plt.ylabel(r'$\Gamma$')
    plt.legend()
    if save:
        plt.savefig('./figures/Gammas.pdf')
    if plotting:
        plt.show()
    plt.close()
#Plot FazimND
def plot_Fazim(plotting, save):
    plt.figure()
    plt.plot(LL_6['points'], LL_6['FazimND'],marker='o', color='r', label='LL:TSR = 6')
    plt.plot(LL_8['points'], LL_8['FazimND'],marker='o', color='g', label='LL:TSR = 8')
    plt.plot(LL_10['points'], LL_10['FazimND'],marker='o', color='b',  label='LL:TSR = 10')

    plt.plot(BEM[:,2,0], BEM[:,4,0]/(0.5*1.*Uinf**2*R), '--',marker='v', color='r', label='BEM:TSR = 6')
    plt.plot(BEM[:,2,1], BEM[:,4,1]/(0.5*1.*Uinf**2*R), '--',marker='v', color='g', label='BEM:TSR = 8')
    plt.plot(BEM[:,2,2], BEM[:,4,2]/(0.5*1.*Uinf**2*R), '--',marker='v', color='b', label='BEM:TSR = 10')

    plt.grid()
    plt.xlabel('R_r')
    plt.ylabel(r'$F_{azim}$')
    plt.legend()
    if save:
        plt.savefig('./figures/Fazim.pdf')
    if plotting:
        plt.show()
    plt.close()
#Plot FaxialND
def plot_Fnorm(plotting, save):
    plt.figure()
    plt.plot(LL_6['points'], LL_6['FaxialND'],marker='o', color='r', label='LL:TSR = 6')
    plt.plot(LL_8['points'], LL_8['FaxialND'],marker='o', color='g', label='LL:TSR = 8')
    plt.plot(LL_10['points'], LL_10['FaxialND'],marker='o', color='b',  label='LL:TSR = 10')

    plt.plot(BEM[:,2,0], BEM[:,3,0]/(0.5*1*Uinf**2*R), '--',marker='v', color='r', label='BEM:TSR = 6')
    plt.plot(BEM[:,2,1], BEM[:,3,1]/(0.5*1*Uinf**2*R), '--',marker='v', color='g', label='BEM:TSR = 8')
    plt.plot(BEM[:,2,2], BEM[:,3,2]/(0.5*1*Uinf**2*R), '--',marker='v', color='b', label='BEM:TSR = 10')

    plt.grid()
    plt.xlabel('R_r')
    plt.ylabel(r'$F_{axial}$')
    plt.legend()
    if save:
        plt.savefig('./figures/Fnorm.pdf')
    if plotting:
        plt.show()
    plt.close()

def plot_ConvSpeed(plotting, save):
    plt.figure()
    plt.plot(LL_U_CS_8['AW'], LL_U_CS_8['Ct'],marker='o', color='r', label='LL:TSR = 8, Ct')
    plt.plot(LL_U_CS_8['AW'], LL_U_CS_8['Cp'], marker='o', color='b', label='LL:TSR = 8, Cp')

    plt.grid()
    plt.xlabel('$a_w$')
    plt.ylabel(r'$C_P, C_T$')
    plt.legend()
    if save:
        plt.savefig('./figures/ConvectionSpeed.pdf')
    if plotting:
        plt.show()
    plt.close()

'''plot_alpha(plotting, save)
plot_inflow(plotting, save)
plot_gammas(plotting, save)
plot_Fazim(plotting, save)
plot_Fnorm(plotting, save)'''
plot_ConvSpeed(plotting, False)