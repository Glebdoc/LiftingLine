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
LL_C_CS_8 = np.genfromtxt('LL_C_CS_8.csv', delimiter=',', skip_header=1, names=['AW', 'Ct', 'Cp'], dtype=None, encoding=None)
LL_C_blade_8 = np.genfromtxt('LL_C_blade_8.csv', delimiter=',', skip_header=1, names=['N', 'Ct', 'Cp'], dtype=None, encoding=None)
LL_U_blade_8 = np.genfromtxt('LL_U_blade_8.csv', delimiter=',', skip_header=1, names=['N', 'Ct', 'Cp'], dtype=None, encoding=None)
LL_C_wake_8 = np.genfromtxt('LL_C_wake_8.csv', delimiter=',', skip_header=1, names=['dt', 'Ct', 'Cp'], dtype=None, encoding=None)
LL_U_wake_8 = np.genfromtxt('LL_U_wake_8.csv', delimiter=',', skip_header=1, names=['dt', 'Ct', 'Cp'], dtype=None, encoding=None)
LL_C_Rotation_8 = np.genfromtxt('LL_C_Rotation_8.csv', delimiter=',', skip_header=1, names=['NRotation', 'Ct', 'Cp'], dtype=None, encoding=None)
LL_U_Rotation_8 = np.genfromtxt('LL_U_Rotation_8.csv', delimiter=',', skip_header=1, names=['NRotation', 'Ct', 'Cp'], dtype=None, encoding=None)
LL_C_Rotation_Err_8 = np.genfromtxt('LL_C_Rotation_Err_8.csv', delimiter=',', skip_header=1, names=['NRotation', 'iter', 'err'], dtype=None, encoding=None)
LL_U_Rotation_Err_8 = np.genfromtxt('LL_U_Rotation_Err_8.csv', delimiter=',', skip_header=1, names=['NRotation', 'iter', 'err'], dtype=None, encoding=None)
LL_C_Rotation_Err_all_8 = np.genfromtxt('LL_C_Rotation_Err_all_8.csv', delimiter=',', skip_header=1, names=['iter1', 'iter2', 'iter3', 'err1', 'err2', 'err3'], dtype=None, encoding=None)
LL_U_Rotation_Err_all_8 = np.genfromtxt('LL_U_Rotation_Err_all_8.csv', delimiter=',', skip_header=1, names=['iter1', 'iter2', 'iter3', 'err1', 'err2', 'err3'], dtype=None, encoding=None)

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
    #N=15, Nrotation=10, dt=10
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))
    axes[0].plot(LL_C_CS_8['AW'][:-1], LL_C_CS_8['Ct'][:-1],marker='o', color='r', linewidth=1, label='LL:TSR = 8, Cosine')
    axes[0].plot(LL_U_CS_8['AW'][:-1], LL_U_CS_8['Ct'][:-1],'--',marker='v', color='b', linewidth=1, label='LL:TSR = 8, Uniform')
    #axes[0].scatter(LL_C_CS_8['AW'][-1], LL_C_CS_8['Ct'][-1], marker='x', color='g', s=200, linewidth=1, label=f'LL:TSR = 8, Cosine, $a_w$ = {np.round(LL_C_CS_8["AW"][-1], 2)}')
    #axes[0].scatter(LL_U_CS_8['AW'][-1], LL_U_CS_8['Ct'][-1], marker='x', color='black', s=200, linewidth=1, label=f'LL:TSR = 8, Uniform, $a_w$ = {np.round(LL_U_CS_8["AW"][-1], 2)}')
    axes[0].grid()
    axes[0].set_xlabel('$a_w$')
    axes[0].set_ylabel(r'$C_T$')
    axes[0].legend()

    axes[1].plot(LL_C_CS_8['AW'][:-1], LL_C_CS_8['Cp'][:-1],marker='o', color='r', linewidth=1, label='LL:TSR = 8, Cosine')
    axes[1].plot(LL_U_CS_8['AW'][:-1], LL_U_CS_8['Cp'][:-1], '--',marker='v', color='b', linewidth=1, label='LL:TSR = 8, Uniform')
    #axes[1].scatter(LL_C_CS_8['AW'][-1], LL_C_CS_8['Cp'][-1], marker='x', color='g', s=200, linewidth=1, label=f'LL:TSR = 8, Cosine, $a_w$ = {np.round(LL_C_CS_8["AW"][-1], 2)}')
    #axes[1].scatter(LL_U_CS_8['AW'][-1], LL_U_CS_8['Cp'][-1],  marker='x', color='black', s=200, linewidth=1, label=f'LL:TSR = 8, Uniform, $a_w$ = {np.round(LL_U_CS_8["AW"][-1], 2)}')
    axes[1].grid()
    axes[1].set_xlabel('$a_w$')
    axes[1].set_ylabel(r'$C_P$')
    axes[1].legend()

    if save:
        plt.savefig('./figures/Uniform_vs_Cosine_CS.pdf')
    if plotting:
        plt.show()
    plt.close()


def plot_uniform_vs_cosine(plotting, save):
    #dt=10, Nrotation=10
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))
    axes[0].plot(LL_C_blade_8['N'], LL_C_blade_8['Ct'] ,marker='o', color='r', label='LL:TSR = 8, Cosine')
    axes[0].plot(LL_U_blade_8['N'], LL_U_blade_8['Ct'] ,'--' ,marker='v', color='b', label='LL:TSR = 8, Uniform')
    axes[0].grid()
    axes[0].set_xlabel('$N$')
    axes[0].set_ylabel(r'$C_T$')
    axes[0].legend()

    axes[1].plot(LL_C_blade_8['N'], LL_C_blade_8['Cp'], marker='o', color='r', label='LL:TSR = 8, Cosine')
    axes[1].plot(LL_U_blade_8['N'], LL_U_blade_8['Cp'], '--', marker='v', color='b', label='LL:TSR = 8, Uniform')
    axes[1].grid()
    axes[1].set_xlabel('$N$')
    axes[1].set_ylabel(r'$C_P$')
    axes[1].legend()

    if save:
        plt.savefig('./figures/Uniform_vs_Cosine.pdf')
    if plotting:
        plt.show()
    plt.close()

def plot_uniform_vs_cosine_wake(plotting, save):
    #N=15, Nrotation=10
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))
    axes[0].plot(LL_C_wake_8['dt'], LL_C_wake_8['Ct'] ,marker='o', color='r', label='LL:TSR = 8, Cosine')
    axes[0].plot(LL_U_wake_8['dt'], LL_U_wake_8['Ct'] ,'--' ,marker='v', color='b', label='LL:TSR = 8, Uniform')
    axes[0].grid()
    axes[0].set_xlabel('Number of wake segments per rotation')
    axes[0].set_ylabel(r'$C_T$')
    axes[0].legend()

    axes[1].plot(LL_C_wake_8['dt'], LL_C_wake_8['Cp'], marker='o', color='r', label='LL:TSR = 8, Cosine')
    axes[1].plot(LL_U_wake_8['dt'], LL_U_wake_8['Cp'], '--', marker='v', color='b', label='LL:TSR = 8, Uniform')
    axes[1].grid()
    axes[1].set_xlabel('Number of wake segments per rotation')
    axes[1].set_ylabel(r'$C_P$')
    axes[1].legend()

    if save:
        plt.savefig('./figures/Uniform_vs_Cosine_wake.pdf')
    if plotting:
        plt.show()
    plt.close()

def plot_uniform_vs_cosine_Rotation(plotting, save):
    #N=15, dt10
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))
    axes[0].plot(LL_C_Rotation_8['NRotation'], LL_C_Rotation_8['Ct'] ,marker='o', color='r', label='LL:TSR = 8, Cosine')
    axes[0].plot(LL_U_Rotation_8['NRotation'], LL_U_Rotation_8['Ct'] ,'--' ,marker='v', color='b', label='LL:TSR = 8, Uniform')
    axes[0].grid()
    axes[0].set_xlabel('Number of Rotations')
    axes[0].set_ylabel(r'$C_T$')
    axes[0].legend()

    axes[1].plot(LL_C_Rotation_8['NRotation'], LL_C_Rotation_8['Cp'], marker='o', color='r', label='LL:TSR = 8, Cosine')
    axes[1].plot(LL_U_Rotation_8['NRotation'], LL_U_Rotation_8['Cp'], '--', marker='v', color='b', label='LL:TSR = 8, Uniform')
    axes[1].grid()
    axes[1].set_xlabel('Number of Rotations')
    axes[1].set_ylabel(r'$C_P$')
    axes[1].legend()

    if save:
        plt.savefig('./figures/Uniform_vs_Cosine_Rotation.pdf')
    if plotting:
        plt.show()
    plt.close()

def plot_uniform_vs_cosine_Rotation_Err(plotting, save):
    #N=15, dt10
    plt.figure()
    plt.plot(LL_C_Rotation_Err_8['NRotation'], LL_C_Rotation_Err_8['iter'] ,marker='o', color='r', label='LL:TSR = 8, Cosine')
    plt.plot(LL_U_Rotation_Err_8['NRotation'], LL_U_Rotation_Err_8['iter'] ,'--' ,marker='v', color='b', label='LL:TSR = 8, Uniform')
    plt.grid()
    plt.xlabel('Number of Rotations')
    plt.ylabel(r'Iteration')
    plt.legend()

    if save:
        plt.savefig('./figures/Uniform_vs_Cosine_Rotation_Err.pdf')
    if plotting:
        plt.show()
    plt.close()

def plot_uniform_vs_cosine_Rotation_Err_all(plotting, save):
    # N=15, dt10
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))
    axes[1].plot(LL_U_Rotation_Err_all_8['iter1'], np.log10(LL_U_Rotation_Err_all_8['err1']), color='g',
                 label='LL:TSR = 8, $N_{Rotation}$ = 4')
    axes[1].plot(LL_U_Rotation_Err_all_8['iter2'], np.log10(LL_U_Rotation_Err_all_8['err2']), color='r',
                 label='LL:TSR = 8, $N_{Rotation}$ = 15')
    axes[1].plot(LL_U_Rotation_Err_all_8['iter3'], np.log10(LL_U_Rotation_Err_all_8['err3']), color='b',
                 label='LL:TSR = 8, $N_{Rotation}$ = 25')
    axes[1].grid()
    axes[1].set_xlabel('Iteration')
    axes[1].set_ylabel('Error')
    axes[1].legend()

    axes[0].plot(LL_C_Rotation_Err_all_8['iter1'], np.log10(LL_C_Rotation_Err_all_8['err1']), color='g',
                 label='LL:TSR = 8, $N_{Rotation}$ = 4')
    axes[0].plot(LL_C_Rotation_Err_all_8['iter2'], np.log10(LL_C_Rotation_Err_all_8['err2']), color='r',
                 label='LL:TSR = 8, $N_{Rotation}$ = 15')
    axes[0].plot(LL_C_Rotation_Err_all_8['iter3'], np.log10(LL_C_Rotation_Err_all_8['err3']), color='b',
                 label='LL:TSR = 8, $N_{Rotation}$ = 25')
    axes[0].grid()
    axes[0].set_xlabel('Iteration')
    axes[0].set_ylabel('Error')
    axes[0].legend()

    if save:
        plt.savefig('./figures/Uniform_vs_Cosine_Rotation_Err_all.pdf')
    if plotting:
        plt.show()
    plt.close()


'''plot_alpha(plotting, save)
plot_inflow(plotting, save)
plot_gammas(plotting, save)
plot_Fazim(plotting, save)
plot_Fnorm(plotting, save)'''
#plot_ConvSpeed(plotting, save)
#plot_uniform_vs_cosine(plotting, save)
#plot_uniform_vs_cosine_wake(plotting, save)
#plot_uniform_vs_cosine_Rotation(plotting, save)
#plot_uniform_vs_cosine_Rotation_Err(plotting, save)
#plot_uniform_vs_cosine_Rotation_Err_all(plotting, save)