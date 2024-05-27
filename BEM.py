# import necessary libraries
import matplotlib.pyplot as plt
import numpy as np  
import pandas as pd


data = pd.read_excel('polars.xlsx')
polar_cl = data['Cl'].to_numpy()
polar_cd = data['Cd'].to_numpy()
polar_alpha = data['Alfa'].to_numpy()
  
    
def ainduction(CT):
    a = np.zeros(np.shape(CT))
    CT1=1.816
    CT2=2*np.sqrt(CT1)-CT1
    a[CT>=CT2] = 1 + (CT[CT>=CT2]-CT1)/(4*(np.sqrt(CT1)-1))
    a[CT<CT2] = 0.5-0.5*np.sqrt(1-CT[CT<CT2])
    return a


def PrandtlTipRootCorrection(r_R, rootradius_R, tipradius_R, TSR, NBlades, axial_induction):
    temp1 = -NBlades/2*(tipradius_R-r_R)/r_R*np.sqrt( 1+ ((TSR*r_R)**2)/((1-axial_induction)**2))
    Ftip = np.array(2/np.pi*np.arccos(np.exp(temp1)))
    Ftip[np.isnan(Ftip)] = 0
    temp1 = NBlades/2*(rootradius_R-r_R)/r_R*np.sqrt( 1+ ((TSR*r_R)**2)/((1-axial_induction)**2))
    Froot = np.array(2/np.pi*np.arccos(np.exp(temp1)))
    Froot[np.isnan(Froot)] = 0
    return Froot*Ftip, Ftip, Froot


# define function to determine load in the blade element
def loadBladeElement(vnorm, vtan, r_R, chord, twist, polar_alpha, polar_cl, polar_cd):
    """
    calculates the load in the blade element
    """
    vmag2 = vnorm**2 + vtan**2
    inflowangle = np.arctan2(vnorm,vtan)
    alpha = twist + inflowangle*180/np.pi
    cl = np.interp(alpha, polar_alpha, polar_cl)
    cd = np.interp(alpha, polar_alpha, polar_cd)
    lift = 0.5*vmag2*cl*chord
    drag = 0.5*vmag2*cd*chord
    fnorm = lift*np.cos(inflowangle)+drag*np.sin(inflowangle)
    ftan = lift*np.sin(inflowangle)-drag*np.cos(inflowangle)
    gamma = 0.5*np.sqrt(vmag2)*cl*chord
    return fnorm , ftan, gamma, inflowangle*180/np.pi, alpha, cl, cd, lift, drag

def solveStreamtube(Uinf, r1_R, r2_R, rootradius_R, tipradius_R , Omega, Radius, NBlades, chord, twist, polar_alpha, polar_cl, polar_cd ):

    Area = np.pi*((r2_R*Radius)**2-(r1_R*Radius)**2) #  area streamtube
    r_R = (r1_R+r2_R)/2 # centroide

    a = 0.0 # axial induction
    aline = 0.0 # tangential induction factor
    sigma = NBlades / (2 * np.pi * r_R)
    
    Niterations = 200
    Erroriterations =0.00001 # error limit for iteration rpocess, in absolute value of induction
    
    for i in range(Niterations):
        Urotor = Uinf*(1-a) # axial velocity at rotor
        Utan = (1+aline)*Omega*r_R*Radius # tangential velocity at rotor
        # calculate loads in blade segment in 2D (N/m)
        fnorm, ftan, gamma, phi, AoA, cl_chord, cd_chord, L_chord, D_chord = loadBladeElement(Urotor, Utan, r_R,chord, twist, polar_alpha, polar_cl, polar_cd)
        load3Daxial =fnorm*Radius*(r2_R-r1_R)*NBlades # 3D force in axial direction

        CT = load3Daxial/(0.5*Area*Uinf**2)

        anew =  ainduction(CT)
        
        # correct new axial induction with Prandtl's correction
        Prandtl, Prandtltip, Prandtlroot = PrandtlTipRootCorrection(r_R, rootradius_R, tipradius_R, Omega*Radius/Uinf, NBlades, anew);
        if (Prandtl < 0.00001):
            Prandtl = 0.00001 # avoid divide by zero
        anew = anew/Prandtl # correct estimate of axial induction
        a = 0.75*a+0.25*anew # for improving convergence, weigh current and previous iteration of axial induction
        if (a>0.95):
            a=0.95

        aline = ftan*NBlades/(2*np.pi*Uinf*(1-a)*Omega*2*(r_R*Radius)**2)
        aline = aline/Prandtl
        
        if (np.abs(a-anew) < Erroriterations): 
            # print("iterations")
            # print(i)
            break

    return [a , aline, r_R, fnorm , ftan, gamma, phi, AoA, Prandtl, Prandtltip, Prandtlroot, cl_chord, cd_chord, L_chord, D_chord, sigma]


def solveStreamtube_wP(Uinf, r1_R, r2_R, rootradius_R, tipradius_R, Omega, Radius, NBlades, chord, twist, polar_alpha,
                    polar_cl, polar_cd):
    Area = np.pi * ((r2_R * Radius) ** 2 - (r1_R * Radius) ** 2)  # area streamtube
    r_R = (r1_R + r2_R) / 2  # centroide

    a = 0.0  # axial induction
    aline = 0.0  # tangential induction factor
    sigma = NBlades / (2 * np.pi * r_R)

    Niterations = 200
    Erroriterations = 0.00001  # error limit for iteration rpocess, in absolute value of induction

    for i in range(Niterations):
        Urotor = Uinf * (1 - a)  # axial velocity at rotor
        Utan = (1 + aline) * Omega * r_R * Radius  # tangential velocity at rotor
        # calculate loads in blade segment in 2D (N/m)
        fnorm, ftan, gamma, phi, AoA, cl_chord, cd_chord, L_chord, D_chord = loadBladeElement(Urotor, Utan, r_R, chord,
                                                                                              twist, polar_alpha,
                                                                                              polar_cl, polar_cd)
        load3Daxial = fnorm * Radius * (r2_R - r1_R) * NBlades  # 3D force in axial direction

        CT = load3Daxial / (0.5 * Area * Uinf ** 2)

        anew = ainduction(CT)

        # correct new axial induction with Prandtl's correction
        Prandtl, Prandtltip, Prandtlroot = PrandtlTipRootCorrection(r_R, rootradius_R, tipradius_R,
                                                                    Omega * Radius / Uinf, NBlades, anew);
        # if (Prandtl < 0.00001):
        #    Prandtl = 0.00001 # avoid divide by zero
        # anew = anew/Prandtl # correct estimate of axial induction
        a = 0.75*a+0.25*anew # for improving convergence, weigh current and previous iteration of axial induction
        if (a > 0.95):
            a = 0.95

        aline = ftan * NBlades / (2 * np.pi * Uinf * (1 - a) * Omega * 2 * (r_R * Radius) ** 2)
        # aline = aline/Prandtl

        if (np.abs(a - anew) < Erroriterations):
            break

    return [a, aline, r_R, fnorm, ftan, gamma, phi, AoA, Prandtl, Prandtltip, Prandtlroot, cl_chord, cd_chord, L_chord,
            D_chord, sigma]


def BEM():
    for j in range(len(TSR)):
        Omega = Uinf*TSR[j]/Radius

        results =np.zeros([len(r_R)-1,16])

        for i in range(len(r_R)-1):
            chord = np.interp((r_R[i]+r_R[i+1])/2, r_R, chord_distribution)
            twist = np.interp((r_R[i]+r_R[i+1])/2, r_R, twist_distribution)
            results[i,:] = solveStreamtube(Uinf, r_R[i], r_R[i+1], RootLocation_R, TipLocation_R , Omega, Radius, NBlades, chord, twist, polar_alpha, polar_cl, polar_cd)
        final_results[:,:,j] = results[:,:]
        dr = (r_R[1:]-r_R[:-1])*Radius
        Fnorm = final_results[:,3,j]
        Ftan = final_results[:,4,j]
        T_total = np.sum(0.5*1.225*Uinf*Uinf*NBlades*chords*dr*Fnorm)
        Q_total = np.sum(0.5*1.225*Uinf*Uinf*NBlades*chords*final_results[:,2,j]*Radius*Ftan*dr)
        Q.append(Q_total)
        T.append(T_total)
        CP.append(np.sum(dr*final_results[:,4, j]*final_results[:,2, j]*NBlades*Radius*Omega/(0.5*Uinf**3*np.pi*Radius**2)))
        CT.append(np.sum(dr*results[:,3]*NBlades/(0.5*Uinf**2*np.pi*Radius**2)))
    np.save(f'BEM_results.npy', final_results)
        # if mycos == False:
        #     CT_convergence.append(np.sum(dr*results[:,3]*NBlades/(0.5*Uinf**2*np.pi*Radius**2)))
        #     NumberOfPoints.append(len(r_R))
        # else:
        #     CT_convergence_cos.append(np.sum(dr*results[:,3]*NBlades/(0.5*Uinf**2*np.pi*Radius**2)))

def BEM_wP():
    for j in range(len(TSR)):
        Omega = Uinf*TSR[j]/Radius

        results_wP =np.zeros([len(r_R)-1,16])

        for i in range(len(r_R)-1):
            chord = np.interp((r_R[i]+r_R[i+1])/2, r_R, chord_distribution)
            twist = np.interp((r_R[i]+r_R[i+1])/2, r_R, twist_distribution)
            results_wP[i,:] = solveStreamtube_wP(Uinf, r_R[i], r_R[i+1], RootLocation_R, TipLocation_R , Omega, Radius, NBlades, chord, twist, polar_alpha, polar_cl, polar_cd)
        final_results_wP[:,:,j] = results_wP[:,:]
   

Uinf = 10 # unperturbed wind speed in m/s
NBlades = 3
TipLocation_R =  1
RootLocation_R =  0.2
Radius = 50
pitch = 2 # degrees

delta_r_R = 0.01
r_R = np.arange(0.2, 1+delta_r_R/2, delta_r_R)
cosine_spacing = np.cos(np.linspace(np.arccos(RootLocation_R), 0, len(r_R)))
chord_distribution = 3*(1-r_R)+1 # meters
chords = np.delete(chord_distribution, -1)
twist_distribution = -14*(1-r_R)+pitch # degrees

TSR = np.array([6,8,10])
CP = []
CP_rR = []
CT = []
T = []
Q = []
final_results = np.zeros([len(r_R)-1,16,len(TSR)])
BEM()
final_results_wP = np.zeros([len(r_R)-1,16,len(TSR)])
BEM_wP()


#################################################################

def pressure(final_results, TSR):
    p_static = 101325
    # Creating subplots
    fig, axs = plt.subplots(2, 2, figsize=(10, 8))
    for i in range(len(TSR)):
        Omega = Uinf*TSR[i]/Radius
        rR = final_results[:,2,i]
        CT_perAnulli = final_results[:,3,i]*NBlades/(0.5*1.225*(Uinf**2)*2*np.pi*rR*Radius)

        a = 0.5*(1-np.sqrt(1-CT_perAnulli))
        Urotor = Uinf*(1-a) # axial velocity at rotor
        VR = np.sqrt(Urotor**2)
        U4rotor = Uinf*(1-2*a)
        V4 = np.sqrt(U4rotor**2)

        x = rR
        y1 = (p_static + 0.5*1.225*Uinf**2)*1e-5
        y2 = (p_static + 0.5*1.225*Uinf**2)*1e-5
        y4 = (p_static + 0.5*1.225*(V4**2))*1e-5
        y3 = y4
                # Plotting data on subplots
        axs[0, 0].plot([0,1], [y1,y1])
        axs[0, 0].set_title('Point 1')
        axs[0, 0].grid(True)  # Add grid
        axs[0, 0].set_ylim(0.98, 1.2)  # Set y-axis limits
        axs[0, 0].set_xlabel(r'$r/R$', fontsize=12)
        axs[0, 0].set_ylabel(r'$P_{stag}$ [bar]', fontsize=12)

        axs[0, 1].plot([0,1], [y2,y2])
        axs[0, 1].set_title('Point 2')
        axs[0, 1].grid(True)  # Add grid
        axs[0, 1].set_ylim(0.98, 1.2)  # Set y-axis limits
        axs[0, 1].set_xlabel(r'$r/R$', fontsize=12)
        axs[0, 1].set_ylabel(r'$P_{stag}$ [bar]', fontsize=12)

        axs[1, 0].plot(x, y3)
        axs[1, 0].set_title('Point 3')
        axs[1, 0].grid(True)  # Add grid
        axs[1, 0].set_ylim(1.0132, 1.014)  # Set y-axis limits
        axs[1, 0].set_xlabel(r'$r/R$', fontsize=12)
        axs[1, 0].set_ylabel(r'$P_{stag}$ [bar]', fontsize=12)

        axs[1, 1].plot(x, y4)
        axs[1, 1].set_title('Point 4')
        axs[1, 1].grid(True)  # Add grid
        axs[1, 1].set_ylim(1.0132, 1.014)  # Set y-axis limits
        axs[1, 1].set_xlabel(r'$r/R$', fontsize=12)
        axs[1, 1].set_ylabel(r'$P_{stag}$ [bar]', fontsize=12)

        plt.tight_layout()
    plt.show()
    # plt.show()
    # plt.cla()    
def plot_alpha_rR(final_results, TSR, save=False):
    for i in range(len(TSR)):
        plt.plot(final_results[:,2,i], final_results[:,7,i], label = 'TSR = '+str(TSR[i]))
    plt.ylabel(r'$\alpha$')
    plt.xlabel('r/R')
    plt.grid()
    plt.legend()
    if save:
        plt.savefig('alpha_rR.png')
        plt.cla()
    else:
        plt.show()
def plot_alpha_rR_Prandtl(final_results, final_results_wP, save=False):
    plt.plot(final_results[:,2,1], final_results[:,7,1], label = 'With tip and loss correction')
    plt.plot(final_results_wP[:, 2, 1], final_results_wP[:, 7, 1], label='Without tip and loss correction')
    plt.ylabel(r'$\alpha$ [deg]')
    plt.xlabel('r/R')
    plt.grid()
    plt.legend()
    if save:
        plt.savefig('alpha_rR_Prandtl.png')
        plt.cla()
    else:
        plt.show()
def plot_phi_rR(final_results, TSR, save=False):
    for i in range(len(TSR)):
        plt.plot(final_results[:,2,i], final_results[:,6,i], label = 'TSR = '+str(TSR[i]))
    plt.ylabel(r'$\phi$')
    plt.xlabel('r/R')
    plt.grid()
    plt.legend()
    if save:
        plt.savefig('phi_rR.png')
        plt.cla()
    else:
        plt.show()
def plot_phi_rR_Prandtl(final_results, final_results_wP, save=False):
    plt.plot(final_results[:,2,1], final_results[:,6,1], label = 'With tip and loss correction')
    plt.plot(final_results_wP[:, 2, 1], final_results_wP[:, 6, 1], label='Without tip and loss correction')
    plt.ylabel(r'$\phi$ [deg]')
    plt.xlabel('r/R')
    plt.grid()
    plt.legend()
    if save:
        plt.savefig('phi_rR_Prandtl.png')
        plt.cla()
    else:
        plt.show()
def plot_a_rR(final_results, TSR, save=False):
    for i in range(len(TSR)):
        plt.plot(final_results[:,2,i], final_results[:,0,i], label = 'TSR = '+str(TSR[i]))
    plt.ylabel(r'$a$')
    plt.xlabel('r/R')
    plt.grid()
    plt.legend()
    if save:
        plt.savefig('a_rR.png')
        plt.cla()
    else:
        plt.show()
def plot_a_rR_Prandtl(final_results, final_results_wP, save=False):
    plt.plot(final_results[:,2,1], final_results[:,0,1], label = 'With tip and loss correction')
    plt.plot(final_results_wP[:, 2, 1], final_results_wP[:, 0, 1], label='Without tip and loss correction')
    plt.ylabel(r'$a$')
    plt.xlabel('r/R')
    plt.grid()
    plt.legend()
    if save:
        plt.savefig('a_rR_Prandtl.png')
        plt.cla()
    else:
        plt.show()
def plot_aprime_rR(final_results, TSR, save=False):
    for i in range(len(TSR)):
        plt.plot(final_results[:,2,i], final_results[:,1,i], label = 'TSR = '+str(TSR[i]))
    plt.ylabel(r'$a\prime$')
    plt.xlabel('r/R')
    plt.grid()
    plt.legend()
    if save:
        plt.savefig('aprime_rR.png')
        plt.cla()
    else:
        plt.show()
def plot_fnorm_rR(final_results, TSR, save=False):
    for i in range(len(TSR)):
        plt.plot(final_results[:,2,i], final_results[:,3,i]/(0.5*Uinf**2*Radius), label = 'TSR = '+str(TSR[i]))
    plt.ylabel(r'$C_{norm}$')
    plt.xlabel('r/R')
    plt.grid()
    plt.legend()
    if save:
        plt.savefig('fnorm_rR.png')
        plt.cla()
    else:
        plt.show()
def plot_ftan_rR(final_results, TSR, save=False):
    for i in range(len(TSR)):
        plt.plot(final_results[:,2,i], final_results[:,4,i]/(0.5*Uinf**2*Radius), label = 'TSR = '+str(TSR[i]))
    plt.ylabel(r'$C_{tan}$')
    plt.xlabel('r/R')
    plt.grid()
    plt.legend()
    if save:
        plt.savefig('ftan_rR.png')
        plt.cla()
    else:
        plt.show()
def plot_ftan_rR_Prandtl(final_results, final_results_wP, save=False):
    plt.plot(final_results[:,2,1], final_results[:,4,1]/(0.5*Uinf**2*Radius), label = 'With tip and loss correction')
    plt.plot(final_results_wP[:, 2, 1], final_results_wP[:, 4, 1] / (0.5 * Uinf ** 2 * Radius), label='Without tip and loss correction')
    plt.ylabel(r'$C_{tan}$')
    plt.xlabel('r/R')
    plt.grid()
    plt.legend()
    if save:
        plt.savefig('ftan_rR_Prandtl.png')
        plt.cla()
    else:
        plt.show()
def plot_tiprootloss(final_results, TSR, save:False):
    for j in range(3):
        for i in range(len(TSR)):
            plt.plot(final_results[:,2,i], final_results[:,8+j,i], label='TSR = '+str(TSR[i]))
        plt.xlabel('r/R')
        plt.grid()
        plt.legend()
        if j == 0:
            plt.ylabel(r'$f(\mu)$ [-]')
            if save:
                plt.savefig('Prandtl.png')
                plt.cla()
            else:
               plt.show()
        if j == 1:
            plt.ylabel(r'$f(\mu)_{tip}$ [-]')
            if save:
                plt.savefig('Prandtl_tip.png')
                plt.cla()
            else:
               plt.show()
        if j == 2:
            plt.ylabel(r'$f(\mu)_{root}$ [-]')
            if save:
                plt.savefig('Prandtl_root.png')
                plt.cla()
            else:
               plt.show()
def plot_cl_chord(final_results, TSR, save=False):
    for i in range(len(TSR)):
        plt.plot(final_results[:,2,i], final_results[:,11,i], label = 'TSR = '+str(TSR[i]))
    plt.ylabel(r'$Cl$')
    plt.xlabel('r/R')
    plt.grid()
    plt.legend()
    if save:
        plt.savefig('Cl_chord.png')
        plt.cla()
    else:
        plt.show()
def plot_cd_chord(final_results, TSR, save=False):
    for i in range(len(TSR)):
        plt.plot(final_results[:,2,i], final_results[:,12,i], label = 'TSR = '+str(TSR[i]))
    plt.ylabel(r'$Cd$')
    plt.xlabel('r/R')
    plt.grid()
    plt.legend()
    if save:
        plt.savefig('Cd_chord.png')
        plt.cla()
    else:
        plt.show()
def plot_L_chord(final_results, TSR, save=False):
    for i in range(len(TSR)):
        plt.plot(final_results[:,2,i], final_results[:,13,i], label = 'TSR = '+str(TSR[i]))
    plt.ylabel(r'$L$')
    plt.xlabel('r/R')
    plt.grid()
    plt.legend()
    if save:
        plt.savefig('L_chord.png')
        plt.cla()
    else:
        plt.show()
def plot_D_chord(final_results, TSR, save=False):
    for i in range(len(TSR)):
        plt.plot(final_results[:,2,i], final_results[:,14,i], label = 'TSR = '+str(TSR[i]))
    plt.ylabel(r'$D$')
    plt.xlabel('r/R')
    plt.grid()
    plt.legend()
    if save:
        plt.savefig('D_chord.png')
        plt.cla()
    else:
        plt.show()
def plt_TQ(T,Q,save=False):
    T=np.array(T)
    Q = np.array(Q)
    plt.plot(TSR, T/1000, label=r'$T$')
    plt.plot(TSR, Q/1000, label=r'$Q$')
    plt.xlabel('TSR')
    plt.ylabel(r'$T[kN],Q [kNm]$')
    plt.legend()
    plt.grid()
    if save:
        plt.savefig('totalTandQ_TSR.png', bbox_inches='tight')
        plt.cla()
    else:
        plt.show()
        plt.cla()
def plt_CP(CT, CP, save=False):
    TSR_max = np.argmax(CP)
    CP_max = CP[TSR_max]
    plt.plot(TSR[TSR_max], CP_max, 'ro')
    plt.plot(TSR, CP, label=r'$C_P$')
    plt.plot(TSR, CT, label=r'$C_T$')
    plt.xlabel(r'$TSR$')
    plt.ylabel(r'$C_P, C_T$')
    # Add textbox next to the plotted point
    plt.text(TSR[TSR_max], CP_max, r'$C_{P_{max}} = $' + f'{CP_max:.3f})', ha='right', va='bottom')
    plt.legend()
    plt.grid()
    if save:
        plt.savefig('CP_CT.png')
        plt.cla()
    else:
        plt.show()
def plot_Cl_vs_alpha(polar_cl, polar_alpha, polar_cd, save=False):
    plt.plot(polar_alpha, polar_cl)
    plt.ylabel(r'$Cl$ [-]')
    plt.xlabel(r'$\alpha$ [deg]')
    plt.grid()
    max_cl_index = np.argmax(polar_cl)
    max_cl_alpha = polar_alpha[max_cl_index]
    max_cl_value = polar_cl[max_cl_index]
    plt.plot(max_cl_alpha, max_cl_value, 'rx', markersize=10, markeredgewidth=1)
    cl_cd_ratio = polar_cl / polar_cd
    max_ratio_index = np.argmax(cl_cd_ratio)
    max_ratio_alpha = polar_alpha[max_ratio_index]
    max_ratio_cl = polar_cl[max_ratio_index]
    plt.plot(max_ratio_alpha, max_ratio_cl, 'kx', markersize=10, markeredgewidth=1)
    if save:
        plt.savefig('Cl_vs_alpha.png')
        plt.cla()
    else:
        plt.show()
def plot_Cl_vs_Cd(polar_cl, polar_cd, save=False):
    plt.plot(polar_cd, polar_cl)
    plt.ylabel(r'$Cl$ [-]')
    plt.xlabel('Cd [-]')
    plt.grid()
    plt.xlim([0,0.1])
    max_cl_index = np.argmax(polar_cl)
    max_cl_cd = polar_cd[max_cl_index]
    max_cl_value = polar_cl[max_cl_index]
    plt.plot(max_cl_cd, max_cl_value, 'rx', markersize=10, markeredgewidth=1)
    cl_cd_ratio = polar_cl / polar_cd
    max_ratio_index = np.argmax(cl_cd_ratio)
    max_ratio_cd = polar_cd[max_ratio_index]
    max_ratio_cl = polar_cl[max_ratio_index]
    plt.plot(max_ratio_cd, max_ratio_cl, 'kx', markersize=10, markeredgewidth=1)
    if save:
        plt.savefig('Cl_vs_Cd.png')
        plt.cla()
    else:
        plt.show()
def plot_Ct_rR(final_results, TSR, save=False):
    for i in range(len(TSR)):
        plt.plot(final_results[:,2,i], final_results[:,3,i]*final_results[:,15,i]/(0.5 * Uinf**2 * Radius), label = 'TSR = '+str(TSR[i]))
    plt.ylabel(r'$C_{T}$')
    plt.xlabel('r/R')
    plt.grid()
    plt.legend()
    if save:
        plt.savefig('Ct_rR.png')
        plt.cla()
    else:
        plt.show()
def plot_Ct_rR_Prandtl(final_results, final_results_wP, save=False):
    plt.plot(final_results[:,2,1], final_results[:,3,1]*final_results[:,15,1]/(0.5 * Uinf**2 * Radius), label = 'With tip and loss correction')
    plt.plot(final_results_wP[:, 2, 1], final_results_wP[:, 3, 1] * final_results_wP[:, 15, 1] / (0.5 * Uinf ** 2 * Radius), label='Without tip and loss correction')
    plt.ylabel(r'$C_{T}$')
    plt.xlabel('r/R')
    plt.grid()
    plt.legend()
    if save:
        plt.savefig('Ct_rR_Prandtl.png')
        plt.cla()
    else:
        plt.show()

save = False
#pressure(final_results, TSR)
# plt_TQ(T,Q,save)
# plt_CP(CT, CP, save=False)

# areas = (r_R[1:]**2-r_R[:-1]**2)*np.pi*Radius**2
# dr = (r_R[1:]-r_R[:-1])*Radius
# print(dr)

# areas = (r_R[1:]**2-r_R[:-1]**2)*np.pi*Radius**2
# dr = (r_R[1:]-r_R[:-1])*Radius
# CT = np.sum(dr*results[:,3]*NBlades/(0.5*Uinf**2*np.pi*Radius**2))
# CP = np.sum(dr*results[:,4]*results[:,2]*NBlades*Radius*Omega/(0.5*Uinf**3*np.pi*Radius**2))


#plot_alpha_rR(final_results, TSR, save)
plot_alpha_rR_Prandtl(final_results, final_results_wP, save)
#plot_phi_rR(final_results, TSR, save)
plot_phi_rR_Prandtl(final_results, final_results_wP, save)
# plot_a_rR(final_results, TSR, save)
plot_a_rR_Prandtl(final_results, final_results_wP, save)
# plot_aprime_rR(final_results, TSR, save)
# plot_fnorm_rR(final_results, TSR, save)
plot_ftan_rR(final_results, TSR, save)
plot_ftan_rR_Prandtl(final_results, final_results_wP, save)
plot_tiprootloss(final_results, TSR, save)
# plot_cl_chord(final_results, TSR, save)
# plot_cd_chord(final_results, TSR, save)
# plot_L_chord(final_results, TSR, save)
# plot_D_chord(final_results, TSR, save)
plot_Cl_vs_alpha(polar_cl, polar_alpha, polar_cd, save)
plot_Cl_vs_Cd(polar_cl, polar_cd, save)
plot_Ct_rR(final_results, TSR, save)
plot_Ct_rR_Prandtl(final_results, final_results_wP, save)