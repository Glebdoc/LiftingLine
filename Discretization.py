import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from utilits import *
from geometry import Point, Vortex, HorseShoe


class Rotor:
    def __init__(self, point, R, TSR, NB, N, r_R, chord_distribution, twist_distribution, aw, t):
        self.x0 = point.x
        self.y0 = point.y
        self.z0 = point.z
        self.R = R
        self.t = t
        self.TSR = TSR
        self.NB = NB
        self.N = N
        self.r_R = r_R
        self.chord_distribution = chord_distribution
        self.twist_distribution = twist_distribution
        self.aw = aw
        self.horseShoes = []
        self.collocations = np.vstack((np.zeros(N), ((r_R[:-1] * R + r_R[1:] * R) / 2), np.zeros(N)))
        self.Lift = np.zeros(N)

    def assemble(self, Uinf):
        x2 = self.r_R * np.sin(-self.twist_distribution * np.pi / 180)
        z2 = self.chord_distribution * np.cos(-self.twist_distribution * np.pi / 180)
        Uwake = Uinf * (1 - self.aw)
        xw = self.t * Uwake
        n_azim = [np.array([0, 0, 1])]
        translation = np.array([self.x0, self.y0, self.z0])
        for blade in range(self.NB):
            shift = blade * 2 * np.pi / self.NB
            for i in range(len(self.r_R) - 1):
                # assemble the left vortex
                vector1 = spatial_transform(shift, np.array([x2[i], self.r_R[i] * self.R, -z2[i]])) + translation
                vector2 = spatial_transform(shift, np.array([0, self.r_R[i] * self.R, 0])) + translation
                left = [Vortex(Point(*vector1), Point(*vector2), 1)]
                # assemble the centre vortex
                vector1 = spatial_transform(shift, np.array([0, self.r_R[i] * self.R, 0])) + translation
                vector2 = spatial_transform(shift, np.array([0, self.r_R[i + 1] * self.R, 0])) + translation
                centre = Vortex(Point(*vector1), Point(*vector2), 1)
                # assemble the right vortex
                vector1 = spatial_transform(shift, np.array([0, self.r_R[i + 1] * self.R, 0])) + translation
                vector2 = spatial_transform(shift,
                                            np.array([x2[i + 1], self.r_R[i + 1] * self.R, -z2[i + 1]])) + translation
                right = [Vortex(Point(*vector1), Point(*vector2), 1)]

                yw_left = self.r_R[i] * self.R * np.cos((Omega) * t)
                zw_left = z2[i] + self.r_R[i] * self.R * np.sin((Omega) * t)
                zw_left = -zw_left

                yw_right = self.r_R[i + 1] * self.R * np.cos((Omega) * t)
                zw_right = z2[i + 1] + self.r_R[i + 1] * self.R * np.sin((Omega) * t)
                zw_right = -zw_right
                left = left + assembleLeftVortex(x2[i] + xw, yw_left, zw_left, shift, translation)
                right = right + assembleRightVortex(x2[i + 1] + xw, yw_right, zw_right, shift, translation)

                self.horseShoes.append(HorseShoe(left, centre, right))

            if blade == 0:
                continue
            n_azim.append(spatial_transform(shift, n_azim[0]))
            self.collocations = np.hstack(
                (self.collocations, spatial_transform(shift, self.collocations[:, :N], flatten=False)))
        collocations = self.collocations + translation.reshape(3, 1)
        collocations = collocations.T
        self.collocations = collocations

        return self.collocations, self.horseShoes, n_azim


############################################
# Read polar data
data = pd.read_excel('polars.xlsx')
polar_cl = data['Cl'].to_numpy()
polar_cd = data['Cd'].to_numpy()
polar_alpha = data['Alfa'].to_numpy()

# Read BEM data
# keys [a, aline, r_R, fnorm, ftan, gamma, phi, AoA, Prandtl, Prandtltip, Prandtlroot, cl_chord, cd_chord, L_chord, D_chord, sigma]
BEM_TSR6 = np.load('BEM_results.npy')  # [:, key, tsr]
AW = [np.mean(BEM_TSR6[:,0,0]), np.mean(BEM_TSR6[:,0,1]), np.mean(BEM_TSR6[:,0,2])]
TSR = [6, 8, 10]
# Wind turbine geometry
R = 50  # [m]
NB = 3  # -
TIP_LOCATION = 1
ROOT_LOCATION = 0.2
PITCH = 2  # degrees
N_new = np.arange(1, 42, 2)
discretization_type = 'cosine'#'uniform'  # 'cosine'

NRotations = 10  # number of full rotations in the wake
dt = 10  # time steps per rotation
# t = np.linspace(0, 30, 150)
# Wind
Uinf = 10  # unperturbed wind speed in m/s


plotting = True

Ct_array = []
Cp_array = []

for k in range(len(N_new)):
    for tsr in range(1, 2):
        N = N_new[k]
        # Discretize geometry
        r_R = spanwise_discretization(N, ROOT_LOCATION, TIP_LOCATION, discretization_type)
        chord_distribution = 3 * (1 - r_R) + 1  # meters
        twist_distribution = -14 * (1 - r_R) + PITCH  # degrees
        aw = AW[tsr]
        localTSR = TSR[tsr]
        Omega = Uinf * localTSR / R
        Tfinal = 2 * np.pi * NRotations / Omega
        t = np.linspace(0, Tfinal, dt * NRotations)
        rotor1 = Rotor(Point(0, 0, 0), R, localTSR, NB, N, r_R, chord_distribution, twist_distribution, aw, t)
        # rotor2 = Rotor(Point(0,0,(2*R)*2), R, localTSR, NB, N, r_R, chord_distribution, twist_distribution, aw, t)
        rotors = [rotor1]
        NR = len(rotors)

        total_colloc, total_horses, total_nazim = [], [], []

        for rotor in rotors:
            colloc, horses, n_azim = rotor.assemble(Uinf)
            total_nazim.append(n_azim)
            total_colloc.append(colloc)
            total_horses.append(horses)

        total_colloc = np.array(total_colloc).reshape(-1, 3)
        total_horses = np.array(total_horses).flatten()
        total_nazim = np.array(total_nazim)

        u_influences = np.zeros((N * NB * NR, N * NB * NR))
        v_influences = np.zeros((N * NB * NR, N * NB * NR))
        w_influences = np.zeros((N * NB * NR, N * NB * NR))

        for i, colloc in enumerate(total_colloc):
            for j, horse in enumerate(total_horses):
                x = colloc[0]
                y = colloc[1]
                z = colloc[2]
                u_vector = horse.velocity(Point(x, y, z))
                u_influences[i, j] = u_vector[0]
                v_influences[i, j] = u_vector[1]
                w_influences[i, j] = u_vector[2]

        Gammas = np.ones(N * NB * NR)

        err = 1.0
        iter = 0

        while (err > 1e-6 and iter < 1200):
            iter += 1
            u = u_influences @ Gammas
            v = v_influences @ Gammas
            w = w_influences @ Gammas

            vazim = np.zeros(N * NB * NR)
            vaxial = u + Uinf
            for j in range(NR):
                n_azim = total_nazim[j]
                for i in range(N * NB):
                    count = i % N
                    k = i // N
                    vrot = np.cross([Omega, 0, 0], total_colloc[i + j * N * NB])
                    vel1 = np.array([Uinf + u[i + j * N * NB] + vrot[0], 0 + vrot[1] + v[i + j * N * NB],
                                     0 + vrot[2] + w[i + j * N * NB]])
                    vazim[i + j * N * NB] = np.dot(n_azim[k], vel1)

            vtotal = np.sqrt(vaxial ** 2 + vazim ** 2)

            inflowangle = np.arctan2(vaxial, vazim)
            twist = np.interp(total_colloc[:N, 1], r_R * R, twist_distribution)
            chord = np.interp(total_colloc[:N, 1], r_R * R, chord_distribution)
            alpha = twist + inflowangle[:N] * 180 / np.pi

            # Compute Cl from alphas
            Cl = np.interp(alpha, polar_alpha, polar_cl)

            Gammas_old = np.copy(Gammas)
            weight = 0.15

            # assemble gammas
            for rotor in range(NR):
                for blade in range(NB):
                    for i in range(N):
                        Gammas[i + blade * N + rotor * NB * N] = weight * Cl[i] * 0.5 * chord[i] * vtotal[
                            i + blade * N + rotor * NB * N] + (1 - weight) * Gammas_old[i + blade * N + rotor * NB * N]
            err = np.linalg.norm(Gammas - Gammas_old)

        Cd = np.interp(alpha, polar_alpha, polar_cd)
        Lift = 0.5 * chord * 1.225 * vtotal[:N] * vtotal[:N] * Cl
        Drag = 0.5 * chord * 1.225 * vtotal[:N] * vtotal[:N] * Cd
        LiftND = Lift / (0.5 * 1.225 * Uinf ** 2 * R)
        DragND = Drag / (0.5 * 1.225 * Uinf ** 2 * R)

        Fazim = Lift * np.sin(inflowangle[:N]) - Drag * np.cos(inflowangle[:N])
        FazimND = Fazim / (0.5 * 1.225 * Uinf ** 2 * R)
        Faxial = Lift * np.cos(inflowangle[:N]) + Drag * np.sin(inflowangle[:N])
        FaxialND = Faxial / (0.5 * 1.225 * Uinf ** 2 * R)
        GammasND = Gammas / (Uinf * Uinf * np.pi / (NB * Omega))
        dr = r_R[1:] * R - r_R[:-1] * R
        Areas = 2 * np.pi * R * dr
        Ct = np.sum(Faxial * NB * dr) / (0.5 * 1.225 * Uinf ** 2 * np.pi * R ** 2)
        Cp = np.sum(Fazim * NB * dr * r_R[:-1] * Omega * R) / (0.5 * 1.225 * Uinf ** 3 * np.pi * R ** 2)
        Ct_array = np.append(Ct_array, Ct)
        Cp_array = np.append(Cp_array, Cp)
        '''results = np.column_stack((GammasND[:N], alpha, inflowangle[:N] * 180 / np.pi, Cl, chord,
                                   rotor1.collocations[:N, 1] / R, FazimND, FaxialND, np.ones((len(alpha))) * Ct,
                                   np.ones((len(alpha))) * Cp))
        np.savetxt(f'LL_{localTSR}.csv', results, delimiter=',',
                   header='GammasND,alpha,inflow,Cl,chord, points, FazimND, FaxialND, Ct, Cp', comments='')'''


        print('TSR:', localTSR, 'iterations:', iter, 'error:', err)
print('Ct = ', Ct_array)
print('Cp', Cp_array)
results = np.column_stack((N_new, Ct_array, Cp_array))
np.savetxt(f'LL_C_blade_{localTSR}.csv', results, delimiter=',',
                   header='N_new, Ct, Cp', comments='')