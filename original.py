import numpy as np
import scipy.linalg as la
import matplotlib.pyplot as plt
import pandas as pd
import pyvista as pv

class Point:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

class Vortex(Point):
    def __init__(self, point1, point2, Gamma):
        self.x1 = point1.x
        self.y1 = point1.y
        self.z1 = point1.z
        self.x2 = point2.x
        self.y2 = point2.y
        self.z2 = point2.z
        self.Gamma = Gamma

    def crossX(self, x, y, z):
        return (y - self.y1)*(z - self.z2) - (z - self.z1)*(y - self.y2)
    
    def crossY(self, x, y, z):
        return (z - self.z1)*(x - self.x2) - (x - self.x1)*(z - self.z2)
    
    def crossZ(self, x, y, z):
        return (x - self.x1)*(y - self.y2) - (y - self.y1)*(x - self.x2)
    
    def crossMag(self, x, y, z):
        return np.sqrt(self.crossX(x, y, z)**2 + self.crossY(x, y, z)**2 + self.crossZ(x, y, z)**2)
    
    def dot1(self, x, y, z):
        return (self.x2 - self.x1)*(x - self.x1) + (self.y2 - self.y1)*(y - self.y1) + (self.z2 - self.z1)*(z - self.z1)
    
    def dot2(self, x, y, z):
        return (self.x2 - self.x1)*(x - self.x2) + (self.y2 - self.y1)*(y - self.y2) + (self.z2 - self.z1)*(z - self.z2)
    
    def r1(self, x, y, z):
        return np.sqrt((x - self.x1)**2 + (y - self.y1)**2 + (z - self.z1)**2)
    
    def r2(self, x, y, z): 
        return np.sqrt((x - self.x2)**2 + (y - self.y2)**2 + (z - self.z2)**2)

    def length(self):
        return np.sqrt((self.x2 - self.x1)**2 + (self.y2 - self.y1)**2)
    
    def velocity(self, point):
        x = point.x
        y = point.y
        z = point.z

        crossMag = self.crossMag(x, y, z)
        if crossMag < 1e-6:
            return np.array([0, 0, 0])
        r1 = self.r1(x, y, z)
        if r1 < 1e-6:
            return np.array([0, 0, 0])
        r2 = self.r2(x, y, z)
        if r2 < 1e-6:
            return np.array([0, 0, 0])
        
        K = (self.Gamma/(4*np.pi*crossMag*crossMag)) * (self.dot1(x, y, z)/r1 - self.dot2(x, y, z)/r2)
        return K*np.array([self.crossX(x, y, z), self.crossY(x, y, z), self.crossZ(x, y, z)])

class HorseShoe(Vortex):
    def __init__(self, left, centre, right):
        self.leftset = left
        self.centre = centre
        self.rightset = right
        

    def velocity(self, point):
        vort1 = 0
        vort3 = 0
        for i in range(len(self.leftset)):
            vort1 += self.leftset[i].velocity(point)
            vort3 += self.rightset[i].velocity(point)
        vort2 = self.centre.velocity(point)
        return vort1 + vort2 + vort3


# Function for spanwise discretization
def spanwise_discretization(N, ROOT_LOCATION, TIP_LOCATION, discretization_type='uniform'):
    if discretization_type == 'uniform':
        return np.linspace(ROOT_LOCATION, TIP_LOCATION, N+1)
    if discretization_type == 'cosine':
        # Generate cosine-spaced points between 0 and π
        theta = np.linspace(np.pi, 0, N+1)
        
        # Cosine-spaced points in the range [-1, 1]
        cos_spaced = np.cos(theta)
        
        # Map these points to the spanwise locations
        half_span = (TIP_LOCATION - ROOT_LOCATION) / 2
        spanwise_locations = ROOT_LOCATION + half_span + (cos_spaced) * half_span
        
        return spanwise_locations


def assembleLeftVortex(xw, yw, zw, shift):
    left = []
    for i in range(len(xw)-1):
        vector1 = np.array([xw[i+1], yw[i+1], zw[i+1]])
        vector2 = np.array([xw[i], yw[i], zw[i]])
        left.append(Vortex(Point(*spatial_transform(shift, vector1)), Point(*spatial_transform(shift, vector2)), 1))
    return left

def assembleRightVortex(xw, yw, zw, shift):
    right = []
    for i in range(len(xw)-1):
        vector1 = np.array([xw[i], yw[i], zw[i]])
        vector2 = np.array([xw[i+1], yw[i+1], zw[i+1]])
        right.append(Vortex(Point(*spatial_transform(shift, vector1)), Point(*spatial_transform(shift, vector2)), 1))
    return right

def spatial_transform(theta, vector, flatten=True):
    vector.transpose()
    Rx = np.array([[1, 0, 0], 
                   [0, np.cos(theta), -np.sin(theta)], 
                   [0, np.sin(theta), np.cos(theta)]])
    if flatten:
        return Rx@vector.flatten()
    else:
        return Rx@vector

############################################
data = pd.read_excel('polars.xlsx')
polar_cl = data['Cl'].to_numpy()
polar_cd = data['Cd'].to_numpy()
polar_alpha = data['Alfa'].to_numpy()
  
# Wind turbine geometry
R = 50 #[m]
NB = 3 #-
TIP_LOCATION =  1
ROOT_LOCATION =  0.2
PITCH = 2 # degrees 
N = 24
aw = 0.33
TSR = 8
discretization_type = 'uniform' # 'cosine'
t = np.linspace(0, 5, 150)
plotting = False

# Wind 
Uinf = 10 # unperturbed wind speed in m/s

# Discretize geometry
r_R= spanwise_discretization(N, ROOT_LOCATION, TIP_LOCATION, discretization_type)
chord_distribution = 3*(1-r_R)+1 # meters
twist_distribution = -14*(1-r_R)+PITCH # degrees

# Define collocation points
collocations = (r_R[:-1]*R + r_R[1:]*R) / 2
collocations = np.vstack((np.zeros(N), collocations, np.zeros(N)))

# Main surface coordinates 
x2 =r_R*np.sin(-twist_distribution*np.pi/180)
z2 = chord_distribution*np.cos(-twist_distribution*np.pi/180)

Omega = Uinf*TSR/R
Uwake = Uinf*(1-aw)

colors = ['r', 'g', 'b']
myHorses = []

xw = t*Uwake
collocations_total = []
n_azim = [np.array([0, 0, 1])]

for blade in range(NB):
    shift = blade*2*np.pi/NB
    for i in range(len(r_R)-1):
        # assemble the left vortex
        vector1 = np.array([x2[i], r_R[i]*R, -z2[i]])
        vector2 = np.array([0, r_R[i]*R, 0])
        left = [Vortex(Point(*spatial_transform(shift, vector1)), Point(*spatial_transform(shift, vector2)), 1)] 
        # assemble the centre vortex
        vector1 = np.array([0, r_R[i]*R, 0])
        vector2 = np.array([0, r_R[i+1]*R, 0])
        centre = Vortex(Point(*spatial_transform(shift, vector1)), Point(*spatial_transform(shift, vector2)), 1) 
        # assemble the right vortex
        vector1 = np.array([0,r_R[i+1]*R, 0])
        vector2 = np.array([x2[i+1], r_R[i+1]*R, -z2[i+1]])
        right = [Vortex(Point(*spatial_transform(shift, vector1)), Point(*spatial_transform(shift, vector2)), 1)] 

        yw_left = r_R[i]*R * np.cos((Omega)*t)
        zw_left = z2[i] + r_R[i]*R * np.sin((Omega)*t)
        zw_left = -zw_left
        
        yw_right = r_R[i+1]*R * np.cos((Omega)*t)
        zw_right = z2[i+1] + r_R[i+1]*R * np.sin((Omega)*t)
        zw_right = -zw_right
        left = left + assembleLeftVortex(x2[i] + xw, yw_left, zw_left, shift) 
        right = right + assembleRightVortex(x2[i+1] + xw, yw_right, zw_right, shift)

        myHorses.append(HorseShoe(left, centre, right))

        #Update collocations for shifts 

    if blade == 0:
        continue
    n_azim.append(spatial_transform(shift, n_azim[0]))
    collocations = np.hstack((collocations, spatial_transform(shift, collocations[:,:N], flatten=False)))
collocations = collocations.T

u_influences = np.zeros((N*NB, N*NB))
v_influences = np.zeros((N*NB, N*NB))
w_influences = np.zeros((N*NB, N*NB))

for i, colloc in enumerate(collocations):
    for j, horse in enumerate(myHorses):
        x = colloc[0]
        y = colloc[1]
        z = colloc[2]
        u_vector = horse.velocity(Point(x,y,z))
        u_influences[i,j] = u_vector[0]
        v_influences[i,j] = u_vector[1]
        w_influences[i,j] = u_vector[2]

Gammas = np.ones(N*NB)

err = 1.0
iter = 0

print('n_azim', n_azim)
while (err > 1e-5 and iter<1000):
    iter+=1
    u = u_influences@Gammas
    v = v_influences@Gammas
    w = w_influences@Gammas

    vazim = np.zeros(N*NB)
    vaxial = u + Uinf

    for i in range(N*NB):
        count = i%N
        j = i//N
        vrot  = np.cross([Omega, 0 , 0], collocations[i])
        vel1 = np.array([Uinf + u[i] + vrot[0], 0 + vrot[1] + v[i], 0 + vrot[2] + w[i]])
        vazim[i] = np.dot(n_azim[j], vel1)

    vtotal = np.sqrt(vaxial**2 + vazim**2)

    inflowangle = np.arctan2(vaxial  ,vazim)
    twist = np.interp(collocations[:N,1],  r_R*R, twist_distribution)
    chord = np.interp(collocations[:N,1],  r_R*R, chord_distribution)
    alpha = twist + inflowangle[:N]*180/np.pi

    # Compute Cl from alphas 
    Cl = np.interp(alpha, polar_alpha, polar_cl)
    Gammas_old = Gammas
    weight = 0.2
    Gammas = weight*Cl*0.5*chord*vtotal[:N] + (1-weight)*Gammas_old[:N]
    Gammas = np.tile(Gammas, NB)
    #print(Gammas[:N])
    err = np.linalg.norm(Gammas - Gammas_old)
    print('iter:', iter, 'err',err)

Cd = np.interp(alpha, polar_alpha, polar_cd)
Lift = 0.5*chord*1.225*vtotal[:N]*vtotal[:N]*Cl
Drag = 0.5*chord*1.225*vtotal[:N]*vtotal[:N]*Cd
LiftND = Lift/(0.5*1.225*Uinf**2*R)
DragND = Drag/(0.5*1.225*Uinf**2*R)

Fazim = Lift*np.sin(inflowangle[:N]) - Drag*np.cos(inflowangle[:N])
FazimND = Fazim/(0.5*1.225*Uinf**2*R)
Faxial = Lift*np.cos(inflowangle[:N]) + Drag*np.sin(inflowangle[:N])
FaxialND = Faxial/(0.5*1.225*Uinf**2*R)

GammasND = Gammas/(Uinf*Uinf*np.pi/(NB*Omega))

#####################

######################
plt.plot(collocations[:N,1]/R, GammasND[:N], label='Gammas ND')
plt.show()
plt.close()
plt.plot(collocations[:N,1]/R, alpha, label='alpha')
plt.show()
plt.close()
plt.plot(collocations[:N,1]/R, -((u[:N]+vrot[0])/Uinf), label='vaxial')
plt.show()
plt.close()
plt.plot(collocations[:N,1]/R, FazimND, label='Fazim')
plt.plot(collocations[:N,1]/R, FaxialND, label='Faxial')
plt.grid()
plt.xlabel('r/R')
plt.legend()
plt.ylabel('F')
plt.show()

if plotting:
    # BIG PLOT
    ax = plt.figure().add_subplot(111,projection='3d')
    ax.scatter(collocations[:,0], collocations[:,1], collocations[:,2], color='black')
    for horse in myHorses:
        for vortex in horse.leftset:
            plt.plot([vortex.x1, vortex.x2], [vortex.y1, vortex.y2], [vortex.z1, vortex.z2],'r', alpha=0.5) 
        for vortex in horse.rightset:
            plt.plot([vortex.x1, vortex.x2], [vortex.y1, vortex.y2], [vortex.z1, vortex.z2],'g', alpha=0.5)
        vortex = horse.centre
        plt.plot([vortex.x1, vortex.x2], [vortex.y1, vortex.y2], [vortex.z1, vortex.z2],'b')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.show()
