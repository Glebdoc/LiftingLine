import numpy as np

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

