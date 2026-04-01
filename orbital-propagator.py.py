import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from poliastro.bodies import Earth
from poliastro.twobody import Orbit
import astropy.units as u
mu = 3.986e14; a=6778000; e=0
class OrbitaPropagatore:
    def __init__(self,a,e):
        self.a=a
        self.e=e
    def periodoorbitale(self):
        T=2*np.pi*(np.sqrt((self.a**3)/mu))
        return T
    def propaga(self):
        T=self.periodoorbitale()
        v0=np.sqrt(mu/self.a)
        c0=[self.a,0,0,v0]
        def f(t,y):
            rx=y[0]
            ry=y[1]
            r=np.sqrt(rx**2+ry**2)
            ax=-(mu/(r**3))*rx
            ay=-(mu/(r**3))*ry
            drxdt=y[2]
            drydt=y[3]
            dvxdt=ax
            dvydt=ay
            func=[drxdt,drydt,dvxdt,dvydt]
            return func
        sol=solve_ivp(f,[0,T],c0,t_eval=np.linspace(0, T, 1000))
        return sol
    def plotta(self):
        sol=self.propaga()
        fig,ax=plt.subplots()
        ax.plot(sol.y[0],sol.y[1])
        ax.set_title("Propagazione orbitale")
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_aspect('equal')
        ax.set_xlim(min(sol.y[0]),max(sol.y[0]))
        ax.set_ylim(min(sol.y[1]),max(sol.y[1]))
        plt.show()
    def salva(self,filename):
        sol=self.propaga()
        dati=np.column_stack((sol.y[0],sol.y[1],sol.t))
        np.savetxt(filename,dati,header="colonne",delimiter=",")
        return
    def confronta(self):
        sol=self.propaga()
        theta = np.linspace(0, 2*np.pi, 1000)
        x_poli = self.a * np.cos(theta)
        y_poli = self.a * np.sin(theta)
        fig,ax=plt.subplots()
        ax.plot(sol.y[0],sol.y[1],color='b',label="Mia orbita")
        ax.plot(x_poli,y_poli,color='r',linestyle='--',label="Orbita di poliastro")
        ax.set_title("Propagazione orbitale mia vs poliastro")
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_aspect('equal')
        ax.set_xlim(min(sol.y[0]),max(sol.y[0]))
        ax.set_ylim(min(sol.y[1]),max(sol.y[1]))
        ax.legend()
        plt.show()

Terra=OrbitaPropagatore(a,e)
T=Terra.periodoorbitale()/60
print(T)
#sol=Terra.propaga()
#print(sol.y.shape)
#Terra.plotta()
#Terra.salva("orbita-csv")
#import os
#print(os.getcwd())
Terra.confronta()