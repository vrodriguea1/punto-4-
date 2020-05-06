
import numpy as np
import matplotlib.pyplot as plt
import math
from math import exp
from scipy import constants

#constantes 

pi= 3.1416
c= 3e8 #velocidad de la luz 
e0 = constants.epsilon_0
u0 = constants.mu_0
e1= 12  #indice del plastico 


#parametros del sistema 

px= 500  #tamaño universo 
pxm=int(px/2) #mitad universo 
Hy= np.zeros(px, dtype= "float64")  
Ez= np.zeros(px, dtype= "float64")
lol0= np.zeros(px, dtype= "float64")
lol1= np.zeros(px, dtype= "float64")
E1= 0.0254/2 
freq= 10.3*10**9 #frecuencia en hz
dx= 0.03 # factor de discretización 
dt= dx/c
dy= dt
Ttot= 100
E0= 100
xmax= 0.0254 #el lmite del cuadrado de plastico en m 
ymax= xmax 
No= math.sqrt(u0/e0) #impedancia intrinseca del aire 
N1= math.sqrt(u0/(e0*e1)) #impedancia intrinseca de la capa de plastico  
index= 1.23 #indice de refracción
lan= c/freq #longuitud de onda
m=1
d=((m-1/2)*lan)/2*index #grosor de las placas antirefelctivas 
e2= (u0*index**2/e0*No) #epsilon de las capas 
N2= math.sqrt(u0/(e0*e2)) #impedacia intrinseca de las cpas antireflectantes 
factor_norm = math.sqrt(e0/u0)
t0 = 40
spread = 12
Npasos = 450
#d= lan/(4*1,24) #distancia las  

#condiciones de frontera 
#
#Fron_D = [0, 0]
#Fron_U = [0, 0]
#Ez[0] = Fron_D.pop(0)
#Fron_D.append(Ez[1])
#Ez[px - 1] = Fron_U.pop(0)
#Fron_U.append(Ez[px - 2])

medio1= int((px*d)/xmax)  
medio2= int(px-medio1-(E1))
medio3= int(E1) 

lol0[:]= u0 #todo el universo 
lol1[0: medio1]= e0 #medio aire
lol1[medio1:medio2]= 12*e0 #medio capas reflectantes
lol1[medio2:]= e0*e1 #medio plastico 
 

#FDTD 1D
for n in range(1, Ttot+1):
    
    #condiciones de frontera 
    Fron_D = [0, 0]
    Fron_U = [0, 0]
    Ez[0] = Fron_D.pop(0)
    Fron_D.append(Ez[1])
    Ez[px - 1] = Fron_U.pop(0)
    Fron_U.append(Ez[px - 2])
    
   #para el medio absorvente
   #Ez[px] = Ez[px - 1]+((c*(dt-dx))/(c*(dt+dx)))*(Ez[px -1] - Ez[px])
    
    #ingresamos la onda
    puls= exp(-0.5*((t0-n)/spread)**2)
    Ez[pxm]= puls 
    
    plt.clf()
    #se calcula Hy
    for y in range(0,(px-1)):
         #Hy[n,y] = Hy[n, y]+(dt/(u0)*(dx))*((Ez[n, y + 1] -Ez[n, y])/dy)
     Hy[y] = Ez[y] + (dt/(lol0[y]*dx))*(Ez[y+1] - Ez[y])
    #para calcular Ez
    for k in range(1, px):
       # Ez[n,k] = Ez[n,k]+(dt/(e0)*(dx))*(Hy[n, k]- Hy[n, k-1]) 
      Ez[k] = Ez[k] + (dt/(lol1[k]*dx))*(Hy[k] - Hy[k-1])


#%%graficas     
    plt.subplot(211)
    plt.plot(Ez, color='b', linewidth=1)
    plt.ylabel('Campo eléctrico [V/m]', fontsize='8')
    plt.xlim(0, px)
    plt.ylim(-(E0+0.2)*factor_norm, (E0+0.2)*factor_norm)
    
   
    plt.subplot(212)
    plt.plot(Hy, color='r', linewidth=1)
    plt.ylabel('Campo magnetico [T]', fontsize='8')
    plt.xlim(0, px)
    plt.ylim(-(E0+0.2)*factor_norm, (E0+0.2)*factor_norm)
    
    plt.pause(0.001)
    plt.show()   