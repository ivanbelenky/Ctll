import math
import numpy as np
from .. import units as u

"""Calculador de los parámetros orbitales para los tiempos de revisita deseados.

 SS := Sun Synchronous. Para mas informacion referirse a 
 https://www.sciencedirect.com/science/article/abs/u.pii/S1270963814002351"""


#funciones utilidad para el cálculo de perturbaciones

def n0(a):
    return (u.mu/(a**3))**(1/2)

def get_i(a,n):
    return math.acos( -(2/3) * ((a/u.R_e)**2) * (u.w_s*86400/(2*u.pi*n*u.J_2)) )

def dm(a,i):
    return u.J_2*((u.R_e/a)**2)*((3/2)-(9/4)*(math.sin(i)**2))

def F(a,i):
    
    return dm(a,i) + ((u.J_2**2)*(u.R_e/a)**4) * ( (27/8) - (15/16)*(math.sin(i)**2) - (411/64)*(math.sin(i)**4) ) 

def G(a,i):
    return (u.J_2*((u.R_e/a)**2))*(3 - (15/4)*(math.sin(i)**2) ) \
+ ((u.J_2**2)*((u.R_e/a)**4))*( 9 - (309/16)*(math.sin(i)**2) + (645/64)*(math.sin(i)**4)) + \
u.J_4*((u.R_e/a)**4)*( (15/2) - (465/16)*(math.sin(i)**2) + (735/32)*(math.sin(i)**4) ) 
    

#Funcion principal para obtención de parámetros orbitales en una órbita circular.

def get_a_i_N(D):
    """
    Calcula los parametros de las posibles orbitas SS dado los dias de revisita 
    Recibe: D dias de revisita 
    Devuele: tres arrays radios, inclinaciones, mean-motions
    """
    
    a_ = [ ]
    i_ = [ ]
    N_ = [ ]
    
    for i in range(7*D,17*D+1):
        
        N_.append(i/D)
        n_0 = i/D
        T_d = u.s_m*D/i
        a_0 = ( u.mu*((T_d/(2*u.pi))**2) )**(1/3)
        i_0 = get_i(a_0,n_0)
        

        count = 0
        error = float('inf')
        i_ant = i_0
        
        
        while (count < u.ITERMAX and error > u.TOL):
            n_ = n_0*(1 + dm(a_0,i_ant))
            i_sig = get_i(a_0,n_)


            d_n = (1/(i/D)) * F(a_0,i_sig)
            d_w = (1/(i/D)) * G(a_0,i_sig)
            n_0 = (i/D)*(1-d_n-d_w)
            a_0 = ( u.mu * ((u.s_m/(2*u.pi*n_0))**2) )**(1/3)

            count += 1
            error = abs(i_sig-i_ant)
            i_ant = i_sig


        i_.append(i_sig)
        a_.append(a_0)
        
    return np.array(a_), np.array(i_), np.array(N_)




def omdot(a=0,e=0,i=0,nd=0):
    return -3/2*nd*u.J_2*np.cos(i)*((u.R_e/a)**2)/((1-e)**2)