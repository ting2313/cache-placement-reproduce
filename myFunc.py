import math
import numpy as np
from sympy import *
import scipy.integrate
from scipy.special import expi

# pico BS r = 0.15km

v = 1
P_0 = 20 # W
P_n = 1 # pico power
noise = -105 # dBm
alpha = 3.76 # origin=4
lamda = 500
euler = 0.577
F = 1000
L = 1 # Mbits
m_0 = 1
A_m = 0.15*0.15*math.pi
M = 2 # number of pico BS
a_total = 0
C = 1000 #Mbits
D = 5 #s
W = 10 #MHz

# q_f
def popularity(f):
    list = []
    for num in range(1,F-1):
        list.append(1/pow(num,v))
    return (1/pow(f,v))/sum(list)

# d_m is the distance between m and BS
# d_m must be less than 0.4821446
def b_m(d_m):
    x = pow(noise,2)/(P_0*pow(d_m,-alpha))
    return -math.log(2)/np.exp(x)/expi(-x)

# A_m : scale of coverage
def k_m():
    x = lamda*A_m
    return np.exp(-x)/(1-np.exp(-x))*(expi(x)-math.log(x)-euler)

def w_f(f):
    ans = 0.0
    for num in range(f, F+1):
        ans += popularity(num)*L
    return ans


def a_m(d_mList, d_elseList):
    r= symbols("r")
    func = 0
    
    for d_m_pico in d_mList:
        user_num = 0
        dist = P_n*d_m_pico**(-alpha)
        inner = math.log(2)*(2**r)*(noise**2)/dist
        
        for d_pico in d_elseList[user_num]:
            y = d_pico**(-alpha)
            inner += math.log(2)*(2**r)*P_n*y/(((2**r)-1)*P_n*y+P_n*y)
        
        inner *= exp(-(2**r)*(noise**2)/dist)

        for d_pico in d_elseList[user_num]:
            y = P_n*d_pico**(-alpha)
            inner *= y/(((2**r)-1)*y+y)
        
        inner *= lamda
        
        user_num += 1
        func += inner

    func *= r
    return k_m()*scipy.integrate.quad(lambda x:func.subs(r,x),0,np.Infinity)[0]

def a_init(userList, userElseList):
    term_2 = 0
    a = 0
    for f in range(1,F+1):
        term_2 += popularity(f)*L
    for pbs in range(0,M+1):
        a += a_m(userList[pbs], userElseList[pbs])
    a_total = pow(term_2*a,0.5)
    return a_total

def u_m(m, s_list, d_pbs_list):
    term_1 = 0
    for pbs in range(1,M+1):
        inner = 0
        if pbs == m:
            continue
        for f in range(1,F+1):
            inner += b_m(d_pbs_list[pbs])*popularity(f)*L*(1-s_list[pbs][f])
        term_1 += pow(inner,0.5)
    
    print(term_1+a_total)

    return term_1+a_total

def a_iteration(m, s_list, d_pbs_list):
    for f in range(1,F+1):
        J = popularity(f)*C-popularity(f)*L*(f-1)-w_f(f)
        b = b_m(d_pbs_list[m])
        u = u_m(m, s_list, d_pbs_list)
        q = popularity(f)

        z_0 = b**2*u*J
        z_1 = b**3*q*L*J+b**2*q**2*L*D*W*J
        z_2 = 2*b**2*u*q*L*J
        z_3 = -2*b**2*q*L*J
        z_4 = b*u*q*L
        z_5 = b*q*L

        coeff = [z_5, z_4, z_3, z_2, z_1, z_0]
        ans = np.roots(coeff)
        ans = ans.real[abs(ans.imag)<1e-5]
        S_root = abs(ans).max
        S_mf = (b*(q*L)*(F-f+1)-max_root**2)/b*q*L
        

u_list = [[0.05],[0.02,0.01],[]]
else_list = [[[0.5,0.5]],[[0.5,0.9],[0.3,0.7]],[]]
s_list = [[0],[0],[0]]
for x in range(0,1000):
    s_list[0].append(0)
    s_list[1].append(0)
    s_list[2].append(0)
a_total = a_init(u_list, else_list)
u_m(1, s_list, [0,0.3,0.4])
u_m(2, s_list, [0,0.3,0.4])