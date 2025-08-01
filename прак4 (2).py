import numpy as np
import matplotlib.pyplot as plt
from sympy import *
from sympy import Symbol, diff, atan, cos, sqrt, sin, atan2
from sympy.parsing.sympy_parser import parse_expr
from sympy import lambdify
from matplotlib.animation import FuncAnimation
from matplotlib.animation import PillowWriter
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import proj3d
from mpl_toolkits.mplot3d import Axes3D

#Данные задачи
T = 30         #Время полета

x0 = 1         #Точка начального положения квадрокоптера
y0 = 1
z0 = 1

x_dif0 = 0      #Значения в 0 производной
y_dif0 = 0
z_dif0 = 0

x_T = -1        #Точка конечного положения квадрокоптера
y_T = -2
z_T = 3

x_difT = -1    #Значения в T производной
y_difT = -2
z_difT = 3

theta0 = 3.1415/2       #Значения углов в начальной точке
fi0 = 3.1415/2

theta_T = 3.1415/2      #Значения углов в конечной точке
fi_T = 3.1415/2

theta_dif0 = 3.1415/2   #Значения производной углов в 0
fi_dif0 = 3.1415/2

theta_difT = 3.1415/2  #Значения в T производной углов
fi_difT = 3.1415/2


l = 0.25        #Константы(взяты из статьи)
m = 2
k = 6.11*10e-8
g = 9.8


h01 = np.array([x0, x_dif0, -g*np.sin(theta0)/2, -g*np.cos(theta0)*theta_dif0/6,
                x_T, x_difT, -g*np.sin(theta_T), -g*np.cos(theta_T)*theta_difT]) 
h02 = np.array([y0, y_dif0, g*np.cos(theta0)*np.sin(fi0)/2, g*(fi_dif0*np.cos(theta0)*np.cos(fi0)-theta_dif0*np.sin(theta0)*np.sin(fi0))/6,
                y_T, y_difT, g*np.cos(theta_T)*np.sin(fi_T), g*(fi_difT*np.cos(theta_T)*np.cos(fi_T)-theta_difT*np.sin(theta_T)*np.sin(fi_T))])
h03 = np.array([z0, z_dif0, g*(np.cos(theta0)*np.cos(fi0)-1)/2, -g*(theta_dif0*np.sin(theta0)*np.cos(fi0)+fi_dif0*np.cos(theta0)*np.sin(fi0))/6,
                z_T, z_difT,g*(np.cos(theta_T)*np.cos(fi_T)-1), -g*(theta_difT*np.sin(theta_T)*np.cos(fi_T)+fi_difT*np.cos(theta_T)*np.sin(fi_T))])


h = "a7*t**7+a6*t**6+a5*t**5+a4*t**4+a3*t**3+a2*t**2+a1*t+a0"

h1 = parse_expr(h, evaluate=True)

dh1dt = diff(h1, Symbol('t'))
dh1dt2 = diff(dh1dt, Symbol('t'))
dh1dt3 = diff(dh1dt2, Symbol('t'))

a0, a1, a2, a3, t = symbols('a0 a1 a2 a3 t') 

h1 = h1.subs([(a0, h01[0]), (a1, h01[1]), (a2, h01[2]), (a3, h01[3]), (t, T)]) 

dh1dt = dh1dt.subs([(a1, h01[1]), (a2, h01[2]), (a3, h01[3]), (t, T)]) 
dh1dt2 = dh1dt2.subs([(a2, h01[2]), (a3, h01[3]), (t, T)]) 
dh1dt3 = dh1dt3.subs([(a3, h01[3]), (t, T)]) 


a4, a5, a6, a7 = symbols('a4 a5 a6 a7') 

res1 = linsolve([h1-h01[4], dh1dt-h01[5], dh1dt2-h01[6], dh1dt3-h01[7]], (a4, a5, a6,a7))


h1 = parse_expr(h, evaluate=True)
h1 = h1.subs([(a0, h01[0]), (a1, h01[1]), (a2, h01[2]), (a3, h01[3]), (a4, res1.args[0][0]), (a5, res1.args[0][1]), (a6, res1.args[0][2]), (a7, res1.args[0][3])]) 
print(h1)

h2 = parse_expr(h, evaluate=True)
dh2dt = diff(h2, Symbol('t'))
dh2dt2 = diff(dh2dt, Symbol('t'))
dh2dt3 = diff(dh2dt2, Symbol('t'))



h2 = h2.subs([(a0, h02[0]), (a1, h02[1]), (a2, h02[2]), (a3, h02[3]), (t, T)]) 

dh2dt = dh2dt.subs([(a1, h02[1]), (a2, h02[2]), (a3, h02[3]), (t, T)]) 

dh2dt2 = dh2dt2.subs([(a2, h02[2]), (a3, h02[3]), (t, T)]) 
dh2dt3 = dh2dt3.subs([(a3, h02[3]), (t, T)]) 


res2 = linsolve([h2-h02[4], dh2dt-h02[5], dh2dt2-h02[6], dh2dt3-h02[7]], (a4, a5, a6,a7))

h2 = parse_expr(h, evaluate=True)
h2 = h2.subs([(a0, h02[0]), (a1, h02[1]), (a2, h02[2]), (a3, h02[3]), (a4, res2.args[0][0]), (a5, res2.args[0][1]), (a6, res2.args[0][2]), (a7, res2.args[0][3])]) 
print(h2)

h3 = parse_expr(h, evaluate=True)
dh3dt = diff(h3, Symbol('t'))
dh3dt2 = diff(dh3dt, Symbol('t'))
dh3dt3 = diff(dh3dt2, Symbol('t'))

h3 = h3.subs([(a0, h03[0]), (a1, h03[1]), (a2, h03[2]), (a3, h03[3]), (t, T)]) 
dh3dt = dh3dt.subs([(a1, h03[1]), (a2, h03[2]), (a3, h03[3]), (t, T)]) 
dh3dt2 = dh3dt2.subs([(a2, h03[2]), (a3, h03[3]), (t, T)]) 
dh3dt3 = dh3dt3.subs([(a3, h03[3]), (t, T)]) 

res3 = linsolve([h3-h03[4], dh3dt-h03[5], dh3dt2-h03[6], dh3dt3-h03[7]], (a4, a5, a6,a7))

h3 = parse_expr(h, evaluate=True)
h3 = h3.subs([(a0, h03[0]), (a1, h03[1]), (a2, h03[2]), (a3, h03[3]), (a4, res3.args[0][0]), (a5, res3.args[0][1]), (a6, res3.args[0][2]), (a7, res3.args[0][3])]) 
print(h3)


dh1dt = diff(h1, Symbol('t'))
dh1dt2 = diff(dh1dt, Symbol('t'))
dh2dt = diff(h2, Symbol('t'))
dh2dt2 = diff(dh2dt, Symbol('t'))
dh3dt = diff(h3, Symbol('t'))
dh3dt2 = diff(dh3dt, Symbol('t'))

fi = atan2(dh2dt2, dh3dt2+g)      #Выводим фи, тета и u из формул статьи
theta = -1*atan(cos(fi)*dh1dt2/(dh3dt2+g))
u = m*sqrt(dh1dt2**2+dh2dt2**2+(dh3dt2+g)**2)

dfi_dt = diff(fi, Symbol('t'))   #Берем их производные 
dfi_d2t = diff(dfi_dt, Symbol('t'))
dtheta_dt = diff(theta, Symbol('t'))
dtheta_d2t = diff(dtheta_dt, Symbol('t'))

tau_psi = -0.02*dfi_d2t*sin(theta) - 0.02*cos(theta)*dtheta_dt*dfi_dt    #Считаем тау по формулам из статьи
tau_theta = 0.02*cos(fi)*dtheta_d2t + 0.02*sin(fi)*cos(theta)*dfi_d2t - 0.02*sin(fi)*dfi_dt*dtheta_dt + 0.02*(cos(fi)*cos(theta)*dfi_dt-sin(fi)*sin(theta)*dtheta_dt)*dfi_dt
tau_fi = -0.04*sin(fi)*dtheta_d2t + 0.04*cos(fi)*cos(theta)*dfi_d2t - 0.04*cos(fi)*dfi_dt*dtheta_dt - 0.04*(sin(fi)*cos(theta)*dfi_dt + cos(fi)*sin(theta)*dtheta_dt)*dfi_dt


f1 = k*(tau_psi/(4*k) - tau_fi/(2*l) + u/(4*k))
f2 = k*(-tau_psi/(4*k) + tau_fi/(2*l) + u/(4*k))
f3 = k*(tau_psi/(4*k) + tau_fi/(2*l) + u/(4*k))
f4 = k*(-tau_psi/(4*k) - tau_fi/(2*l) + u/(4*k))

f1 = lambdify(Symbol('t'),f1)
f2 = lambdify(Symbol('t'),f2)
f3 = lambdify(Symbol('t'),f3)
f4 = lambdify(Symbol('t'),f4)


times = np.linspace(0, T, 100)    #Графики сил тяги (все значения больше 0)

fig = plt.figure()
ax = fig.add_subplot()
ax.set_xlim(0,30)
ax.set_ylim(0,40)
ax.plot(times, f1(times),'b')
plt.show()
fig = plt.figure()
ax = fig.add_subplot()
ax.set_xlim(0,30)
ax.set_ylim(0,40)
ax.plot(times, f2(times),'r')
plt.show()
fig = plt.figure()
ax = fig.add_subplot()
ax.set_xlim(0,30)
ax.set_ylim(0,40)
ax.plot(times, f3(times),'g')
plt.show()
fig = plt.figure()
ax = fig.add_subplot()
ax.set_xlim(0,30)
ax.set_ylim(0,40)
ax.plot(times, f4(times),'b')
plt.show()


x = lambdify(Symbol('t'), h1)
y = lambdify(Symbol('t'), h2)
z = lambdify(Symbol('t'), h3)
fi_fun = lambdify(Symbol('t'), fi)
theta_fun = lambdify(Symbol('t'), theta)


fig = plt.figure()      #Высчитываем все точки траектории
ax = fig.add_subplot(projection='3d')
ax.scatter(x(0), y(0), z(0), 'g')
ax.set_xlabel('X')
ax.scatter(x(T), y(T), z(T), 'r')
ax.plot(x(times), y(times), z(times))
plt.show()


x1 = x(times)+np.cos(theta_fun(times))*l
x2 = x(times)-np.cos(theta_fun(times))*l
x3, x4 = x(times), x(times)
y1 = y(times)+np.sin(fi_fun(times))*np.sin(theta_fun(times))*l
y2 = y(times)-np.sin(fi_fun(times))*np.sin(theta_fun(times))*l
y3 = y(times)+l*np.cos(fi_fun(times))
y4 = y(times)-l*np.cos(fi_fun(times))
z1 = z(times)+l*np.sin(theta_fun(times))*np.cos(fi_fun(times))
z2 = z(times)-l*np.sin(theta_fun(times))*np.cos(fi_fun(times))
z3 = z(times)-l*np.sin(fi_fun(times))
z4 = z(times)+l*np.sin(fi_fun(times))


fig = plt.figure(figsize=(8,8), dpi=100)    #Анимация
ax = fig.add_subplot(projection='3d')
ax.grid()
x_pl, y_pl, z_pl = [], [], []

plx, = ax.plot3D([], [], [], color='k')
ply, = ax.plot3D([], [], [], color='k')
traj, = ax.plot3D([], [], [], color='red')
xdots, = ax.plot3D([], [], [], color='b',marker='o',linestyle='')
ydots, = ax.plot3D([], [], [], color='b',marker='o',linestyle='')
cdots, = ax.plot3D([], [], [], color='k',marker='o',linestyle='')

start_p, = ax.plot3D([],[],[], color='k',marker='o')
fin_p, = ax.plot3D([],[],[], color='k',marker='o')

ax.text(x0, y0, z0+0.5, 'Start') 
ax.text(x_T, y_T, z_T+0.5, 'Finish')     

ax.set_xlim3d([x0-4*l, x0+4*l])
ax.set_ylim3d([y0-4*l, y0+4*l])
ax.set_zlim3d([z0-4*l, z0+4*l])

ax.view_init(elev=5, azim=-175)

xs, ys, zs = x(times), y(times), z(times)
fi_angles, theta_angles = fi_fun(times), theta_fun(times)

def init():
    plx.set_data_3d([], [], [])
    ply.set_data_3d([], [], [])
    traj.set_data_3d([], [], [])
    xdots.set_data_3d([], [], [])
    ydots.set_data_3d([], [], [])
    cdots.set_data_3d([], [], [])
    start_p.set_data_3d([x0],[y0],[z0])
    fin_p.set_data_3d([x_T],[y_T],[z_T])
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.set_title("Квадрокоптер в полете")
    return traj, plx, ply, xdots, ydots, cdots,start_p, fin_p

def update(i):
    x_pl.append(xs[i])
    y_pl.append(ys[i])
    z_pl.append(zs[i])
    
    traj.set_data_3d(x_pl, y_pl, z_pl)
    
    plx.set_data_3d([x1[i],x2[i]], [y1[i], y2[i]], [z1[i], z2[i]])
    ply.set_data_3d([x3[i],x4[i]], [y3[i], y4[i]], [z3[i], z4[i]])
    xdots.set_data_3d([x1[i], x2[i]], [y1[i], y2[i]], [z1[i], z2[i]])
    ydots.set_data_3d([x3[i], x4[i]], [y3[i], y4[i]], [z3[i], z4[i]])

    cdots.set_data_3d([xs[i]], [ys[i]], [zs[i]])
    
    ax.set_xlim3d([xs[i]-4*l, xs[i]+4*l])
    ax.set_ylim3d([ys[i]-4*l, ys[i]+4*l])
    ax.set_zlim3d([zs[i]-4*l, zs[i]+4*l])
    
    ax.set_title(f'Квадрокоптер в полете\nt = {np.round(times[i],2)}')
   
    return traj, plx, ply, xdots, ydots, cdots,


ani = FuncAnimation(fig, update, frames=1000, init_func=init, blit=False, interval=1, repeat=False)
#writer = PillowWriter(fps=30)
#ani.save('прак42.gif', writer=writer)
plt.show()
