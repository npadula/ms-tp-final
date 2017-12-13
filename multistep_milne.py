from math import sqrt
import numpy
import matplotlib.pyplot as plt
 
#Runge-Kutta de orden 4
def rk4(f, x0, y0, xfinal, n):
    xs = [0] * (n + 1)
    ys = [0] * (n + 1)
    h = (xfinal - x0) / float(n)
    xs[0] = x = x0
    ys[0] = y = y0
    for i in range(1, n + 1):
        k1 = h * f(x, y)
        k2 = h * f(x + 0.5 * h, y + 0.5 * k1)
        k3 = h * f(x + 0.5 * h, y + 0.5 * k2)
        k4 = h * f(x + h, y + k3)
        xs[i] = x = x0 + i * h
        ys[i] = y = y + (k1 + k2 + k2 + k3 + k3 + k4) / 6
    return xs, ys

#Calcula un paso de Milne para un h dado, f es dy/dt, xs, ys y fs son listas de valores de x, y e y' respectivamente
def milne_step(f,xs,ys,fs,h):
    size = len(fs) - 1
    xnuevo = xs[size] + h
    yp = ys[size-3] + ((4/3)*h) * (2*fs[size-2] - fs[size-1] + 2*fs[size])
    fyp = f(xnuevo,yp)
    yc = ys[size-1] + (h/3) * (fs[size-1] + 4*fs[size] + fyp)
    
    #Estimamos error (Richardson)
    eyc = (yc-yp)/-29
    #Mejoramos el valor
    yce = yc + eyc
    return xnuevo,yce,f(xnuevo,yce) #Devolvemos los datos para completar una nueva columna de la tabla

def milne(f,xs,ys,fs,xfinal):
    size = len(xs) - 1
    h = xs[size] - xs[size-1]
    while(xs[size] < xfinal):
        xnuevo,ynuevo,fnuevo = milne_step(f,xs,ys,fs,h)
        xs.append(xnuevo)
        ys.append(ynuevo)
        fs.append(fnuevo)
        size+=1


def f(x,y):
    return (x*y)/10.0

#Solucion analitica del PVI
def sol_analitica(x):
    return 0.1902*numpy.exp((numpy.power(x,2))/20)


#Tomamos como punto final x=10, con 1000000 puntos en total
#Calculamos los valores con el R-K
xs, ys = rk4(f, 1, 0.2, 10, 9)

#Tomamos los tres primeros valores para inicializar el P-C de Milne
xsm = xs[0:3]
ysm = ys[0:3]

#Calculamos los valores de las derivadas en los puntos
fs = []
for x,y in zip(xsm,ysm):
    fs.append(f(x,y))
    
#Resolvemos por Milne
milne(f,xsm,ysm,fs,10)


#Solucion exacta
yexacta = []
err_RK = []
err_milne = []
for (xrk,xmilne,yrk,ymilne) in zip(xs,xsm,ys,ysm):
    yex = sol_analitica(xrk)
    
    yexacta.append(yex)
    erk = numpy.abs(yex - yrk )
    emilne = numpy.abs(yex - ymilne)
    err_RK.append(erk)
    err_milne.append(emilne)
    
    
    
#Grafico de los datos obtenidos
#y(x)
figure, axes = plt.subplots(2,1)
axes[0].set_xlabel('x')
axes[0].set_ylabel('y(x)')
#axes[0].set_yscale('log',basey=10)
axes[0].plot(xs,ys, color='blue', label='R-K4')
axes[0].plot(xsm,ysm, color='red', label='Milne P-C')
axes[0].plot(xs,yexacta, color='yellow', label='Solucion Analitica')


#|error|
axes[1].set_xlabel('x')
axes[1].set_ylabel('|Error|')
axes[1].plot(xs,err_RK, color='blue', label='R-K4')
axes[1].plot(xsm,err_milne,color='red',label='Milne P-C')
axes[1].plot([],[],color='yellow',label='Solucion Analítica')


plt.legend(bbox_to_anchor=(0.5, 2.45), loc=2, borderaxespad=0.)
plt.show()


