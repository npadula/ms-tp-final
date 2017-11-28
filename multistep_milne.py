from math import sqrt
 
def rk4(f, x0, y0, x1, n):
    xs = [0] * (n + 1)
    ys = [0] * (n + 1)
    h = (x1 - x0) / float(n)
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
 
def f(x, y):
    return x * sqrt(y)

def fej(x,y):
    return (x*y)/10.0
 
xs, ys = rk4(fej, 1, 0.2, 5, 4)


for x,y in zip(xs,ys):
    print(fej(x,y))