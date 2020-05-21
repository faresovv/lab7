from numpy import array, sin, cos, exp, sqrt
from scipy import integrate
import matplotlib.pyplot as plt

data = []
with open("ans.dat") as file:
    for line in file:
        data.append([float(line.split(',')[0]), float(line.split(',')[1]), float(line.split(',')[2])])

data = array(data)

x = data[:, 2]

y = (609/1202 - 733*sqrt(5)/3606)*exp(4*x*(-2 - sqrt(5))) + (609/1202 + 733*sqrt(5)/3606)*exp(4*x*(-2 + sqrt(5))) - 8*cos(4*x)*exp(x)/601 - 5*exp(x)*sin(4*x)/1803
dy = -44*cos(4*x)*exp(x)/1803 + 91*exp(x)*sin(4*x)/1803 + (-8 - 4*sqrt(5))*(609/1202 - 733*sqrt(5)/3606)*exp((4*x)*(-2 - sqrt(5))) + (-8 + 4*sqrt(5))*(609/1202 + (733*sqrt(5))/3606)*exp((4*x)*(-2 + sqrt(5)))

def sp(y, x):
    return [y[1], sin(4*x)*exp(x) - 16*y[1] + 16*y[0]]

ys = integrate.odeint(sp, [1, 0], x)[:, 0]
dys = integrate.odeint(sp, [1, 0], x)[:, 1]

plt.figure(figsize=(16, 9))
plt.title("Решение", fontsize=14)
plt.plot(x, data[:, 0], label="$y(x)$")
plt.plot(x, data[:, 1], label="$y'(x)$")
plt.legend(fontsize=14)
plt.minorticks_on()
plt.grid(which="both")
plt.show()

plt.figure(figsize=(16, 9))
plt.title("Фазовый портрет", fontsize=14)
plt.plot(data[:, 0], data[:, 1])
plt.xlabel("$y(x)$", fontsize=14)
plt.ylabel("$y'(x)$", fontsize=14)
plt.minorticks_on()
plt.grid(which="both")
plt.show()

plt.figure(figsize=(16, 9))
plt.title("Отклонение от точного решения", fontsize=14)
plt.plot(x, abs(y - data[:, 0]), label="$ \Delta y(x)$")
plt.plot(x, abs(dy - data[:, 1]), label="$ \Delta y'(x)$")
plt.legend(fontsize=14)
plt.minorticks_on()
plt.grid(which="both")
plt.show()

plt.figure(figsize=(16, 9))
plt.title("Отклонение от решения scipy.integrate.odeint", fontsize=14)
plt.plot(x, abs(ys - data[:, 0]), label="$ \Delta y(x)$")
plt.plot(x, abs(dys - data[:, 1]), label="$ \Delta y'(x)$")
plt.legend(fontsize=14)
plt.minorticks_on()
plt.grid(which="both")
plt.show()
