import matplotlib.pyplot as plt
import numpy as np

def modulus(x): return np.sqrt(x.real**2 + x.imag**2)
def f(x, y): return modulus(1 + x + y*1j + 0.5 * (x + y*1j)**2 + (x + y*1j)**3 / 6)

x = np.linspace(-5, 5, 1000)
y = np.linspace(-5, 5, 1000)

X, Y = np.meshgrid(x, y)
Z = f(X, Y)

cs = plt.contour(X, Y, Z, [1], colors='black')

plt.show()

p = cs.collections[0].get_paths()[0]

f = open("stab_rk3.plt", "w")
for i in p.vertices:
  f.write(str(i[0]) + "\t" + str(i[1]) + "\n")
f.close()
