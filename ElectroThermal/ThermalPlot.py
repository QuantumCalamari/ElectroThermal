#matplotlib inline
import matplotlib.pyplot as plt
plt.style.use('seaborn-white')
import numpy as np

def f(x, y):
    return x*y

#x_fileName = "x.txt"
#t_fileName = "t.txt"

x_fileName = "../../x.txt"
t_fileName = "../../t.txt"

x_lineList = [line.rstrip('\n') for line in open(x_fileName)]
t_lineList = [line.rstrip('\n') for line in open(t_fileName)]


with open("../../u.txt") as textFile:
#with open("u.txt") as textFile:
    lines = [line.split() for line in textFile]

x = list(map(float, x_lineList))
y = list(map(float, t_lineList))

X, Y = np.meshgrid(x, y)
Z = lines

plt.contourf(X, Y, Z, 20, cmap='rainbow');
plt.colorbar();
#plt.show();

fig = plt.gcf()
fig.set_size_inches(18.5, 10.5)
fig.savefig('thermal.png', dpi=300)
fig.savefig('thermal.pdf', dpi=300)
