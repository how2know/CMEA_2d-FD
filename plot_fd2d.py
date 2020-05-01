import mpl_toolkits.mplot3d
from numpy import *
from pylab import *

u = loadtxt("u_fd.txt")
N = int(sqrt(u.shape[0]))
print("N=%d" % N)
u.shape = (N,N)
u = u.T; # correct transpose

fig = figure()
ax = fig.gca(projection="3d")
x,y = mgrid[0:1:N*1j, 0:1:N*1j]
ax.plot_surface(x,y,u,cmap=cm.viridis)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('u')
show()

#fig2 = figure()
#ax = fig2.gca(projection="3d")
#ax.plot_surface(x,y, sin(2*pi*x)*cos(2*pi*y),cmap=cm.viridis)
#ax.set_xlabel('x')
#ax.set_ylabel('y')
#ax.set_zlabel('u')
#show()


pcolormesh(x,y,abs(u - sin(2*pi*x)*cos(2*pi*y)))
colorbar()
title("Approximation error")
show()

n = int(N/2);
plot(x[:, n], u[:, n], label='u')
plot(x[:,n], sin(2*pi*x[:,n])*cos(2*pi*y[:,n]), label='exact')
title("Cut at y = " + str(round(1.*n/(N+1),4)))
legend()
show()

