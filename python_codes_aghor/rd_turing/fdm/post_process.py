#######################################
# spin-up Hyun and Park JFM, 1992;
# A compact FDM on staggered grid for Navier-Stokes flows, Zhang et. al., Int. J. Numer. Meth. Fluids. 2006;
# Author: Pratik Aghor
import numpy as np
import interpolate
from interpolate import *
import params
from params import *

#Import plotting functions:
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import animation
# from matplotlib.animation import FuncAnimation
#######################################
def make_movie(Nt, nsave, x_full, y_full, Lx, Ly, Nx_full, Ny_full, varname="u", name="movie", show=1):
    X, Y = np.meshgrid(x_full, y_full)
    color_map="coolwarm"
    file_count = int(Nt/nsave)
    framerate = 20
    intvl = 20
    fig = plt.figure()
    ax = plt.axes(xlim=(0, Lx), ylim=(0, Ly))
    ax.set_xlabel(r'$x$')  # Set x label
    ax.set_ylabel(r'$y$')  # Set y label

    ims = []
    for nt in range(0, Nt+1, nsave):
        u = np.loadtxt('data/'+varname+str(nt)+'.asc')
        im = ax.contourf(X, Y, u, cmap=color_map)
        add_arts = im.collections
        ims.append(add_arts)

    anim = animation.ArtistAnimation(fig, ims, interval=50, blit=True,
                                repeat_delay=1000)


    # anim = animation.FuncAnimation(
    #     fig, animate, init_func=init, frames=file_count, interval=intvl, blit=True)

    anim.save(str(name)+".mp4", fps=framerate, extra_args=['-vcodec', 'libx264'])

    if show==1:
        plt.show()

    plt.close()
#######################################
# X, Y = np.meshgrid(x_full, y_full)
# color_map="coolwarm"
#
# for nt in range(0, Nt+1):
#     if(nt % nsave == 0):
#         u = np.loadtxt('data/u'+str(nt)+'.asc')
#         # v = np.loadtxt('data/v'+str(nt)+'.asc')
#
#         fig = plt.figure(1)  # Create a figure instance
#         ax = fig.gca()  # projection='3d' to Get current axes in 3D projection
#         ax.contourf(X, Y, u, cmap=color_map)
#         ax.set_xlabel(r'$x$')  # Set x label
#         ax.set_ylabel(r'$y$')  # Set y label
#         plt.savefig('data/u'+str(nt)+'_contours.png')
#
#         # fig = plt.figure(2)  # Create a figure instance
#         # ax = fig.gca()  # projection='3d' to Get current axes in 3D projection
#         # ax.contourf(Y, X, v, cmap=color_map)
#         # ax.set_xlabel(r'$x$')  # Set x label
#         # ax.set_ylabel(r'$y$')  # Set y label
#         # plt.savefig('data/v'+str(nt)+'_contours.png')

#######################################
make_movie(Nt, nsave, x_full, y_full, Lx, Ly, nx+1, ny+1, varname="u", name="u-movie", show=1)
make_movie(Nt, nsave, x_full, y_full, Lx, Ly, nx+1, ny+1, varname="v", name="v-movie", show=0)
#######################################
