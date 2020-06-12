# post processing
import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.ioff()
import os
from matplotlib import animation

from dedalus import public as de

#######################################
'''
make_movie(): modified from Pulkit Dubey's function
inputs:
    xrange - array containing min and max x values
    yrange - array containing min and max y values
    name - a string containing the name and format, e.g. "movie.mp4"
    show - 0 on 1 depending on whether a popup window for the movie is desired

output:
    An mp4 / gif format video saved as <name>.mp4 in the root directory
    A popup window of the animation if show is 1 (default)

example call:
    # all default options
    make_movie()
    # all specified options
    make_movie(10, 2, Lx, Nx, [0, 5], varname="myvarname", name="mymovie.mp4", show=1)
'''
#######################################
def make_dedalus_movie(Nt, nsave, Lx, Nx, yr=[-1, 1], varname="u", name="movie", show=1):

    x_basis = de.Fourier('x', Nx, interval=(-Lx/2, Lx/2), dealias=3/2)
    domain = de.Domain([x_basis], np.float64)

    x = domain.grid(0)

    framerate = 20
    intvl = 20
    list = os.listdir("data/")
    file_count = len(list)

    fig = plt.figure()
    ax = plt.axes(xlim=(-Lx/2, Lx/2), ylim=(yr[0], yr[1]))
    ax.set_xlabel(r'$x$')  # Set x label
    ax.set_ylabel(r'$c$')  # Set y label
    plt.title(name)

    line, = ax.plot([], [], lw=2)

    def init():
        line.set_data([], [])
        return line,


    def animate(j):
        n = int(j*nsave)
        if(n <= Nt):
            n_str = str(n)
            filename = str(varname) + n_str
            hf = h5py.File('data/' + filename + '.h5', 'r')
            ut = hf.get(str(varname))
            u = np.array(ut)

            line.set_data(x, u)

        return line,

    anim = animation.FuncAnimation(
        fig, animate, init_func=init, frames=file_count, interval=intvl, blit=True)

    anim.save(str(name)+".mp4", fps=framerate, extra_args=['-vcodec', 'libx264'])

    if show==1:
        plt.show()

    plt.close()
