import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import scipy.integrate as ode
import time

def initialParticleState(nparticles):
    state = np.zeros(nparticles*4)
    particlepositions = np.array([  0, 0,
                                    1, 0,
                                    -1, 0
                                    ])
    particlevelocities = np.array([ 0, 0,
                                    0, 1,
                                    0, -1
                                    ])
    particlevelocities = particlevelocities * 2
    state[::2] = particlepositions
    state[1::2] = particlevelocities
    particlemasses = np.array([10, 1, 2])

    if len(particlepositions)/2 != nparticles or len(particlevelocities)/2 != nparticles or len(particlemasses) != nparticles:
        print("Error - initial states size does not match particle size")
        pass
    return state, particlemasses

def printParticleState(state, n):
    print("particle {} \t x {}  {} \t v {}  {}"\
    .format(n, state[4*n], state[4*n + 2], state[4*n+1], state[4*n+3]))

def printHistory(variable):
    imax = variable.shape[0]
    for i in range(imax):
        currentstate = variable[i,:]
        for k in currentstate:
            print(k, end = '\t')
        print()

def gnuplotOutput(t,x,y):
    t = np.around(t, 3)
    x = np.around(x, 3)
    y = np.around(y, 3)
    imax = len(t)
    for i in range(imax):
        print(t[i], end = '\t')
        for k in range(nparticles):
            print(x[i,k], "\t", y[i,k], end = '\t')
        print()

def particlePositionsandVelocities(state, plotaswell):
    # Convert state to column vectors
    x = state[::4]
    y = state[2::4]
    vx = state[1::4]
    vy = state[3::4]

    if plotaswell:
        k = 1
        plt.quiver(x,y,vx,vy)
        plt.grid()
        plt.xlim(k * (np.min(x) - np.max(np.abs(vx))), k * (np.max(x) + np.max(np.abs(vx))))
        plt.ylim(k * (np.min(y) - np.max(np.abs(vy))), k * (np.max(y) + np.max(np.abs(vy))))
        plt.show()
    return x, y, vx, vy

def distance(state, n1, n2):
    # returns the distance between two particles
    return np.sqrt((state[4*n1] - state[4*n2])**2 + (state[4*n1 + 2] - state[4*n2 + 2])**2)

def acceleration(state):
    # returns the 2nd derivative of all particles as a result of the interactions with all other particles.
    G = -1

    nparticles = int(len(state)/4)
    y2dot = np.zeros(nparticles * 2)

    for i in range(nparticles): # active particle
        for k in range(nparticles): # other particles
            if k != i:
                # GM/r^3 we reject the other mass, since are only concerned with acceleration
                temp = G * masses[k] * distance(state, i, k)**(-3)
                # y2dotx
                y2dot[2*i] += temp * (state[4*i] - state[4*k])
                # y2doty
                y2dot[2*i + 1] += temp * (state[4*i + 2] - state[4*k + 2])
                # print("y2dot \t {}".format(y2dot))
    return y2dot

def statederivative(t, y):
    """
    Returns the derivative of the space state variables.
    Inputs:     t - A scalar time value
                y - An array of the state space of the system, length 4xnparticles since each particle is
                    uniquely defined by 4 variables.
    Outputs:    ydot - A derivative of y for the pendulum problem
    """
    ydot = np.empty(len(y))
    ydot[::2] = y[1::2] # set the derivate of the position of the particles equal to their velocity
    particleacceleration = acceleration(y)
    ydot[1::2] = particleacceleration
    return ydot

def solve(initialstate, dt):
    t0 = time.time()
    tmax = 100
    solution = ode.RK45(statederivative, t0 = 0.0, y0 = state, t_bound = tmax, max_step = dt)

    N = 800 # number of steps
    t = np.zeros(N)
    x = np.zeros([N, nparticles])
    y = np.zeros([N, nparticles])
    vx = np.zeros([N, nparticles])
    vy = np.zeros([N, nparticles])

    # generate data
    for i in range(N):
        t[i] = solution.t
        x[i,:] = solution.y[::4]
        y[i,:] = solution.y[2::4]
        vx[i,:] = solution.y[1::4]
        vy[i,:] = solution.y[3::4]
        solution.step()

    print("Data assembled in {}s".format(time.time() - t0))
    return t, x, y

def constanttimeintervals(t,x,y, dt):
    N = len(t)

    tmax = t[-1]
    Nnew = int(tmax/dt)
    tnew = np.zeros(Nnew)
    xnew = np.zeros([Nnew, nparticles])
    ynew = np.zeros([Nnew, nparticles])

    i = 0 # write to new array
    ilower = 0 # lower index in old array
    while ilower < N-1:
        timegap = 0
        iupper = ilower # generate upper index in old array
        while timegap < dt * 0.99: # timegap starts as zero
            iupper += 1
            if iupper == N:
                break
            timegap = t[iupper] - t[ilower]
        tnew[i] = t[ilower]
        xnew[i,:] = x[ilower,:]
        ynew[i,:] = y[ilower,:]
        i += 1
        ilower = iupper
    tnew = tnew[:-1] # for some reason the last tdiff is negative, so we remove.
    xnew = xnew[:-1]
    ynew = ynew[:-1]
    return tnew, xnew, ynew

def plotanimation(x,y, saveasmp4):
    t0 = time.time()
    fig, ax = plt.subplots()
    ntail = 30
    axs = [] # list of images to be plotted
    for i in range(len(x)):
        imin = 0
        if i >= ntail:
            imin = i - ntail
        imax = i

        ax = plt.plot(x[imin: i + 1, 0], y[imin: i + 1, 0], 'b-')
        ax += plt.plot(x[i, 0], y[i, 0], 'ro')

        for k in range(1, nparticles):
            ax += plt.plot(x[imin: i + 1, k], y[imin: i + 1, k], 'b-')
            ax += plt.plot(x[i, k], y[i, k], 'ro')
        axs.append(ax)

    plt.xlim(-2,2)
    plt.ylim(-2,2)

    im_ani = animation.ArtistAnimation(fig, axs, interval = 50, repeat_delay = 3000, blit = True)

    if saveasmp4:
        Writer = animation.writers['ffmpeg']
        writer = Writer(fps=15, metadata=dict(artist='Oliver Normand'), bitrate=1800)
        im_ani.save('im1.mp4', writer = writer)
        print("Animation save time {}".format(time.time() - t0))
    if not saveasmp4:
        print("Animation generated in {}s".format(time.time() - t0))
        plt.show()


nparticles = 3
dt = 0.005
state, masses = initialParticleState(nparticles) # x1, vx1, y1, vy1, x2, vx2, y2, vy2 ...
t, x, y = solve(state, dt)
t, x, y = constanttimeintervals(t,x,y,dt)

plotanimation(x,y, True)
