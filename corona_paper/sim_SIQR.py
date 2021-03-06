import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib import animation
from itertools import combinations

class Particle:
    """A class representing a two-dimensional particle."""

    def __init__(self, x, y, vx, vy, type, radius=0.01, styles=None):
        """Initialize the particle's position, velocity, and radius.

        Any key-value pairs passed in the styles dictionary will be passed
        as arguments to Matplotlib's Circle patch constructor.

        """
        self.type=type
        self.r = np.array((x, y))
        self.v = np.array((vx, vy))
        self.radius = radius

        if(self.type == "I"):
            self.sickstart = 0
        else:
            self.sickstart=None
        if(self.type == "Q"):
            self.sickvel = np.array((0,0))
        else:
            self.sickvel = None

        self.setStyle()

    def quarantine(self, r):
        self.r = np.array(r)
        self.sickvel = self.v
        self.v = np.array((0,0))
    def setStyle(self):
        if self.type=="S":
            self.styles = {'color': 'b',  'fill': True}
        elif self.type == "I":
            self.styles = {'color': 'orange', 'fill': True}
        elif self.type == "Q":
            self.styles = {'color': 'red', 'fill': True}
        elif self.type == "R":
            self.styles = {'color': 'g', 'fill': True}
    # For convenience, map the components of the particle's position and
    # velocity vector onto the attributes x, y, vx and vy.
    @property
    def x(self):
        return self.r[0]
    @x.setter
    def x(self, value):
        self.r[0] = value
    @property
    def y(self):
        return self.r[1]
    @y.setter
    def y(self, value):
        self.r[1] = value
    @property
    def vx(self):
        return self.v[0]
    @vx.setter
    def vx(self, value):
        self.v[0] = value
    @property
    def vy(self):
        return self.v[1]
    @vy.setter
    def vy(self, value):
        self.v[1] = value

    def overlaps(self, other):
        """Does the circle of this Particle overlap that of other?"""

        return np.hypot(*(self.r - other.r)) < self.radius + other.radius

    def draw(self, ax):
        """Add this Particle's Circle patch to the Matplotlib Axes ax."""

        circle = Circle(xy=self.r, radius=self.radius, **self.styles)
        ax.add_patch(circle)
        return circle
    def setType(self, type, time):
        if self.type == type:
            pass
        else:
            self.type=type
            if type == "I":
                self.sickstart=time
            self.setStyle()

    def advance(self, dt):
        """Advance the Particle's position forward in time by dt."""

        self.r += self.v * dt

        # Make the Particles bounce off the walls
        if self.type != "Q":
            if self.x - self.radius < 0:
                self.x = self.radius
                self.vx = -self.vx
            if self.x + self.radius > 1:
                self.x = 1-self.radius
                self.vx = -self.vx
            if self.y - self.radius < 0:
                self.y = self.radius
                self.vy = -self.vy
            if self.y + self.radius > 1:
                self.y = 1-self.radius
                self.vy = -self.vy

class Simulation:
    """A class for a simple hard-circle molecular dynamics simulation.

    The simulation is carried out on a square domain: 0 <= x < 1, 0 <= y < 1.

    """

    def __init__(self, n, vr, p, types, dt, sim_time, sicknessperiod, incubationperiod, radius=0.01,  styles=None):
        self.nrofcollisions=0
        self.nrofcollisionsFor=0
        self.nrofcollisionsMid=0
        self.nrofcollisionsBack=0
        self.time=dt
        self.daytime=dt
        self.p = p
        self.Iarr=[]
        self.Sarr=[]
        self.Qarr=[]
        self.Tarr=[]
        self.Rarr=[]
        self.dt = dt
        self.sicknessperiod = sicknessperiod
        self.sim_time = sim_time
        self.nrofparticles = n
        self.incubationperiod = incubationperiod
        self.q=0
        """Initialize the simulation with n Particles with radii radius.

        radius can be a single value or a sequence with n values.

        Any key-value pairs passed in the styles dictionary will be passed
        as arguments to Matplotlib's Circle patch constructor when drawing
        the Particles.

        """

        self.init_particles(n,vr, types,radius, styles)

    def init_particles(self, n, vr, types, radius, styles=None):
        """Initialize the n Particles of the simulation.

        Positions and velocities are chosen randomly; radius can be a single
        value or a sequence with n values.

        """
        """
        try:
            iterator = iter(radius)
            assert n == len(radius)
        except TypeError:
            # r isn't iterable: turn it into a generator that returns the
            # same value n times.
            def r_gen(n, radius):
                for i in range(n):
                    yield radius
            radius = r_gen(n, radius)
        """
        rad = radius
        self.n = n
        self.particles = []
        for type in types:
            for i in range(0,type["num"]):
        #for i, rad in enumerate(radius):
            # Try to find a random initial position for this particle.
                while True:
                    # Choose x, y so that the Particle is entirely inside the
                    # domain of the simulation.
                    x, y = rad + (1 - 2*rad) * np.random.random(2)
                    # Choose a random velocity (within some reasonable range of
                    # values) for the Particle.
                    #vr = 0.1 * np.random.random() + 0.05
                    vphi = 2*np.pi * np.random.random()
                    vx, vy = vr * np.cos(vphi), vr * np.sin(vphi)
                    particle = Particle(x, y, vx, vy, type["type"], rad, styles)
                    # Check that the Particle doesn't overlap one that's already
                    # been placed.
                    for p2 in self.particles:
                        if p2.overlaps(particle):
                            break
                    else:
                        self.particles.append(particle)
                        break
    def unquarantine(self, p):
        while True:
            # Choose x, y so that the Particle is entirely inside the
            # domain of the simulation.
            x, y = p.radius + (1 - 2*p.radius) * np.random.random(2)
            # Choose a random velocity (within some reasonable range of
            # values) for the Particle.
            #vr = 0.1 * np.random.random() + 0.05
            vphi = 2*np.pi * np.random.random()
            particle = Particle(x, y, 0, 0, "R", p.radius)
            # Check that the Particle doesn't overlap one that's already
            # been placed.
            for p2 in self.particles:
                if p2.overlaps(particle):
                    break
            else:
                p.r = np.array((x,y))
                p.v = p.sickvel
                break
    def get_quarantine_space(self):

        cols = 3
        rows = 10
        if self.q >= cols*rows:
            self.q=0
        self.q +=1#states = self.countStates()
        q = self.q
        nra=q
        nr = nra-1
        cols = 3
        col = nr%3
        row = int((nr-col)/cols)

        x=1+0.4/(cols + 1)*(col+1)
        y=0+1/(rows+1)*(row+1)
        return (x,y)
    def countStates(self):
        states = [0,0,0,0]
        for p in self.particles:
            if p.type=="S":
                states[0]+=1
            elif p.type == "I":
                states[1]+=1
            elif p.type == "Q":
                states[2]+=1
            elif p.type == "R":
                states[3]+=1
        return states

    def handle_collisions(self):
        """Detect and handle any collisions between the Particles.

        When two Particles collide, they do so elastically: their velocities
        change such that both energy and momentum are conserved.

        """

        def change_velocities(p1, p2):
            self.nrofcollisions+=1
            """
            Particles p1 and p2 have collided elastically: update their
            velocities.

            """

            m1, m2 = p1.radius**2, p2.radius**2
            M = m1 + m2
            r1, r2 = p1.r, p2.r
            d = np.linalg.norm(r1 - r2)**2
            v1, v2 = p1.v, p2.v
            u1 = v1 - 2*m2 / M * np.dot(v1-v2, r1-r2) / d * (r1 - r2)
            u2 = v2 - 2*m1 / M * np.dot(v2-v1, r2-r1) / d * (r2 - r1)
            p1.v = u1
            p2.v = u2

        # We're going to need a sequence of all of the pairs of particles when
        # we are detecting collisions. combinations generates pairs of indexes
        # into the self.particles list of Particles on the fly.
        pairs = combinations(range(self.n), 2)
        for i,j in pairs:
            if self.particles[i].overlaps(self.particles[j]):
                if np.random.binomial(1,self.p) and ( (self.particles[i].type == "S" and self.particles[j].type=="I") or (self.particles[j].type == "S" and self.particles[i].type=="I") ):
                    self.particles[j].setType("I", self.time)
                    self.particles[i].setType("I",self.time)
                change_velocities(self.particles[i], self.particles[j])

    def advance_animation(self, dt):
        """Advance the animation by dt, returning the updated Circles list."""

        for i, p in enumerate(self.particles):
            p.advance(dt)
            self.circles[i].center = p.r
        self.handle_collisions()
        return self.circles

    def advance(self, dt):
        """Advance the animation by dt."""
        for i, p in enumerate(self.particles):
            p.advance(dt)
        self.handle_collisions()



    def init(self):
        """Initialize the Matplotlib animation."""

        self.circles = []
        for particle in self.particles:
            self.circles.append(particle.draw(self.ax))
        return self.circles
    def init2(self):#graphs
        #self.circles.append(self.ax2.plot(self.Sarr))
        #self.line.set_data([], [])
        #return self.line,
        self.line.set_data([], [])
        self.line2.set_data([], [])
        self.line3.set_data([], [])
        self.line4.set_data([], [])
        return self.line,self.line2, self.line3, self.line4,

    def animate2(self,i):
        #self.line.set_data(self.Tarr, self.Sarr)
        #return self.line,
        x = np.array(self.Tarr)
        y = np.array(self.Sarr)
        self.line.set_data(x, y)
        x2 = np.array(self.Tarr)
        y2 = np.array(self.Iarr)
        self.line2.set_data(x2, y2)
        x3 = np.array(self.Tarr)
        y3 = np.array(self.Qarr)
        self.line4.set_data(x3, y3)
        x4 = np.array(self.Tarr)
        y4 = np.array(self.Rarr)
        self.line3.set_data(x4, y4)
        return self.line, self.line2, self.line3, self.line4,

    def animate(self, i):
        """The function passed to Matplotlib's FuncAnimation routine."""
        #print("nr of collisions per time" , self.nrofcollisions/self.time)
        #meanvelocity = [(p.v[0]**2+p.v[1]**2)**0.5 for p in self.particles]
        #meanvelocity=sum(meanvelocity)/len(meanvelocity)
        #print("mean velocity", meanvelocity)
        #print("calculated col per time", 2*0.015*30*meanvelocity/1)
        if self.daytime>1:#one day elapsed
            self.nrofcollisionsBack=self.nrofcollisionsMid # dy/dx = y(x+h)-y(x-h)/2h
            self.nrofcollisionsMid=self.nrofcollisionsFor
            self.nrofcollisionsFor = self.nrofcollisions
            print("nr of collisions (1 day lag):", (self.nrofcollisionsFor-self.nrofcollisionsBack)/2)
            print(self.time)
            self.daytime = 0
            states = self.countStates()#tab the next lines to reduce to plotting daily
            self.Sarr.append(states[0])
            self.Iarr.append(states[1])
            self.Qarr.append(states[2])
            self.Rarr.append(states[3])
            self.Tarr.append(self.time)

        if self.time > self.sim_time:
            print("Sarr=",self.Sarr)
            print("Iarr=",self.Iarr)
            print("Qarr=", self.Qarr)
            print("Rarr=",self.Rarr)
            input()
            exit()
        #print("time ", self.time)
        #if self.time>5:
        #    print(self.Iarr)
        for p in self.particles:
            if p.type=="I" and p.sickstart+self.incubationperiod<self.time:
                p.setType("Q", self.time)
                p.quarantine(self.get_quarantine_space())
        for p  in self.particles:
            if p.type == "Q" and p.sickstart+self.incubationperiod+self.sicknessperiod<self.time:
                p.setType("R", self.time)
                self.unquarantine(p)
        self.advance_animation(self.dt)
        self.time +=self.dt
        self.daytime +=self.dt
        return self.circles#, self.myplot

    def do_animation(self, save=False):
        """Set up and carry out the animation of the molecular dynamics.

        To save the animation as a MP4 movie, set save=True.
        """

        fig, (self.ax,self.ax2) = plt.subplots(1,2)
        for s in ['top','bottom','left','right']:
            self.ax.spines[s].set_linewidth(2)
        self.ax.set_aspect('equal', 'box')
        self.ax.set_xlim(0, 1.4)#0.4 for quarantine
        self.ax.set_ylim(0, 1)
        self.ax.xaxis.set_ticks([])
        self.ax.yaxis.set_ticks([])
        self.ax.plot([1,1],[0,1], color='black')
        #fig2, self.ax2 = plt.subplots()
        #fig2 = plt.figure()
        #self.ax2 = plt.axes(xlim=(0, self.sim_time), ylim=(0, self.nrofparticles))
        self.ax2.set_xlim([0, self.sim_time])
        self.ax2.set_ylim([0, self.nrofparticles])
        self.line, = self.ax2.plot([], [], lw=3)
        self.line2, = self.ax2.plot([], [], lw=3)
        self.line3, = self.ax2.plot([], [], lw=3)
        self.line4, = self.ax2.plot([], [], lw=3)
        anim = animation.FuncAnimation(fig, self.animate, init_func=self.init,
                               frames=1, interval=1, blit=True)

        anim2 = animation.FuncAnimation(fig, self.animate2, init_func=self.init2,
                               frames=1, interval=1, blit=True)
        if save:
            Writer = animation.writers['ffmpeg']
            writer = Writer(fps=100, bitrate=1800)
            anim.save('collision.mp4', writer=writer)
        else:
            plt.show()


if __name__ == '__main__':
    radii = 0.02
    initvel = 0.15
    #stylesS = {'edgecolor': 'C0', 'linewidth': 2, 'fill': None}
    types = [{"type":"S", "num":49}, {"type":"I", "num":1},{"type":"Q", "num":0}, {"type":"R", "num":0}]
    nparticles = sum([x["num"] for x in types])
    #stylesI = {'edgecolor': '', 'linewidth': 2, 'fill': None}
    #stylesR = {'edgecolor': 'C0', 'linewidth': 2, 'fill': None}
    dt=0.05
    p=0.5
    incubationperiod = 5
    sicknessperiod = 10
    sim_time=50
    sim = Simulation(nparticles, initvel, p, types, dt, sim_time, sicknessperiod, incubationperiod, radii)
    sim.do_animation(save=False)
