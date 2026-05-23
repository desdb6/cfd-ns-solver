# Class to simulate an incompressible Newtonian fluid

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import time

class NewtonianFluid:
    """A class simulating Newtonian fluids."""
    def __init__(self, Lx, Ly, dx, dy, dt, mu=1, rho=1, U=0, Fx=0, Fy=0): 
        """Initialisation of dimensions of the system."""
        ### Typechecks ###
        typeCheck([Lx, Ly, dx, dy, dt, mu, rho, U, Fx, Fy], float)
        posCheck([Lx, Ly, dx, dy, dt, mu, rho])
        try:
            if not U>=0:
                raise ValueError
        except:
            raise ValueError('Error, U must be a numerical value grater than zero!')

        self._Lx=Lx
        self._Ly=Ly
        self._dx=dx
        self._dy=dy
        self._dt=dt
        self._nx=int(np.floor(Lx/dx))
        self._ny=int(np.floor(Ly/dy))

        self._p=np.zeros([self.nx, self.ny])
        self._u=np.zeros([self.nx, self.ny])
        self._v=np.zeros([self.nx, self.ny])

        self._mu=mu
        self._rho=rho
        self._U=U
        self._reynolds_number=np.max([self.Lx, self.Ly])*self.U*self.rho/self.mu
        self._Fx=Fx
        self._Fy=Fy

        self._obstacles=[]
        self._iterations=0

        self.set_boundary_conditions_p()
        self.set_boundary_conditions_u()
        self.set_boundary_conditions_v()

    ### Properties ###

    @property
    def Lx(self):
        """ Returns the x-dimention of the system. """
        return self._Lx  
        
    @property
    def Ly(self):
        """ Returns the y-dimention of the system. """
        return self._Ly
     
    @property
    def dx(self):
        """ Returns the x-spacing of the grid. """
        return self._dx
  
    @property
    def dy(self):
        """ Returns the y-spacing of the grid. """
        return self._dy
      
    @property
    def dt(self):
        """ Returns the timestep of the system. """
        return self._dt
      
    @property
    def mu(self):
        """ Returns the viscosity of the fluid. """
        return self._mu

    @property
    def rho(self):
        """ Returns the density of the fluid. """
        return self._rho
  
    @property
    def U(self):
        """ Returns the characteristic speed of the system. """
        return self._U
      
    @property
    def Fx(self):
        """ Returns the horizontal force on the fluid. """
        return self._Fx

    @property
    def Fy(self):
        """ Returns the vertical force on the fluid. """
        return self._Fy

    @property
    def nx(self):
        """ Returns the horizontal grid size of the system. """
        return self._nx
    
    @property
    def ny(self):
        """ Returns the vertical grid size of the system. """
        return self._ny
    
    @property
    def p(self):
        """ Returns the pressure field. """
        return self._p
    
    @property
    def u(self):
        """ Returns the horizontal speed field. """
        return self._u
    
    @property
    def v(self):
        """ Returns the vertical speed field. """
        return self._v

    @property
    def zeta(self):
        """ Returns the pressure field. """
        return self._zeta

    @property
    def reynolds_number(self):
        """ Returns the Reynolds Number. """
        return self._reynolds_number
    
    @property
    def obstacles(self):
        """ Returns the obstacles of the system. """
        return self._obstacles
    
    @property
    def p_BVC(self):
        """ Returns the pressure bounadry conditions of the system. """
        return self._p_BVC
    
    @property
    def u_BVC(self):
        """ Returns the horizontal speed boundary conditions of the system. """
        return self._u_BVC
    
    @property
    def v_BVC(self):
        """ Returns the vertical speed boundary conditions of the system. """
        return self._v_BVC
    
    @property
    def iterations(self):
        """ Returns the iterations of the system. """
        return self._iterations

### Initializing Methods ###

    def set_parameters(self, dx=None, dy=None, dt=None, mu=None, rho=None, U=None, Fx=None, Fy=None):
        """Sets parameters for the simulation."""
        ### Set default values ###
        if dx == None:
            dx=self.dx
        if dy == None:
            dy=self.dy
        if dt == None:
            dt=self.dt
        if mu == None:
            mu=self.mu
        if rho == None:
            rho=self.rho
        if U == None:
            U=self.U
        if Fx == None:
            Fx=self.Fx
        if Fy == None:
            Fy=self.Fy

        ### Typechecks ###
        typeCheck([dx, dy, dt, mu, rho, U, Fx, Fy], float)
        posCheck([dx, dy, dt, mu, rho])
        try:
            if not U>=0:
                raise ValueError
        except:
            raise ValueError('Error, U must be a numerical value grater than zero!')

        self._dx=dx
        self._dy=dy
        self._dt=dt
        self._mu=mu
        self._rho=rho
        self._U=U
        self._reynolds_number=self.Lx*self.U*self.rho/self.mu
        self._Fx=Fx
        self._Fy=Fy

        self._nx=int(np.floor(self.Lx/self.dx))
        self._ny=int(np.floor(self.Ly/self.dy))

        for obstacle in self.obstacles:
            obstacle._update_grid(self.Lx, self.Ly, self.dx, self.dy)

        print(f'Parameters set at: dx={self.dx}, dy={self.dy}, dt={self.dt}, mu={self.mu}, rho={self.rho}, U={self.U}, Fx={self.Fx}, Fy={self.Fy}, Reynolds number={self.reynolds_number}.\n')

    def set_boundary_conditions_p(self, BVC_up="Neumann", BVC_down="Neumann", BVC_left="Neumann", BVC_right="Neumann"):
        """Sets boundary conditions for the pressure."""
        ### Typechecks ###
        bvcCheck([BVC_up, BVC_down, BVC_left, BVC_right])
        self._p_BVC={'up':BVC_up, 'down':BVC_down, 'left':BVC_left, 'right':BVC_right}

    def set_boundary_conditions_u(self, BVC_up=0, BVC_down=0, BVC_left=0, BVC_right=0):
        """Sets boundary conditions for the horizontal speed."""
        ### Typechecks ###
        bvcCheck([BVC_up, BVC_down, BVC_left, BVC_right])
        self._u_BVC={'up':BVC_up, 'down':BVC_down, 'left':BVC_left, 'right':BVC_right}

    def set_boundary_conditions_v(self, BVC_up=0, BVC_down=0, BVC_left=0, BVC_right=0):
        """Sets boundary conditions for the vertical speed."""
        ### Typechecks ###
        bvcCheck([BVC_up, BVC_down, BVC_left, BVC_right])
        self._v_BVC={'up':BVC_up, 'down':BVC_down, 'left':BVC_left, 'right':BVC_right}

    def add_obstacle(self):
        """Adds an obstacle to the system."""
        obstacle=NoSlipObstacle(self)
        self._obstacles.append(obstacle)
        return obstacle

    def check_stability(self):
        stability= (self.dt<=(np.min([self.dx, self.dy]))**2/(2*self.mu+self.U*np.min([self.dx, self.dy])))
        if stability==True:
            print(f"Stability expected, dt={self.dt} <= {(np.min([self.dx, self.dy]))**2/(2*self.mu+self.U*np.min([self.dx, self.dy]))}.\n")
            return True
        if stability==False:
            print(f"Stability not expected, dt={self.dt} > {(np.min([self.dx, self.dy]))**2/(2*self.mu+self.U*np.min([self.dx, self.dy]))}.\n")
            return False

    def update_dt(self, factor=10):
        ### Typechecks ###
        typeCheck([factor], float)
        posCheck([factor])

        ### Update time step ###
        self._dt=1/factor*(np.min([self.dx, self.dy]))**2/(2*self.mu+self.U*np.min([self.dx, self.dy]))
        return
    

### Solving Methods ###

    def simulate_system(self, tol=5*10**(-4), max_iteratrions_initial=500, max_iterations=10000):
        """Calculate solution to Navier-Stokes equations"""
        ### Typechecks ###
        typeCheck([max_iterations, max_iteratrions_initial], int)
        typeCheck([tol], float)
        posCheck([tol, max_iterations, max_iteratrions_initial])

        print(f'Simulating system')

        ### Reset ###
        self._p=np.zeros([self.nx, self.ny])
        self._u=np.zeros([self.nx, self.ny])
        self._v=np.zeros([self.nx, self.ny])
        self._apply_all_boundary_conditions()
        for obstacle in self.obstacles:
            self._apply_all_obstacle_boundary_conditions(obstacle)

        begin_time=time.time()

        try:

            ### Find initial speed solution ###
            converged=False
            iterations=0
            while (converged==False and iterations<max_iteratrions_initial):
                ### Apply discretisation scheme ###
                u_temp=np.copy(self.u)
                self._update_u()

                v_temp=np.copy(self.v)
                self._update_v()


                ### BVC ###
                self._apply_all_boundary_conditions()
                for obstacle in self.obstacles:
                    self._apply_all_obstacle_boundary_conditions(obstacle)

                ### Check convergence ###
                u_diff=np.max(abs((self.u-u_temp))/abs((self.u + tol)))
                v_diff=np.max(abs((self.v-v_temp))/abs((self.v + tol)))
                if np.all([u_diff<tol, v_diff<tol]):
                    converged=True

                iterations+=1

            ### Solve system ###
            converged=False
            while (converged==False and iterations<max_iterations):
                ### Apply discretisation scheme ###
                self._update_p()
                
                u_temp=np.copy(self.u)
                self._update_u()
                
                v_temp=np.copy(self.v)
                self._update_v()
                
                ### BVC ###
                self._apply_all_boundary_conditions()
                for obstacle in self.obstacles:
                    self._apply_all_obstacle_boundary_conditions(obstacle)

                ### Check convergence ###
                u_diff=np.max(abs((self.u-u_temp))/abs((self.u + tol)))
                v_diff=np.max(abs((self.v-v_temp))/abs((self.v + tol)))
                if np.all([u_diff<tol, v_diff<tol]):
                    converged=True
                    print(f'Converged after {iterations} iterations, computation time: {time.time()-begin_time}s.\n')

                iterations+=1
            
            self._iterations=iterations-1

            

            if iterations==max_iterations:
                print(f'No convergence reached after {iterations} iterations, largest change equals {np.max([u_diff, v_diff])}, computation time: {time.time()-begin_time}s.\n')
        except:
            self._iterations=iterations-1
            print(f'Error, divergence occured after {self._iterations} iterations!\n')
            self._p=np.zeros([self.nx, self.ny])
            self._u=np.zeros([self.nx, self.ny])
            self._v=np.zeros([self.nx, self.ny])
        ### Calculate vorticity ###
        self._calculate_vorticity()

    def update_reynolds_number(self):
        U=np.average(np.sqrt(self.u**2+self.v**2))
        self._reynolds_number=np.max([self.Lx, self.Ly])*U*self.rho/self.mu
        print(f"Reynolds number updated, Re={self.reynolds_number}.\n")

### Plotting Methods ###
       
    def plot_quiver_speed(self):
        """Makes a quiver plot of the system."""
        ### Create meshgrids ###
        x = np.linspace(0, self.Lx, self.nx)
        y = np.linspace(0, self.Ly, self.ny)
        X, Y = np.meshgrid(x, y)

        fig, ax = plt.subplots() # Create figure

        contour_filled = ax.contourf(X, Y, np.sqrt(self.u**2+self.v**2).T, levels=50, cmap='rainbow_r') # Plot filled contours
        ax.quiver(X, Y, self.u.T, self.v.T) # Quiver plot
        for obstacle in self.obstacles:  # Draw obstacles
            rect = patches.Rectangle((obstacle.x0, obstacle.y0), (obstacle.x1-obstacle.x0), (obstacle.y1-obstacle.y0), linewidth=1, edgecolor='black', facecolor='white')
            ax.add_patch(rect)  
        # Layout
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_xlim([0, self.Lx])
        ax.set_ylim([0, self.Ly])
        ax.set_aspect('equal')
        ax.set_title('Quiver plot')
        clb=fig.colorbar(contour_filled, ax=ax, shrink=0.6) # Colorbar
        clb.set_label('Speed', rotation=270, labelpad=15)

        plt.show()

    def plot_quiver_pressure(self):
        """Makes a quiver plot of the system."""
        ### Create meshgrids ###
        x = np.linspace(0, self.Lx, self.nx)
        y = np.linspace(0, self.Ly, self.ny)
        X, Y = np.meshgrid(x, y)

        fig, ax = plt.subplots() # Create figure

        contour_filled = ax.contourf(X, Y, self.p.T, levels=50, cmap='rainbow_r') # Plot filled contours
        ax.quiver(X, Y, self.u.T, self.v.T) # Quiver plot
        for obstacle in self.obstacles:  # Draw obstacles
            rect = patches.Rectangle((obstacle.x0, obstacle.y0), (obstacle.x1-obstacle.x0), (obstacle.y1-obstacle.y0), linewidth=1, edgecolor='black', facecolor='white')
            ax.add_patch(rect)  
        # Layout
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_xlim([0, self.Lx])
        ax.set_ylim([0, self.Ly])
        ax.set_aspect('equal')
        ax.set_title('Quiver plot')
        clb=fig.colorbar(contour_filled, ax=ax, shrink=0.6) # Colorbar
        clb.set_label('Pressure', rotation=270)

        plt.show()

    def plot_streamline(self, linewidth=None):
        """Makes a streamline plot of the system."""
        ### Typechecks ###
        if not linewidth==None:
            typeCheck([linewidth], float)

        ### Set default value ###
        if linewidth==None:
            linewidth=2*np.sqrt(self.u.T**2 + self.v.T**2)/np.max(np.sqrt(self.u.T**2 + self.v.T**2))

        ### Create meshgrids ###
        x = np.linspace(0, self.Lx, self.nx)
        y = np.linspace(0, self.Ly, self.ny)
        X, Y = np.meshgrid(x, y)

        fig, ax = plt.subplots() # Create figure

        contour_filled = ax.contourf(X, Y, self.zeta.T, levels=50, cmap='rainbow_r') # Plot filled contours
        ax.streamplot(X, Y, self.u.T, self.v.T, color='black', linewidth=linewidth) # Plot streamlines
        for obstacle in self.obstacles: # Draw obstacles
            rect = patches.Rectangle((obstacle.x0, obstacle.y0), (obstacle.x1-obstacle.x0), (obstacle.y1-obstacle.y0), linewidth=1, edgecolor='black', facecolor='white')
            ax.add_patch(rect) 
        # Layout
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_xlim([0, self.Lx])
        ax.set_ylim([0, self.Ly])
        ax.set_aspect('equal')
        ax.set_title('Streamline plot')
        clb=fig.colorbar(contour_filled, ax=ax, shrink=0.6) # Colorbar
        clb.set_label('Vorticity', rotation=270, labelpad=15)

        plt.show()

    def plot_vorticity(self):
        """Makes a vorticity plot of the system."""
        ### Create meshgrids ###
        x = np.linspace(0, self.Lx, self.nx)
        y = np.linspace(0, self.Ly, self.ny)
        X, Y = np.meshgrid(x, y)

        fig, ax = plt.subplots() # Create figure

        plot = ax.contourf(X, Y, self.zeta.T, cmap='rainbow_r', alpha=0.5, levels=50) # Plot vorticity
        for obstacle in self.obstacles: # Draw obstacles   
            rect = patches.Rectangle((obstacle.x0, obstacle.y0), (obstacle.x1-obstacle.x0), (obstacle.y1-obstacle.y0), linewidth=1, edgecolor='black', facecolor='white')
            ax.add_patch(rect) 
        # Layout
        ax.set_title('Vorticity plot')
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_xlim([0, self.Lx])
        ax.set_ylim([0, self.Ly])
        ax.set_aspect('equal')
        clb=fig.colorbar(plot, ax=ax, shrink=0.6) # Colorbar
        clb.set_label('Vorticity', rotation=270, labelpad=15)

        plt.show()
    
    def plot_vorticity_symlog(self):
        """Makes a symmetrical logarithmic voriticty plot of the system."""
        ### Create meshgrids ###
        x = np.linspace(0, self.Lx, self.nx)
        y = np.linspace(0, self.Ly, self.ny)
        X, Y = np.meshgrid(x, y)

        fig, ax = plt.subplots() # Create figure

        logzeta=np.sign(self.zeta)*np.log(np.abs(self.zeta)+1) # Calculate symlog of voriticty

        
        plot = ax.contourf(X, Y, logzeta.T, cmap='rainbow_r', alpha=0.5, levels=50) # Plot voriticity
        for obstacle in self.obstacles:  # Draw obstacles
            rect = patches.Rectangle((obstacle.x0, obstacle.y0), (obstacle.x1-obstacle.x0), (obstacle.y1-obstacle.y0), linewidth=1, edgecolor='black', facecolor='white')
            ax.add_patch(rect) 
        # Layout
        ax.set_title('Logarithmic vorticity plot')
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_xlim([0, self.Lx])
        ax.set_ylim([0, self.Ly])
        ax.set_aspect('equal')
        clb=fig.colorbar(plot, ax=ax, shrink=0.6) # Colorbar
        clb.set_label('Symlog Vorticity', rotation=270, labelpad=15)

        plt.show()

    def plot_subplots(self):
        """ Creates a figure containing all general plots. """
        ### Create meshgrids ###
        x = np.linspace(0, self.Lx, self.nx)
        y = np.linspace(0, self.Ly, self.ny)
        X, Y = np.meshgrid(x, y)

        fig, ax = plt.subplots(2, 2, figsize=(7, 7), constrained_layout=False)  # Create figure

        ### Quiver plot speed ###
        contour_filled_quiver_speed = ax[0, 0].contourf(
            X, Y, np.sqrt(self.u**2 + self.v**2).T, levels=50, cmap='rainbow_r'
        )  # Plot filled contours
        ax[0, 0].quiver(X, Y, self.u.T, self.v.T)  # Quiver plot
        for obstacle in self.obstacles:  # Draw obstacles
            rect = patches.Rectangle(
                (obstacle.x0, obstacle.y0),
                (obstacle.x1 - obstacle.x0),
                (obstacle.y1 - obstacle.y0),
                linewidth=1,
                edgecolor='black',
                facecolor='white'
            )
            ax[0, 0].add_patch(rect)
        # Layout
        ax[0, 0].set_xlabel('X')
        ax[0, 0].set_ylabel('Y')
        ax[0, 0].set_xlim([0, self.Lx])
        ax[0, 0].set_ylim([0, self.Ly])
        ax[0, 0].set_aspect('equal')
        ax[0, 0].set_title('Quiver plot')
        clb = fig.colorbar(contour_filled_quiver_speed, ax=ax[0, 0], shrink=0.6)
        clb.set_label('Speed', rotation=270, labelpad=10)

        ### Quiver plot pressure ###
        contour_filled_quiver_pressure = ax[1, 0].contourf(
            X, Y, self.p.T, levels=50, cmap='bwr'
        )  # Plot filled contours
        ax[1, 0].quiver(X, Y, self.u.T, self.v.T)  # Quiver plot
        for obstacle in self.obstacles:  # Draw obstacles
            rect = patches.Rectangle(
                (obstacle.x0, obstacle.y0),
                (obstacle.x1 - obstacle.x0),
                (obstacle.y1 - obstacle.y0),
                linewidth=1,
                edgecolor='black',
                facecolor='white'
            )
            ax[1, 0].add_patch(rect)
        # Layout
        ax[1, 0].set_xlabel('X')
        ax[1, 0].set_ylabel('Y')
        ax[1, 0].set_xlim([0, self.Lx])
        ax[1, 0].set_ylim([0, self.Ly])
        ax[1, 0].set_aspect('equal')
        ax[1, 0].set_title('Quiver plot')
        clb = fig.colorbar(contour_filled_quiver_pressure, ax=ax[1, 0], shrink=0.6)
        clb.set_label('Pressure', rotation=270, labelpad=10)

        ### Streamline plot vorticity ###
        contour_filled_streamline_vorticity = ax[0, 1].contourf(
            X, Y, self.zeta.T, levels=50, cmap='rainbow_r'
        )  # Plot filled contours
        ax[0, 1].streamplot(X, Y, self.u.T, self.v.T, color='black', linewidth=1)
        for obstacle in self.obstacles:  # Draw obstacles
            rect = patches.Rectangle(
                (obstacle.x0, obstacle.y0),
                (obstacle.x1 - obstacle.x0),
                (obstacle.y1 - obstacle.y0),
                linewidth=1,
                edgecolor='black',
                facecolor='white'
            )
            ax[0, 1].add_patch(rect)  # Note: follows original code
        # Layout
        ax[0, 1].set_xlabel('X')
        ax[0, 1].set_ylabel('Y')
        ax[0, 1].set_xlim([0, self.Lx])
        ax[0, 1].set_ylim([0, self.Ly])
        ax[0, 1].set_aspect('equal')
        ax[0, 1].set_title('Streamline plot')
        clb = fig.colorbar(contour_filled_streamline_vorticity, ax=ax[0, 1], shrink=0.6)
        clb.set_label('Vorticity', rotation=270, labelpad=10)

        ### Symlog vorticity plot ###
        logzeta = np.sign(self.zeta) * np.log(np.abs(self.zeta) + 1)
        symlog_vorticity = ax[1, 1].contourf(
            X, Y, logzeta.T, levels=50, cmap='rainbow_r', alpha=0.5
        )  # Plot vorticity
        for obstacle in self.obstacles:  # Draw obstacles
            rect = patches.Rectangle(
                (obstacle.x0, obstacle.y0),
                (obstacle.x1 - obstacle.x0),
                (obstacle.y1 - obstacle.y0),
                linewidth=1,
                edgecolor='black',
                facecolor='white'
            )
            ax[1, 1].add_patch(rect)
        # Layout
        ax[1, 1].set_title('Logarithmic vorticity plot')
        ax[1, 1].set_xlabel('X')
        ax[1, 1].set_ylabel('Y')
        ax[1, 1].set_xlim([0, self.Lx])
        ax[1, 1].set_ylim([0, self.Ly])
        ax[1, 1].set_aspect('equal')
        clb = fig.colorbar(symlog_vorticity, ax=ax[1, 1], shrink=0.6)
        clb.set_label('Symlog Vorticity', rotation=270, labelpad=10)

        fig.tight_layout()
        plt.show()

    def plot_u_profile(self, x=None):
        """Plots a horizontal speed profile of the system."""
        ### Set default value ###
        if x==None:
            x=self.Lx/2

        ### Typechecks ###
        typeCheck([x], float)
        boundsCheck([x], 0, self.Lx)

        ### Extract data ###
        x_coor=int(np.floor(x/self.Lx*self.nx))
        u_profile=self.u[x_coor, :]

        fig, ax = plt.subplots() # Create figure

        ax.plot(u_profile, self.Ly*np.arange(0, self.ny)/self.ny, color='blue', linestyle='dashed', linewidth=1, label="Horizontal speed") # Plot u-profile
        # Layout
        ax.set_xticks(ticks=np.arange(np.min(u_profile), np.max(u_profile), (np.max(u_profile)-np.min(u_profile))/10))
        ax.set_yticks(ticks=np.arange(0, self.Ly+0.1, self.Ly/10))
        ax.tick_params('both', direction='in')
        ax.set_title(f'Horizontal speed profile at x={x}', fontsize=14)
        ax.set_xlabel('$V_{x}(y)$', fontsize=14)
        ax.set_ylabel('$y$', fontsize=14)
        ax.set_box_aspect(1)
        ax.grid()
        ax.legend()

        plt.show()

    def plot_v_profile(self, y=None):
        """Plots a vertical speed profile of the system."""
         ### Set default value ###
        if y==None:
            y=self.Ly/2

        ### Typechecks ###
        typeCheck([y], float)
        boundsCheck([y], 0, self.Ly)

        ### Extract data ###
        y_coor=int(np.floor(y/self.Ly*self.ny))
        v_profile=self.v[:, y_coor]

        fig, ax = plt.subplots()  # Create figure

        ax.plot(self.Lx*np.arange(0, self.nx)/self.nx, v_profile, color='red', linestyle='dashed', linewidth=1, label="Vertical speed") # Plot v-profile
        #Layout
        ax.set_xticks(ticks=np.arange(0, self.Lx+0.1, self.Lx/10))
        ax.set_yticks(ticks=np.arange(np.min(v_profile), np.max(v_profile), (np.max(v_profile)-np.min(v_profile))/10))
        ax.tick_params('both', direction='in')
        ax.set_title(f'Vertical speed profile at y={y}', fontsize=14)
        ax.set_xlabel('$x$', fontsize=14)
        ax.set_ylabel('$V_{y}(x)$', fontsize=14)
        ax.set_box_aspect(1)
        ax.grid()
        ax.legend()
        plt.show()    

    def plot_profiles_adami(self, x=None, y=None):
        """Plots horizontal and vertical speed profiles of the system, in the style of the S. Adami paper (see BlackBoard)."""
        if x==None:
            x=self.Lx/2
        if y==None:
            y=self.Ly/2

        ### Typechecks ###
        typeCheck([x, y], float)
        boundsCheck([x], 0, self.Lx)
        boundsCheck([y], 0, self.Ly)

        ### Extract data ###
        x_coor_u = int(np.floor(x / self.Lx * self.nx))
        y_coor_u = self.Ly * np.arange(0, self.ny) / self.ny
        u_profile = self.u[x_coor_u, :]
        y_coor_v = int(np.floor(y / self.Ly * self.ny))
        x_coor_v = self.Lx * np.arange(0, self.nx) / self.nx
        v_profile = self.v[:, y_coor_v]
        
        fig, ax = plt.subplots() # Create figure
        
        ax.plot(u_profile, y_coor_u, color='blue', linewidth=1, linestyle='dashed', label='Horizontal speed') # Plot u-profile
        # Layout
        ax.set_xlabel('$V_x(y)$', fontsize=14)
        ax.set_ylabel('$y$', fontsize=14)
        ax.set_title(f'Profiles at x = {x} and y = {y}', fontsize=14)
        ax.set_yticks(np.arange(0, self.Ly + 0.1, self.Ly / 10))
        ax.tick_params(direction='in')
        ax.grid()
        ax.set_box_aspect(1)
        
        ax_top = ax.twiny() # Copy axes
        twin_y = ax_top.twinx()

        twin_y.plot(x_coor_v, v_profile, color='red', linewidth=1, linestyle='dashed', label='Vertical speed') # Plot v-profile
        # Layout
        ax_top.set_xticks(np.arange(0, self.Lx + 0.1, self.Lx / 10))
        ax_top.set_xlabel('$x$', fontsize=14)
        ax_top.tick_params(direction='in')
        ax_top.xaxis.set_label_position('top')
        ax_top.xaxis.tick_top()    
        twin_y.set_ylabel('$V_y(x)$', fontsize=14)
        twin_y.tick_params(axis='both', direction='in')
        twin_y.set_box_aspect(1)
        
        plt.show()

    def plot_profiles(self, x=None, y=None):
        """Plots horizontal and vertical speed profiles of the system."""
        ### Set default value ###
        if x==None:
            x=self.Lx/2
        if y==None:
            y=self.Ly/2

        ### Typechecks ###
        typeCheck([x, y], float)
        boundsCheck([x], 0, self.Lx)
        boundsCheck([y], 0, self.Ly)

        ### Exctract data ###
        x_coor=int(np.floor(x/self.Lx*self.nx))
        u_profile=self.u[x_coor, :]
        y_coor=int(np.floor(y/self.Ly*self.ny))
        v_profile=self.v[:, y_coor]

        fig, ax = plt.subplots(1, 2) # Create figure

        ax[0].plot(u_profile, self.Ly*np.arange(0, self.ny)/self.ny, color='blue', linewidth=1, linestyle='dashed', label='Horizontal speed') # Plot u-profile
        #Layout
        ax[0].set_title(f'Horizontal speed profile at x={x}', fontsize=14)
        ax[0].set_xlabel('$V_{x}(y)$', fontsize=14)
        ax[0].set_ylabel('$y$', fontsize=14)
        ax[0].set_box_aspect(1)
        ax[0].grid()
        ax[0].legend(fontsize=14)

        ax[1].plot(self.Lx*np.arange(0, self.nx)/self.nx, v_profile,  color='red', linewidth=1, linestyle='dashed', label='Vertical speed') # Plot v-profile
        #Layout
        ax[1].set_title(f'Vertical speed profile at y={y}', fontsize=14)
        ax[1].set_xlabel('$x$', fontsize=14)
        ax[1].set_ylabel('$V_{y}(x)$', fontsize=14)
        ax[1].set_box_aspect(1)
        ax[1].grid()
        ax[1].legend(fontsize=14)

        plt.show()

### Hidden Methods ###

    def _update_p(self):
        """Applies the relaxation scheme to the pressure field."""
        self._p[1:self.nx-1, 1:self.ny-1]=1/(2*(self.dx**2+self.dy**2))*((self.p[0:self.nx-2, 1:self.ny-1]+self.p[2:self.nx, 1:self.ny-1])*self.dy**2+(self.p[1:self.nx-1, 0:self.ny-2]+self.p[1:self.nx-1, 2:self.ny])*self.dx**2) \
                -(self.rho*self.dx**2*self.dy**2)/(2*(self.dx**2+self.dy**2))* \
                (1/self.dt*((self.u[2:self.nx, 1:self.ny-1]-self.u[0:self.nx-2, 1:self.ny-1])/(2*self.dx)+(self.v[1:self.nx-1, 2:self.ny]-self.v[1:self.nx-1, 0:self.ny-2])/(2*self.dy))-((self.u[2:self.nx, 1:self.ny-1]-self.u[0:self.nx-2, 1:self.ny-1])/(2*self.dx))**2-2*((self.u[1:self.nx-1, 2:self.ny]-self.u[1:self.nx-1, 0:self.ny-2])/(2*self.dy))*((self.v[2:self.nx, 1:self.ny-1]-self.v[0:self.nx-2, 1:self.ny-1])/(2*self.dx))-((self.v[1:self.nx-1, 2:self.ny]-self.v[1:self.nx-1, 0:self.ny-2])/(2*self.dy))**2)

    def _update_u(self):
        """Applies the relaxation scheme to the horizontal speed field."""
        self._u[1:self.nx-1, 1:self.ny-1]=self.u[1:self.nx-1, 1:self.ny-1]-self.dt/self.dx*np.multiply(self.u[1:self.nx-1, 1:self.ny-1], self.u[1:self.nx-1, 1:self.ny-1]-self.u[0:self.nx-2, 1:self.ny-1]) \
                -self.dt/self.dy*np.multiply(self.v[1:self.nx-1, 1:self.ny-1], self.u[1:self.nx-1, 1:self.ny-1]-self.u[1:self.nx-1, 0:self.ny-2]) \
                -self.dt/(2*self.rho*self.dx)*(self.p[2:self.nx, 1:self.ny-1]-self.p[0:self.nx-2, 1:self.ny-1]) \
                +self.mu*self.dt*(1/self.dx**2*(self.u[2:self.nx, 1:self.ny-1]+self.u[0:self.nx-2, 1:self.ny-1]-2*self.u[1:self.nx-1, 1:self.ny-1])+1/self.dy**2*(self.u[1:self.nx-1, 0:self.ny-2]+self.u[1:self.nx-1, 2:self.ny]-2*self.u[1:self.nx-1, 1:self.ny-1])) \
                +self.dt*self.Fx
        
    def _update_v(self):
        """Applies the relaxation scheme to the vertical speed field."""
        self._v[1:self.nx-1, 1:self.ny-1]=self.v[1:self.nx-1, 1:self.ny-1]-self.dt/self.dx*np.multiply(self.u[1:self.nx-1, 1:self.ny-1], self.v[1:self.nx-1, 1:self.ny-1]-self.v[0:self.nx-2, 1:self.ny-1]) \
                -self.dt/self.dy*np.multiply(self.v[1:self.nx-1, 1:self.ny-1], self.v[1:self.nx-1, 1:self.ny-1]-self.v[1:self.nx-1, 0:self.ny-2]) \
                -self.dt/(2*self.rho*self.dy)*(self.p[1:self.nx-1, 2:self.ny]-self.p[1:self.nx-1, 0:self.ny-2]) \
                +self.mu*self.dt*(1/self.dx**2*(self.v[2:self.nx, 1:self.ny-1]+self.v[0:self.nx-2, 1:self.ny-1]-2*self.v[1:self.nx-1, 1:self.ny-1])+1/self.dy**2*(self.v[1:self.nx-1, 0:self.ny-2]+self.v[1:self.nx-1, 2:self.ny]-2*self.v[1:self.nx-1, 1:self.ny-1])) \
                +self.dt*self.Fy
    
    def _calculate_vorticity(self):
        """Calculates the vorticity of the system."""
        ### Internal ###
        self._zeta=np.zeros([self.nx, self.ny])
        self._zeta[1:self.nx-1, 1:self.ny-1]=-(self.u[1:self.nx-1, 2:self.ny]-self.u[1:self.nx-1, 0:self.ny-2])/(2*self.dy)+(self.v[2:self.nx, 1:self.ny-1]-self.v[0:self.nx-2, 1:self.ny-1])/(2*self.dx)

        ### Edges ###
        self._zeta[0, 1:self.ny-1]=-(self.u[0, 2:self.ny]-self.u[0, 0:self.ny-2])/(2*self.dy)+(self.v[1, 1:self.ny-1]-self.v[0, 1:self.ny-1])/(self.dx)
        self._zeta[self.nx-1, 1:self.ny-1]=-(self.u[self.nx-1, 2:self.ny]-self.u[self.nx-1, 0:self.ny-2])/(2*self.dy)+(self.v[self.nx-1, 1:self.ny-1]-self.v[self.nx-2, 1:self.ny-1])/(self.dx)
        self._zeta[1:self.nx-1, 0]=-(self.u[1:self.nx-1, 1]-self.u[1:self.nx-1, 0])/(self.dy)+(self.v[2:self.nx, 0]-self.v[0:self.nx-2, 0])/(2*self.dx)
        self._zeta[1:self.nx-1, self.ny-1]=-(self.u[1:self.nx-1, self.ny-1]-self.u[1:self.nx-1, self.ny-2])/(self.dy)+(self.v[2:self.nx, 0]-self.v[0:self.nx-2, 0])/(2*self.dx)

        ### Corners ###                                                                            
        self._zeta[0, 0]=-(self.u[0, 1]-self.u[0, 0])/(self.dx)+(self.v[1, 0]-self.v[0, 0])/(self.dy)
        self._zeta[0, self.ny-1]=-(self.u[0, self.ny-1]-self.u[0, self.ny-2])/(self.dx)+(self.v[1, self.ny-1]-self.v[0, self.ny-1])/(self.dy)
        self._zeta[self.nx-1, 0]=-(self.u[self.nx-1, 1]-self.u[self.nx-1, 0])/(self.dx)+(self.v[self.nx-1, 0]-self.v[self.nx-2, 0])/(self.dy)
        self._zeta[self.nx-1, self.ny-1]=-(self.u[self.nx-1, self.ny-1]-self.u[self.nx-1, self.ny-2])/(self.dx)+(self.v[self.nx-1, self.ny-1]-self.v[self.nx-2, self.ny-1])/(self.dy)

    def _apply_boundary_conditions(self, matrix, BVC):
        """Applies boundary conditions to a field."""
        # Up
        if BVC["up"]=="Neumann":
            matrix[1:self.nx-1, self.ny-1]=matrix[1:self.nx-1, self.ny-2]
        elif BVC["up"]=="Periodic":
            matrix[1:self.nx-1, self.ny-1]=(matrix[1:self.nx-1, self.ny-2]+matrix[1:self.nx-1, 1])/2
        else:
            matrix[1:self.nx-1, self.ny-1]=BVC["up"]

        # Down
        if BVC["down"]=="Neumann":
            matrix[1:self.nx-1, 0]=matrix[1:self.nx-1, 1]
        elif BVC["down"]=="Periodic":
            matrix[1:self.nx-1, 0]=(matrix[1:self.nx-1, self.ny-2]+matrix[1:self.nx-1, 1])/2
        else:
            matrix[1:self.nx-1, 0]=BVC["down"]

        # Left
        if BVC["left"]=="Neumann":
            matrix[0, 1:self.ny-1]=matrix[1, 1:self.ny-1]
        elif BVC["left"]=="Periodic":
            matrix[0, 1:self.ny-1]=(matrix[1, 1:self.ny-1]+matrix[self.nx-2, 1:self.ny-1])/2
        else:
            matrix[0, 1:self.ny-1]=BVC["left"]

        # Right
        if BVC["right"]=="Neumann":
            matrix[self.nx-1, 1:self.ny-1]=matrix[self.nx-2, 1:self.ny-1]
        elif BVC["left"]=="Periodic":
            matrix[self.nx-1, 1:self.ny-1]=(matrix[1, 1:self.ny-1]+matrix[self.nx-2, 1:self.ny-1])/2
        else:
            matrix[self.nx-1, 1:self.ny-1]=BVC["right"]

        return matrix

    def _apply_all_boundary_conditions(self):
        """Applies boundary conditions to all three fields."""
        # Pressure
        self._p=self._apply_boundary_conditions(self.p, self.p_BVC)
        # Horizontal
        self._u=self._apply_boundary_conditions(self.u, self.u_BVC)
        # Pressure
        self._v=self._apply_boundary_conditions(self.v, self.v_BVC)

    def _apply_obstacle_boundary_conditions(self, matrix, obstacle, BVC):
        """Applies boundary conditions for obstacles to a field."""
        # Up
        if not obstacle.touching_bounds["up"]:
            if BVC["up"] == "Neumann":
                matrix[obstacle.x_coors, obstacle.y1_coor] = matrix[obstacle.x_coors, obstacle.y1_coor+1]
            else:
                matrix[obstacle.x_coors, obstacle.y1_coor] = BVC["up"]

        # Down
        if not obstacle.touching_bounds["down"]:
            if BVC["down"] == "Neumann":
                matrix[obstacle.x_coors, obstacle.y0_coor] = matrix[obstacle.x_coors, obstacle.y0_coor-1]
            else:
                matrix[obstacle.x_coors, obstacle.y0_coor] = BVC["down"]

        # Left
        if not obstacle.touching_bounds["left"]:
            if BVC["left"] == "Neumann":
                matrix[obstacle.x0_coor, obstacle.y_coors] = matrix[obstacle.x0_coor-1, obstacle.y_coors]
            else:
                matrix[obstacle.x0_coor, obstacle.y_coors] = BVC["left"]

        # Right
        if not obstacle.touching_bounds["right"]:
            if BVC["right"] == "Neumann":
                matrix[obstacle.x1_coor, obstacle.y_coors] = matrix[obstacle.x1_coor+1, obstacle.y_coors]
            else:
                matrix[obstacle.x1_coor, obstacle.y_coors] = BVC["right"]

        #Interior
        matrix[obstacle.x0_coor+1: obstacle.x1_coor, obstacle.y0_coor+1: obstacle.y1_coor]=0

        return matrix
    
    def _apply_all_obstacle_boundary_conditions(self, obstacle):
        """Applies boundary conditions for obstacles to all fields."""
        # Pressure
        self._p=self._apply_obstacle_boundary_conditions(self.p, obstacle, obstacle.p_BVC)
        # Horizontal speed
        self._u=self._apply_obstacle_boundary_conditions(self.u, obstacle, obstacle.u_BVC)
        # Vertical speed
        self._v=self._apply_obstacle_boundary_conditions(self.v, obstacle, obstacle.v_BVC)

class NoSlipObstacle:
    """An obstacle with no-slip boundary conditions."""
    def __init__(self, fluid):
        """Initializes no-slip object"""
        ### Typechecks ###
        typeCheck([fluid], NewtonianFluid)

        self._Lx=fluid.Lx
        self._Ly=fluid.Ly
        self._dx=fluid.dx
        self._dy=fluid.dy
        self._nx=fluid.nx
        self._ny=fluid.ny

        self._p_BVC={'up':'Neumann', 'down':'Neumann', 'left':'Neumann', 'right':'Neumann'}
        self._u_BVC={'up':0, 'down':0, 'left':0, 'right':0}
        self._v_BVC={'up':0, 'down':0, 'left':0, 'right':0}

        self.fluid=fluid
        self.index=len(self.fluid._obstacles)

    ### Properties ###

    @property
    def Lx(self):
        """ Returns the x-dimention of the system. """
        return self._Lx  
        
    @property
    def Ly(self):
        """ Returns the y-dimention of the system. """
        return self._Ly
     
    @property
    def dx(self):
        """ Returns the x-spacing of the grid. """
        return self._dx
  
    @property
    def dy(self):
        """ Returns the y-spacing of the grid. """
        return self._dy

    @property
    def nx(self):
        """ Returns the vertical grid size of the system. """
        return self._nx

    @property
    def ny(self):
        """ Returns the vertical grid size of the system. """
        return self._ny

    @property
    def x0(self):
        """ Returns the left x-position of the obstacle. """
        return self._x0  
    
    @property
    def x1(self):
        """ Returns the right x-position of the obstacle. """
        return self._x1 

    @property
    def y0(self):
        """ Returns the bottom y-position of the obstacle. """
        return self._y0 
    
    @property
    def y1(self):
        """ Returns the top y-position of the obstacle. """
        return self._y1 

    @property
    def x0_coor(self):
        """ Returns the left x-coordinate of the obstacle. """
        return self._x0_coor  
    
    @property
    def x1_coor(self):
        """ Returns the right x-coordinate of the obstacle. """
        return self._x1_coor

    @property
    def y0_coor(self):
        """ Returns the bottom y-coordinate of the obstacle. """
        return self._y0_coor
    
    @property
    def y1_coor(self):
        """ Returns the top y-coordinate of the obstacle. """
        return self._y1_coor
    
    @property
    def x_coors(self):
        """ Returns the x-coordinates of the horizontal walls of the obstacle. """
        return self._x_coors
    
    @property
    def y_coors(self):
        """ Returns the y-coordinates of the vertical walls of the obstacle. """
        return self._y_coors  
    
    @property
    def x_dim(self):
        """ Returns the x-dimention of the system. """
        return self._x_dim
    
    @property
    def y_dim(self):
        """ Returns the y-dimention of the system. """
        return self._y_dim
    
    @property
    def touching_bounds(self):
        """ Returns dictionary containing booleans (if the walls of the obstacle touch the walls of the object). """
        return self._touching_bounds
    
    @property
    def p_BVC(self):
        """ Returns the pressure boundary conditions of the obstacle. """
        return self._p_BVC
    
    @property
    def u_BVC(self):
        """ Returns the horizontal speed boundary conditions of the obstacle. """
        return self._u_BVC

    @property
    def v_BVC(self):
        """ Returns the vertical speed boundary conditions of the obstacle. """
        return self._v_BVC

    ### Initialization methods ###

    def set_dimensions(self, x0, x1, y0, y1):
        """Sets dimensions of the obstacle."""
        ### Typechecks ###
        typeCheck([x0, x1, y0, y1], float)

        boundsCheck([x0], 0, x1)
        boundsCheck([x1], x0, self.Lx)
        boundsCheck([y0], 0, y1)
        boundsCheck([y1], y0, self.Ly)

        self._x0=x0
        self._x1=x1
        self._y0=y0
        self._y1=y1
        

        self._x0_coor=int(np.floor(x0/self.Lx*self.nx))
        self._x1_coor=int(np.floor(x1/self.Lx*self.nx))
        self._y0_coor=int(np.floor(y0/self.Ly*self.ny))
        self._y1_coor=int(np.floor(y1/self.Ly*self.ny))

        self._x_coors=np.arange(self.x0_coor+1, self.x1_coor)
        self._y_coors=np.arange(self.y0_coor+1, self.y1_coor)

        self._x_dim=np.size(self.x_coors)
        self._y_dim=np.size(self.y_coors)

        self._touching_bounds={'up':self.y1_coor==self.ny, 'down':self.y0_coor==0, 'left':self.x0_coor==0, 'right':self.x1_coor==self.nx}

    def set_boundary_conditions_p(self, BVC_up="Neumann", BVC_down="Neumann", BVC_left="Neumann", BVC_right="Neumann"):
        """Sets custom boundary conditions for the pressure."""
        ### Typechecks ###
        bvc_obstacleCheck([BVC_up, BVC_down, BVC_left, BVC_right])
        self._p_BVC={'up':BVC_up, 'down':BVC_down, 'left':BVC_left, 'right':BVC_right}

    def set_boundary_conditions_u(self, BVC_up=0, BVC_down=0, BVC_left=0, BVC_right=0):
        """Sets custom boundary conditions for the horizontal speed."""
        ### Typechecks ###
        bvc_obstacleCheck([BVC_up, BVC_down, BVC_left, BVC_right])
        self._u_BVC={'up':BVC_up, 'down':BVC_down, 'left':BVC_left, 'right':BVC_right}

    def set_boundary_conditions_v(self, BVC_up=0, BVC_down=0, BVC_left=0, BVC_right=0):
        """Sets custom boundary conditions for the vertical speed."""
        ### Typechecks ###
        bvc_obstacleCheck([BVC_up, BVC_down, BVC_left, BVC_right])
        self._v_BVC={'up':BVC_up, 'down':BVC_down, 'left':BVC_left, 'right':BVC_right}

    def delete(self):
        """Removes the obstacle from the system."""
        self.fluid.obstacles.pop(self.index)
        print(f"Deleted {self} from {self.fluid}")
        self.__del__()

    ### Hidden methods ###

    def _update_grid(self, Lx, Ly, dx, dy):
        """Updates object attributes in case grid size changes."""
        self._Lx=Lx
        self._Ly=Ly
        self._dx=dx
        self._dy=dy
        self._nx=int(np.floor(Lx/dx))
        self._ny=int(np.floor(Ly/dy))

        self._x0_coor=int(np.floor(self.x0/self.Lx*self.nx))
        self._x1_coor=int(np.floor(self.x1/self.Lx*self.nx))
        self._y0_coor=int(np.floor(self.y0/self.Ly*self.ny))
        self._y1_coor=int(np.floor(self.y1/self.Ly*self.ny))

        self._x_coors=np.arange(self.x0_coor+1, self.x1_coor)
        self._y_coors=np.arange(self.y0_coor+1, self.y1_coor)

        self._x_dim=np.size(self.x_coors)
        self._y_dim=np.size(self.y_coors)

    def __del__(self):
        self=None

def typeCheck(objects, check_type):
    """Typechecks a list of objects"""
    for object in objects:
        try: 
            if check_type in [float, list]:
                object_to_type=check_type(object)
            elif not isinstance(object, check_type):
                raise TypeError               
        except:
            raise TypeError(f'Error, the variable {object} has to be of type {check_type}!')

def bvcCheck(BVCs):
    """Checks the validity of boundary conditions."""
    for BVC in BVCs:
        if not np.any([BVC=="Neumann", BVC=="Periodic", isinstance(BVC, float), isinstance(BVC, int)]):
            raise ValueError('Error, boundary value conditions must be "Neumann", "Periodic", or a numeric value!')
    if (BVCs[0]=="Periodic" and not BVCs[1]=="Periodic") or (BVCs[1]=="Periodic" and not BVCs[0]=="Periodic") or (BVCs[2]=="Periodic" and not BVCs[3]=="Periodic") or (BVCs[3]=="Periodic" and not BVCs[2]=="Periodic"):
        raise ValueError('Error, opposite boundary conditions must both, or neither be periodic!')
             

def bvc_obstacleCheck(BVCs):
    """Checks the validity of boundary conditions for obstacles."""
    for BVC in BVCs:
        if not np.any([BVC=="Neumann", isinstance(BVC, float), isinstance(BVC, int)]):
            raise ValueError('Error, boundary value conditions for obstacles must be "Neumann", or a numeric value!')  

def posCheck(numbers):
    """Checks if list of numbers is positive."""
    for number in numbers:
        if number <=0:
            raise ValueError(f'Error, the variable {number} has to be greater than 0!')
        
def boundsCheck(numbers, lb, ub):
    """Checks if list of numbers is in a given interval."""
    for number in numbers:
        if np.any([number<lb, number>ub]):

            raise ValueError(f'Error, the variable {number} has to be between {lb} and {ub}!')
