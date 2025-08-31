
The code integrates the motion of a system of N planets, provided as input, for which the mass, initial position, and velocity are known.

In particular, the following objects can be initialized:

1. Point(list):
   A point with 3 coordinates. `list` must contain 2 or 3 numbers. In case of two numbers, the third dimension is set to 0.

2. Planet( position: list | Point,  velocity: list | Point, mass: float =0., name: str ='', color: str =''):
   `position` and `velocity` can be either `Point` objects or lists of 2 or 3 dimensional coordinates.
   `color` is a parameter in case you want to plot the planet’s motion using a `Plot` object.
   `T_energy()`: method to calculate the kinetic energy.

3. Integrator(planetlist:list | str, Ndim: int =None, terminal_print: bool =False):
   ATTRIBUTES:
   - `planetlist`: can be either
     - a list of `Planet` objects
     - a file with 9 columns (3×position, 3×velocity, 1×mass, planet name, color for plotting)
       Example: `0.1 0.1 0.1 0 0 0 100 planet1 blue`
   - `Ndim`: specifies the dimension of the `Point` objects, can be 2 or 3.
     If `Ndim=None` and all planets have a zero coordinate in positions and velocities, the default `Ndim=2` is used.
   - `terminal_print`: control flag to print positions at each integration step.

   METHODS:
   3.1 `reboot()`: restores the initial conditions.
   
   3.2 `printpos()`: prints the position vectors of the list of planets.
       `print_list(method='')`: prints positions, velocities, and names of all planets in a list; `method` is a string printed before the planets.

   3.3 `leapfrog(self, dt: numbers.Real, N: int)`: integration using leapfrog with step `dt` and `N` steps.
   
   3.4 `verlet(self, dt: numbers.Real, N: int)`: integration using Verlet with step `dt` and `N` steps.
       `step_verlet(dt, planetlist=None)`: single integration step with Verlet.
       By default, `planetlist` is the list used to initialize the Integrator, but a different list can be provided.

   3.5 `RK4(self, dt: numbers.Real, N: int)`: integration using RK4 with step `dt` and `N` steps.
       `step_RK4(dt, planetlist=None)`: single integration step with RK4.
       By default, `planetlist` is the list used to initialize the Integrator, but a different list can be provided.

5. Plot(planetlist:list | str, Ndim: int=None, figsize: tuple[int, int]=(10,10), dtperpoint=None, blackbackground:bool =True):
   Integrator object that produces a 2D/3D plot at the end of the integration.

   ATTRIBUTES:
   - `planetlist`: same as in Integrator.
   - `Ndim`: determines whether the plots are 2D or 3D.
   - `figsize`: sets the figure size.
   - `dtperpoint`: number of integration steps executed before plotting a point.
     Useful if there are many steps to avoid slowing down the plot.
     Default: if the number of steps < 50000, `dtperpoint=1`; otherwise `dtperpoint=Number of steps / 50000`.
   - `blackbackground`: flag to control the figure background color.

   METHODS:
   
   4.1 All Integrator methods, except `leapfrog`, `verlet`, and `RK4`, which are redefined.

   4.2 `show()`: method to plot at the end of integration.

   4.3 `leapfrog(self,  dt: numbers.Real, N: int, position=111, planetlist: list=None)`:
       `dt`, `N`, and `planetlist` as in Integrator;
       `position` specifies the subplot or figure position.
       Similarly for `verlet(dt, N, position=111, planetlist=None)` and `RK4(dt, N, position=111, planetlist=None)`.

   4.4 `all_methods(self, dt: numbers.Real, N: int)`: integrates and plots the system using all three methods; the `planetlist` attribute of Graph is updated at step N with the RK4 method.

   4.5 `reboot()`: restores initial conditions.

   4.6 `clear_figure()`: resets/clears the figure.


7. Interactive_plot(planetlist:list | str, Ndim:int =None, figsize: tuple[int, int]=(10,10), delay: numbers.Real =10, fixed_axes: bool=False, Nmaxintegration: int =10000000, blackbackground: bool=True):
   Integrator object that produces a 2D/3D plot evolving over time.

   ATTRIBUTES:
   - `planetlist`, `Ndim`, `figsize`, `blackbackground`: same as in Graph.
   - `delay`: time in milliseconds between frames (updating all planet positions).
     It is adjusted if the computational integration time is high.
   - `fixed_axes`: flag to keep axes fixed.
   - `Nmaxintegration`: maximum number of plot updates; stops integration when reached.

   METHODS:
   
   5.1 All methods from Integrator and Graph, except `leapfrog`, `verlet`, and `RK4`, which are redefined.

   5.2 `leapfrog(self, dt: numbers.Real, N: int, position=111)`:
       `dt`: integration step,
       `N`: number of steps executed BEFORE plotting a point; delay may be adjusted based on computation time (maximum between delay and the 1st integration time),
       `position` and `planetlist` same as in Graph.
       Similarly for `verlet(dt, N, position=111, planetlist=None)` and `RK4(dt, N, position=111, planetlist=None)`.
       `all_methods(self, dt: numbers.Real, N: int)`: integrates and plots the system with all three methods; `planetlist` is updated using RK4.
   5.3 `reboot()`: restores initial conditions and recreates the previously closed interactive window if needed.
