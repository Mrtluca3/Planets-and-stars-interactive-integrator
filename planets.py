import warnings
import code #for interactive mode
import matplotlib.pyplot as plt
from copy import deepcopy
from math import *
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg as FCTA
import tkinter as tk
import random
from time import perf_counter
import argparse
import numbers


G=6.67430* 10**(-11) #SI units: m^3 kg^-1 s^-2

#Classes
class Point():
    def __init__(self, position: list): 
        if isinstance(position, list):
          if all(isinstance(pos, numbers.Real) for pos in position): #checking the numbers in the list
           if len(position)==2: 
            self.x=position[0]  
            self.y=position[1]
            #default z=0 if list has just 2 elements
            self.z=0
           elif len(position)==3:
            self.x=position[0]  
            self.y=position[1]         
            self.z=position[2]
           else:
            raise ValueError("Definition not allowed: too many arguments, must be 2 or 3")              
          else:
           raise ValueError("Definition not allowed: must be float or int")          
        else:
          raise ValueError("Definition not allowed: position must be a list")
           
    def __getitem__(self, index):
        if index==0:
          return self.x
        elif index==1:
          return self.y
        elif index==2:
          return self.z
        else:
          raise ValueError("Index not allowed")
    def __setitem__(self, index, value):
        if index == 0:
            self.x = value
        elif index == 1:
            self.y = value
        elif index == 2:
            self.z = value
        else:
            raise ValueError("Index not allowed")
    #OPERATIONS
    def mod(self): #method: return the module of a Point object
        return (self.x*self.x+self.y*self.y+self.z*self.z)**(0.5)
    def __sub__(self, point2):
        #checking the object type
        if point2.__class__==self.__class__:
          return Point([self.x-point2.x, self.y-point2.y, self.z-point2.z]) 
        else:
          raise ValueError("Subtraction not allowed")
    def __add__(self, point2):
        #checking the object type
        if point2.__class__==self.__class__:
          return Point([self.x+point2.x, self.y+point2.y, self.z+point2.z]) 
        else:
          raise ValueError("Addition not allowed")
    def __iadd__(self, point2):
        #checking the object type
        if point2.__class__==self.__class__:
          return Point([self.x+point2.x, self.y+point2.y, self.z+point2.z]) 
        else:
          raise ValueError("Addition not allowed")
    def __mul__(self, number): 
        if isinstance(number, (int, float)):
          return Point([self.x*number, self.y*number, self.z*number]) 
        else:
          raise ValueError("Multiplication not allowed")
    def __truediv__(self, number): 
        if isinstance(number, (int, float)):
          inv_number=1./number
          return self*inv_number
        else:
          raise ValueError("Division not allowed")
    def __eq__(self, other):
      #Method: check if the objects types are the same
      if self.__class__==other.__class__:
        if all([self[i]==other[i] for i in range(3)]):
          return True
        else:
          return False
      else: raise ValueError("Comparison not allowed")
    def __str__(self): #print coordinates instead of memory address
        return str(self.x)+" "+str(self.y)+" "+str(self.z)

class Material_point(): 
    def __init__(self, position: list | Point, velocity: list | Point, mass: float): 
        if isinstance(position, list):  
           self.pos=Point(position)
        elif isinstance(position, Point):
           self.pos=position
        else:
           raise ValueError("Definition not allowed: position must be a list or Point object")
      
        if isinstance(velocity, list):
           self.vel=Point(velocity)
        elif isinstance(velocity, Point):
           self.vel=velocity
        else:
           raise ValueError("Definition not allowed: velocity must be a list or Point object")
        self.mass=mass
    def T_energy(self): #Method: calculate kinetic energy
        return 0.5*self.mass * self.vel.mod()*self.vel.mod()

class Planet(Material_point):
    def __init__(self, position: list | Point,  velocity: list | Point, mass: float =0., name: str ='', color: str =''):
        #inheritance
        super().__init__(position, velocity, mass)
        self.name=name
        if color!='':
          self.color=color
        else:
          self.color='blue'
        self.size=8 #marker size for plotting
    def __str__(self): 
        #Method: calling 'print(Planet)', return Planet name
        return self.name
    def print(self):
        print("Name", self.name, " mass", self.mass, "\nPosition: ", self.pos, "\nVelocity: ", self.vel)
    def __eq__(self, other):
        # Method: check if two Planets have the same position or velocity
        if other.__class__==self.__class__:
          if self.pos==other.pos and self.vel==other.vel:
             return True
          else:
             return False          
        else:
          raise ValueError("Comparison not allowed")
    
#Classes to calculate acceleration
class Force():
    def phi(planet1: Planet, planet2: Planet, Ndim: int):
        #acceleration due to a single gravitational field 
        pos1=planet1.pos 
        pos2=planet2.pos
        mass2=planet2.mass
        #pos1: position of the planet that receives acceleration
        #pos2: position of the planet that gives acceleration
        phivect=Point([0, 0, 0]) 
        for d in range(Ndim):
          phivect[d]= G*mass2*(pos2[d]-pos1[d])/((pos2-pos1).mod()**(3)) 
        return phivect
    def total_phi(planet1: Planet, planetlist: list, Ndim: int):  
        #Method: sum all the accelerations due to each planet in the list
        #Note: there is a checking that planet1 is not included
        phivect=Point([0, 0, 0])
        for planet in planetlist:
          if not planet1==planet:
              phivect+=Force.phi(planet1, planet, Ndim)
        return phivect

#Integrator Class
class Integrator():
    def __init__(self, planetlist:list | str, Ndim: int =None, terminal_print: bool =False):
      if isinstance(planetlist, list):
        self.planetlist=planetlist
      elif isinstance(planetlist, str): #file reading
        try: 
         file=open(planetlist, 'r')     
        except:
          raise ImportError(planetlist+" file doesn't exist")
        l=[]
        for line in file:
          line=line.strip()
          if not line.startswith("#"):
            data = line.split()
            l.append(Planet([float(num) for num in data[0:3]], #position
                           [float(num) for num in data[3:6]], #velocity
                           mass=float(data[6]),
                           name=data[7],
                           color=data[8]
                           )
                   )
        self.planetlist=l
        print('Initial conditions successfully read from file')
      else:
        raise ValueError('planetlist must be a list or a file.txt with 9 columns')
      
      if Ndim==2 or Ndim==3: #user can select the number of dimensions Ndim
        self.Ndim=Ndim
      else: #automatic selection Ndim=3 if at least one planet has a dimension !=0
        if Ndim!=2 and Ndim!=3 and Ndim!=None:
          print("Error, Ndim has to be 2 or 3 or None. The dimension will be automatically based on your system.")
        if any([p.pos[2]!=0 or p.vel[2]!=0 for p in self.planetlist]):
         self.Ndim=3
        else:
         self.Ndim=2 
      self.terminal_print=terminal_print #checking the terminal output
      self.copy_reboot=deepcopy(self.planetlist) #copy used in method 'reboot'
    def reboot(self): 
      #Method: reboot the initial conditions
      self.planetlist=deepcopy(self.copy_reboot)
      print("Robooting the initial conditions")

    def printpos(self): #print positions of a planets list
      for planet in self.planetlist:
        for d in range(self.Ndim):
         print(planet.pos[d], end=' ')
        print(end='  ')
      print('')
    def print_list(self, method: str='', planetlist: list=None): #Method: print all the planets of a list
      if planetlist==None:
       planetlist=self.planetlist
      print(method) #printing integration method
      for p in planetlist:
         p.print()
    ##########
    #Leapfrop#
    ##########

    def step0leap(self, dt: numbers.Real, planetlist: list=None):
        if planetlist==None:
          planetlist=self.planetlist #they are pointers: modifying planetlist will also modify self.planetlist
        # Note: the option to pass 'planetlist' is provided to allow using this method with
        #       the 'Graph' class, in its 'all' method, where it can plot all three methods,
        #       each evolving independently. This is why the step should not be executed
        #       exclusively on the instance attribute 'self.planetlist'.

        #First step of leapfrog integration using Eulero's approximation
        #In input
        #planetlist: r0, v_(0)
        #step1:
        #v1/2=v0+ phi(x0)* dt/2
        #x1=x0+v(1/2)*dt
        #Output
        #planetlist: r0, v_(0)
        #planetsfrog: r1, v_(1/2)
        v12=[deepcopy(planetlist[i].vel + Force.total_phi(planetlist[i], planetlist, self.Ndim)*(dt*0.5)) 
             for i in range(len(planetlist))]

        planetsfrog=[Planet(planetlist[i].pos + v12[i]*dt, v12[i], mass=planetlist[i].mass)
                      for i in range(len(planetlist))]
        return planetsfrog
    
    def stepleap(self, dt: numbers.Real, planetlist=None, planetsfrog=None):
        if planetlist==None:
           planetlist=self.planetlist 
        if planetsfrog==None:   
           planetsfrog=self.planetsfrog
        #Input
        #planetlist: r(n-1), v(n-1)
        #planetsfrog: rn e v(n-1/2)

        #LEAP FROG STEP:        
        #x_(n+1)=x_n + v_(n+1/2) * dt
        #v_(n+1/2)=v_(n-1/2) + phi(x_n) * dt

        #output
        #planetlist: r(n), v(n)
        #planetsfrog: r(n+1) e v(n+1/2)
 
      
        vlist=[deepcopy(frog.vel) for frog in planetsfrog] # list with velocities (n-1/2)-th, used to calculate vn mean 
        #v_(n+1/2):
        v12=[deepcopy(planetsfrog[i].vel + 
                           Force.total_phi(planetsfrog[i], planetsfrog, self.Ndim) * dt) 
                           for i in range(len(planetlist))] 
        # i=0 planetlist ha poizioni r_1, v_1, i=1: r_2, v_2, ..., i=N-1: r_N, v_N 
        # Note: planetlist contains the n-th positions (x_n) and velocities (v_n), updated in each cycle, 
        #       taken from planetsfrog.
        # In fact, Leapfrog performs an initial step ('step0') that advances the planetsfrog list by one step.
        # In the loop over N steps, planetlist updates its coordinates exactly like in the other integrators:
        # i=0: planetlist has positions r_1, velocities v_1
        # i=1: r_2, v_2
        # ...
        # i=N-1: r_N, v_N

        for i in range(len(planetlist)):
          planetlist[i].pos=planetsfrog[i].pos 
          planetlist[i].vel=(v12[i]+vlist[i])*0.5
          # update planetsfrog according to the algorithm (x(n+1), v(n+1/2))
          planetsfrog[i]=Planet(planetsfrog[i].pos + v12[i]*dt, v12[i], mass=planetlist[i].mass) 

    def leapfrog(self, dt: numbers.Real, N: int):
      if self.terminal_print:
        # Terminal Output: integration method and planets names
        # For each step, positions will be updated
        print("leapfrog: ", end='')
        for planet in self.planetlist:
         print(planet, end=', ') #print of the names: column headers 
        print('')
      
      t_ex=perf_counter() #starts the time of excution
      self.planetsfrog=self.step0leap(dt)  
      for n in range(N):
        self.stepleap(dt)
        if self.terminal_print:
         self.printpos()
      if not self.terminal_print: 
        #Terminal output: just final position at time Tmax
        self.print_list(method='\nLeapfrog') 
      t_ex-=perf_counter() #time ends
      print("Execution time: ", -t_ex, "[s]")
    ########
    #verlet#
    ########
    def stepverlet(self, dt: numbers.Real, planetlist=None):
       if planetlist==None:
         planetlist=self.planetlist 
       #Input: xn e vn
       #Output: xn+1 e vn+1
       #xn+1= xn+ vn*dt+ phi(xn)* dt^2/2
       #vn+1=vn+ (phi(xn)+phi(xn+1)) * dt/2

       #saving n-th positions by copying n-th list
       planetlistcopy=deepcopy(planetlist)
       #update position xn+1
       for i in range(len(planetlist)):
          #calculating xn+1
          planetlist[i].pos=planetlistcopy[i].pos + planetlistcopy[i].vel*dt+ Force.total_phi(planetlistcopy[i], planetlistcopy, self.Ndim) *(dt*dt*0.5)
   
       for i in range(len(planetlist)):     
          #calculating vn+1
          planetlist[i].vel=planetlistcopy[i].vel + (Force.total_phi(planetlistcopy[i], planetlistcopy, self.Ndim) + Force.total_phi(planetlist[i], planetlist, self.Ndim)) *(dt*0.5)        
          # Note: it is strictly necessary to perform two loops, since phi(x_{n+1}) depends on all
          #       the positions updated at step n+1
   
    def verlet(self, dt: numbers.Real, N: int):
      if self.terminal_print:
        print("Verlet: ", end='')
        for planet in self.planetlist:
         print(planet, end=', ')  
        print('')
      t_ex=perf_counter()    
      for n in range(N):
        self.stepverlet(dt)
        if self.terminal_print:
          self.printpos()

      if not self.terminal_print: 
  
          self.print_list(method='\nVerlet') 
      t_ex-=perf_counter()
      print("Execution time: ", -t_ex, "[s]")
    #####
    #RK4#
    #####
    def F(planetlist, dt: numbers.Real, Ndim):
        #useful function for RK4
        #Output: list of Material Point and updating of a step dt 
        return [Material_point(planetlist[i].vel*dt, 
                                Force.total_phi(planetlist[i], planetlist, Ndim)*dt,  
                                mass=0.) 
                                for i in range(len(planetlist))]
        #Note: they are massless for the moment

    def stepRK4(self, dt: numbers.Real, planetlist=None):
        if planetlist==None:
           planetlist=self.planetlist 
        Ndim=self.Ndim
        #using k1 as a list of Material point
        #(with velocity * dt in position)
        #(with acceleration * dt in velocità)
        #mass is null
        k1=Integrator.F(planetlist, dt, Ndim)
        #
        liststep=[Planet(planet.pos+ k.pos*0.5, 
                           planet.vel+ k.vel*0.5, 
                           planet.mass)
                           for planet,k in zip(planetlist, k1)]
        k2=Integrator.F(liststep, dt, Ndim)
        #
        liststep=[Planet(planet.pos+ k.pos*0.5, 
                           planet.vel+ k.vel*0.5, 
                           planet.mass)
                           for planet,k in zip(planetlist, k2)]
        k3=Integrator.F(liststep, dt, Ndim)
        #
        liststep=[Planet(planet.pos+k.pos, 
                           planet.vel+k.vel, 
                           planet.mass)
                           for planet,k in zip(planetlist, k3)]
        k4=Integrator.F(liststep, dt, Ndim)
        for i in range(len(planetlist)):
         planetlist[i].pos=planetlist[i].pos+ (k1[i].pos+k2[i].pos*2+k3[i].pos*2+k4[i].pos)/6 
         planetlist[i].vel=planetlist[i].vel+ (k1[i].vel+k2[i].vel*2+k3[i].vel*2+k4[i].vel)/6 
   
    def RK4(self, dt: numbers.Real, N: int):
      if self.terminal_print:
        print("Rk4: ", end='')
        for planet in self.planetlist:
         print(planet, end=', ') 
        print('')
      t_ex=perf_counter() 
      for n in range(N):
         self.stepRK4(dt)
         if self.terminal_print:
          self.printpos()
      
      if not self.terminal_print: 
          self.print_list( method='\nRk4') 
      t_ex-=perf_counter() 
      print("Execution time: ", -t_ex, "[s]")





#Classes for plotting
class Plot(Integrator):
   def __init__(self, planetlist:list | str, Ndim: int=None, figsize: tuple[int, int]=(10,10), dtperpoint=None, blackbackground:bool =True):
    super().__init__(planetlist, Ndim=Ndim)
    self.fig=plt.figure(figsize=figsize)
    #Flag controlling the 1st plot
    self.step0=True
    self.figsize=figsize
    #Flag controlling the background color
    self.blackbackground=blackbackground
    if blackbackground:
     self.fig.patch.set_facecolor('black')
    self.dtperpoint=dtperpoint #number of integration steps before plotting a point in the Plot
    if self.Ndim==2:
      self.stringdim=None
    else:
      self.stringdim='3d'
   def clear_figure(self):
     self.fig.clear()
     self.__init__(self.planetlist, self.Ndim, self.figsize, self.dtperpoint, self.blackbackground)
     self.step0=True
   def show(self):
     self.figreboot=deepcopy(self.fig)
     self.fig.show()
     self.fig=self.figreboot  
      

   def initialpointplot(self, planetlist, ax, methodname):
     #Point Size of the planet scale with planets dimensions:
     maxmass=max([planet.mass for planet in planetlist])
     minmass=min([planet.mass for planet in planetlist])
     #calculating planet.size in the next cycle
     if self.Ndim==2:
      for planet in planetlist:
        if maxmass!=minmass:
           planet.size=(planet.mass-minmass)/(maxmass-minmass)*10+5 #markersize: between 5 to 15

        if (len(planetlist)<15): # label included only if they are less than 15 planets
          l=planet.name
        else:
          l=None
        ax.plot(planet.pos.x, planet.pos.y, marker='.',  label=l, color=planet.color, markersize=planet.size)
      ax.set_title('Planets with '+methodname)
      ax.set_xlabel("X")
      ax.set_ylabel("Y")
      if self.blackbackground:
       ax.tick_params(axis='x', colors='white')  # tick color ax X
       ax.tick_params(axis='y', colors='white')  # tick color ax Y
       ax.xaxis.label.set_color('white')  # label color ax X
       ax.yaxis.label.set_color('white')  # label color ax Y
       ax.set_facecolor('black')

     else:  #Ndim=3
      for planet in planetlist:
        if maxmass!=minmass:
           planet.size=(planet.mass-minmass)/(maxmass-minmass)*10+5 #markersize: between 5 to 15
        if (len(planetlist)<15): 
          l=planet.name
        else:
          l=None
        ax.plot(planet.pos.x, planet.pos.y, planet.pos.z, marker='.',  label=planet.name, color=planet.color, markersize=planet.size)
      ax.set_title('Planets with '+methodname)
      ax.set_xlabel("X")
      ax.set_ylabel("Y")
      ax.set_ylabel("Z")
      if self.blackbackground:
       ax.set_facecolor('black')
       ax.tick_params(axis='x', colors='white')  
       ax.tick_params(axis='y', colors='white')  
       ax.tick_params(axis='z', colors='white')  
       ax.xaxis.label.set_color('white')  
       ax.yaxis.label.set_color('white')  
       ax.zaxis.label.set_color('white')  
      
   def pointplot(self, planetlist, ax): #No printing of labels
     if self.Ndim==2:
      for planet in planetlist:
       ax.plot(planet.pos.x, planet.pos.y, marker='.', color=planet.color, markersize=planet.size)
     else:
      for planet in planetlist:
       ax.plot(planet.pos.x, planet.pos.y, planet.pos.z, marker='.', color=planet.color, markersize=planet.size)

   #Integration with plotting
   def leapfrog(self,  dt: numbers.Real, N: int, position=111, planetlist: list=None):  
     if planetlist==None:
        planetlist=self.planetlist 
     if self.dtperpoint==None:
        if N>50000: 
          #  Number of plotted points is limited, in order to avoid overloading the plot when N is too large
          self.dtperpoint=int(N/50000)
        else:
          self.dtperpoint=1
     if self.step0:
       self.ax=self.fig.add_subplot(position, projection=self.stringdim)
       #initial point plotting
       self.initialpointplot(planetlist, self.ax, 'leapfrog') 
       self.step0=False
     #timer starts
     t_ex=perf_counter()
     #step0 of leapfrog
     self.planetsfrog=self.step0leap(dt, planetlist=planetlist)
     for n in range(N):
        self.stepleap(dt, planetlist=planetlist)
        if (n%self.dtperpoint==0):
          self.pointplot(planetlist, self.ax) 
     t_ex-=perf_counter() 
     self.print_list('\nleap frog', planetlist=planetlist)
     print("Execution time: ", -t_ex, "[s]")
     self.ax.grid(), self.ax.legend(), plt.tight_layout()
 
   def verlet(self,dt: numbers.Real, N: int, position=111, planetlist: list =None):
     if planetlist==None:
        planetlist=self.planetlist
     if self.dtperpoint==None:
        if N>50000:
          self.dtperpoint=int(N/50000)
        else:
          self.dtperpoint=1
     if self.step0:
       self.ax=self.fig.add_subplot(position, projection=self.stringdim)
       self.initialpointplot(planetlist, self.ax, 'verlet')
       self.step0=False
     t_ex=perf_counter()
     for n in range(N):
       self.stepverlet(dt, planetlist=planetlist)
       if (n % self.dtperpoint==0):
        self.pointplot(planetlist, self.ax) 
     t_ex-=perf_counter()
     self.print_list(method='\nVerlet', planetlist=planetlist)
     print("Execution time: ", -t_ex, "[s]")
     self.ax.grid(), self.ax.legend(), plt.tight_layout()     
  
   def RK4(self, dt: numbers.Real, N: int, position=111, planetlist: list =None):
      if planetlist==None:
        planetlist=self.planetlist 
      if self.dtperpoint==None:
        if N>50000:
          self.dtperpoint=int(N/50000)
        else:
          self.dtperpoint=1
      if self.step0:
       self.ax=self.fig.add_subplot(position, projection=self.stringdim)  
       self.initialpointplot(planetlist, self.ax, 'RK4')
       self.step0=False
      t_ex=perf_counter()
      for n in range(N):
         self.stepRK4(dt, planetlist=planetlist)
         if (n % self.dtperpoint==0):
          self.pointplot(planetlist, self.ax) 
      t_ex-=perf_counter()
      self.print_list(method='\nRK4', planetlist=planetlist)
      print("Execution time: ", -t_ex, "[s]")
      self.ax.grid(), self.ax.legend(), plt.tight_layout()
   
   def all_methods(self, dt: numbers.Real, N: int):
      self.leapfrog(dt, N, position=221, planetlist=deepcopy(self.planetlist))  
      self.step0=True  
        
      self.verlet(dt, N, position=222, planetlist=deepcopy(self.planetlist))
      self.step0=True
      self.RK4(dt, N, position=223) 
      self.step0=True
      plt.tight_layout()

class Interactive_plot(Plot):
  def __init__(self,planetlist:list | str, Ndim:int =None, figsize: tuple[int, int]=(10,10), delay: numbers.Real =10, fixed_axes: bool=False, Nmaxintegration: int =10000000, blackbackground: bool=True): 
    super().__init__(planetlist, Ndim=Ndim, figsize=figsize, blackbackground=blackbackground)
    #Opening tk window
    self.root=tk.Tk()
    self.root.title("Planets plot")
    self.fig=plt.figure(figsize=(6,6))
    self.widget=FCTA(self.fig, master=self.root)
    self.widget.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

    #Button to save the plot
    self.buttonsave=tk.Button(master=self.root, command=self.save_pdf, text='Save in pdf')
    self.buttonsave.pack(side=tk.BOTTOM)

    # counter of number of updating to stop the integration if self.counter>self.Nmaxintegration
    self.counter=0
    self.Nmaxintegration=Nmaxintegration
    #Updating time of the Figure is given by 'delay', set at 0.1sec. if not specified
    #It is dynamically adjusted according to the computational time of integration
    self.delay=delay 
    #If fixed_axes==true, axes limits are fixed to the initial ones
    self.fixed_axes=fixed_axes
    #Flag to check window opening, needed for rebooting
    self.openroot=True
    #String Flag describing the id of a 'tk.after' programmed, needed for rebooting 
    self.id=None
  def save_pdf(self):
    self.fig.savefig('plot.pdf', format='pdf')
    print('Figure saved in pdf as "plot.pdf" ')
  def mainloop(self):
    self.root.mainloop()
  def reboot(self):
    super().reboot()
    if not self.id==None:
      self.root.after_cancel(self.id)
      self.root.quit()
      # reinitialize
      self.__init__(self.planetlist, Ndim=self.Ndim, figsize=self.figsize, delay=self.delay, fixed_axes=self.fixed_axes, Nmaxintegration=self.Nmaxintegration, blackbackground=self.blackbackground)
  def leapfrog(self, dt: numbers.Real, N: int, position=111):
    planetlist=self.planetlist
    self.counter+=1
    #step0
    if self.step0: 
     self.ax=self.fig.add_subplot(position, projection=self.stringdim)
     self.planetsfrog=self.step0leap(dt)
     self.initialpointplot(planetlist, self.ax, 'Leapfrog')
     if self.fixed_axes:
      self.ax.relim(), self.ax.autoscale_view(), self.ax.set_autoscale_on(False)
     self.ax.grid(), self.ax.legend(), plt.tight_layout()
     self.step0=False

    # Note: in interactive plots, unlike in a standard plot, the positions at the end of 
    #       N integration steps with timestep dt are plotted at each frame (determined by 'self.delay').
    #       Furthermore, the window update time 'self.delay' is dynamically adjusted
    #       based on the computational time required to execute all_methods for the N steps.

    t_ex=perf_counter()
    for _ in range(N): #N steps excuted without plotting
      self.stepleap(dt)
    t_ex-=perf_counter()
    #delay updating
    self.delay=max(int(-t_ex*1000), self.delay)
    #plot
    self.pointplot(planetlist, self.ax)
    self.widget.draw()
    if self.counter<self.Nmaxintegration:
      self.id=self.root.after(self.delay, self.leapfrog, dt, N)
    if self.openroot: 
      # Check if this is the first integration and open the window.
      # In subsequent cycles called by root.after, self.mainloop() will not be called again.
      self.openroot=False 
      self.mainloop()  

  def verlet(self, dt: numbers.Real, N: int, position=111):
    planetlist=self.planetlist
    self.counter+=1
    ##step0
    if self.step0: 
     self.ax=self.fig.add_subplot(position, projection=self.stringdim)
     self.initialpointplot(planetlist, self.ax, 'Verlet')
     self.ax.grid(), self.ax.legend(), plt.tight_layout()
     if self.fixed_axes:
        self.ax.relim(), self.ax.autoscale_view(), self.ax.set_autoscale_on(False)
     self.step0=False

    t_ex=perf_counter() 
    for _ in range(N): 
      self.stepverlet(dt)
    t_ex-=perf_counter() 
    self.delay=max(int(-t_ex*1000), self.delay) 
    self.pointplot(planetlist, self.ax)
    self.widget.draw()
    if self.counter<self.Nmaxintegration:
      self.id=self.root.after(self.delay, self.verlet,  dt, N)
    if self.openroot: 
      self.openroot=False 
      self.mainloop()  
 
  def RK4(self, dt: numbers.Real, N: int, position=111):
    planetlist=self.planetlist
    self.counter+=1
    ##step0
    if self.step0: 
     self.ax=self.fig.add_subplot(position, projection=self.stringdim)
     self.initialpointplot(planetlist, self.ax, 'RK4')
     self.ax.grid(), self.ax.legend(), plt.tight_layout()
     if self.fixed_axes:
        self.ax.relim(), self.ax.autoscale_view(), self.ax.set_autoscale_on(False)
     self.step0=False
    t_ex=perf_counter()
    for _ in range(N): 
      self.stepRK4(dt)
    t_ex-=perf_counter()
    self.delay=max(int(-t_ex*1000), self.delay) 
    self.pointplot(planetlist, self.ax)
    self.widget.draw()
    if self.counter<self.Nmaxintegration:
      self.id=self.root.after(self.delay, self.RK4, dt, N)
    if self.openroot: 
      self.openroot=False 
      self.mainloop()  

  def all_methods(self, dt: numbers.Real, N: int):
    planetlist=self.planetlist
    #step0#
    self.counter+=1
    if self.step0: 
     # Note: creating a list of copies of the system, so the 3 integrations are update simultaneously 
     self.lists=[deepcopy(planetlist), deepcopy(planetlist), planetlist] 
     self.axes=[self.fig.add_subplot(2, 2, i, projection=self.stringdim) for i in range(1, 4)]

     self.planetsfrog=self.step0leap(dt, planetlist=planetlist)
     for j,ax,lab in zip(self.lists, self.axes, ['Leapfrog', 'Verlet','Rk4']):     
        self.initialpointplot(j, ax, lab)
        ax.grid(), ax.legend(), plt.tight_layout()
        if self.fixed_axes:
          ax.relim(), ax.autoscale_view(), ax.set_autoscale_on(False)
     self.step0=False

    t_ex=perf_counter()
    for _ in range(N): 
      self.stepleap(dt, planetlist=self.lists[0])
      self.stepverlet(dt, planetlist=self.lists[1])
      self.stepRK4(dt, planetlist=self.lists[2])  
      #Oss: planetlists are updated
    t_ex-=perf_counter()
    self.delay=max(int(-t_ex*1000), self.delay)

    for lis, ax in zip(self.lists, self.axes):     
        self.pointplot(lis, ax)
    self.widget.draw()
    if self.counter<self.Nmaxintegration:
      self.id=self.root.after(self.delay, self.all_methods, dt, N) #check
    if self.openroot: 
      self.openroot=False 
      self.mainloop()  
    

# Command handled via the terminal
def main():
    parser=argparse.ArgumentParser(usage="Integrator of a system of planets/celestial bodies, using Verlet, RK4, or leapfrog method. \nSearch the link below to read the Readme file.\nhttps://drive.google.com/file/d/1Ijzu6Y3HtT49FA5Z2jiv9xNgvc7CqHMN/view?usp=sharing ")                                   
    # Adding nargs
    parser.add_argument('file', 
                        nargs="?", #number not fixed, so interactive mode is also available 
                        help="file.txt with initial conditions of planets in rows and 9 columns with:"+
                        "3D position, 3D velocity, mass, name, color of each planet, no delimiter needed, comments start with: "+"#"+".\nInitial conditions will be integrated with Leapfrog (dt=1000[s], Number of steps=1000).\nOnline file with initial conditions is available at https://drive.google.com/file/d/1S782eHvvw5998bK1rPyfQRtaujBv6WoD/view?usp=sharing")

    parser.add_argument("-i", 
                        action="store_true", #to store i as true if specified '-i' 
                        help='Interactive mode')

    # Parsing of arguments 
    arguments = parser.parse_args()
        
    #SI units
    # Sun starting velocity=[-1.12474e+01, 7.54876e+00, 2.68723e-01]
    # Sun
    Sun = Planet(
        [0.0, 0.0, 0.0],            # position x, y, z [m]
        [0.0, 0.0, 0.0],            # velocity vx, vy, vz [m/s]
        mass=1.98854e+30,
        name='Sun',
        color='yellow'
    )

    # Mercury
    Mercury = Planet(
        [-56939499000.0, -28343030000.0, 2907607800.0],
        [11649.7, -41479.3, -4459.52],
        mass=3.302e+23,
        name='Mercury',
        color='dimgray'
    )

    # Venus
    Venus = Planet(
        [42666101000.0, 99089170000.0, -1102842200.0],
        [-32293.0, 13696.0, 2050.91],
        mass=4.8685e+24,
        name='Venus',
        color='khaki'
    )

    # Earth
    Earth = Planet(
        [-143959899000.0, -40990530000.0, 1990300.0],
        [7651.51, -28751.4, 2.08354],
        mass=5.97219e+24,
        name='Earth',
        color='royalblue'
    )

    # Mars
    Mars = Planet(
        [-114927899000.0, -197277830000.0, -1313202200.0],
        [21836.9, -10113.2, -747.957],
        mass=6.4185e+23,
        name='Mars',
        color='firebrick'
    )

    # Jupiter
    Jupiter = Planet(
        [-567080899000.0, -578478830000.0, 15091377800.0],
        [9167.93, -8532.44, -169.767],
        mass=1.89813e+27,
        name='Jupiter',
        color='sandybrown'
    )

    # Saturn
    Saturn = Planet(
        [81869401000.0, -1503393830000.0, 22872377800.0],
        [9113.12, 496.372, -371.643],
        mass=5.68319e+26,
        name='Saturn',
        color='goldenrod'
    )

    # Uranus
    Uranus = Planet(
        [2624878101000.0, 1401746170000.0, -28782322200.0],
        [-3259.37, 5688.78, 63.2569],
        mass=8.68103e+25,
        name='Uranus',
        color='cadetblue'
    )

    # Neptune
    Neptune = Planet(
        [4302818101000.0, -1243213830000.0, -73569822200.0],
        [1471.32, 5253.63, -142.701],
        mass=1.0241e+26,
        name='Neptune',
        color='mediumblue'
    )

    # Pluto
    Pluto = Planet(
        [1655358101000.0, -4736013830000.0, 27812077800.0],
        [5245.41, 638.51, -1607.09],
        mass=1.307e+22,
        name='Pluto',
        color='lightgray'
    )

    # List collecting all planets and the Sun
    solar_system = [Sun, Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, Neptune, Pluto]
    print("Created list of planets 'solar_system'")
 
 
    # Actions for each arguments
    if arguments.i:
        print("Interactive mode is on.")
        # Interactive mode
        code.interact(local={**globals(), **locals()})#access to global and local variables
        

    elif not arguments.file==None:
       if isinstance(arguments.file, str):
        plot=Plot(arguments.file)
        plot.leapfrog(1000, 1000)
        plt.show()
    else:
        print("Invalid command, try -i for interactive mode, provide the file name with initial conditions, or try -h for help")


if __name__ == "__main__":
    # Name is '__main__' if the program is launched from the terminal,
    # otherwise it is the file name if imported as a library
    main()
