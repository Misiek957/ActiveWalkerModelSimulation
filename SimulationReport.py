#%%#Authors: Eli Gumble, Peter Brommer, Harry Brown. Incline and Elevation modifications created by Michal Nahlik

#Initialisation
import matplotlib.pyplot as plt
import numpy as np
# from scipy.integrate import simps
from scipy import signal as sg
from scipy.interpolate import RectBivariateSpline as ReBiSpline
# from numpy import ma
# from matplotlib import colors, ticker, cm
# from random import choice
import scipy.ndimage.filters as filters
import scipy.ndimage.morphology as morphology
import timeit
import math
import cv2
import random
from PIL import Image
import elevation_import as elev
%matplotlib inline

#%% User input Variables

areaName = 'Rhyd Ddu Path, Snowdon, Wales'
locSelect = [52.061386, -4.085166]
entrySelect = [[53.060900, -4.087302],
               [53.061653, -4.084131]]
areaWidth = 500

p_alpha = 1.6 # value of persistence
deg_max = 8 # maximum permissible ascent angle (30 degrees) # 2.860 = desirable, 3.732 = acceptable, 4.760 = min safety
theta_max = np.deg2rad(deg_max)

#%%

elevData = elev.Elevation(areaName, locSelect[0], locSelect[1], areaWidth, entrySelect)  # Extract the elevation data at selected location
elevValue = np.transpose(elevData.elev_interpolation())  # Store the interpolated elevation values
plt.imshow(elevValue)
plt.colorbar()

#%%
# Read grids from image
im = Image.open("images/flat500x500.bmp")
im = im.resize((areaWidth, areaWidth), Image.ANTIALIAS)

#%%

Base = np.array(im)
# Define internal quantities and variables
g_max = None
g_height = route = []
scale = 1  # m per pixel
Nx = Base[:,0,0].size  # N appears to be resolution
Ny = Base[0,:,0].size  # Nx,Ny is size, Nz is RGB level
xmin=-scale*0.5*(Nx-1)
xmax=scale*0.5*(Nx-1)
ymin=-scale*0.5*(Ny-1)
ymax=scale*0.5*(Ny-1)
x = np.linspace(xmin, xmax, Nx) # This is defining the axes and full space
y = np.linspace(ymin, ymax, Ny)
Y, X= np.meshgrid(y, x)
TrailPotential = np.zeros((Nx,Ny))
DestinationPotential=np.zeros((Nx,Ny))
Weight=np.zeros((Nx,Ny))  # Create gradient to sit on Nx, Ny
intens=np.zeros((Nx,Ny))
q_alpha=np.zeros((Nx,Ny))
expdist=np.zeros((2*Nx-1,2*Ny-1))
dest=np.zeros(2)
start=np.zeros(2)
grad=np.zeros((2,Nx,Ny))
vel=np.asarray([0.,0.])
pos=np.asarray([0.,0.])
#desdirx=ReBiSpline(x,y,grad[0,:,:],s=2)
#desdiry=ReBiSpline(x,y,grad[1,:,:],s=2)
intens[:]=0.

# Parameters
t_track=50. # Track decay time - after 50 walkers ignore a trail, it decays by 1/e
dt=0.1  # dt per time step, continuous markings every dt metres
dvel=3 # desired walker velocity in m/s
tau=5.
isigma=1./2. # trail potential
conv_thresh=10.e-4
precision=1.**2 #distance to target.
eps=0.025 #random motion contribution, same for all
max_steps = 25000 # maximum amount of steps taken by walker, must be greater than 500

#%%

class Slope:

    prevDir = 0 # record of previous angle of direction

    def __init__(self, posx, posy):
        try:
            xdir = 0.8*g_slope[0][posx][posy] # horizontal slope component at pos, scaled
            ydir = 0.8*g_slope[1][posx][posy] # vertical slope component at pos
        except IndexError:
            raise Exception("Walker travelled out of bounds..."
                  "\nInstruction: Decrease the walker persistence, or decrease the maximum permissible angle")
        self.xdir = xdir
        self.ydir = ydir
        self.mod = math.sqrt(xdir**2 + ydir**2)  # modulus of the gradient
        self.angle = np.arctan2(self.ydir,self.xdir)  # direction of the gradient vector  # np.arctan2() - gets result in correct quadriant
    # end __init__()

    def forbidden_angle(self, curDir):
        # Calculate forbidden angle (alpha)
        thetaSlope = math.atan(self.mod)  # angle of the slope relative from ground to the height
        if ((self.mod or thetaSlope) == 0 ) or (thetaSlope < theta_max):
            # print("Continue (CurDir = " + str(curDir) + ") (grad = " + str(self.angle) + ") slope = ",thetaSlope)
            return curDir
        alpha = (math.pi/2 - math.asin(math.tan(theta_max)/math.tan(thetaSlope)))  # Forbidden angle from the gradient value
        # Distinguish the deviation of the current direction from the closest forbidden angle and assign the new current direction accordingly
        curDir = curDir % (2*np.pi) # convert to 2pi
        # Difference between current direction and gradient, zipped around so it is always between +pi and -pi
        deviation = (curDir-self.angle + 0.5*np.pi) % np.pi - 0.5*np.pi  # deviation between pi/2 and -pi/2, absolute value
        # print('deviation = ', deviation)
        if np.abs(deviation)>alpha:
            # print('Final value unchanged \n (curDir = ' + str(curDir) + " grad = " + str(self.angle) + ") (alpha =" + str(alpha)+ ") slope = ",thetaSlope,")")
            return curDir # Current direction remains unchanged as it is outside of the forbidden zone
        if deviation < 0:
            # print("push right (CurDir =" + str(curDir) + ") (grad = " + str(self.angle) + ") (alpha =" + str(alpha)+ ") slope = ",thetaSlope,")")
            curDir = curDir-deviation-alpha  # curDir - deviation = grad (or opposite), curDir not subjected to mod
        else:
            # print("push left (CurDir =" + str(curDir) + ") (grad = " + str(self.angle) + ") (alpha =" + str(alpha)+ ") slope = ",thetaSlope,")")
            curDir = curDir-deviation+alpha
        # Convert back to [-pi,pi]
        if (0 <= curDir <= np.pi) or (curDir > 2 * np.pi):
            curDir = (curDir % (2*np.pi))
            return curDir
        elif (np.pi < curDir < 2 * np.pi) or (curDir < 0): # Place curDir within [-pi:0]
            curDir = (curDir % (2*np.pi)) - 2 * np.pi
            return curDir
    # end forbidden_angle()

    @staticmethod
    def persistence(curDir):
        global prevDir
        # Apply persistence of direction formula (Gilks equation 6) - weighted average of the movement
        gammax= (p_alpha*np.cos(prevDir) + (1 - p_alpha)* np.cos(curDir))
        gammay= (p_alpha*np.sin(prevDir) + (1 - p_alpha)* np.sin(curDir))
        gamma = np.arctan2(gammay,gammax)  # Previous direction as estimation
        return gamma

#%%

class ActiveWalkerModel:

    @staticmethod
    def init():
        global z, g_max, g_nat, g_grad ,g_height, route, Base
        # Set up map, Create blank arrays for map
        z = np.zeros((Nx,Ny))
        g_max=np.zeros((Nx,Ny)) # empty matrix
        g_nat=np.zeros((Nx,Ny))
        g_grad=np.zeros((Nx,Ny))
        g_nat=np.maximum(np.ones_like(g_nat),np.float64(Base[:,:,0])) # red channel, np.ones_like() sets minimum value to 1
        g_max=np.maximum(np.ones_like(g_max),np.float64(Base[:,:,1])) # green channel
        # Assign the requested elevation data as the height matrix for the model
        g_height = np.transpose(elevData.elev_interpolation())
        # g_height=np.fromfunction(lambda i, j: 0.27*(i ), (Nx, Ny), dtype=float) # Set custom gradient
        # g_height = cv2.GaussianBlur(g_height,(5,5),2) # Apply 2D convolution using a gaussian kernel
        z=g_nat

    @staticmethod
    def setup_weights():
        global Weight
        #Setup weight matrix, here trapezoid rule.
        Weight[:,:]=1
        Weight[1:-1,:]=2
        Weight[:,1:-1]=2
        Weight[1:-1,1:-1]=4
        Weight*=0.25*((x[-1]-x[0])/(Nx-1))*((y[-1]-y[0])/(Ny-1))

    @staticmethod
    def setup_distance_matrix():
        # Setup distance matrix
        for xi in range(1,Nx+1):
            for yi in range(1,Ny+1):
                expdist[xi-1,yi-1]=np.exp(-isigma*np.sqrt((x[Nx-xi]-xmin)**2+(y[Ny-yi]-ymin)**2))
                expdist[-xi,-yi]  = expdist[xi-1,yi-1]
                expdist[-xi,yi-1] = expdist[xi-1,yi-1]
                expdist[xi-1,-yi] = expdist[xi-1,yi-1]
        # find index range > conv_thresh
        subexpdist=expdist[(expdist>conv_thresh).any(1)]
        subexpdist=subexpdist[:, np.any(subexpdist>conv_thresh, axis=0)]
        #subexpdist=subexpdist[:,np.any(subexpdist>conv_thresh, axis=0)]
        #expdist[subexpdist]=0.
        #expdist
        #subexpdist
        return subexpdist

    @staticmethod
    def calc_tr_new():
        global TrailPotential, z, Weight, subexpdist
        TrailPotential[:,:]=sg.convolve2d(z[:,:]*Weight[:,:],subexpdist[:,:],mode="same")  # 2D convolution

    @staticmethod
    def set_up_walker(route_id): # input route_id commented for one route
        #set up walker
        global vel,pos,track,intens,dest,start,route
        start_id = route_id
        end =  len(elevData.entryMet)
        r = np.concatenate((np.array(range(route_id)),np.array(range(route_id+1,end))))
        end_id = int(random.choice(r))
        start=np.array(elevData.entryMet[start_id])  # commented for simplicity
        dest=np.array(elevData.entryMet[end_id]) # commented for simplicity
        vel=np.array([0.,0.])
        pos=np.array(start)
        track=np.zeros((max_steps,2)) # Set up for 5000 steps maximum
        #track[0,:]=pos[:]

    @staticmethod
    def setup_potentials():
        #Calculate gradients (eq 19 Helbing), Trail gradient
        global grad,desdirx,desdiry,dest
        grad=0.003*np.array(np.gradient(TrailPotential))
        #Destination potential
        DestinationPotential=-np.sqrt((dest[0]-x[:,None])**2+(dest[1]-y[None,:])**2)
        #Combine gradients
        grad+=np.array(np.gradient(DestinationPotential)[:])
        #Normalise
        #grad[:,:,:]/=(np.sqrt(grad[0,:,:]**2+grad[1,:,:]**2))
        desdirx=ReBiSpline(x,y,grad[0,:,:],s=2) # gradeint plus magnitude, Spline approximation over a rectangular mesh
        desdiry=ReBiSpline(x,y,grad[1,:,:],s=2)

    @staticmethod
    def calc_path():
        global pos,vel,intens,track,dest,dvel,tau, prevDir, desdirx, desdiry
        iterator=0
        hist=10
        samp=10
        avpos=np.zeros((2,hist))
        #Setup While loop to run until either the walker reaches the destination or the walker has passed max_steps movement cycles to
        #attempt to get there
        while np.dot(pos-dest,pos-dest)>precision and iterator<max_steps: # takes 5000 steps maximum
            #set the postiion of the walker on its first then subsequent cycles
            #conditional logic saying to update the average position of the walker every 10 iterations
            #if (iterator%samp==0): avpos[:,(i%hist)//samp]=pos[:] #ORIGINAL
            if iterator%samp==0:
                avpos[:,(iterator%(hist*samp))//samp]=pos[:]
            #print((iterator%hist)//samp)
            # print(avpos)
            gradmagnitude=max(0.0001,np.sqrt(desdirx(pos[0],pos[1])**2+desdiry(pos[0],pos[1])**2))
            xi=np.array(np.random.normal(0,1,2))
            # Equation 6 in Helbing, differential in position, eliminised velocity decay components
            # gradmagnitude makes sure it is normalised, desdir not normalised
            # pos[0]+= dt *(dvel * desdirx(pos[0],pos[1])/gradmagnitude +np.sqrt(2.*eps/tau)*xi[0])  # x-position vector component
            # pos[1]+= dt *(dvel * desdiry(pos[0],pos[1])/gradmagnitude +np.sqrt(2.*eps/tau)*xi[1])  # y-position vector component
            # posGrad = math.degree(math.atan(pos[0]/pos[1]) # future position
            curDir = math.atan2(desdiry(pos[0],pos[1]),desdirx(pos[0],pos[1])) # atan2 not to flip walker around
            posx = int(((pos[0]-xmin)*(Nx-1))/(xmax-xmin)) # index value of position x
            posy = int(((pos[1]-ymin)*(Ny-1))/(ymax-ymin)) # index value of position y
            # print('x =' + str(posx) + '  y = ' + str(posy)) # Debugging analysis
            # print('Curdir = ',curDir,' Persistent = ', persistence(curDir))
            if iterator>0: # Exception for first step
                curDir = Slope.persistence(curDir) # Apply persistence of movement
            curGradient = Slope(posx,posy)  # set the current gradient to the Slope class
            curDir = curGradient.forbidden_angle(curDir) # apply forbidden angle rule
            # print("CurDir at " + str(i) + "= " + str(curDir))
            pos[0] += dt * (dvel * math.cos(curDir))
            pos[1] += dt * (dvel * math.sin(curDir))
            prevDir = curDir # applying a back trace of the previous direction# Calculate current facing direction
            # Original
            # pos+=dt*vel
            #vel[0]+=-1/tau*vel[0] + (dvel/tau)*desdirx(pos[0],pos[1])/gradmagnitude+np.sqrt(2.*eps/tau)*xi[0]   # Eqiation 5 in Helbing, differential in velocity
            #vel[1]+=-1/tau*vel[1] + (dvel/tau)*desdiry(pos[0],pos[1])/gradmagnitude+np.sqrt(2.*eps/tau)*xi[1]
            #Set the current position of the walker into the track array for the current iteration
            track[iterator,:]=pos[:]
            # (pos[0]-xmin)*(Nx-1)/(xmax-xmin)  - scale to resolution [Nx-1] * distance to left edge / width of whole field
            try:
                intens[int((pos[0]-xmin)*(Nx-1)/(xmax-xmin)),int((pos[1]-ymin)*(Ny-1)/(ymax-ymin))]+=1.
            except IndexError:
                raise Exception("Walker travelled out of bounds..."
                      "\nInstruction: Decrease the walker persistence, or decrease the maximum permissible angle")
            iterator+=1
            if iterator%(hist*samp)==0:
                meanpos=np.mean(avpos,axis=1)
                if np.dot(pos-meanpos,pos-meanpos)<precision:
                    print ("Stalled progress ",pos,meanpos,vel, dest)
                    break
        if iterator==(max_steps - 500): print ("Missed goal\n",dest,pos)
        return iterator

    @staticmethod
    def update_ground():
        # Calculate Q_alpha (strength of markings) eq 15
        global q_alpha,intens,z,g_max,t_track,g_nat
        q_alpha=intens*(1.-z/g_max)
        # Time evolution of ground potential
        #zdiff=(1./t_track)*(g_nat-z)+q_alpha
        z+=(1./t_track)*(g_nat-z)+q_alpha
        #cs = plt.contourf(X, Y, zdiff, cmap=cm.PuBu_r)
        #cbar = plt.colorbar()
        #plt.show
        #z[140:160,45:75]

#%%

class Plot:
    @staticmethod
    def path_elev():
        # plt.contourf(Y, X, g_height, levels=np.linspace(g_height.min(),g_height.max(),1000),cmap='PuBu_r')
        plt.contour(Y, X, g_height, levels=np.linspace(np.amin(g_height),np.amax(g_height),16),inline=1)
        plt.colorbar()
        plt.scatter(track[0:max_steps-1,1],track[0:max_steps-1,0],1,c='y')
        plt.show(block=False)

    @staticmethod
    def path():
        # plt.contourf(X, Y, z, levels=np.linspace(z.min(),z.max(),1000),cmap='PuBu_r')
        plt.contourf(Y, X, z, levels=np.linspace(z.min(),z.max(),1000),cmap='PuBu_r')
        plt.colorbar()
        # plt.scatter(track[0:1999,0],track[0:1999,1],1,c='y')
        plt.show(block=False)

    @staticmethod
    def directions():
        #Plot the direction
        scgrad=np.arctan2(grad[1],grad[0])
        levels = np.linspace(-np.pi, np.pi, 360)
        cs = plt.contourf(Y, X,scgrad, levels=levels,cmap='hsv')
        # cbar = plt.colorbar()
        plt.colorbar()
        plt.scatter(track[0:1999,1],track[0:1999,0])
        #plt.scatter(start, dest)
        print(start)
        print(dest)
        plt.show()

    @staticmethod
    def potentials():
        global dest
        TotPot = np.zeros((Nx,Ny))
        TotPot =- np.sqrt((dest[0]-x[:,None])**2+(dest[1]-y[None,:])**2)
        TotPot += 0.003*TrailPotential
        maxima=ActiveWalkerModel.detect_local_maxima(TotPot)
        cs = plt.contourf(Y, X, TotPot, levels=np.linspace(TotPot.min(),TotPot.max(),1000),cmap='PuBu_r')
        # cbar = plt.colorbar()
        plt.colorbar()
        print(maxima)
        plt.scatter(y[maxima[1]],x[maxima[0]])
        plt.show()
        # commit test

#%% SETUP

elevData = elev.Elevation(areaName, locSelect[0], locSelect[1], areaWidth, entrySelect)  # Extract the elevation data at selected location
elevValue = elevData.elev_interpolation()  # Store the interpolated elevation values
elevValue = np.transpose(elevValue)
plt.contourf(elevValue)
plt.colorbar()
ActiveWalkerModel.init()
ActiveWalkerModel.setup_weights()
p_alpha = p_alpha/2 # scaling factor
subexpdist = ActiveWalkerModel.setup_distance_matrix()
# subexpdist.shape
timeit.timeit(ActiveWalkerModel.calc_tr_new,number=1)

#%% RUN SIMULATION

g_slope = np.gradient(g_height,1)
grady = g_slope[0] # vertical slope compontent
gradx = g_slope[1] # horizontal slope component
for i in range(0,50):
    ActiveWalkerModel.calc_tr_new()
    intens[:]=0.
    for j in range(0,2):
        ActiveWalkerModel.set_up_walker(np.random.randint(0,len(elevData.entryMet))) # np.random.randint(0,len(elevData.entryMet))
        ActiveWalkerModel.setup_potentials()
        # ActiveWalkerModel.calc_path()
        print(i, start," -> ", dest, pos, ActiveWalkerModel.calc_path())
    ActiveWalkerModel.update_ground()
    #plot_path()
print("End of simulation")

#%% Plot results images

elevData.map_image_get()
print("Entry = " + str(elevData.entryMet[0]) + '\n Dest'+ str(elevData.entryMet[1]))

Plot.path_elev()
Plot.path()
Plot.directions()
Plot.potentials()