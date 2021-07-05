# -*- coding: utf-8 -*-
"""
Created on Thu Feb 28 09:47:08 2019

@author: jrhillae
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Jan 25 13:40:18 2019

@author: jrhillae
@author: dbonte - reading input marram from neutral landscapes (28/11/19)
                - changes in lateral growth
                - initialisation phase of 3years (1095 days): allowing sand profile around existing vegetation
                - growth for another 15 y 
                - disturbance 2y before end of 20y simulations
                - storm event: some diaspores remain for subsequent vegetation development


"""
from mpl_toolkits.mplot3d import axes3d
import numpy as np
import random as rnd
import matplotlib.animation as animation
from matplotlib import cm
import matplotlib.pyplot as plt
import math
import pandas as pd
from nlmpy import nlmpy 
import time



                                                     
class Sand:
    
    #intialisatie
    def __init__(self,
                 max_x,
                 max_y,
                 sandflux,
                 sand_bulk_density,
                 length_cell):
        """Initializing class sand"""
        #numbers
        #boundary condition
        self.sandflux=sandflux
        #sand_bulk_density
        self.sand_bd = sand_bulk_density 
        
        #cumulative outflux of sand per day
        self.outflux=0
        #tracing cumulative outflux of sand per year at the landside (S)- sea is positioned North (N)
        self.outflux_landside=0
        #trace most recent average outflux (per cell) at lateral side=> input flux at a lateral side in the future, initial influx 0
        self.outflux_lateral=0  
        #arrays
        #initial height of sand is 5 meter at the back of the landscape
        initial_value = 5*sand_bulk_density*length_cell**2
        #total volume of sand, sand volume equals 0 in first row, and gradually increases with slope of 34° to 5m  #Adjusted to 10°
        self.Mass = np.zeros((max_x,max_y))
        for x in range(0,max_x):
            for y in range (0,max_y):
                self.Mass[x,y]=min(initial_value, x*math.tan(math.radians(10))*length_cell*\
                         sand_bulk_density*length_cell**2)
    
        
        #array with max q (saturated sand flux) per day per location
        self.qs=np.zeros((max_x,max_y))  
        #array with shear veloctiy per location (depending on vegetation)
        self.shear_velocity=np.zeros((max_x,max_y))
        #array with shear stress per location (depending on vegetation)
        self.shear_stress=np.zeros((max_x,max_y))
        #array with q per location. (important output for debugging)
        self.q=np.zeros((max_x,max_y))  
        #array with change in sand per cell per timestep
        self.change_sand=np.zeros((max_x,max_y))
        #array with growth speed of sand flux per location
        self.growthspeed=np.zeros((max_x,max_y))  
        #tracing total time cells have been without depostion
        self.time_no_deposition=np.zeros((max_x, max_y))
        #array with shelter
        self.shelter=np.zeros((max_x, max_y))
        #tracing total time cells have too much deposition
        self.time_too_much_deposition=np.zeros((max_x,max_y))
        #tracing total amount of sand deposited during winter per cell
        self.change_sand_winter=np.zeros((max_x, max_y))
        
                
class Ammophila:
    "Initializing class marram grass"
    def __init__(self, max_x, max_y, max_height, P, H, n):        
        self.Density = np.zeros((max_x,max_y))
        #Each cell has a chance P to be occupied with marram grass at the start of the simulation
        # clustered environment with 25% having resources
        
        #procedure to get the P correct
        n1=1  #at least one tussock of marram
        n0=n1*(1-P)/P

        nlm = nlmpy.mpd(nRow=max_y, nCol=max_x, h=H)
        ##print('nlm',nlm)
        self.Density = 20 * (nlmpy.classifyArray(nlm, [n0,n1]))   ##number zeros on ones; 
        
        nlmpy.exportASCIIGrid("grid"+str(P)+"_"+str(H)+"_"+str(n)+".txt", self.Density/20)
        
 #save initial density for joint count stats cfr field work       
        '''title='P'+str(P*100)+'H'+str(H*100)+'.txt'
        file=open(title, 'w')   
        file.write(self.Density)
        file.close() '''
        
        
        ##print('marram', self.Density)
        '''for x in range(max_x):
           for y in range (max_y):
                if rnd.random()<P:        
                    self.Density[x,y]+=20'''
        #define array to calculate growth speed per location 
        self.Growthspeed=np.zeros((max_x,max_y))
        #array necessary for lateral growth    
        self.new=np.zeros((max_x,max_y))
        #maximum height of marram grass
        self.max_height= max_height
        #density of marram grass before winter
        self.density_before_winter=np.zeros((max_x, max_y))
  
class Dune_dynamics:
    
    def __init__(self,
                 max_x,
                 max_y,
                 sandflux,
                 dist_winds,
                 mean_wind_velocity,
                 deviation_wind_velocity,
                 max_height,
                 sand_bulk_density,
                 length_cell,
                 P,
                 H,
                 n,
                 stormtime,
                 storm,
                 rain,
                 timesteps,
                 output_name,movie):
        '''Initialization and regulation of dune dynamics'''             
        print('start')
        #attributes of class dune dynamics
        self.max_x=max_x
        self.max_y=max_y
        self.sand = Sand(max_x,
                         max_y,
                         sandflux,
                         sand_bulk_density,
                         length_cell)
        self.ammophila = Ammophila(max_x,
                                   max_y,
                                   max_height,
                                   P,
                                   H, n)
        self.direction_wind='N'
        self.distribution_winds=dist_winds
        self.length_cell=length_cell
        self.movie=movie
        #print(self.sand.Mass)
        
        '''first=time.time()
        
        #defining output
        #2D output
        ims=[]
        fig1 = plt.figure()        
        ax1=fig1.add_subplot(1,2,1)
        ax2=fig1.add_subplot(1,2,2)
        
        #3D output
        ims2=[]
        fig2 = plt.figure()
        ax = fig2.gca(projection='3d')
        ax.set_zlim(0, 10)
        X = np.arange(0, max_x, 1)
        Y = np.arange(0, max_y, 1)
        X,Y= np.meshgrid(X,Y)'''
        
        #Debugging file
        title='Debugging'+'output_name'+'.txt'
        output=open(title, 'w')     
        
        #Initiating data collection for pandas dataframe
        '''data={"P":[], "H":[], "Replicate":[], "Time":[],
                  "Max_hoogte_zand":[],
                  "Tot_density_marram":[],
                  "Cumulative_outflux_per_year":[],
                  "Total_dune_volume":[]}'''
            
        #Running the simulation
        for t in range(timesteps):
            print('Time', t)
            
            #1: determine wind velocity
            self.determine_wind_velocity(mean_wind_velocity,
                                         deviation_wind_velocity,
                                         t)
         
            #2:In case a storm occurs, it happens after 60% of a simulation
            if storm==True and t==(timesteps*stormtime):
                self.Storm()
                
            #3:Rain events might occur with a chance of 20% from the middle of the simulation onwards
            if rain==True and t>timesteps/2 and t%365==0:
                if rnd.random()<0.20:
                    self.ammophila.Density+=1
            
            #4: determine wind direction
            self.wind_direction()
            
            #5: calculate shelters
            self.define_shelter(output)
            
            #6: implement sand dynamics
            self.Sand_dynamics(t, output)
            
            #7: save outflux of sand dynamics
            #only when wind direction is N (the sea), outflux is added to self.sand.outflux_landside
            if self.direction_wind=='N': 
                self.sand.outflux_landside+=self.sand.outflux
            elif self.direction_wind=='E' or self.direction_wind=='W':
                self.sand.outflux_lateral=self.sand.outflux/self.max_x
            
            #8: implement avalanches based on gravity
            self.gravity(output)
            
            #9: Marram grass dynamics
            # first day of spring: sprouting and set array change_sand_winter back to 0
            if t%365==0 and t>0:
                #self.SproutingSpring()
                self.sand.change_sand_winter.fill(0)
            #the first day of winter: save density of ammophila
            elif t%365==153:
                #important that this array is a copy and is stored separately in memory!!
                self.ammophila.density_before_winter=self.ammophila.Density.copy()
            #spring and summer: growth    
            if t%365<153:                  
                self.LocalGrowth()
            
            #NEW VERSION
            #hereunder if statement to allow lateral growth only after t*01,i.e. after 3y initialisation 
                if t>timesteps*0.15:
                    self.Growth_across_cells(output, t)
                
            #autumn and winter: submersion
            else:self.Submerging()
                
            #10: put landscape back in original position, wind north
            self.wind_north()
            
            #11: create beach
            self.First_row_zero()
           
            #At the end of a day, regulate output 
            #1) movies: in total each movie consists of 100 frames    
            
            
            '''##2D frames
            if t%(timesteps/100)==0:
                ima=ax1.imshow(self.ammophila.Density, animated=True,
                               cmap='Greens', interpolation='none', origin="upper")
                imb=ax2.imshow(self.sand.Mass, animated=True,
                               cmap='YlOrBr', interpolation='none', origin="upper")
                ims.append([ima, imb])
            ##3D frames
            if t%(timesteps/100)==0:
                
                imc=ax.plot_surface(X, Y, self.sand.Mass/(self.sand.sand_bd*(self.length_cell**2)), cmap=cm.YlOrBr,
                       linewidth=0, antialiased=False)
                #imc=ax.plot_wireframe(X, Y, self.sand.Mass/(self.sand_bd*(self.length_cell**2)), rstride=5, cstride=5)
                ims2.append([imc])
            
            '''
            #2) arrays and output parameters are saved each year
            if t%365==0:
                #print('yes')
                data["P"].append(P)
                data["H"].append(H)
                data["Replicate"].append(n)
                data["Time"].append(t)
                data["Max_hoogte_zand"].append(np.max(self.sand.Mass))
                data["Tot_density_marram"].append(np.sum(self.ammophila.Density))
                data["Cumulative_outflux_per_year"].append(self.sand.outflux_landside)
                data["Total_dune_volume"].append(np.sum(self.sand.Mass/(self.sand.sand_bd*self.length_cell**2)))#kg
                '''np.save('Topography'+str(t)+'.npy', self.sand.Mass/(self.sand.sand_bd*self.length_cell**2))
                np.save('Marram'+str(t)+'.npy', self.ammophila.Density) '''
                
                #print(data)
                #set self.sand.outflux_landside back to zero before the start of the next year
                self.sand.outflux_landside=0
            
        #save and export pandas dataframe
        df=pd.DataFrame.from_dict(data)
        df.to_csv('dataframe'+str(output_name)+'.csv')
        #exporting array as dataframe:
        #pd.DataFrame.from_records(x)
        #now exporting array as array
        '''np.save('Sand'+str(output_name)+'.npy', self.sand.Mass)
        np.save('Marram'+str(output_name)+'.npy', self.ammophila.Density)
        #create and save 2D and 3D animation  
        
        
        if movie:  
            ani1 = animation.ArtistAnimation(fig2, ims2, interval=200, blit=False, repeat_delay=1000) 
            ani2 = animation.ArtistAnimation(fig1, ims, interval=200, blit=False, repeat_delay=1000)
            ani1.save('3D'+str(output_name)+'.mp4',  writer='ffmpeg', dpi=200)
            ani2.save('2D'+str(output_name)+'.mp4',  writer='ffmpeg', dpi=200)
        output.close()
        delta=time.time()-first
        print('time of simulation:',delta)
        
        plt.show()'''
        
    def determine_wind_velocity(self,
                                mean_wind_velocity,
                                deviation_wind_velocity,
                                t):
        '''Sample wind velocity from normal distribution based on average and standard deviation of each month'''
        wv_index=int(t/30.42%12)
        self.wind_velocity=np.random.normal(mean_wind_velocity[wv_index],deviation_wind_velocity[wv_index])
        while self.wind_velocity<0 or self.wind_velocity>9:
            self.wind_velocity=np.random.normal(mean_wind_velocity[wv_index], deviation_wind_velocity[wv_index])
            
    def Storm(self):
        '''Implementing storm event within the first quarter rows of the landscape'''
        for y in range (self.max_y):
            for x in range (int(self.max_x/4)):
                if self.ammophila.Density[x,y]==0: self.sand.Mass[x,y]=0  #***without vegation all sand is eroded
                
                else: 
                #***new version to reduce ammophila density; roots remain present 
                #if storm then 5% of sand remains and 5% of vegetaion
                    self.ammophila.Density[x,y]=self.ammophila.Density[x,y]*0.01
                    self.sand.Mass[x,y]=self.sand.Mass[x,y]*0.05
                
    def define_slope(self, Tot_height):
        '''Calculate average slope of Tot_height'''
        #array with length of slope per location
        length_of_slope=np.zeros((self.max_x, self.max_y))
        slope=np.zeros((self.max_x, self.max_y))
         
        for y in range(self.max_y):
            
            #x==0, determine the start condition
            Start=0
            Flat, Increase, Decrease= False, False, False
            #determine start condition:
            if Tot_height[0,y]==Tot_height[1,y]:
                Flat= True
                
            elif Tot_height[0,y]< Tot_height[1,y]:
                Increase=True
                
            elif Tot_height[0,y]>Tot_height[1,y]:
                Decrease=True 
            
            for x in range (1, self.max_x-1):
                                  
                if Tot_height[x,y]==Tot_height[x+1,y]:
                    if Flat==False: 
                        Flat, Increase, Decrease= True, False, False
                        Slope= (Tot_height[x,y]-Tot_height[Start,y])/((x-Start)*self.length_cell)
                        for a in range (Start+1, x+2):
                            slope[a,y]=Slope  
                            length_of_slope[a,y]=x-Start             
                        Start=x
                    
                elif Tot_height[x,y]<Tot_height[x+1,y]:
                    if Increase==False:
                        Flat, Increase, Decrease= False, True, False
                        Slope= (Tot_height[x,y]-Tot_height[Start,y])/((x-Start)*self.length_cell)
                        for a in range (Start+1, x+2):
                            slope[a,y]=Slope 
                            length_of_slope[a,y]=x-Start
                        Start=x              
                
                elif Tot_height[x,y]>Tot_height[x+1,y]:
                    if Decrease==False:
                        Flat, Increase, Decrease= False, False, True
                        Slope= (Tot_height[x,y]-Tot_height[Start,y])/((x-Start)*self.length_cell)
                        for a in range (Start+1, x+2):
                            slope[a,y]=Slope
                            length_of_slope[a,y]=x-Start 
                        Start=x
            #determine end condition:
            x=self.max_x-1
            Slope= (Tot_height[x,y]-Tot_height[Start,y])/((x-Start)*self.length_cell)
            for a in range (Start+1, x+1):
                slope[a,y]=Slope 
                length_of_slope[a,y]=x-Start 
                
        return slope, length_of_slope
        
    def First_row_zero(self):
        '''Turn sand mass of first row to zero to create beach'''
        for y in range(self.max_y):
            self.sand.Mass[0,y]=0
        
    def wind_direction(self):
        '''Determine wind direction and rotate landscape accordingly'''
        #define wind direction:
        chance=rnd.random()
        #defining_Future Wind
        if chance < self.distribution_winds[0]:
            Future_wind='N'
        elif chance < self.distribution_winds[0]+self.distribution_winds[1]:
            Future_wind='W'
        elif chance < 1- self.distribution_winds[3]:
            Future_wind='S'
        else: 
            Future_wind='E'
            
        Winds=['N', 'W', 'S', 'E']
        
        
        #index future wind
        ind_fw=Winds.index(Future_wind)
        
        if ind_fw>0:
            self.sand.Mass=np.rot90(self.sand.Mass, ind_fw)
            self.ammophila.Density=np.rot90(self.ammophila.Density, ind_fw)
            self.sand.time_no_deposition=np.rot90(self.sand.time_no_deposition, ind_fw)
            self.sand.time_too_much_deposition=np.rot90(self.sand.time_too_much_deposition, ind_fw)
            self.sand.change_sand_winter=np.rot90(self.sand.change_sand_winter, ind_fw)
            self.ammophila.density_before_winter=np.rot90(self.ammophila.density_before_winter, ind_fw)
            
        self.direction_wind=Future_wind   
        
    def wind_north(self):
        '''Put landscape back in original position with wind direction North'''
        #determine index of current wind        
        Winds=['N', 'W', 'S', 'E']
        ind_cw=Winds.index(self.direction_wind)
        
        if self.direction_wind!='N':
            
            self.sand.Mass=np.rot90(self.sand.Mass, 4 - ind_cw)
            self.ammophila.Density=np.rot90(self.ammophila.Density,4-ind_cw)
            self.sand.time_no_deposition=np.rot90(self.sand.time_no_deposition,4-ind_cw)
            self.sand.time_too_much_deposition=np.rot90(self.sand.time_too_much_deposition,4-ind_cw)
            self.sand.change_sand_winter=np.rot90(self.sand.change_sand_winter,4-ind_cw)
            self.ammophila.density_before_winter=np.rot90(self.ammophila.density_before_winter,4-ind_cw)
            
            self.direction_wind='N'
            
    def define_shelter(self, output):
        '''Determine position of shelters in landscape
        Updating shelter array, shelter: 1 , no shelter: 0 '''
        #step 1: fill the array with zeros
        self.sand.shelter.fill(0)   
        
        #step 2: define height per cell (based on vegetation and sand)
        height_sand=(self.sand.Mass/(self.sand.sand_bd*self.length_cell**2)) #m
        height_vegetation= ((np.sqrt(self.ammophila.Density/100))*self.ammophila.max_height)#m
        Tot_height= height_sand + height_vegetation #m
        
        ##print('Tot_height', Tot_height)
        #step 3, determine slope per location and total length of slope
        slope, length_of_slope= self.define_slope(Tot_height)
        
        #step 4, determine shelter areas
        for y in range(self.max_y):
            #per negative slope, the shelter is only determined once (in the first cell of the slope)
            Neg_slope=False
            for x in range(self.max_x):
                #only slopes more negative than -14 degrees slope have a shelter
                if (slope[x,y]< math.tan(math.radians(-14))) and Neg_slope==False:
                    #calculate lenght of the shelter
                    shelter_length=((slope[x,y]*(length_of_slope[x,y]*self.length_cell))/math.tan(math.radians(-14)))/self.length_cell
                    #print(shelter_length)
                    #loop over each cell in the shelter and change shelter array to 1 
                    #loop stops when the height of a cell exceeds the height of the shelter (considering the shelter is a triangle)
                    q=x
                    #print(int(x+length_of_slope[x,y])-1)
                    #print( q, x+shelter_length, math.tan(math.radians(14))*((shelter_length-(q-x))*0.20))
                    #shelter_length is float, therefore ceil is applied
                    while (q<self.max_x) and (q<math.ceil(x+shelter_length)) and \
                    (height_sand[q,y]<(math.tan(math.radians(14))*((shelter_length-(q-x))*self.length_cell)) +\
                     Tot_height[int(x+length_of_slope[x,y]-1), y]):
                        self.sand.shelter[q,y]=1
                        q+=1
                    Neg_slope=True  
                
                if (slope[x,y]>0) or (slope[x,y]==0):
                    Neg_slope=False  
                 
        
    def Sand_dynamics(self, 
                      t,
                      output):
        '''Implement sand dynamics'''
        self.sand.change_sand.fill(0)
        self.sand.outflux=0
        #a wind velocity of 3.87 m/s at a height of 10m is the threshold for sand dynamics
        if self.wind_velocity>3.87:
            #First, calculate shear velocity, shear stress and qs per location in grid
            #general rule: shear_velocity=sqrt(shear_stress/sigma_a)
            shear_stress_max=((self.wind_velocity*0.058)**2) *1.25
            
            #updating shear stress based on distribution of vegetation
            self.sand.shear_stress=shear_stress_max/(1+(16*(self.ammophila.Density/100)))

                
            #including ventura effect
            for x in range (self.max_x):  
                #occupied: marram was present in previous cell
                #start: x coordinate of first cell with marram in continuous row
                #sum_reduction: cumulative reduction of shear stress in continuous row of marram grass
                Occupied, start, Sum_reduction = False, False, 0            
                for y in range (self.max_y):
                    
                    if self.ammophila.Density[x,y]>0:
                        if Occupied==False:
                            Occupied, start, Sum_reduction= True, y, (shear_stress_max - self.sand.shear_stress[x,y])
                            
                        else: 
                            Sum_reduction += (shear_stress_max - self.sand.shear_stress[x,y])
                    
                    else:
                        if Occupied == True:
                            #avoid wrapped boundaries
                            if start-1>=0:
                                self.sand.shear_stress[x, start-1] = self.sand.shear_stress[x, start-1]+(Sum_reduction*0.25)
                            self.sand.shear_stress[x, y] = self.sand.shear_stress[x, y]+(Sum_reduction*0.25)
                            Occupied, start, Sum_reduction = False, False, 0         
                #checking last cell
                if Occupied == True:
                        #avoid wrapped boundaries
                        if start-1>0:
                            self.sand.shear_stress[x, start-1] = self.sand.shear_stress[x, start-1]+(Sum_reduction*0.25)
            
            output.write('shear stress after ventura' + '\n')
            for y in range(10):
                for x in range (10):
                    output.write(str(round(self.sand.Mass[x,y],4))+ '\t')
                output.write('\n')                
            
            #updating array with shear velocity
            self.sand.shear_velocity= np.sqrt(self.sand.shear_stress/1.25)
    
            #calculate qs per location per day based on prevalent shear velocity
            #Bagnold formula
            #qs= C*sigma_a/ g*sqrt(dn/Dn)*(u-ut)**3
            #C= 1.8, sigma_a= 1.25 kg/m^3, g=9.81 m/s^2, alpha=0.058, sigma_p= 2650 kg/m^3,
            #dn= 335 micrometer, Dn= 250 micrometer, z=10 m
            #original formula expresses qs in kg/m/s => correct for lenght of one cell
            self.sand.qs=(1.8* (1.25/9.81)*math.sqrt(335/250)*\
                          (self.sand.shear_velocity-(3.87*0.058))**3)*3600*24*self.length_cell
                          
            self.sand.qs[self.sand.qs<0]=0
            
            #define boundary condition at windward side,
            #expressed relative to qmax (percentage determined by parameter self.sand.sandflux)
            
            if self.direction_wind=='N':
                flux_firstcell= min(10, self.sand.sandflux*(1.8* (1.25/9.81)*math.sqrt(335/250)*\
                                                ((self.wind_velocity-3.87)*0.058)**3)*3600*24*self.length_cell)
            elif self.direction_wind=='E' or self.direction_wind=='W':
                flux_firstcell= self.sand.outflux_lateral
            else: flux_firstcell=0
            
            #Second, calculate growth speed of sand dynamics (i.e. 1/lsat) per location
            self.sand.growthspeed=0.2*self.sand.shear_velocity
                    
            for y in range(self.max_y):            
                Flux= flux_firstcell
                for x in range(self.max_x):
                    #save flux that enters the cell
                    self.sand.q[x,y]=Flux
                    #erosion
                    if Flux< self.sand.qs[x,y] and self.sand.shelter[x,y]==0:
                        change_flux=Flux*self.sand.growthspeed[x,y]*(1-(Flux/self.sand.qs[x,y]))
                        if self.sand.Mass[x,y]<change_flux:
                            change_flux=(self.sand.Mass[x,y])
                    #deposition
                    elif Flux> self.sand.qs[x,y]:
                        change_flux= -(Flux-self.sand.qs[x,y])/5                  
                    #Flux*0.2*(1-(Flux/self.sand.qs[x,y]))
                    #no change
                    else: change_flux=0
                    
                    #trace sand dynamics
                    Flux+=change_flux                
                    self.sand.change_sand[x,y]=-change_flux
                #add remaining sand at end of landscape to outflux    
                self.sand.outflux+=Flux
            
            #update array with total amount of sand
            self.sand.Mass+=self.sand.change_sand      
                        
    def gravity(self, output):
        '''Implement avalanches by gravity'''
        avalanches=np.zeros((self.max_x, self.max_y))
        
    
        directions=[(1,0),(0,1),(1,1),(-1,0),(0,-1),(-1,1),(1,-1),(-1,-1)]
        #Daily, shuffle items witin directions
        rnd.shuffle(directions)
        #Max_slope becomes steeper with increasing ammophila density
        max_slope=34+0.034*self.ammophila.Density
        
        for x in range (self.max_x):
            for y in range (self.max_y):
                #checking for too steep slopes
                #initially assume there is no avalanche
                            
                avalanche=False
                direction_avalanche=(0,0)  
                #trace direction of steepest slope > 0.67 or 34 degrees
                for k,l in directions:
                    #all tested cells are positioned within the landscape
                    if 0<=x+k<self.max_x and 0<=y+l<self.max_y:
                        #determine slope in that direction
                        slope=(((self.sand.Mass[x,y]/(self.sand.sand_bd*self.length_cell**2))-\
                                (self.sand.Mass[x+k, y+l]/(self.sand.sand_bd*self.length_cell**2)))/(math.sqrt((k*self.length_cell)**2+(l*self.length_cell)**2)))
                        #avalanche when slope steeper than 0.67
                        if slope>math.tan(math.radians(max_slope[x,y])):
                            if avalanche:
                                #comparing slope with other too steep slopes starting from same cell
                                #avalanche occurs in direction of steepest slope
                                if slope>avalanche:
                                    direction_avalanche=k,l
                                    avalanche= slope
                            #previously, no steep slope detected
                            else:
                                direction_avalanche=k,l
                                avalanche=slope
                                
                #implementing movement of sand by avalanche
                if avalanche != False:
                    v,w=direction_avalanche
                    amount_sand=((self.sand.Mass[x,y]/(self.sand.sand_bd*self.length_cell**2))-\
                                  ((self.sand.Mass[x+v, y+w]/(self.sand.sand_bd*self.length_cell**2)+\
                                    math.tan(math.radians(max_slope[x,y]))*math.sqrt((v*self.length_cell)**2+(w*self.length_cell)**2))))*\
                                    (self.sand.sand_bd*self.length_cell**2)
                    
                    #avalanche
                    avalanches[x,y]-=amount_sand
                    avalanches[x+v, y+w]+=amount_sand
                  
        
        #after implementing sand and gravity dynamics        
        #update array time_no_deposition, + 1 when no deposition, otherwise 0
        self.sand.Mass+= avalanches
        self.sand.change_sand+= avalanches
        
        self.sand.time_no_deposition[self.sand.change_sand<0.000001]+=1
        self.sand.time_no_deposition[self.sand.change_sand>=0.000001]=0
        #maximum 1 meter of sand during growth season
        self.sand.time_too_much_deposition[self.sand.change_sand>1/153]+=1
        self.sand.time_too_much_deposition[self.sand.change_sand<=1/153]=0
        
    def SproutingSpring(self):
        '''First day of spring, a large sprouting event occurs depending on
        density of marram grass before winter and amount of sand deposited during winter'''
        #calculate density of marram grass first day of spring based on density before winter and sand deposition
        #according to this linear relationship, per 1 unit of density, its density stays 1 if no sand was deposited, 
        #but density becomes 0 if 1 m of sand was deposited during winter
        self.ammophila.Density= 5*(self.ammophila.density_before_winter*\
            (-(self.sand.change_sand_winter/(self.sand.sand_bd*self.length_cell**2))+1))
        
        self.ammophila.Density=np.clip(self.ammophila.Density,0,99.99)
        #if change in sand larger than 1m, density of marram grass becomes 0
        self.ammophila.Density= np.where((self.sand.change_sand_winter/(self.sand.sand_bd*self.length_cell**2))>1,
                                         0, self.ammophila.Density)
        #if change in sand negative during winter, density of marram grass does not change
        self.ammophila.Density= np.where((self.sand.change_sand_winter/(self.sand.sand_bd*self.length_cell**2))<0,
                                         self.ammophila.density_before_winter, self.ammophila.Density)
             
    def LocalGrowth(self):
        '''Regulate local growth of marram grass'''
        #determining local growth speed based on sand deposition
        self.ammophila.Growthspeed=-462.08*(self.sand.change_sand/(self.sand.sand_bd*self.length_cell**2)\
                                            -(0.5/152))**2+ 0.005
        #print('growthspeed ammophila', self.ammophila.Growthspeed)
        self.ammophila.Growthspeed[self.ammophila.Growthspeed<0.000000000000001]=0
        #negative growth in case of too much or too few sand deposition
        self.ammophila.Growthspeed[self.sand.time_no_deposition>10]=-0.05
        self.ammophila.Growthspeed[self.sand.time_too_much_deposition>10]=-0.02
        #adding randomness to growthspeed
        self.ammophila.Growthspeed+= np.random.normal(0, 0.001, (self.max_x, self.max_y))
        #logistic growth
        self.ammophila.Density += self.ammophila.Growthspeed*\
        self.ammophila.Density*(100-self.ammophila.Density)/100
        #maximum should be 99.99 not 100 (when 100, value remains 100 and does not decrease when growth is negative).
        self.ammophila.Density=self.ammophila.Density.clip(0,99.99)
        #0.00001 is the cut off value for presence of ammophila
        self.ammophila.Density[self.ammophila.Density<0.00001]=0
        
      
    def Growth_across_cells(self, output, t):
        '''Regulate lateral growth'''
        self.ammophila.new.fill(0)
        #calculate chance of lateral growth, depending on sand deposition
        #***new version with chance of lateral growth indepednent of sand flux
                        
        chance_lateral_growth=-14440*(self.sand.change_sand/(self.sand.sand_bd*self.length_cell**2)\
                                       -(0.4/152))**2+0.1
        
                                       
        randomness=np.random.normal(0, 0.005, (self.max_x, self.max_y))
        
        
        chance_lateral_growth= np.where(self.ammophila.Density>0, chance_lateral_growth + randomness, 0)
        
                
        #print('chance_lateral_growth', chance_lateral_growth)
        for x in range(self.max_x):
            for y in range(self.max_y):
                #determining whether cell will sprout
                
                #*** NEW VERSION make a 0.25% chance anyway if density is high (set to 80* fill)
                if (rnd.random()<0.0025) and (self.ammophila.Density[x,y]>80): chance_lateral_growth[x,y]=1
                
                if rnd.random()<chance_lateral_growth[x,y]: 
                    #assert self.ammophila.Density[x,y]>0, print(chance_lateral_growth[x,y], self.ammophila.Density[x,y])
                    #Where will cell sprout?: 10% chance of sprouting 3 cells away  - levy walk
                    if rnd.random()<0.1:
                        a,b=rnd.choice([(w,v) for w in [-3,3] for v in range(-3,4)] +\
                                        [(v,w) for w in [-3,3] for v in range(-2,3)])
                    else:
                    #90% chance of sprouting in neighbouring cell
                        a,b= rnd.choice([(-1,-1), (0,-1), (1,-1), (-1,0), (1,0), (-1,1), (0,1), (1,1)])
                        
                    #checking wheter new cell is positioned in the landscape
                    if (0<=x+a<self.max_x) and (0<=y+b<self.max_y):
                        #print("yes")
                        if self.ammophila.Density[x+a,y+b] < 90:
                            self.ammophila.new [x+a,y+b] += 1
                    
        self.ammophila.Density+=self.ammophila.new

            
    def Submerging (self):
        '''In winter, marram grass is submerged'''
        #trace netto change in sand deposition or erosion per cell during autumn and winter
        self.sand.change_sand_winter+= self.sand.change_sand
        #updating density of marram grass based on level of submerging during winter
        #print('Density', self.ammophila.Density)
        #h1: height before submerging
        h1=(np.sqrt(self.ammophila.Density/100))*self.ammophila.max_height
        ##print('h1', h1)
        #h2[h1>0] efficienter maar werkte niet, kijken hoe dit komt
        h3=np.zeros((self.max_x, self.max_y))
        #submersion only occurs where ammophila is present h3: 1: ammophila present, 0: no ammophila present
        h3[h1>0]=1
        #h2: height after submerging
        h2= np.where(self.sand.change_sand>0, (h1-(self.sand.change_sand/(self.sand.sand_bd*(self.length_cell**2))))*h3, h1)    
        #print('h2', h2)        
        # plant completely submerged=> height becomes 0
        h2[h2<0]=0
        #calculate new local densities based on h2
        #self.ammophila.Density=((h2/self.ammophila.max_height)**2)*100
        #print('Density', self.ammophila.Density)

#start of simulation

##♣open data for pandas
data={"P":[], "H":[], "Replicate":[], "Time":[],
                  "Max_hoogte_zand":[],
                  "Tot_density_marram":[],
                  "Cumulative_outflux_per_year":[],
                  "Total_dune_volume":[]}


for rep in range(1,6):
    
        
        for SpatCor in range(0,11):
            Hirsh=SpatCor/10
                       
            for PropMarram in range(1,10):
                Prop=PropMarram/10
     
            
                print("P=",Prop," H=",Hirsh," Repnr",rep)
                Dune_dynamics(
                    max_x= 100, #number of cells
                    max_y= 100, #number of cells
                    sandflux= 0.2,#relative percentage of qmax
                    dist_winds=[0.7,0.15,0.15,0], #distribution of wind directions[north, west, south, east]
                    #average wind strenght and its standard deviation per month, first month April, last month March 
                    #first day of a simulation is April 1st (first day of marram growth in spring)
                    mean_wind_velocity=[5.59, 5.35, 6.83, 5.31, 5.77, 6.21, 6.29, 5.80, 7.56, 8.29, 7.71, 6.23],#m/s
                    deviation_wind_velocity=[3.33, 2.74, 2.42, 2.59, 2.95, 3.07, 3.58, 2.79, 4.02, 4.75, 3.21, 3.41],
                    max_height=0.90, #m [Huiskes 1979]
                    sand_bulk_density= 1500,#kg/m^3
                    length_cell=0.2, #length of one cell #m
                    P=Prop,#Percentage of marram grass in landscape at initialization
                    H=Hirsh, #spatial correlation, Hirsh factir at initialization
                    n=rep, #for having replicates in the loop
                    stormtime=0.9,  #relative time when storm should appear - 0.9 is last two years
                    storm=False,  #storm at timesteps*stormtime
                    rain=False,
                    timesteps=20*365+1,#total runtime
                    output_name='_P10_100x100_per_rep_all',
                    movie=False) 


