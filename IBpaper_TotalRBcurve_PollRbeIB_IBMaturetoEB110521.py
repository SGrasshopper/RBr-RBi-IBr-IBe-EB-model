import random
from CellModeller.Regulation.ModuleRegulator import ModuleRegulator
from CellModeller.Biophysics.BacterialModels.CLBacterium import CLBacterium
import numpy
import math

from CellModeller.Integration.CLCrankNicIntegrator import CLCrankNicIntegrator 
from CellModeller.Integration.CLEulerSigIntegrator import CLEulerSigIntegrator
from CellModeller.Signalling.GridDiffusion import GridDiffusion 





#Launch Code: python Scripts/CellModellerGUI.py
#Launch Code Batch: python Scripts/batch.py
#for openGL:  conda install -c conda-forge pyopencl
#pip install pyopencl

max_cells = 2**15

#Specify parameter for solving diffusion dynamics #Add
grid_size = (4, 4, 4) # grid size
grid_dim = (64, 8, 12) # dimension of diffusion space, unit = number of grid
grid_orig = (-128, -14, -8) # where to place the diffusion space onto simulation space


def setup(sim):

    # Set biophysics, signalling, and regulation models
    # jitter turns on 3d
    # gamma controls growth inhibition from neighbors
    biophys = CLBacterium(sim, jitter_z=False, gamma = 200000)

    # add the planes to set physical  boundaries of cell growth
    #biophys.addPlane((0,-16,0), (0,1,0), 1)
    #biophys.addPlane((0,16,0), (0,-1,0), 1)
    sig = GridDiffusion(sim, 1, grid_dim, grid_size, grid_orig, [10.0])

    # Here we set up the numerical integration:
    # Crank-Nicholson method:
    integ = CLCrankNicIntegrator(sim, 1, 4, max_cells, sig, boundcond='reflect')
    # Alternative is to use the simple forward Euler method:
    #integ = CLEulerSigIntegrator(sim, 1, 2, max_cells, sig, boundcond='reflect')

    # use this file for reg too
    regul = ModuleRegulator(sim, sim.moduleName)    
    # Only biophys and regulation
    sim.init(biophys, regul, sig, integ)


    # Specify the initial cell and its location in the simulation
    sim.addCell(cellType=0, pos=(0,0,0)) 
    if sim.is_gui:

        # Add some objects to draw the models
        from CellModeller.GUI import Renderers
        therenderer = Renderers.GLBacteriumRenderer(sim)
        sim.addRenderer(therenderer)
        #sigrend = Renderers.GLGridRenderer(sig, integ)
        #sim.addRenderer(sigrend) #Add

    sim.pickleSteps = 10
    




def init(cell):


    # Specify mean and distribution of initial cell size
    #cell.targetVol = 1.0 + random.uniform(-0.01,0.01) # +/- (start, end)
    cell.targetVol = 2 #* numpy.random.normal(1, 0.05) #for normally distributed variance (middle, std)
    
    # Specify growth rate of cells
    cell.growthRate = 1.0 + random.uniform(-0.05,0.05) #remove either targetVol or GrowthRate variability
    cell.parentGrowth = [0] #progenitor cell logged growthRate
    
    cell.color = [2.0, 0.5, 1.5]#specify color of cell
    #cell.color = [1.0,0.0,0.0] #red
    #cell.color = [0.0,1.0,0.0] #green
    #cell.color = [0.0,0.0,1.0] #blue
    
    #RNA and protein 
    cell.rnaamt = [0,0,0,0] # RNA levels, used, in part, to drive geneamt levels
    cell.geneamt = [0.0, 0.0, 0.0, 0.0]   #[0]= hctA-GFP, [1]=Euo, [2]=HctA, [3]=HctB
    
    #EB to RB germination time
    cell.germTime = [(100 + random.uniform(-20,20))] #based on livecell and single cell expansion data: need to measure actually germ time variation and fit to dist

    #RBr > RBe conversion percent
    cell.percentchance = [0,0] #curve that drives RBr > RBe conversion


    #Specify initial concentration of chemical 
    cell.species[:] = [0.0] #species is concentration, I checked TC, cell growth causes steady state, cell division increases concentration
    cell.signals[:] = [0.0]

def specRateCL(): # Signal adds at rate k1
    return '''
    const float k0 = 0.0f;
    const float d0 = 0.3f;
    float x0 = species[0];
    rates[0] = k0 - d0*x0;
    '''
    # k1 = production rate of x0
    # d1 = degradation rate of x0

def sigRateCL(): #Add
    return '''
    const float k1 = 1.0f;
    float x0 = signals[0];
    rates[0] = k1;
    ''' 

time = 0
def update(cells):
    global time
    global n0
    time += 1
    time2 = (time/10)
    print('time = ' + str(time))
    print('time2 = ' + str(time2)) 

    #Iterate through each cell and flag cells that reach target size for division
    
    # Celltypes: 0=germ_EB, 1=RBr, 2=RBe, 3=IB, 4=pre_EB, 5=EB, 6=non-dividing RBs

    for (id, cell) in cells.items():
    
        if time >= cell.germTime[0]:

            #cell.percentchance[0] = (105/(1 + numpy.exp((1.76663094e+01-time2)*3.77251916e-01)) - 5) #on singlecell LVA counts
            cell.percentchance[0] = (97.81/(1 + numpy.exp((2.15841312e+01-((time2*cell.growthRate)))*6.77630536e-01)) + 2.19) #on livecell data and early RBe counts, percent chance of RB conversion
            #cell.percentchance[0] = (96.45042921/(1 + numpy.exp((13.60222209-time2)*1.4553212)) + 1.68390956)#early RBe counts
        
        
        #pr = RNA production rate
        #nr = RNA degradation rate
        pr0 = 0.04
        pr1 = 0.02
        pr2 = 0.04
        pr3 = 0.06
        nr0 = 0.01
        nr1 = 0.02
        nr2 = 0.01
        nr3 = 0.024        
                

        #p = protein production rate
        #n = protein degradation rate
        
        p0 = 1.0
        p1 = 0.5
        p2 = 1.0
        p3 = 0.5
        n0 = 0.05
        n1 = 0.08
        n2 = 0.05
        n3 = 0.01
        
        
        if cell.volume > cell.targetVol:
            cell.divideFlag = True
            
            
        if cell.cellType == 0: #germinating EB>RB
            cell.divideFlag = False
            cell.growthRate = 0.0
            if time >= cell.germTime[0]: #Become RBr if time is reached
                cell.cellType = 1 #RBr
                cell.growthRate = 1.0 + random.uniform(-0.05,0.05)
                cell.parentGrowth[0] = cell.growthRate
                
                

        if cell.cellType == 1: #RBr
            #cell.rnaamt[0] = cell.rnaamt[0] + (pr0 * cell.growthRate) - (nr0 * cell.rnaamt[0] * cell.growthRate) #extra RNA
            #cell.geneamt[0] = cell.geneamt[0] + (p0 * cell.growthRate * cell.rnaamt[0]) - (n0 * cell.growthRate * cell.geneamt[0]) #extra protein
            
            cell.rnaamt[1] = cell.rnaamt[1] + (pr1 * cell.growthRate) - (nr1 * cell.rnaamt[1] * cell.growthRate) #Euo RNA
            cell.geneamt[1] = cell.geneamt[1] + (p1 * cell.growthRate * cell.rnaamt[1]) - (n1 * cell.growthRate * cell.geneamt[1]) #Euo
            cell.geneamt[0] = 0 # GFP
            cell.geneamt[2] = 0 # HctA
            cell.geneamt[3] = 0 # HctB
            cell.color = [[1/cell.geneamt[1], 1, 1/cell.geneamt[1]]]
            #print('growthRate = ' + str(cell.growthRate))
            #print('percentchance = ' + str(cell.percentchance[0]))
            print('RBe Trigger Value = ' + str(cell.percentchance[0]))
            print('time2 = ' + str(time2))
            
            
            if time2.is_integer() and random.uniform(0,100) <= cell.percentchance[0]: # time2.is_integer() and  RBr to RBe. species[0] = magic Rbr>RBe signal
                #print('time2 is int' + str(time2))
                print('im an RBe')
                cell.cellType = 2 #RBe conversion
                cell.rnaamt[1] = cell.rnaamt[1] + (pr1 * cell.growthRate) - (nr1 * cell.rnaamt[1] * cell.growthRate) #Euo RNA
                cell.geneamt[1] = cell.geneamt[1] + (p1 * cell.growthRate * cell.rnaamt[1]) - (n1 * cell.growthRate * cell.geneamt[1]) #Euo
                cell.color = [[1/cell.geneamt[1], 1, 1/cell.geneamt[1]]] # color magic, fix
                #cell.geneamt[2] = 0 # HctA
                #cell.geneamt[3] = 0 # HctB
            
        #if cell.cellType == 1 and time2 >= 30: #RBr stop divide for CDI model
        #    cell.cellType = 6
        #    cell.color = [[1/cell.geneamt[1], 1, 1/cell.geneamt[1]]]
        #    cell.growthRate = 0    
        #    cell.geneamt[1] = cell.geneamt[1] - (n1 * cell.geneamt[1]) #Euo
        #    cell.geneamt[2] = 0 # HctA
        #    cell.geneamt[3] = 0 # HctB 


        if cell.cellType == 2: #need to keep if for dead cell
            cell.rnaamt[1] = cell.rnaamt[1] + (pr1 * cell.growthRate) - (nr1 * cell.rnaamt[1] * cell.growthRate) #Euo RNA
            cell.geneamt[1] = cell.geneamt[1] + (p1 * cell.growthRate * cell.rnaamt[1]) - (n1 * cell.growthRate * cell.geneamt[1]) #Euo 
            #cell.geneamt[2] = 0 # HctA
            #cell.geneamt[3] = 0 # HctB
            
            #if time2 >= 30: #RBe stop divide for CDI model
            #    cell.cellType = 6
            #    cell.color = [[1/cell.geneamt[1], 1, 1/cell.geneamt[1]]]
            #    cell.growthRate = 0
            #    cell.geneamt[1] = cell.geneamt[1] - (n1 * cell.geneamt[1]) #Euo
            #    cell.geneamt[2] = 0 # HctA
            #    cell.geneamt[3] = 0 # HctB
                
        if cell.cellType == 6: #AB
            cell.color = [[1/cell.geneamt[1], 1, 1/cell.geneamt[1]]] 
            cell.growthRate = 0
            cell.geneamt[1] = cell.geneamt[1] - (n1 * cell.geneamt[1]) #Euo

        if cell.cellType == 3: #IB
            cell.geneamt[1] = cell.geneamt[1] - (n1 * cell.parentGrowth[0] * cell.geneamt[1]) # Euo
            cell.rnaamt[2] = cell.rnaamt[2] + (pr2 * cell.parentGrowth[0])  - (nr2 * cell.rnaamt[2]) #hctA RNA
            cell.geneamt[2] = cell.geneamt[2] + (p2 * cell.parentGrowth[0] * cell.rnaamt[2]) - (n2 * cell.parentGrowth[0] * cell.geneamt[2]) #hctA
            cell.rnaamt[0] = cell.rnaamt[0] + (pr2 * cell.parentGrowth[0])  - (nr2 * cell.rnaamt[0]) #GFP RNA
            cell.geneamt[0] = cell.geneamt[0] + (p2 * cell.parentGrowth[0] * cell.rnaamt[0]) #GFP no degradation
            #cell.geneamt[3] = 0 # HctB
            cell.color = [[0, 0, cell.geneamt[2]*5]] #blue fast
            if cell.geneamt[2] >= 3.5:
                cell.cellType = 4 #pre_EB

                
       
        if cell.cellType == 4: #pre_EB
            cell.geneamt[0] = cell.geneamt[0] #GFP
            cell.geneamt[2] = cell.geneamt[2] - (n2 * cell.parentGrowth[0] * cell.geneamt[2]) #hctA no deg for mEos
            cell.rnaamt[3] = cell.rnaamt[3] + (pr3 * cell.parentGrowth[0])  - (nr3 * cell.rnaamt[3]) #hctB RNA
            cell.geneamt[3] = cell.geneamt[3] + (p3 * cell.parentGrowth[0] * cell.rnaamt[3])  - (n3 * cell.parentGrowth[0] * cell.geneamt[3]) #hctB
            cell.growthRate = 0
            cell.color = [[cell.geneamt[3]/10, 0, cell.geneamt[3]/40 + cell.geneamt[2]/7]] #blue to black to pink
            if cell.geneamt[3] >= 20: #hctB
                cell.cellType = 5 #infectious EB
                    

        if cell.cellType == 5:  #infectious EB
            cell.geneamt[0] = cell.geneamt[0] #GFP no degradation
            cell.geneamt[2] = cell.geneamt[2] - (n2 * cell.parentGrowth[0] * cell.geneamt[2]) #hctA deg
            cell.rnaamt[3] = cell.rnaamt[3] + (pr3 * cell.parentGrowth[0])  - (nr3 * cell.rnaamt[3]) #hctB RNA
            cell.geneamt[3] = cell.geneamt[3] + (p3 * cell.parentGrowth[0] * cell.rnaamt[3])  - (n3 * cell.parentGrowth[0] * cell.geneamt[3]) #hctB
            cell.color = [2.0, 0.0, 0.5] #pink
            

def divide(parent, d1, d2):
    # Specify target cell size that triggers cell division
    # Celltype1=RBr, Celltype2=RBe, Celltype3=IB, Celltype4=immature EB, Celltype5=mature EB
    
    
    if parent.cellType == 1: # If RBr: make 2RBrs
        d1.cellType = 1
        d1.targetVol = 2
        d1.growthRate = parent.parentGrowth[0] * numpy.random.normal(1, 0.05)
       
        d2.cellType = 1
        d2.targetVol = 2 #* numpy.random.normal(1, 0.05)
        d2.growthRate = parent.parentGrowth[0] * numpy.random.normal(1, 0.05)
        
        d1.geneamt[0] = parent.geneamt[0]/2
        d2.geneamt[0] = parent.geneamt[0]/2
        d1.geneamt[1] = parent.geneamt[1]/2
        d2.geneamt[1] = parent.geneamt[1]/2 


    if parent.cellType == 2: # If RBe: make 1RBe, 1IB
        d1.cellType = 2
        d1.targetVol = 2 
        d1.growthRate = parent.parentGrowth[0] * numpy.random.normal(1, 0.05)
        
        d2.cellType = 3
        d2.growthRate = 0
        d2.parentGrowth[0] = parent.growthRate
        d2.targetVol = 10
        
        d1.geneamt[0] = parent.geneamt[0]/2
        d2.geneamt[0] = parent.geneamt[0]/2
        d1.geneamt[1] = parent.geneamt[1]/2
        d2.geneamt[1] = parent.geneamt[1]/2 
        
               


#this model: GermEB lavendar, Rbr green, Rbe green, IB blue>black>red, EB  hot pink


# Rbr matures into Rbe based on percentchance curve from empirical data


# Rbe divides > Rbe and IB


# IB starts with hctA == 0, matures to EB(celltype == 5), based on HctB accumulation                 


        
# ISSUES WITH CURRENT MODEL




        
