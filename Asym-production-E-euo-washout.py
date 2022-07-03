import random
import numpy
import math
from CellModeller.Regulation.ModuleRegulator import ModuleRegulator
from CellModeller.Biophysics.BacterialModels.CLBacterium import CLBacterium
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
    biophys = CLBacterium(sim, jitter_z=False, gamma = 2000)

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
    cell.rnaamt = [0.0, 0.0, 0.0, 0.0] # RNA levels, used, in part, to drive geneamt levels
    cell.geneamt = [0.0, 0.0, 0.0, 0.0]   #[0]= Ectopic protein, [1]=Euo, [2]=HctA, [3]=HctB
    
    #EB to RB germination time
    cell.germTime = [(100 + random.uniform(-20,20))] #based on livecell and single cell expansion data: need to measure actually germ time variation and fit to dist

    #RBr > RBe conversion percent
    cell.percentchance = [0,0] #curve that drives RBr > RBe conversion

    #Specify initial concentration of chemical 
    cell.species[:] = [0.0] #species is concentration, normal per cell = * volume
    cell.signals[:] = [0.0]

def specRateCL(): # Signal adds at rate k0
    return '''
    const float k0 = 20.0f;
    const float d0 = 0.0f;
    float x0 = species[0];
    rates[0] = k0 - d0*x0;
    '''
    # k0 = production rate of x0
    # 10.0f adds 0.5 every update step
    # 20.0f adds 1 every update step
    # 100.0f adds 5 every update step
    # 200.0f adds 10 every update step
    # d0 = degradation rate of x0

def sigRateCL(): #Add
    return '''
    const float k1 = 1.0f;
    float x0 = signals[0];
    rates[0] = k1;
    ''' 
    
time = 0
def update(cells):
    global time
    time += 1
    #print('time hours = ' + str(time/10))
    
    # Celltypes: 0=germ_EB, 1=RBr, 2=RBe, 3=IB, 4=pre_EB, 5=EB, 6=non-dividing RBs
    # Iterate through each cell
    for (id, cell) in cells.items():
        print('cell sp = ' + str(cell.species[0]))
        print('cell sp norm = ' + str(cell.species[0]*cell.volume))
        print('time = '+ str(time))
    
        if time >= cell.germTime[0]:

            #cell.percentchance[0] = (105/(1 + numpy.exp((1.76663094e+01-(time/10))*3.77251916e-01)) - 5) #on singlecell LVA counts
            cell.percentchance[0] = (97.81/(1 + numpy.exp((2.15841312e+01-(((time/10)*cell.growthRate)))*6.77630536e-01)) + 2.19) #on livecell data and early RBe counts, percent chance of RB conversion
            #cell.percentchance[0] = (96.45042921/(1 + numpy.exp((13.60222209-(time/10))*1.4553212)) + 1.68390956)#early RBe counts
       
        #pr = RNA production rate
        #dr = RNA degradation rate
        pr0 = 0.04
        pr1 = 0.02
        pr2 = 0.04
        pr3 = 0.06
        dr0 = 0.01
        dr1 = 0.02
        dr2 = 0.01
        dr3 = 0.024        
        
        #pp = protein production rate
        #dp = protein degradation rate
        pp0 = 1.0
        pp1 = 0.5
        pp2 = 1.0
        pp3 = 0.5
        dp0 = 0.05
        dp1 = 0.08
        dp2 = 0.05
        dp3 = 0.01
            
        #flag cells that reach target size for division
        if cell.volume > cell.targetVol:
            cell.divideFlag = True
            
        if cell.cellType == 0: #germinating EB>RB
            cell.divideFlag = False
            cell.growthRate = 0.0
            if time >= cell.germTime[0]: #Become RBr if time is reached
                cell.cellType = 1 #RBr
                cell.growthRate = 1.0 + random.uniform(-0.05,0.05)
                cell.parentGrowth[0] = cell.growthRate #why is this here?
                
        if cell.cellType == 1: #RBr
            
            #cell.rnaamt[0] = cell.rnaamt[0] + (pr0 * cell.growthRate) - (nr0 * cell.rnaamt[0] * cell.growthRate) #extra RNA
            #cell.geneamt[0] = cell.geneamt[0] + (p0 * cell.growthRate * cell.rnaamt[0]) - (n0 * cell.growthRate * cell.geneamt[0]) #extra protein
            cell.geneamt[0] = 0 # Germination could tie it to germ time. ie a count up or down color
            cell.rnaamt[1] = cell.rnaamt[1] + (pr1 * cell.growthRate) - (dr1 * cell.rnaamt[1] * cell.growthRate) #Euo RNA
            cell.geneamt[1] = cell.geneamt[1] + (pp1 * cell.growthRate * cell.rnaamt[1]) - (dp1 * cell.growthRate * cell.geneamt[1]) #Euo
            cell.geneamt[2] = 0 # HctA
            cell.geneamt[3] = 0 # HctB
            cell.color = [[1/cell.geneamt[1], 1, 1/cell.geneamt[1]]]
            #print('growthRate = ' + str(cell.growthRate))
            #print('percentchance = ' + str(cell.percentchance[0]))
            #print('RBe Trigger Value = ' + str(cell.percentchance[0]))
            
            if (time/10).is_integer() and random.uniform(0,100) <= cell.percentchance[0]:
                print('im an RBe')
                cell.cellType = 2 #RBe conversion
                cell.rnaamt[1] = cell.rnaamt[1] + (pr1 * cell.growthRate) - (dr1 * cell.rnaamt[1] * cell.growthRate) #Euo RNA
                cell.geneamt[1] = cell.geneamt[1] + (pp1 * cell.growthRate * cell.rnaamt[1]) - (dp1 * cell.growthRate * cell.geneamt[1]) #Euo
                cell.color = [[1/cell.geneamt[1], 1, 1/cell.geneamt[1]]] # color magic, fix
            #if (time/10) >= 30: #RBr stop divide for CDI model
            #   cell.cellType = 6
            #   cell.growthRate = 0    
            #   cell.geneamt[1] = cell.geneamt[1] - (n1 * cell.geneamt[1]) #Euo
            #   cell.color = [[1/cell.geneamt[1], 1, 1/cell.geneamt[1]]]

        if cell.cellType == 2: #RBe
            cell.rnaamt[1] = cell.rnaamt[1] + (pr1 * cell.growthRate) - (dr1 * cell.rnaamt[1] * cell.growthRate) #Euo RNA
            cell.geneamt[1] = cell.geneamt[1] + (pp1 * cell.growthRate * cell.rnaamt[1]) - (dp1 * cell.growthRate * cell.geneamt[1]) #Euo 
            cell.color = [[1/cell.geneamt[1], 1, 1/cell.geneamt[1]]]
            #if (time/10) >= 30: #RBe stop divide for CDI model
            #    cell.cellType = 6
            #    cell.growthRate = 0
            #    cell.geneamt[1] = cell.geneamt[1] - (n1 * cell.geneamt[1]) #Euo
            #    cell.color = [[1/cell.geneamt[1], 1, 1/cell.geneamt[1]]]
                
        #if cell.cellType == 6: #AB
        #    cell.color = [[1/cell.geneamt[1], 1, 1/cell.geneamt[1]]] 
        #    cell.growthRate = 0
        #    cell.geneamt[1] = cell.geneamt[1] - (n1 * cell.geneamt[1]) #Euo

        if cell.cellType == 3: #IB
            dr1 = 0.08
            dp1 = 0.08
            cell.rnaamt[1] = cell.rnaamt[1] - (dr1 * cell.parentGrowth[0] * cell.rnaamt[1]) # Euo
            cell.geneamt[1] = cell.geneamt[1] - (dp1 * cell.parentGrowth[0] * cell.geneamt[1]) # Euo
            cell.rnaamt[2] = cell.rnaamt[2] + (pr2 * cell.parentGrowth[0])  - (dr2 * cell.rnaamt[2]) #hctA RNA
            cell.geneamt[2] = cell.geneamt[2] + (pp2 * cell.parentGrowth[0] * cell.rnaamt[2]) - (dp2 * cell.parentGrowth[0] * cell.geneamt[2]) #hctA
            print('Euo= '+str(cell.geneamt[1]))
            print('HctA= '+str(cell.geneamt[2]))
            if (time/10) < 1: #theo washout with E-Euo-flag
                cell.geneamt[1] = 1.5 #keeps Euo high to repress ctcB and hctA
            cell.color = [[0, cell.geneamt[1], cell.geneamt[2]]] #blue fast
            #if cell.geneamt[2] >= 3.5: #switch based on HctA amount need to switch over to euo amount or have HctA driven by Euo amount
            if cell.geneamt[1] <= 0.6: #switch based on Euo amount. When it drops IB to pre_EB occors
                cell.cellType = 4 #pre_EB      
       
        if cell.cellType == 4: #pre_EB.  Drive this with HctA->ctcB?
            pr3 = 0.08
            dr3 = 0.024
            pp3 = 0.5
            dp3 = 0.001
            cell.geneamt[2] = cell.geneamt[2] - (dp2 * cell.parentGrowth[0] * cell.geneamt[2]) #hctA 
            cell.rnaamt[3] = cell.rnaamt[3] + (pr3 * cell.parentGrowth[0])  - (dr3 * cell.rnaamt[3]) #hctB RNA
            cell.geneamt[3] = cell.geneamt[3] + (pp3 * cell.parentGrowth[0] * cell.rnaamt[3])  - (dp3 * cell.parentGrowth[0] * cell.geneamt[3]) #hctB
            cell.growthRate = 0
            cell.color = [[cell.geneamt[3]/20, cell.geneamt[1], cell.geneamt[2]]] #blue to black to pink
            if cell.geneamt[3] >= 70: #hctB use this number to delay converstion 
                cell.cellType = 5 #infectious EB
                    
        if cell.cellType == 5:  #infectious EB
            cell.geneamt[2] = cell.geneamt[2] - (dp2 * cell.parentGrowth[0] * cell.geneamt[2]) #hctA 
            cell.rnaamt[3] = cell.rnaamt[3] + (pr3 * cell.parentGrowth[0])  - (dr3 * cell.rnaamt[3]) #hctB RNA
            cell.geneamt[3] = cell.geneamt[3] + (pp3 * cell.parentGrowth[0] * cell.rnaamt[3])  - (dp3 * cell.parentGrowth[0] * cell.geneamt[3]) #hctB
            cell.color = [[cell.geneamt[3], 0.0, cell.geneamt[2]]] #pink
            
def divide(parent, d1, d2):
    # Specify target cell size that triggers cell division
    # Celltype1=RBr, Celltype2=RBe, Celltype3=IB, Celltype4=immature EB, Celltype5=mature EB
    
    if parent.cellType == 1: # If RBr: make 2RBrs
        print('p sp norm = ' + str(parent.species[0]*parent.volume))
        print('p vol = ' + str(parent.volume))
        print('time = '+ str(time))
        d1.cellType = 1
        d1.targetVol = 2
        d1.growthRate = parent.parentGrowth[0] * numpy.random.normal(1, 0.05)
       
        d2.cellType = 1
        d2.targetVol = 2 #* numpy.random.normal(1, 0.05)
        d2.growthRate = parent.parentGrowth[0] * numpy.random.normal(1, 0.05)
        print('d sp norm = ' + str(d2.species[0]*d2.volume))
        print('d vol = ' + str(d2.volume))
        print('time = '+ str(time))
        
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
# Curently gene expression controls cell type and this is good.  
# But would like to make gene expression controlling gene expression 
# [Euo] should inversely control [HctA] and [CtcB] and use [hctA] feedback to turn it off
# [ctcB] should control [hctB] and [hctB] should repress everthing.
# [HctB] threshold should triger infectivity or formation of the True EB.


        
