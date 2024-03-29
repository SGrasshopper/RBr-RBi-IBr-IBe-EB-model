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

max_cells = 1500

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
    cell.parentGrowth = [0] #progenitor cell logged growthRate. "Metabolism rate" as these cells don't realy grow
    
    cell.color = [2.0, 0.5, 1.5]#specify color of cell
    #cell.color = [1.0,0.0,0.0] #red
    #cell.color = [0.0,1.0,0.0] #green
    #cell.color = [0.0,0.0,1.0] #blue
    
    #RNA and protein 
    cell.rnaamt = [0.0,0.0,0.0,0.0,0.0] # RNA levels, used, in part, to drive geneamt levels #[0]= ectExp, [1]=Euo, [2]=HctA, [3]=CtcB, [4]=HctB
    cell.geneamt = [0.0,0.0,0.0,0.0,0.0]   #[0]= ectExp, [1]=Euo, [2]=HctA, [3]=CtcB, [4]=HctB
    
    #EB to RB germination time
    cell.germTime = [(120 + random.uniform(-40,40))] #based on livecell and single cell expansion data: need to measure actually germ time variation and fit to dist

    #RBr > RBi conversion percent
    cell.percentchance = [0,0] #curve that drives RBr > RBi conversion

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
   
def sigRateCL(): #Add

    return '''
    const float k1 = 1.0f;
    float x0 = signals[0];
    rates[0] = k1;
    ''' 
    
time = 0
def update(cells): #Iterate through each cell update and flag cells that reach target size for division
    global time
    #global n0 #whats this and why global?
    time += 1
    time2 = (time/10) 
    print('time = ' + str(time))
    print('time2 = ' + str(time2)) 
    
    # Celltypes: 0=germ_EB, 1=RBr, 2=RBi, 3=IBr, 4=IBe, 5=EB

    for (id, cell) in cells.items():
    
        if time >= cell.germTime[0]:
            #cell.percentchance[0] = (105/(1 + numpy.exp((1.76663094e+01-time2)*3.77251916e-01)) - 5) #on singlecell LVA counts
            cell.percentchance[0] = (97.81/(1 + numpy.exp((2.15841312e+01-((time2*cell.growthRate)))*6.77630536e-01)) + 2.19) #on livecell data and early RBe counts, percent chance of RB conversion
            #cell.percentchance[0] = (96.45042921/(1 + numpy.exp((13.60222209-time2)*1.4553212)) + 1.68390956)#early RBi counts
        #print('growthRate = ' + str(cell.growthRate))
        #print('percentchance = ' + str(cell.percentchance[0]))

        if cell.volume > cell.targetVol:
            cell.divideFlag = True
            
        if cell.cellType == 0: #germinating EB>RB
            cell.divideFlag = False
            cell.growthRate = 0.0
            if time >= cell.germTime[0]: #Become RBr if time is reached
                cell.cellType = 1 #RBr
                cell.growthRate = 1.0 + random.uniform(-0.05,0.05)
                
        if cell.cellType == 1: #RBr
            pr0 = 0.0 #RNA production rate of ectExp
            nr0 = 0.0 #RNA degredation rate of ectExp
            p0 = 0.0  #Protein production rate of ectExp
            n0 = 0.0  #Protein degredation rate of ectExp
            
            pr1 = 0.02 #RNA production rate of Euo
            nr1 = 0.02 #RNA degredation rate of Euo
            p1 = 0.5  #protein production rate of Euo
            n1 = 0.08 #protein degredation rate of Euo
        	
        	#I guess you are using RNA->protein to do your math? pP*pRNA-nP*nRNA?
            cell.rnaamt[0] = cell.rnaamt[0] + (pr0 * cell.growthRate) - (nr0 * cell.rnaamt[0]) # RNA of ectopic expressed protein
            cell.geneamt[0] = cell.geneamt[0] + (p0 * cell.rnaamt[0] * cell.growthRate) - (n0 * cell.geneamt[0]) #Ectopic expresed protein
            
            cell.rnaamt[1] = cell.rnaamt[1] + (pr1 * cell.growthRate) - (nr1 * cell.rnaamt[1]) #Euo RNA
            cell.geneamt[1] = cell.geneamt[1] + (p1 * cell.rnaamt[1] * cell.growthRate) - (n1 * cell.geneamt[1]) #Euo protein
            
            cell.color = [[1/cell.geneamt[1], 1, 1/cell.geneamt[1]]]
            #print('growthRate = ' + str(cell.growthRate))
            #print('percentchance = ' + str(cell.percentchance[0]))
            print('RBi Trigger Value = ' + str(cell.percentchance[1]))
            print('time2 = ' + str(time2))
            
            if time2.is_integer() and random.uniform(0,100) <= cell.percentchance[0]: # time2.is_integer() and  RBr to RBe. species[0] = magic Rbr>RBe signal
                #print('time2 is int' + str(time2))
                print('im an RBi')
                cell.cellType = 2 #RBi conversion
                cell.color = [[1/cell.geneamt[1], 1, 1/cell.geneamt[1]]] # color magic, fix
        
        if  cell.cellType == 2: #RBi 
            pr0 = 0.0 #RNA production rate of ectExp
            nr0 = 0.0 #RNA degredation rate of ectExp
            p0 = 0.0  #Protein production rate of ectExp
            n0 = 0.0  #Protein degredation rate of ectExp
            
            pr1 = 0.02 #RNA production rate of Euo
            nr1 = 0.02 #RNA degredation rate of Euo
            p1 = 0.5  #protein production rate of Euo
            n1 = 0.08 #protein degredation rate of Euo
            
            cell.rnaamt[0] = cell.rnaamt[0] + (pr0 * cell.growthRate) - (nr0 * cell.rnaamt[0]) # RNA of ectopic expressed protein
            cell.geneamt[0] = cell.geneamt[0] + (p0 * cell.rnaamt[0] * cell.growthRate) - (n0 * cell.geneamt[0]) #Ectopic expresed protein
            
            cell.rnaamt[1] = cell.rnaamt[1] + (pr1 * cell.growthRate) - (nr1 * cell.rnaamt[1]) #Euo RNA
            cell.geneamt[1] = cell.geneamt[1] + (p1 * cell.rnaamt[1] * cell.growthRate) - (n1 * cell.geneamt[1]) #Euo protein 
            
        if  cell.cellType == 3: #IBr   
            #need to add Euo inhibition of HctA. High Euo blocks HctA production. Inherit high Euo and have it degrade fast. Turn on HctA at low levels of Euo.
            pr0 = 0.0 #RNA production rate of ectExp
            nr0 = 0.0 #RNA degredation rate of ectExp
            p0 = 0.0  #Protein production rate of ectExp
            n0 = 0.0  #Protein degredation rate of ectExp
            
            pr1 = 0.02 #RNA production rate of Euo
            nr1 = 0.02 #RNA degredation rate of Euo
            p1 = 0.5  #protein production rate of Euo
            n1 = 0.08 #protein degredation rate of Euo
            
            pr2 = 0.02 #RNA production rate of HctA
            nr2 = 0.01 #RNA degredation rate of HctA
            p2 = 0.5  #protein production rate of HctA
            n2 = 0.05 #protein degredation rate of HctA
            
            pr3 = 0.06 #RNA production rate of CtcB
            nr3 = 0.024 #RNA degredation rate of CtcB
            p3 = 0.5  #protein production rate of CtcB
            n3 = 0.01 #protein degredation rate of CtcB
            
            if time < 100000:
                pr0 = 0.0 #RNA production rate of ectExp
                nr0 = 0.0 #RNA degredation rate of ectExp
                p0 = 0.0  #Protein production rate of ectExp
                n0 = 0.0  #Protein degredation rate of ectExp
            
            if cell.geneamt[0] > 10:
                pr1 = 0.02 #RNA production rate of Euo
                nr1 = 0.02 #RNA degredation rate of Euo
                p1 = 0.5  #protein production rate of Euo
                n1 = 0.08 #protein degredation rate of Euo
            
            cell.rnaamt[0] = cell.rnaamt[0] + (pr0 * cell.parentGrowth[0]) - (nr0 * cell.rnaamt[0]) # RNA of ectopic expressed protein
            cell.geneamt[0] = cell.geneamt[0] + (p0 * cell.rnaamt[0] * cell.parentGrowth[0]) - (n0 * cell.geneamt[0]) #Ectopic expresed protein

            cell.rnaamt[1] = cell.rnaamt[1] + (pr1 * cell.parentGrowth[0]) - (nr1 * cell.rnaamt[1]) #Euo RNA
            cell.geneamt[1] = cell.geneamt[1] + (p1 * cell.rnaamt[1] * cell.parentGrowth[0]) - (n1 * cell.geneamt[1]) #Euo protein 
            
            if cell.geneamt[1] > 3: # High Euo blocks expression of HctA and CtcB
                pr2 = 0.0 #RNA production rate of HctA
                nr2 = 0.0 #RNA degredation rate of HctA
                p2 = 0.0  #protein production rate of HctA
                n2 = 0.0 #protein degredation rate of HctA
            
                pr3 = 0.0 #RNA production rate of CtcB
                nr3 = 0.0 #RNA degredation rate of CtcB
                p3 = 0.0  #protein production rate of CtcB
                n3 = 0.0 #protein degredation rate of CtcB
            
            cell.rnaamt[2] = cell.rnaamt[2] + (pr2 * cell.parentGrowth[0]) - (nr2 * cell.rnaamt[2]) #hctA RNA
            cell.geneamt[2] = cell.geneamt[2] + (p2 * cell.rnaamt[2] * cell.parentGrowth[0]) - (n2 * cell.geneamt[2]) #hctA link this to Euo levels
            
            cell.rnaamt[3] = cell.rnaamt[3] + (pr3 * cell.parentGrowth[0]) - (nr3 * cell.rnaamt[3]) #ctcB RNA
            cell.geneamt[3] = cell.geneamt[3] + (p3 * cell.rnaamt[3] * cell.parentGrowth[0]) - (n3 * cell.geneamt[3] ) #CtcB protein
    
            cell.color = [[0, 0, cell.geneamt[2]/100]] #need to fix color
            if cell.geneamt[3] >= 4: 
                cell.cellType = 4 #IBe
       
        if cell.cellType == 4: #IBe
            pr2 = 0.02 #RNA production rate of HctA
            nr2 = 0.01 #RNA degredation rate of HctA
            p2 = 0.5  #protein production rate of HctA
            n2 = 0.05 #protein degredation rate of HctA
            
            pr3 = 0.06 #RNA production rate of CtcB
            nr3 = 0.024 #RNA degredation rate of CtcB
            p3 = 0.5  #protein production rate of CtcB
            n3 = 0.01 #protein degredation rate of CtcB
            
            pr4 = 0.06 #RNA production rate of HctB      
            nr4 = 0.024 #RNA degredation rate of HctB  
            p4 = 0.5  #protein production rate of HctB           
            n4 = 0.01 #protein degredation rate of HctB
            
            cell.geneamt[2] = cell.geneamt[2] - (n2 * cell.geneamt[2]) #HctA protein deg
            
            cell.geneamt[3] = cell.geneamt[3] - (n3 * cell.geneamt[3]) #CtcB protein deg
            
            cell.rnaamt[4] = cell.rnaamt[4] + (pr4 * cell.parentGrowth[0]) - (nr4 * cell.rnaamt[4]) #hctB RNA
            cell.geneamt[4] = cell.geneamt[4] + (p4 * cell.rnaamt[4] * cell.parentGrowth[0]) - (n4 * cell.geneamt[4]) #HctB protein
            
            cell.color = [[0, 0, cell.geneamt[4]/100]] #need to fix color
            if cell.geneamt[4] >= 15: #hctB protein
                cell.cellType = 5 #EB
                    
        if cell.cellType == 5: #EB
            pr4 = 0.06 #RNA production rate of HctB      
            nr4 = 0.024 #RNA degredation rate of HctB  
            p4 = 0.5  #protein production rate of HctB           
            n4 = 0.01 #protein degredation rate of HctB
            
            cell.rnaamt[4] = cell.rnaamt[4] + (pr4 * cell.parentGrowth[0]) - (nr4 * cell.rnaamt[4]) #hctB RNA
            cell.geneamt[4] = cell.geneamt[4] + (p4 * cell.rnaamt[4] * cell.parentGrowth[0]) - (n4 * cell.geneamt[4]) #hctB protein
            
            cell.color = [2.0, 0.0, 0.5]
            
def divide(parent, d1, d2):
    # Specify target cell size that triggers cell division
    # Celltype1=RBr, Celltype2=RBi, Celltype3=IBr, Celltype4=IBe, Celltype5=mature EB
    
    if parent.cellType == 1: # If RBr: make 2RBrs
        d1.cellType = 1
        d1.targetVol = 2 #* numpy.random.normal(1, 0.05)
        d1.growthRate = parent.growthRate * numpy.random.normal(1, 0.05)
        
        d2.cellType = 1
        d2.targetVol = 2 #* numpy.random.normal(1, 0.05)
        d2.growthRate = parent.growthRate * numpy.random.normal(1, 0.05)
        
        d1.geneamt[0] = parent.geneamt[0]/2
        d2.geneamt[0] = parent.geneamt[0]/2
        d1.geneamt[1] = parent.geneamt[1]/2
        d2.geneamt[1] = parent.geneamt[1]/2 

    if parent.cellType == 2: # If RBi: make 1RBi, 1IBe
        d1.cellType = 2
        d1.targetVol = 2 #* numpy.random.normal(1, 0.05)
        d1.growthRate = parent.growthRate * numpy.random.normal(1, 0.05)
        
        d2.cellType = 3
        d2.growthRate = 0
        d2.parentGrowth[0] = parent.growthRate
        d2.targetVol = 10
        
        d1.geneamt[0] = parent.geneamt[0]/2
        d2.geneamt[0] = parent.geneamt[0]/2
        d1.geneamt[1] = parent.geneamt[1]/2
        d2.geneamt[1] = parent.geneamt[1]/2 
        
    #print('d1.percentchance[1] =' + str(d1.percentchance[1]))
    #print('d2.percentchance[1] =' + str(d2.percentchance[1]))         

#this model: GermEB lavendar, Rbr green, Rbe green, IB black, EB  hot pink

# RBr matures into RBi based on percentchance curve from empirical data

# RBi divides -> RBi and IBr
# IBr -> IBe is controled by the degredation of Euo which then induces HctA and CtcB. HctA represses Euo while CtcB induces HctB and the sigma 54 regulon
# IBe -> EB is tigered by high levels of HctB and the sigma 54 regulon
