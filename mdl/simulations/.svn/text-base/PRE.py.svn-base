# Simulation of alanine dipeptide in vacuum
# Using Replica Exchange

from MDL import *
import Constants
import numpy as np
import numpy.random as nr
#import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import time
from mpi4py import MPI
import sys
import pprint
import cPickle
import re

#globald = []
X=[];Y=[]
debug = True

XYZ = True
DCD = True

XYZFile = "ww_e"
DCDFile = "ww_e"

#Global Files
#pdbFile = "data/alanDipeptideBlock/alanC7axial.pdb"
#psfFile = "data/alanDipeptideBlock/blockdialanine.psf"
#parFile = "data/alanDipeptideBlock/par_all27_prot_lipid.inp"
pdbFile = "data/wwd/ww_exteq_nowater1.pdb"
psfFile = "data/wwd/ww_exteq_nowater1.psf"
parFile = "data/wwd/par_all27_prot_lipid.inp"
computeTime = 0.0
ioTime = 0.0
TMP = 300
if (len(sys.argv) > 1):
    workingDir = sys.argv[1]
  #  print workingDir
else:
    workingDir = "."
if not os.path.exists(workingDir):
   # print workingDir
    os.makedirs(workingDir)
MonteCarloSteps = 10
BoundaryConditions = "Vacuum"
minTemp = 300.0
maxTemp = 340.0
MDsteps = 1000
OutputFreq = 100



def metropolis( u_i, u_j, t_i, t_j ):
    # Metropolis for replicas i with potential energy u_i
    #                         temperature t_i
    #                    and  j with potential energy u_j
    #                         temperature t_j
    K_b = Constants.boltzmann()
    deltaE = (1 / (K_b * t_i) - 1/ (K_b * t_j) ) - (u_j - u_i)
    if( deltaE < 0 ):
        return True
    acceptProb = np.exp( -deltaE )
    randNum = nr.random()
    if( randNum < acceptProb ):
        return True
    else:
        return False

def replica_exchange(rank, phys, minT, maxT, nreplicas, mdsteps, mcsteps, outfreq):
    """
    replica exchange method
    molecular system:                     phys
    propagator object:                    prop
    forces to be computed (forcefield):   ff
    temperature range:                    between minT and maxT
    number of replicas:                   nreplicas
    number of md steps between exchanges: mdsteps
    number of exchange attempts:          mcsteps
    """

    accProb = 0.0
    pp = pprint.PrettyPrinter(indent=4)
    if rank == 0:
        allTemps = np.zeros((mcsteps,nreplicas))
        # Assign temperatures to replicas
        incT = (maxT - minT) / (nreplicas-1)
        temperatures = np.arange(minT, maxT+.1, incT)
        if debug: pp.pprint( temperatures)
    else:
        temperatures = None


#    print "processor %d says size of system is %d\n " % (rank, phys.myPDB.size())
    myReplica = Replica(phys,mdsteps,outputfreq=outfreq)

    for j in range(mcsteps):
        if debug and rank==0: print "in MC iteration %d\n" % j

        if rank==0:
            myTemp = np.zeros(1)
            allTemps[j] = temperatures
            myTemp[0] = temperatures[0]
            for i in range(1,nreplicas):
                comm.Send(temperatures[i], dest = i, tag=13)
        else:
            myTemp = np.zeros(1)
            comm.Recv(myTemp, source=0, tag=13)

        u = np.zeros(1)
        #    if debug:
        #        print( "\nprocessor %d says path is \n " % rank )
        #        pp.pprint (sys.path)
        print myTemp[0]
        print rank
        print j
        u = myReplica.compute(T=myTemp[0],myRank=rank,step=j)
        #    if debug: print "processor %d says potential is %f\n " % (rank, u)
        uresult = comm.gather(u, 0)
        if debug and rank==0: print "processor %d says potentials are %s\n " % (rank, uresult)
        if rank == 0:
            minR = mlab.find(temperatures==minT)
            if debug:print "processor %d says lowest temp replica is %d\n" % (rank,minR)        

        #  Attempt exchange of k and k+1 or k and k-1
        if rank==0:
            print "URESULT:"
            print uresult
            k = nr.randint(nreplicas)
            if debug: print "selected replica %d\n " % k
            if (k==nreplicas-1):
                if debug: print "attempt exchange of replicas %d and %d\n" % (k,k-1)
                if metropolis(uresult[k], uresult[k-1], temperatures[k], temperatures[k-1]):
                    if debug: print "exchanged replicas %d and %d \n " % (k,k-1)
                    accProb = accProb + 1.0
                    temp=temperatures[k]
                    temperatures[k]=temperatures[k-1]
                    temperatures[k-1]=temp
            else:
                if debug: print "attemp exchange of  replicas %d and %d \n " % (k,k+1)
                if metropolis(uresult[k], uresult[k+1], temperatures[k], temperatures[k+1]):
                    if debug: print "exchanged replicas %d and %d \n " % (k,k+1)
                    accProb = accProb + 1.0
                    temp=temperatures[k]
                    temperatures[k]=temperatures[k+1]
                    temperatures[k+1]=temp
            if debug: print temperatures
    if rank == 0: return (accProb / mcsteps, allTemps)
    else: return (accProb/mcsteps)
        
class Replica(object):
    def __init__(self,phys, mdsteps, outputfreq=10):
        self.phys = phys
        self.forces = Forces()
        self.ff = self.forces.makeForceField(phys,"charmm")
        self.io = IO()
        self.mdsteps = mdsteps
        self.outfreq = outputfreq # output interval

    def compute(self,T,myRank,step):
        """ perform MD move at given temperature """
        self.phys.temperature = T
        self.io.reset()
	self.io.build()
        ####Here is DCD IO File Specification#####
        if DCD:
            dcdFname = "%s/%s.%d-%d.dcd" % (workingDir,DCDFile,myRank,step)
            self.io.files['dcdtrajpos']=(dcdFname,self.outfreq)
            print "Printing DCD %s" % dcdFname
        if XYZ:
            xyzFname = "%s/%s.%d-%d.xyz" % (workingDir,XYZFile,myRank,step)
            self.io.files['xyztrajpos']=(xyzFname,self.outfreq)
            if debug: print "Printing XYZ %s" % xyzFname
      

        self.prop = Propagator(self.phys, self.forces, self.io)

        nsamples = self.mdsteps / self.outfreq
        for i in range(nsamples):
            self.prop.propagate(scheme="LangevinLeapfrog", \
                                steps=self.outfreq, \
                                dt=2.0, \
                                forcefield=self.ff, \
                                params={'temp':self.phys.temperature,'gamma':10})

        u = self.forces.energies.potentialEnergy(self.phys)
        if XYZ: io.writeXYZPos(phys,xyzFname)
  #      return (u, d)
        return u
        
    
if __name__ == '__main__':
    comm = MPI.COMM_WORLD
    size = MPI.COMM_WORLD.Get_size()
    name = MPI.Get_processor_name()
    rank = MPI.COMM_WORLD.Get_rank()
    computeTime = 0.0
    ioTime = 0.0
    # PHYSICAL SYSTEM
    # Everyone reads the same file
    # We use afs distributed file system
    phys = Physical()
    io = IO()
    io.readPDBPos(phys, pdbFile)
    io.readPSF(phys, psfFile)
    io.readPAR(phys, parFile)
    phys.bc = BoundaryConditions
    phys.temperature = TMP
    phys.exclude = "scaled1-4"

    finalTemps = np.zeros((MonteCarloSteps,size))
    
    if rank==0:
        init_time = time.time()
        (accProb,finalTemps) = replica_exchange(rank, phys, minT=minTemp, maxT=maxTemp, nreplicas=size, \
                           mdsteps=MDsteps, mcsteps=MonteCarloSteps, outfreq=OutputFreq)

        fin_time = time.time()
        computeTime += (fin_time - init_time)
#        print finalTemps
    else:
        init_time = time.time()
        (accProb) = replica_exchange(rank, phys, minT=minTemp, maxT=maxTemp, nreplicas=size, \
                               mdsteps=MDsteps, mcsteps=MonteCarloSteps, outfreq=OutputFreq)
        
        fin_time = time.time()
        computeTime += (fin_time - init_time)

    if rank==0: 
        print "Acceptance Probability: %f" % accProb 

#        if debug: print "dihedrals are %s\n" % globald

#        for i,d in enumerate(globald):
#            if i%2 == 0: X.append( d )
#            else: Y.append( d )
#        outputFile = "%s/rem-%dmd-%dmc.dat" % (workingDir,MDsteps,MonteCarloSteps)
#        with open(outputFile,'w') as fd:
#            cPickle.dump(globald, fd)

        #plt.hexbin(X, Y, bins = 'log')
        # if debug: plt.show()
        patFile = re.compile(r'.*(\d+)\-(\d+).*xyz')
        cmd = "ls -l %s | grep ^- | awk '{print $9}'" % workingDir
        files = []
        init_time = time.time()
        files = os.popen(cmd).read().split()
        for f in files:
            fGroups = patFile.match(f)
            if fGroups:
                frank = fGroups.group(1)
       #         print frank
                fstep = fGroups.group(2)
       #         print fstep
                ftemp = finalTemps[fstep,frank]
                fname = "%s/%s/%s_step%s-replica%s.xyz" % (workingDir,ftemp,XYZFile,fstep,frank)
                d = os.path.dirname(fname)
                if not os.path.exists(d):
                    os.makedirs(d)
                # Move each file to the directory for its corresponding temperature
                cmd2 = "mv %s/%s %s" % (workingDir,f,fname)
        #        os.system(cmd2)
        fin_time = time.time()
        ioTime += fin_time - init_time
        print "IO Time: %f" % ioTime
        print "Compute Time: %f" % computeTime
    
        
    
