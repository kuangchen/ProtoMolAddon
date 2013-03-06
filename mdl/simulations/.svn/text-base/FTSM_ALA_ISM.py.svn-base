from MDL import *

import FTSM

io = IO()
PHI_DIHEDRAL = 11
PSI_DIHEDRAL = 18

# DEFINE THE NUMBER OF POINTS ON THE STRING 
numpoints = 8


# PHYSICAL SYSTEMS
physarray = []
forcearray = []
proparray = []
for i in range(0, numpoints):
    physarray.append([Physical(), Physical()])
    forcearray.append([Forces(), Forces()])
    
    for j in range(0, 2):
        io.readPDBPos(physarray[i][j], "data/alanSolStates/alanC7axial_wb5_min_eq.pdb")
        io.readPSF(physarray[i][j], "data/alanSolStates/alan_wb5.psf")
        io.readPAR(physarray[i][j], "data/diAlanine/par_all27_prot_lipid.inp")
        physarray[i][j].bc = "Vacuum"
        physarray[i][j].temperature = 300
        physarray[i][j].exclude = "scaled1-4"
        physarray[i][j].seed = 1234
        
    proparray.append([Propagator(physarray[i][0],forcearray[i][0],io), Propagator(physarray[i][1],forcearray[i][1],io)])
    
#phys = Physical()
#io.readPDBPos(phys, "examples/alanSolStates/alanC7axial_wb5_min_eq.pdb")
#io.readPSF(phys, "examples/alanSolStates/alan_wb5.psf")
#io.readPAR(phys, "examples/diAlanine/par_all27_prot_lipid.inp")
#phys.bc = "Vacuum"
#phys.temperature = 300
#phys.exclude = "scaled1-4"
#phys.seed = 1234

#forces = Forces()
ff = forcearray[0][0].makeForceField(physarray[0][0])
ff.bondedForces("badihh")
#print ff.forcearray
ff.nonbondedForces("le")
#ff.nonbondedForces("lc")
#ff.setAlgorithm("LennardJonesCoulomb", "SimpleFull")

ff.params['LennardJones'] = {'algorithm':'Cutoff',
                             'switching':'C2',
                             'cutoff':20.0}

ff.params['CoulombDiElec'] = {'algorithm':'Cutoff',
                              'switching':'C2',
                              'cutoff':20.0}


#prop = Propagator()


# SECOND PHYSICAL SYSTEM
#phys2 = Physical()
#io.readPDBPos(phys2, "examples/alanSolStates/alanC7axial_wb5_min_eq.pdb")
#io.readPSF(phys2, "examples/alanSolStates/alan_wb5.psf")
#io.readPAR(phys2, "examples/diAlanine/par_all27_prot_lipid.inp")
#phys2.bc = "Vacuum"
#phys2.temperature = 300
#phys2.exclude = "scaled1-4"
#phys2.seed = 1234

#forces2 = Forces()
#prop2 = Propagator()



###################################################################
# STEP 1: OBTAIN THE INITIAL STRING
S = [] # INITIALIZE TO EMPTY
# PROPAGATE ONCE SO OUR CONSTRAINING FORCE CAN GIVE US THE INITIAL
# PHI, PSI
#phiI = physarray[0][0].phi(PHI_DIHEDRAL)
#psiI = physarray[0][0].phi(PSI_DIHEDRAL)
phiI = -40*numpy.pi/180
psiI = 130*numpy.pi/180
print phiI, " ", psiI

# READ THE SECOND PDB
# PROPAGATE ONCE TO GET THE FINAL PHI,PSI
#io.readPDBPos(physarray[0][0], "examples/alanSolStates/alanC7equat_wb5_min_eq.pdb")
#phiF = physarray[0][0].phi(PHI_DIHEDRAL)
#psiF = physarray[0][0].phi(PSI_DIHEDRAL)
phiF = -40*numpy.pi/180
psiF = -45*numpy.pi/180
print phiF, " ", psiF

# SEPARATE INTO EQUIDISTANT PHI AND PSI
dphi = (phiF - phiI) / (numpoints-1)
dpsi = (psiF - psiI) / (numpoints-1)
S.append([phiI, psiI])
newphi = phiI
newpsi = psiI
for ii in range(1, numpoints-1):
    newphi += dphi
    newpsi += dpsi
    S.append([newphi, newpsi])
S.append([phiF, psiF])


kappa = 40.0
gamma = 20000.0

ff.params['HarmonicDihedral'] = {'kbias':[kappa, kappa],
                                 'dihedralnum':[PHI_DIHEDRAL-1, PSI_DIHEDRAL-1],
                                 'angle':[S[0][0], S[0][1]]}

#io.pause=1
stringgraph=io.newGraph('Phi', 'Psi')

# PRINT INITIAL STRING I0
#print "\nI0: ",
#print S
#io.plotVector(proparray[0][0],stringgraph,S, rangex=[-numpy.pi, numpy.pi], rangey=[-numpy.pi, numpy.pi])

# SET THE STATE BACK TO THE ORIGINAL PDB
#io.readPDBPos(physarray[0][0], "examples/alanSolStates/alanC7axial_wb5_min_eq.pdb")




print "\nI0: ",
print S
io.plotVector(proparray[0][0],stringgraph,S, rangex=[-numpy.pi, 0], rangey=[-100*numpy.pi/180, numpy.pi])
# SET THE STATE BACK TO THE ORIGINAL PDB
#io.readPDBPos(physarray[0][0], "examples/alanSolStates/alanC7axial_wb5_min_eq.pdb")



dt = 1.0
for iter in range(0, 100000): # NUMBER OF FTSM ITERATIONS
    for workpt in range(0, numpoints): # LOOPING OVER POINTS
        if (iter >= 10000 and iter <= 100000):
            kappa += (100.-40.)/90000.
        # UPDATE FREE SPACE
        # USE FIRST SYSTEM TO GET M
        # USE SECOND SYSTEM TO OBTAIN PHI AND PSI DIFFERENCES
        # FROM TARGETS
        M = FTSM.M(physarray[workpt][0], PHI_DIHEDRAL, PSI_DIHEDRAL)
        #print "DIFFERENCE: ", S[workpt][0]-physarray[workpt][1].phi(PHI_DIHEDRAL), " " , S[workpt][1]-physarray[workpt][1].phi(PSI_DIHEDRAL)
        S[workpt][0] -= (kappa/gamma)*M*dt*(S[workpt][0]-physarray[workpt][1].phi(PHI_DIHEDRAL))
        S[workpt][1] -= (kappa/gamma)*M*dt*(S[workpt][1]-physarray[workpt][1].phi(PSI_DIHEDRAL))
        # UPDATE CARTESIAN
        # Dr. Izaguirre: I have checked and this constraint
        # is correct.  The energy is harmonic, but the force (the gradient)
        # is not harmonic.  In fact it is exactly what is in the paper.

        proparray[workpt][0].propagate(scheme="LangevinImpulse", steps=1, dt=dt, forcefield=ff)

        proparray[workpt][1].propagate(scheme="LangevinImpulse", steps=1, dt=dt, forcefield=ff)

        # My own function which sets phi and psi for individual force objects
        # Saves performance since I only change 'angle', I don't want to define
        # all new force objects by changing params.
        if (workpt == numpoints-1):
           FTSM.setPhiPsi(ff, S[0][0], S[0][1])
        else:
           FTSM.setPhiPsi(ff, S[workpt+1][0], S[workpt+1][1])


    #print "\nI"+str(iter+1)+" BEFORE SMOOTH: ",S
    # REPARAMETERIZE S
    #S.reverse()
    #FTSM.switchPhiPsi(S)
    S = FTSM.arclenChris(FTSM.extractPhi(S), FTSM.extractPsi(S))
    #FTSM.switchPhiPsi(S)
    #S.reverse()
    print "\nI"+str(iter+1)+": ",S
    io.plotVector(proparray[0][0],stringgraph,S, rangex=[-numpy.pi, 0], rangey=[-100*numpy.pi/180, numpy.pi])

