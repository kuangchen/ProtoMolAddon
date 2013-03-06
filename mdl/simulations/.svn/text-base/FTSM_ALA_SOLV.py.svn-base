from MDL import *

import FTSM

io = IO()
PHI = 11
PSI = 18

# DEFINE THE NUMBER OF POINTS ON THE STRING 
numpoints = 8
numsteps = 2000

# PHYSICAL SYSTEMS
x = []
y = []
force = []
prop = []
for i in range(0, numpoints):
    x.append(Physical())
    y.append(Physical())
    force.append([Forces(), Forces()])
    io.readPDBPos(x[i], "data/alanDipeptideSol/solvate_eq.pdb")
    io.readPSF(x[i], "data/alanDipeptideSol/solvate.psf")
    io.readPAR(x[i], "data/alanDipeptideSol/par_all27_prot_lipid.inp")
    x[i].bc = "Periodic"
    x[i].temperature = 300
    x[i].exclude = "scaled1-4"
    x[i].seed = 1234
    x[i].copy(y[i])

    prop.append([Propagator(x[i], force[i][0], io), Propagator(y[i], force[i][1], io)])


kappa = float(40)
gamma = float(2000)

ff = force[0][0].makeForceField(x[0])

ff.bondedForces("badihh")
ff.nonbondedForces("lc")


###################################################################
# STEP 1: OBTAIN THE INITIAL STRING
z = [] # INITIALIZE TO EMPTY
# PROPAGATE ONCE SO OUR CONSTRAINING FORCE CAN GIVE US THE INITIAL
# PHI, PSI
#phiI = physarray[0][0].phi(PHI_DIHEDRAL)
#psiI = physarray[0][0].phi(PSI_DIHEDRAL)
phiI = -40*numpy.pi/180
psiI = 130*numpy.pi/180
print phiI, " ", psiI

# READ THE SECOND PDB
# PROPAGATE ONCE TO GET THE FINAL PHI,PSI
phiF = -40*numpy.pi/180
psiF = -45*numpy.pi/180
print phiF, " ", psiF

# SEPARATE INTO EQUIDISTANT PHI AND PSI
dphi = (phiF - phiI) / (numpoints-1)
dpsi = (psiF - psiI) / (numpoints-1)
z.append([phiI, psiI])
newphi = phiI
newpsi = psiI
for ii in range(1, numpoints-1):
    newphi += dphi
    newpsi += dpsi
    z.append([newphi, newpsi])
z.append([phiF, psiF])



ff.params['LennardJones'] = {'algorithm':'Cutoff',
                             'switching':'Cutoff',
                             'cutoff':9}

ff.params['Coulomb'] = {'algorithm':'Cutoff',
                        'switching':'Cutoff',
                        'cutoff':9}

ff.params['HarmonicDihedral'] = {'kbias':[kappa, kappa],
                                 'dihedralnum':[PHI-1, PSI-1],
                                 'angle':[z[0][0], z[0][1]]}

#io.initializePlot('string')
#io.pause=1
#stringgraph=io.newGraph('Phi', 'Psi')

# PRINT INITIAL STRING I0
print "\nI0: ",
print z
#io.plotVector(prop[0][0],stringgraph,z, rangex=[-numpy.pi, numpy.pi], rangey=[-numpy.pi, numpy.pi])


dt = 1.0
for iter in range(0, numsteps): # NUMBER OF FTSM ITERATIONS
    for p in range(0, numpoints): # LOOPING OVER POINTS
        if (iter >= 10000 and iter <= 100000):
            kappa += (100.-40.)/90000.

        if (p != 0 or iter != 0):
           FTSM.setConstraint(PHI, PSI, phi=z[p][0], psi=z[p][1], kappa=kappa, forcefield=ff)
        
        # UPDATE FREE SPACE
        # USE FIRST SYSTEM TO GET M
        # USE SECOND SYSTEM TO OBTAIN PHI AND PSI DIFFERENCES
        # FROM TARGETS
        zp0 = z[p][0]
        z[p][0] -= (kappa/gamma)*dt*(FTSM.M(x[p], PHI, PHI)*(z[p][0]-y[p].dihedral(PHI)) + FTSM.M(x[p], PHI, PSI)*(z[p][1] - y[p].dihedral(PSI)))
        z[p][1] -= (kappa/gamma)*dt*(FTSM.M(x[p], PSI, PHI)*(zp0-y[p].dihedral(PHI)) + FTSM.M(x[p], PSI, PSI)*(z[p][1] - y[p].dihedral(PSI)))
        
        
        # UPDATE CARTESIAN
        prop[p][0].propagate(scheme="velocityscale", steps=1, dt=dt, forcefield=ff, params={'T0':300})

        prop[p][1].propagate(scheme="velocityscale", steps=1, dt=dt, forcefield=ff, params={'T0':300})


        # My own function which sets phi and psi for individual force objects
        # Saves performance since I only change 'angle', I don't want to define
        # all new force objects by changing params.
        #FTSM.setConstraint(phi=z[(p+1)%numpoints][0], psi=z[(p+1)%numpoints][1], kappa=kappa, forcefield=ff)

    z = FTSM.reparam(z)

    print "\nI"+str(iter+1)+": ",z
    #io.plotVector(prop[0][0],stringgraph,z, rangex=[-numpy.pi, numpy.pi], rangey=[-numpy.pi, numpy.pi])

