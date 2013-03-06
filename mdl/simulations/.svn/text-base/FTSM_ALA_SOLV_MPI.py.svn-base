from MDL import *

import FTSM
import mpi
import math

io = IO()
PHI = 11
PSI = 18

# DEFINE THE NUMBER OF POINTS ON THE STRING 
numpoints = 8


# PHYSICAL SYSTEM
x = Physical()
y = Physical()
force = []

io.readPDBPos(x, "data/alanDipeptideSol/solvate_eq.pdb")
io.readPSF(x, "data/alanDipeptideSol/solvate.psf")
io.readPAR(x, "data/alanDipeptideSol/par_all27_prot_lipid.inp")
x.bc = "Periodic"
x.temperature = 300
x.exclude = "scaled1-4"
x.seed = 1234
temperature = 300
#temperature = x.getTemperature()
y = x.copy()

force = [Forces(), Forces()]
prop = [Propagator(x, force[0], io), Propagator(y, force[1], io)]

kappa = 40.0
gamma = 2000.0


ff = force[0].makeForceField(x)
ff2 = force[0].makeForceField(x)
ff.bondedForces("badihh")
ff.nonbondedForces("lc")

ff.params['LennardJones'] = {'algorithm':'Cutoff',
                             'switching':'Cutoff',
                             'cutoff':9}

ff.params['Coulomb'] = {'algorithm':'Cutoff',
                        'switching':'Cutoff',
                        'cutoff':9}


###################################################################
# STEP 1: OBTAIN THE INITIAL STRING
if (mpi.rank == 0):
   z = [] # INITIALIZE TO EMPTY
   avg = []
   for f in range(0, 8):
       avg.append([0,0])
   # PROPAGATE ONCE SO OUR CONSTRAINING FORCE CAN GIVE US THE INITIAL
   # PHI, PSI
   phiI = -40*numpy.pi/180
   psiI = 130*numpy.pi/180
   #phiI = -110*numpy.pi/180
   #psiI = 135*numpy.pi/180

   # READ THE SECOND PDB
   # PROPAGATE ONCE TO GET THE FINAL PHI,PSI

   phiF = -40*numpy.pi/180
   psiF = -45*numpy.pi/180
   #phiF = -105*numpy.pi/180
   #psiF = -35*numpy.pi/180

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
   print "\nI0: ",
   print z
else:
   z = []

#if (mpi.rank == 0):
#   z = mpi.bcast(z)
#else:
#   z = mpi.bcast()
#   
#z_p = z[mpi.rank]
z_p = mpi.scatter(z)[0]
print z_p
import sys
sys.exit(0)
ff.params['HarmonicDihedral'] = {'kbias':[kappa, kappa],
                                 'dihedralnum':[PHI-1, PSI-1],
                                 'angle':[z_p[0], z_p[1]]}


dt = 1.0
kappaincr = (1000.-40.)/100000.
for iter in range(0, 2000): # NUMBER OF FTSM ITERATIONS
    #for workpt in range(0, numpoints): # LOOPING OVER POINTS
        if (iter >= 10001 and iter <= 110000):
            kappa += kappaincr


        if (iter != 0):
           FTSM.setConstraint(PHI, PSI, phi=z_p[0], psi=z_p[1], kappa=kappa, forcefield=ff)       
        # UPDATE FREE SPACE
        # USE FIRST SYSTEM TO GET M
        # USE SECOND SYSTEM TO OBTAIN PHI AND PSI DIFFERENCES
        # FROM TARGETS
        zp0 = z_p[0]
	z_p[0] -= (kappa/gamma)*dt*(FTSM.M(x, PHI, PHI)*(z_p[0]-y.angle(PHI)) + FTSM.M(x, PHI, PSI)*(z_p[1] - y.angle(PSI)))
        z_p[1] -= (kappa/gamma)*dt*(FTSM.M(x, PSI, PHI)*(zp0-y.angle(PHI)) + FTSM.M(x, PSI, PSI)*(z_p[1] - y.angle(PSI)))
        
        # UPDATE CARTESIAN
        # Dr. Izaguirre: I have checked and this constraint
        # is correct.  The energy is harmonic, but the force (the gradient)
        # is not harmonic.  In fact it is exactly what is in the paper.        
        prop[0].propagate(scheme="velocityscale", steps=1, dt=dt, forcefield=ff, params={'T0':300})
        prop[1].propagate(scheme="velocityscale", steps=1, dt=dt, forcefield=ff, params={'T0':300})

        mpi.barrier()
        z = mpi.gather([z_p])

        if (mpi.rank == 0):
            z = FTSM.reparamTrevor(z)
            if (iter >= 110000):
                for g in range(0, 8):
                    for h in range(0, 2):
                        avg[g][h] += z[g][h]
            print "\nI"+str(iter+1)+": ",z
            if (iter+1 == 200000):
                for g in range(0, 8):
                    for h in range(0, 2):
                        avg[g][h] /= 90000.
                print "\nAVG: ", avg
                    


        z_p = mpi.scatter(z)[0]
        mpi.barrier()


