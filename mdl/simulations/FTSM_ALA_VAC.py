from MDL import *

import FTSM


        
io = IO()
PHI = 11
PSI = 18

# DEFINE THE NUMBER OF POINTS ON THE STRING 
numpoints = 20
#numsteps = 200000
numsteps = 50000

# HELPER FUNCTIONS
def copy(string1):
    retval = []
    for i in range(0, len(string1)):
       retval.append([string1[i][0], string1[i][1]])
    return retval

def rmsd(string1, string2):
    sum = 0
    for i in range(0, len(string1)):        
        sum += (numpy.sqrt((string1[i][0]-string2[i][0])**2+
                           (string1[i][1]-string2[i][1])**2))**2
    result = numpy.sqrt(sum)/numpoints
    return result


# PHYSICAL SYSTEMS
x = []
y = []
force = []
prop = []
for i in range(0, numpoints):
    x.append(Physical())
    y.append(Physical())
    force.append([Forces(), Forces()])
    io.readPDBPos(x[i], "data/alanDipeptideVac/minC7eq.pdb")
    io.readPSF(x[i], "data/alanDipeptideVac/alan_mineq.psf")
    io.readPAR(x[i], "data/alanDipeptideVac/par_all27_prot_lipid.inp")
    x[i].bc = "Vacuum"
    x[i].temperature = 300
    x[i].exclude = "scaled1-4"
    x[i].seed = 1234
    x[i].copy(y[i])

    prop.append([Propagator(x[i], force[i][0], io), Propagator(y[i], force[i][1], io)])


kappa = float(40)
gamma = float(2000)

ff = force[0][0].makeForceField(x[0])

#ff.bondedForces("hh")
ff.bondedForces("badihh")
ff.nonbondedForces("lc")


###################################################################
# STEP 1: OBTAIN THE INITIAL STRING
z = [] # INITIALIZE TO EMPTY
# PROPAGATE ONCE SO OUR CONSTRAINING FORCE CAN GIVE US THE INITIAL
# PHI, PSI
#phiI = physarray[0][0].phi(PHI_DIHEDRAL)
#psiI = physarray[0][0].phi(PSI_DIHEDRAL)
#phiI = -40*numpy.pi/180
#psiI = 130*numpy.pi/180
# TAKEN FROM ALANINE DIPEPTIDE WITH OLD STRING METHOD
phiI = -1.8746023933876219
psiI = 2.4420694033344503

# READ THE SECOND PDB
# PROPAGATE ONCE TO GET THE FINAL PHI,PSI
phiF = 1.2745412016180144
psiF = -1.2217429834113103
#phiF = -40*numpy.pi/180
#psiF = -45*numpy.pi/180

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

ff.params['HarmonicDihedral'] = {'kbias':[kappa, kappa],
                                 'dihedralnum':[PHI-1, PSI-1],
                                 'angle':[z[0][0], z[0][1]]}

#rmsd_vals = []
#rmsdgraph = io.newGraph('step', 'rmsd')
#io.initializePlot('string')
#io.pause=1
#stringgraph=io.newGraph('Phi', 'Psi')
# PRINT INITIAL STRING I0
#print "\nI0: ",
#print z
#io.plotVector(prop[0][0],stringgraph,z, rangex=[-numpy.pi, numpy.pi], rangey=[-numpy.pi, numpy.pi])


dt = 1.0
dphi = 0
dpsi = 0

avg = []
for f in range(0, numpoints):
   avg.append([0,0])

#converged = False
#d = 40000
#initial = copy(z)

convergedPhis = []
convergedPsis = []
matlabfile = open("formatlab.txt", "w")

for iter in range(0, numsteps): # NUMBER OF FTSM ITERATIONS
    for p in range(0, numpoints): # LOOPING OVER POINTS

        if (iter >= 15000):# and iter <= 100000):
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
        # My own function which sets phi and psi for individual force objects
        # Saves performance since I only change 'angle', I don't want to define
        # all new force objects by changing params.
        #if (iter != 0 or p != 0):
        #  FTSM.setConstraint(phi=z[p][0], psi=z[p][1], kappa=kappa, forcefield=ff)
        
        # UPDATE CARTESIAN
        # Dr. Izaguirre: I have checked and this constraint
        # is correct.  The energy is harmonic, but the force (the gradient)
        # is not harmonic.  In fact it is exactly what is in the paper.
        prop[p][0].propagate(scheme="velocityscale", steps=1, dt=dt, forcefield=ff, params={'T0':300})
        prop[p][1].propagate(scheme="velocityscale", steps=1, dt=dt, forcefield=ff, params={'T0':300})

        # My own function which sets phi and psi for individual force objects
        # Saves performance since I only change 'angle', I don't want to define
        # all new force objects by changing params.
        #FTSM.setConstraint(phi=z[(p+1)%numpoints][0], psi=z[(p+1)%numpoints][1], kappa=kappa, forcefield=ff)

    z = FTSM.reparam(z)


    #if (iter % 10 == 0):
    #   rmsd_vals.append([iter, rmsd(initial,z)])
    
    #print rmsd(initial, z)
    if (iter >= 15000):
       # Write to the file of points to plot in Matlab
       for bb in range(0, len(z)):
          convergedPhis.append(z[bb][0])
          convergedPsis.append(z[bb][1])
        
    ### AVERAGING CODE ###########################
    #if (iter >= 110000):
    if (iter >= 15000):
       for g in range(0, numpoints):
          for h in range(0, 2):
             avg[g][h] += z[g][h]
    #print "\nI"+str(iter+1)+": ",z
    #if (iter+1 == 200000):
    if (iter+1 == 50000):
       for g in range(0, numpoints):
          for h in range(0, 2):
             #avg[g][h] /= 90000.
             avg[g][h] /= 35000.
       print "\nAVG: ", avg
    ##############################################

    if (iter % 100 == 0):
       print "\nI"+str(iter+1)+": ",z
    #if (iter % 10 == 0):
    #   io.plotVector(prop[0][0],rmsdgraph,rmsd_vals)
    #io.plotVector(prop[0][0],stringgraph,z, rangex=[-numpy.pi, numpy.pi], rangey=[-numpy.pi, numpy.pi])

matlabfile.write(str(convergedPhis))
matlabfile.write('\n')
matlabfile.write(str(convergedPsis))
