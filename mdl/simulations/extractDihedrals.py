from MDL import *

import sys
myPhys=Physical()
myIO=IO()
myIO.readPDBPos(myPhys,"data/alanDipeptideVac/minC7eq.pdb")
myIO.readPSF(myPhys,"data/alanDipeptideVac/alan_mineq.psf")
myIO.readPAR(myPhys,"data/alanDipeptideVac/par_all27_prot_lipid.inp")
myPhys.bc="Vacuum"
myPhys.temperature=300
myPhys.exclude="scaled1-4"
myPhys.seed=1234

PHI = 11
PSI = 18

import numpy

prefix = sys.argv[1]
OUTFILE = open(prefix+'/allsample.dihed.30.protomol', 'w')

for ii in range (0, 30):
    for jj in range (0, 30):
        print "***** BOX ",ii,jj," *****"
        dcdfile = prefix+'/sampleatrc.'+str(ii)+'.'+str(jj)+'.dcd'
        for kk in range(0, 5000):
           myIO.readDCDTrajectoryPos(myPhys, dcdfile)
           phi = myPhys.angle(PHI)
           psi = myPhys.angle(PSI)
           OUTFILE.writelines('\t'+str(phi)+'\t'+str(psi)+'\n')
            
