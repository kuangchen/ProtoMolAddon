from MDL import *
 
 
import sys
import os
 
def setPhiPsi(phi,psi,kp,forcefield):
    flag = False
    for i in range(0, forcefield.forcetypes.__len__()):
        if (forcefield.forcetypes[i] == 'h'):
            if (not flag):
               forcefield.forcearray[i].setReferenceDihedral(phi)
               flag = True
            else:
               forcefield.forcearray[i].setReferenceDihedral(psi)
            forcefield.forcearray[i].setKappa(kp)
 
 
 
numgridsquares = int(sys.argv[1])

xlim=([-numpy.pi,numpy.pi])
ylim=([-numpy.pi,numpy.pi])
 
print xlim[0]
print xlim[1]
print ylim
 
midpts = numpy.zeros((numgridsquares,numgridsquares,2))

print midpts[1,1]
print "TREV"
xdiff = (xlim[1]-xlim[0])/numgridsquares
 
print xdiff
 
for i in range(0,numgridsquares):
    tmpx = xlim[0]+xdiff*(i+0.5)
    for j in range(0,numgridsquares):
        tmpy = ylim[0]+xdiff*(j+0.5)
        print tmpx," ",tmpy
        midpts[i][j][0] = tmpx
        midpts[i][j][1] = tmpy

#sys.exit(1)

#print midpts[0][2][0]," ",midpts[0][2][1]


myPhys=Physical()
myForces=Forces()
myIO=IO()
myIO.readPDBPos(myPhys,"data/alanDipeptideVac/minC7eq.pdb")
myIO.readPSF(myPhys,"data/alanDipeptideVac/alan_mineq.psf")
myIO.readPAR(myPhys,"data/alanDipeptideVac/par_all27_prot_lipid.inp")
myPhys.bc="Vacuum"
myPhys.temperature=300
myPhys.exclude="scaled1-4"
myPhys.seed=1234
myProp = Propagator(myPhys,myForces,myIO)

#print "BFA"
ff = myForces.makeForceField(myPhys)
ff.bondedForces("badihh")
ff.nonbondedForces("le")

kappa=20.0
kappa1=float(sys.argv[2])
PHI=11
PSI=18

outdir=sys.argv[3]

ff.params['HarmonicDihedral'] = {'kbias':[kappa, kappa],
                                 'dihedralnum':[PHI-1, PSI-1],
                                 'angle':[midpts[0][0][0],midpts[0][0][1]]}
import FTSM
#print midpts
ff.params['LennardJones'] = {'algorithm':'Cutoff',
                             'cutoff':12,
                             'switchon':10,
                             'switching':'C2'}

ff.params['CoulombDiElec'] = {'algorithm':'Cutoff',
                              'switching':'Cn',
                              'switchon':10,
                              'cutoff':12,
                              'order':2}


for ii in range(0,numgridsquares):
    for jj in range(0,numgridsquares):
        print "BOX: ", ii, jj
        if (ii !=0):
            FTSM.setConstraint(PHI, PSI, midpts[ii][jj][0],midpts[ii][jj][1],kappa,ff)
        elif (jj != 0) :
            FTSM.setConstraint(PHI, PSI, midpts[ii][jj][0],midpts[ii][jj][1],kappa,ff)

        pts = ii*numgridsquares + jj

        for i in range(0,20):
            myProp.propagate(scheme="LangevinImpulse",steps=100,dt=1.0,forcefield=ff,params={'temp':300,'gamma':91})
            
        FTSM.setConstraint(PHI, PSI, midpts[ii][jj][0],midpts[ii][jj][1],kappa1,ff)

        #Equilibrate at this point for a bit
        myProp.propagate(scheme="LangevinImpulse",steps=50,dt=1.0,forcefield=ff,params={'temp':300,'gamma':91})
        
        ef = outdir+'/sampleatrc.energies.'+str(ii)+'.'+str(jj)
        df = outdir+'/sampleatrc.'+str(ii)+'.'+str(jj)+'.dcd'
        print ef," ",df
        myIO.reset()
        myIO.files={'energies':(ef,10),'dcdtrajpos':(df,10)}
        myIO.build()
        myProp.propagate(scheme="LangevinImpulse",steps=50000,dt=1.0,forcefield=ff,params={'temp':300,'gamma':91})
        sys.stdout.flush()
