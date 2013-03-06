# A DRAFT OF A SIMULATION OF BLOCKED ALANINE DIPEPTIDE
# USING THE NEW STRUCTURE
from MDL import *
from numpy import * # to generate random numbers
#from pylab import semilogy,clf # to use matplotlib

def converge(xf,xc,no,i,n):
	#STRONG CONVERGENCE CRITERION TEST - TRAJECTORY DIFFERENCE
	nn = []
	for i in range(i,n+1):
		nn.append(linalg.norm(xf[i].positions-xc[i].positions))
	no.append(nn)

# DEFINE THE NUMBER OF COARSE POINTS
numpoints = 16
# DEFINE THE NUMBER OF ITERATIONS
numiter = 4

finestep = 0.025          # in fs
coarse2fine = 4          # ratio
coarsestep = finestep * coarse2fine #in fs
fixcoarse = 21
fixfine = 10
gammafine = 80
gammacoarse = gammafine
fdof = 0


# PHYSICAL SYSTEMS
io = IO()
normary = [] 
phif = []

xc0 = []
xf = []
dx = []
xc = []
forceary = []
propary = []
# 2 forces - coarse and fine
# 2 propagators - coarse and fine
forceary.append(Forces())
forceary.append(Forces())

io = IO()

for i in range(0, numpoints+1):
    xc0.append(Physical())
    io.readPDBPos(xc0[i], "data/alanDipeptideVac/minC7eq.pdb")
    io.readPSF(xc0[i], "data/alanDipeptideVac/alan_mineq.psf")
    io.readPAR(xc0[i], "data/alanDipeptideVac/par_all27_prot_lipid.inp")
    io.readEigenvectors(xc0[i], "data/alanDipeptideVac/eigVmC7eq")
    xc0[i].bc = "Vacuum"
    xc0[i].temperature = 300
    xc0[i].exclude = "scaled1-4"
    xc0[i].seed = 1234
    xf.append(xc0[i].copy())
    dx.append(xc0[i].copy())
    xc.append(xc0[i].copy())

dof = xc0[0].numAtoms()*3
temp = numpy.ndarray(dof)
tempv = numpy.ndarray(dof)
temp2 =numpy.ndarray(dof)
temp2v =numpy.ndarray(dof)

ff = forceary[0].makeForceField(xc0[0])
ff.bondedForces("badi")
ff.nonbondedForces("le")
ff.params['LennardJones'] = {'algorithm':'Cutoff',
			     'switching':'C2',
			     'cutoff':12,
			     'switchon':9}
ff.params['CoulombDiElec'] = ff.params['LennardJones'].copy()


ff2 = forceary[0].makeForceField(xc0[0])
ff2.bondedForces("badi")
ff2.nonbondedForces("le")
ff2.params['LennardJones'] = {'algorithm':'Cutoff',
			      'switching':'C2',
			      'cutoff':5.5,
			      'switchon':4.5}
ff2.params['CoulombDiElec'] = ff2.params['LennardJones'].copy()


ff3 = forceary[0].makeForceField(xc0[0])
ff3.bondedForces("badi")
ff3.nonbondedForces("le")
ff3.params['LennardJones'] = {'algorithm':'Cutoff',
			      'switching':'C2',
			      'cutoff':5.5,
			      'switchon':4.5}
ff3.params['CoulombDiElec'] = ff3.params['LennardJones'].copy()


finef = forceary[1].makeForceField(dx[0])
finef.bondedForces("badi")
finef.nonbondedForces("le")
finef.params = ff.params.copy()

finef2 = forceary[1].makeForceField(dx[0])
finef2.bondedForces("badi")
finef2.nonbondedForces("le")
finef2.params = ff2.params.copy()

finef3 = forceary[1].makeForceField(dx[0])
finef3.bondedForces("badi")
finef3.nonbondedForces("le")
finef3.params = ff2.params.copy()


io.screen = 1	

#one step of PIT
#DEBUG
#for i in range(0,len(xc0)):
#	print "norm xc0[%d] %f" % (i,linalg.norm(xc0[i].positions))
#	print "norm xc[%d] %f" % (i,linalg.norm(xc[i].positions))
#	print "norm dx[%d] %f" % (i,linalg.norm(dx[i].positions))
#END DEBUG


#io.files = {'gui':('MDL',1)}

#equilibrate
xf[i].seed=int(round(numpy.random.uniform(0,100000)))
forceary[1].forcevec.clear()
propary.append(Propagator(xc0[i], forceary[0], io))
propary.append(Propagator(xf[i], forceary[1], io))



propary[1].propagate(scheme=['NormalModeDiagonalize', 'NormalModeLangevin', 'NormalModeMinimizer'],
                       steps=coarse2fine*100,
                       cyclelength=[1,1],
                       dt=finestep,
                       forcefield=[finef, finef2, finef3],
                       params={'NormalModeDiagonalize':{'reDiagFrequency':100,
                                                        'minSteps':20,
                                                        'minLim':0.1,
                                                        'removeRand':1},
                               'NormalModeLangevin':{'firstmode':1,
                                                     'numbermodes':fixfine,
                                                     'gamma':gammafine,
                                                     'seed':1234,
                                                     'gencompnoise':0,
                                                     'temperature':300},
                               'NormalModeMinimizer':{'minimlim':0.5,
                                                      'rediag':0,
                                                      'randforce':1,
                                                      'simplemin':1,
                                                      'euFactor':0.5}})

for ii in range(0,dof):
	xc0[0].positions[ii]=xf[0].positions[ii]
	xc0[0].velocities[ii]=xf[0].velocities[ii]
	xc[0].positions[ii]=xf[0].positions[ii]
	xc[0].velocities[ii]=xf[0].velocities[ii]

#create numpoints seeds for coarse and fine propagators
seeds=numpy.random.uniform(0,numpoints*2,numpoints)
for i in range(0,len(seeds)):
	seeds[i]=int(round(seeds[i]))
#initial step: line 1 in pseudocode
for i in range(0,numpoints):
#	print "past line 1: step %d" % (i)
#coarse steps saved into xc0
#line 2 in pseudocode
	for ii in range(0,dof):
		temp[ii]=xc0[i].positions[ii]
		tempv[ii]=xc0[i].velocities[ii]
	forceary[0].forcevec.clear()
	propary[0].phys = xc0[i]
	propary[0].propagate(scheme=['NormalModeDiagonalize', 'NormalModeLangevin', 'NormalModeMinimizer'],
                       steps=1,
                       cyclelength=[1,1],
                       dt=coarsestep,
                       forcefield=[ff, ff2, ff3],
                       params={'NormalModeDiagonalize':{'reDiagFrequency':100,
                                                        'minSteps':20,
                                                        'minLim':0.1,
                                                        'removeRand':1},
                               'NormalModeLangevin':{'firstmode':1,
                                                     'numbermodes':fixfine,
                                                     'gamma':gammafine,
                                                     'seed':1234,
                                                     'gencompnoise':0,
                                                     'temperature':300},
                               'NormalModeMinimizer':{'minimlim':0.5,
                                                      'rediag':0,
                                                      'randforce':1,
                                                      'simplemin':1,
                                                      'euFactor':0.5}})

	for ii in range(0, dof):
		xc0[i+1].positions[ii] = xc0[i].positions[ii]
		xc0[i+1].velocities[ii] = xc0[i].velocities[ii]
		xc0[i].positions[ii] = temp[ii]
		xc0[i].velocities[ii] = tempv[ii]
		temp2[ii]=xf[i].positions[ii]
		temp2v[ii]=xf[i].velocities[ii]
		xf[i].positions[ii] = temp[ii]
		xf[i].velocities[ii] = tempv[ii]
#line 5 in pseudocode
	forceary[1].forcevec.clear()
	propary[1].phys = xf[i]
	propary[1].propagate(scheme=['NormalModeDiagonalize', 'NormalModeLangevin', 'NormalModeMinimizer'],
                       steps=coarse2fine*100,
                       cyclelength=[1,1],
                       dt=finestep,
                       forcefield=[finef, finef2, finef3],
                       params={'NormalModeDiagonalize':{'reDiagFrequency':100,
                                                        'minSteps':20,
                                                        'minLim':0.1,
                                                        'removeRand':1},
                               'NormalModeLangevin':{'firstmode':1,
                                                     'numbermodes':fixfine,
                                                     'gamma':gammafine,
                                                     'seed':1234,
                                                     'gencompnoise':0,
                                                     'temperature':300},
                               'NormalModeMinimizer':{'minimlim':0.5,
                                                      'rediag':0,
                                                      'randforce':1,
                                                      'simplemin':1,
                                                      'euFactor':0.5}})



   	for ii in range(0, dof):
		xf[i+1].positions[ii] = xf[i].positions[ii]
		xf[i+1].velocities[ii] = xf[i].velocities[ii]
        	dx[i+1].positions[ii] = xf[i+1].positions[ii]-xc0[i+1].positions[ii]
		dx[i+1].velocities[ii] = xf[i+1].velocities[ii]-xc0[i+1].velocities[ii]
		xf[i].positions[ii] = temp2[ii]
		xf[i].velocities[ii] = temp2v[ii]


#	print "past line 5: step %d" % (i)
#print "finished initialization"
#DEBUG
#for i in range(0,len(xc0)):
#	print "norm xc0[%d] %f" % (i,linalg.norm(xc0[i].positions))
#	print "norm xc[%d] %f" % (i,linalg.norm(xc[i].positions))
#	print "norm dx[%d] %f" % (i,linalg.norm(dx[i].positions))
#END DEBUG
converge(xf,xc0,normary,1,numpoints)
#main loop
for j in range(numiter):
#create numpoints seeds for coarse and fine propagators
#	seeds=numpy.random.uniform(0,numpoints*2,numpoints)
#	for i in range(0,len(seeds)):
#		seeds[i]=int(round(seeds[i]))
	#line 9 in pseudocode
	for i in range(j,numpoints):
#		print "past line 9: step %d" % (i)
		#line 10 in pseudocode
		for ii in range(0,dof):
			temp[ii]=xc[i].positions[ii]
			tempv[ii]=xc[i].velocities[ii]
		forceary[0].forcevec.clear()
		propary[0].phys = xc[i]
		propary[0].propagate(scheme=['NormalModeDiagonalize', 'NormalModeLangevin', 'NormalModeMinimizer'],
                       steps=1,
                       cyclelength=[1,1],
                       dt=coarsestep,
                       forcefield=[ff, ff2, ff3],
                       params={'NormalModeDiagonalize':{'reDiagFrequency':100,
                                                        'minSteps':20,
                                                        'minLim':0.1,
                                                        'removeRand':1},
                               'NormalModeLangevin':{'firstmode':1,
                                                     'numbermodes':fixfine,
                                                     'gamma':gammafine,
                                                     'seed':1234,
                                                     'gencompnoise':0,
                                                     'temperature':300},
                               'NormalModeMinimizer':{'minimlim':0.5,
                                                      'rediag':0,
                                                      'randforce':1,
                                                      'simplemin':1,
                                                      'euFactor':0.5}})



		#propary[0].propagate("NormModeCoarse", xc[i], forceary[0], io, 1, 1, ff, ('fixmodes', fixcoarse), ('gamma', gammacoarse), ('fdof', fdof), coarsestep, ff2, ('avModeMass',3.0), ('seed',seeds[i]),('dtratio',coarse2fine))
		for ii in range(0, dof):
			xc0[i+1].positions[ii] = xc[i].positions[ii]
			xc0[i+1].velocities[ii] = xc[i].velocities[ii]
			xc[i+1].positions[ii] = xc0[i+1].positions[ii]+dx[i+1].positions[ii]
			xc[i+1].velocities[ii] = xc0[i+1].velocities[ii]+dx[i+1].velocities[ii]
			xc[i].positions[ii] = temp[ii]
			xc[i].velocities[ii] = tempv[ii]
			temp2[ii]=xf[i].positions[ii]
			temp2v[ii]=xf[i].velocities[ii]
			xf[i].positions[ii] = temp[ii]
			xf[i].velocities[ii] = tempv[ii]
		forceary[1].forcevec.clear()
		propary[1].phys = xf[i]
		propary[1].propagate(scheme=['NormalModeDiagonalize', 'NormalModeLangevin', 'NormalModeMinimizer'],
                       steps=coarse2fine*100,
                       cyclelength=[1,1],
                       dt=finestep,
                       forcefield=[finef, finef2, finef3],
                       params={'NormalModeDiagonalize':{'reDiagFrequency':100,
                                                        'minSteps':20,
                                                        'minLim':0.1,
                                                        'removeRand':1},
                               'NormalModeLangevin':{'firstmode':1,
                                                     'numbermodes':fixfine,
                                                     'gamma':gammafine,
                                                     'seed':1234,
                                                     'gencompnoise':0,
                                                     'temperature':300},
                               'NormalModeMinimizer':{'minimlim':0.5,
                                                      'rediag':0,
                                                      'randforce':1,
                                                      'simplemin':1,
                                                      'euFactor':0.5}})




		#propary[1].propagate("NormModeCoarse", xf[i], forceary[1], io, coarse2fine, 1, ff, ('fixmodes', fixfine), ('gamma', gammafine), ('fdof', fdof), finestep, ff2, ('avModeMass', 3.0),('seed',seeds[i]),('dtratio',1))
   		for ii in range(0, dof):
			xf[i+1].positions[ii]=xf[i].positions[ii]
			xf[i+1].velocities[ii]=xf[i].velocities[ii]
			dx[i+1].positions[ii]=xf[i+1].positions[ii]-xc0[i+1].positions[ii]
			dx[i+1].velocities[ii]=xf[i+1].velocities[ii]-xc0[i+1].velocities[ii]
			xf[i].positions[ii]=temp2[ii]
			xf[i].velocities[ii]=temp2v[ii]
#		print "past line 14: step %d" % (i)
	converge(xf,xc,normary,1,numpoints)		

		
#analysis of error
#pylab.clf()
x=range(len(normary[0]))
n=normary
#pylab.semilogy(x,n[0],x,n[1],x,n[2])
#pylab.xlabel('trajectory slice')
#pylab.ylabel('norm of position error')
#pylab.legend(('0th it','1st it','2nd it'),loc='lower right')
#pylab.title('Coarse & fine NML DT=2dt GammaC=GammaF/2')

