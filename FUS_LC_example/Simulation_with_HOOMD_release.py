#######################################################################################################################
# Python script to construct and run coexistence simulations of intrinsically disordered proteins using the HPS model[1] 
# [1] Sequence determinants of protein phase behavior from a coarse-grained model 
#     GL Dignon, W Zheng, YC Kim, RB Best, J Mittal - PLoS computational biology, 2018
#
# Authors: Roshan M Regy, Wenwei Zheng, Gregory L. Dignon, Jeetain Mittal 
#######################################################################################################################
# Prerequisities: This code was tested using HOOMD v2.9.0 with azplugins installed and GSD v2.2.0 
#
# Usage: python Simulation_with_HOOMD.py YOUR_FASTA_PROTEIN_SEQUENCE_FILE 
#
# This code has the following steps:
#
# (1) Create a single chain initial configuration as a linear polymer chain --> generates start.gsd 
#
# (2) Replicate the single chain initial configuration using the replicate() function in HOOMD to 
#   user defined (nx*ny*nz) system size and then reduce the dimensions of the replicated configuration 
#   to form a cubic box with user defined dimensions using hoomd.update.box_resize() while 
#   running an NVT simulation at user defined temperature and time step  --> generates resize.gsd, resize.dcd, resize.log
#
# (3) Unwrap molecules over the periodic boundary in z and extend the box in the z direction to create a slab configuration 
#   --> generates box2slab_extend.gsd 
#
# (4) Run a production NVT simulation at user defined temperature and time step and for user defined number of steps
#   --> generates production.dcd, restart_tmp1.gsd, restart_tmp2.gsd, production.log, stress.log  
#######################################################################################################################
import sys,os,numpy as np
import hoomd, hoomd.md as md
import hoomd.deprecated as old
from hoomd import azplugins
import gsd, gsd.hoomd, gsd.pygsd 

# SYSTEM SIZE
nx=5 # Replicate system in x direction 
ny=5 # Replicate system in y direction 
nz=4 # Replicate system in z direction --> Total system size would be nx*ny*nz
# Please keep magnitudes of nx, ny and nz as close to each other 
# as possible to avoid large replicated system box sizes in a specific direction 

# SLAB DIMENSIONS
boxsize=15.0 # The x and y dimensions of the slab configuration in nanometers
slab_z_length = 280.0 # The z dimension of the slab configuration in nanometers

# RESIZING RUN PARAMETERS
resize_T=300 # Temperature for resizing run in Kelvin
resize_steps=500000 # Number of steps (>500) to be used for resizing the box to boxsize set previously
if resize_steps<500:
    resize_steps=501
resize_dt=0.01 # Time step in picoseconds for box resizing 

# PRODUCTION RUN PARAMETERS
production_dt=0.01 # Time step for production run in picoseconds
production_steps=100000000 # Total number of steps 
production_T=300 # Temperature for production run in Kelvin



seq={'R':'ARG','H':'HIS','K':'LYS','D':'ASP','E':'GLU',
     'S':'SER','T':'THR','N':'ASN','Q':'GLN','C':'CYS',
     'U':'SEC','G':'GLY','P':'PRO','A':'ALA','V':'VAL',
     'I':'ILE','L':'LEU','M':'MET','F':'PHE','Y':'TYR',
     'W':'TRP'}

# ##### Read one letter amino acid sequence from file
filein = sys.argv[1]
fileout='%s_seq3.dat'%(filein)
nline=1
count=0
fout=open(fileout,'w')
with open(filein,'r') as fid:
    for i in fid:
        if i[0]!='#':
            for j in i:
                if j in seq:
                    fout.write(' %s'%seq[j])
                    count+=1
                    if count==nline:
                        fout.write('\n')
                        count=0
fout.close()


# #### 1.2 Read sequence and force field parameters
# ##### Input parameters for all the amino acids (force field)
ff_para = 'stats_module.dat'
aalist={}
with open(ff_para,'r') as fid:
    for i in fid:
        if i[0]!='#':
            tmp=i.rsplit()
            aalist[tmp[0]]=np.loadtxt(tmp[1:],dtype=float)
aakeys=list(aalist.keys())
# This translates each amino acid type into a number, which will be used in HOOMD
# For example, GLY is with an ID of 10
aamass=[]
aacharge=[]
aaradius=[]
aahps=[]
print ('aakeys')
print (aakeys)
print ('aalist[i][1]')
print (aalist[aakeys[1]][1])
for i in aakeys:
    aamass.append(aalist[i][0])
    aacharge.append(aalist[i][1])
    aaradius.append(aalist[i][2])
    aahps.append(aalist[i][3])
# Now we can translate the entire sequence into a number code according to the order in 'aakeys'
chain_id=[]
chain_mass=[]
chain_charge=[]
with open(fileout,'r') as fid:
    for i in fid:
        iname=i.rsplit()[0]
        chain_id.append(aakeys.index(iname))
        chain_mass.append(aalist[iname][0])
        chain_charge.append(aalist[iname][1])

bond_length=0.38
chain_length=len(chain_id)
box_length=bond_length*chain_length+10
# #### 1.3 Now we can build HOOMD data structure for one single frame
s=gsd.hoomd.Snapshot()
s.particles.N = chain_length
s.particles.types = aakeys
s.particles.typeid = chain_id
s.particles.mass = chain_mass
s.particles.charge = chain_charge

#  Build initial position as a linear chain
pos=[]
for i in range(len(chain_id)):
    # Change the z-coordinate to have a linear chain
    pos.append((0,0,(i-int(len(chain_id)/2))*bond_length))
pos=np.array(pos)
s.particles.position= pos

# Initialize bond
nbonds=chain_length-1
s.bonds.N = nbonds
s.bonds.types = ['AA_bond']
s.bonds.typeid = [0]*(nbonds)
bond_pairs=np.zeros((nbonds,2),dtype=int)
for i in range(nbonds):
    bond_pairs[i,:]=np.array([i,i+1])
s.bonds.group = bond_pairs

# Box size
s.configuration.dimensions=3
s.configuration.box=[box_length,box_length,box_length,0,0,0]
s.configuration.step=0

# #### 1.4 Write intial singe chain gsd file
f = gsd.hoomd.open(name='start.gsd', mode='wb')
f.append(s)
f.close()
#################################################################################################
# ## start.gsd contains one single chain of the given protein 
# ----------------------------------------------------------------------------------------------
# ## 2.0 Replicate single chain system to given number of chains for a slab 
################################################################################################
hoomd.context.initialize()
system = hoomd.init.read_gsd('start.gsd')


# ### Replicate the single chain here. Remember total number of chains = nx*ny*nz
system.replicate(nx=nx,ny=ny,nz=nz) 

# #### Bonds
harmonic=hoomd.md.bond.harmonic()
harmonic.bond_coeff.set('AA_bond',k=8368,r0=bond_length)

# #### Neighborlist and exclusions
nl = hoomd.md.nlist.cell()
nl.reset_exclusions(exclusions=['1-2', 'body'])

# #### Pairwise interactions
#nb = hoomd.md.pair.lj_lambda(r_cut=0, nlist=nl) 
nb = azplugins.pair.ashbaugh(r_cut=0, nlist=nl)
for i in aakeys:
    for j in aakeys:
        nb.pair_coeff.set(i,j,lam=(aalist[i][3]+aalist[j][3])/2.,
                          epsilon=0.8368, sigma=(aalist[i][2]+aalist[j][2])/10./2.,r_cut=2.0)    
        
# #### Electrostatics
yukawa = hoomd.md.pair.yukawa(r_cut=0.0, nlist=nl)
for i,atom1 in enumerate(aakeys):
    for j,atom2 in enumerate(aakeys):
        yukawa.pair_coeff.set(atom1,atom2,epsilon=aalist[atom1][1]*aalist[atom2][1]*1.73136, kappa=1.0, r_cut=3.5) 
# #### Group Particles
all = hoomd.group.all()

# #### Set up integrator
hoomd.md.integrate.mode_standard(dt=resize_dt) # Time units in ps
kTinput=resize_T * 8.3144598/1000.
integrator = hoomd.md.integrate.langevin(group=all, kT=kTinput, seed=63535)
# #### Resize the box after replication to 15x15x15nm  
hoomd.update.box_resize(L=hoomd.variant.linear_interp([(0,system.box.Lx),(resize_steps-500,boxsize)]),scale_particles=True)
for cnt,i in enumerate(aakeys):
    integrator.set_gamma(i,gamma=aamass[cnt]/1000.0)

# #### Output log file with box dimensions and restart file after box resizing
hoomd.analyze.log(filename='resize.log', quantities=['potential_energy','kinetic_energy','temperature','pressure_xx','pressure_yy','pressure_zz','lx','ly','lz'], period=10000, overwrite=True, header_prefix='#')
hoomd.dump.gsd('resize.gsd', period=10000, group=all, truncate=True)
hoomd.dump.dcd('resize.dcd', period=10000, group=all, overwrite=True)
# #### Run resizing simulation
hoomd.run(tsteps=resize_steps)
#################################################################################################
# ## resize.gsd contains replicated system with 100 chains and box size of 15x15x15 nm 
# ----------------------------------------------------------------------------------------------
# ## 3.0 Extend 15nm cube to 15x15x280nm slab 
################################################################################################
def extend(s):
    boxdim = s.configuration.box[:3]
    zmin,zmax,dz = -boxdim[2]/2., boxdim[2]/2., boxdim[2]
    pos1 =  s.particles.position
    pos = pos1.copy()
    skip=0
    ncomp=1
    for k in range(ncomp):
        nchain = int(s.particles.N/chain_length)
        nres = chain_length
        for i in range(nchain):
            mol_coord = pos[i*nres+skip:(i+1)*nres+skip,2]
            for j in range(1,nres):
                dist2 = (mol_coord[j] - mol_coord[j-1])**2
                if dist2 > 8:
                    excess = np.sign(mol_coord[j] - mol_coord[j-1])*dz
                    mol_coord[j] = mol_coord[j] - excess 
                com = np.mean(mol_coord)
                if com < zmin:
                    mol_coord += dz
                elif com > zmax:
                    mol_coord -= dz
            pos[i*nres+skip:(i+1)*nres+skip,2] = mol_coord
        skip += nchain*nres
    return pos

f = gsd.pygsd.GSDFile(open('resize.gsd','rb'))
t = gsd.hoomd.HOOMDTrajectory(f)
s1 = t[0]
s = gsd.hoomd.Snapshot()
s.particles.N = s1.particles.N
s.particles.types = s1.particles.types 
s.particles.typeid = s1.particles.typeid 
s.particles.mass = s1.particles.mass
s.particles.charge = s1.particles.charge
s.particles.position = extend(s1)
s.bonds.N = s1.bonds.N
s.bonds.types = s1.bonds.types
s.bonds.typeid = s1.bonds.typeid
s.bonds.group = s1.bonds.group
s.configuration.box = s1.configuration.box
s.configuration.dimensions=3
s.configuration.box = [s1.configuration.box[0],s1.configuration.box[1],slab_z_length,0,0,0] 
s.configuration.step = 0
outfile = gsd.hoomd.open('box2slab_extend.gsd','wb')
outfile.append(s)
outfile.close()
#################################################################################################
# ### Minimized slab formed and saved in minimize.gsd
#---------------------------------------------------------------------------------------------
# ## 4.0. Run a production slab simulation using minimize.gsd from previous step
################################################################################################
hoomd.context.initialize()
system = hoomd.init.read_gsd('box2slab_extend.gsd')

n_steps = production_steps # 1 microseconds

fileroot = 'Production'
nl = hoomd.md.nlist.cell()

## Bonds
harmonic = hoomd.md.bond.harmonic()
harmonic.bond_coeff.set('AA_bond', k=8360, r0=0.381)
## Nonbonded
nl.reset_exclusions(exclusions=['1-2', 'body'])
nb = azplugins.pair.ashbaugh(r_cut=0, nlist=nl)
for i in aakeys:
    for j in aakeys:
        nb.pair_coeff.set(i,j,lam=(aalist[i][3]+aalist[j][3])/2.,
                          epsilon=0.8368, sigma=(aalist[i][2]+aalist[j][2])/10./2.,r_cut=2.0)    

## Electrostatics
yukawa = hoomd.md.pair.yukawa(r_cut=0.0, nlist=nl)
for i,atom1 in enumerate(aakeys):
    for j,atom2 in enumerate(aakeys):
        yukawa.pair_coeff.set(atom1,atom2,epsilon=aalist[atom1][1]*aalist[atom2][1]*1.73136, kappa=1.0, r_cut=3.5) 

## Group Particles
all = hoomd.group.all()

## Set up integrator
hoomd.md.integrate.mode_standard(dt=production_dt) # Time units in ps
temp = production_T*0.00831446
integrator = hoomd.md.integrate.langevin(group=all, kT=temp, seed=399991) # Temp is kT/0.00831446
for cnt,i in enumerate(aakeys):
    integrator.set_gamma(i,gamma=aamass[cnt]/1000.0)
## Outputs
hoomd.analyze.log(filename=fileroot+'.log', quantities=['potential_energy', 'pressure_xx', 'pressure_yy', 'pressure_zz', 'temperature','lx','ly','lz'], period=100000, overwrite=False, header_prefix='#')
hoomd.analyze.log(filename='stress.log', quantities=['pressure_xy', 'pressure_xz', 'pressure_yz'], period=100000, overwrite=False, header_prefix='#') # Output stress tensor
hoomd.dump.gsd('restart_tmp1.gsd', period=1000000, group=all, truncate=True)
hoomd.dump.gsd('restart_tmp2.gsd', period=1000000, group=all, truncate=True, phase=500000)
hoomd.dump.dcd(fileroot+'_dump.dcd', period=100000, group=all, overwrite=False)

## Run simulation
hoomd.run_upto(production_steps, limit_hours=48)
########################################################################################################
