from __future__ import print_function, division
from espressomd import System, interactions, lb, polymer, lbboundaries, shapes
from espressomd.observables import ComPosition
from espressomd.accumulators import Correlator

from numpy import savetxt, zeros
import numpy as np
import sys
import espressomd.lb

	
import os

# AVB: Create an output directory for this to store the output files
outdir = "./Noelle/output-SampleThree"
if not os.path.exists(outdir):
    os.makedirs(outdir)


# Setup constant
time_step = 0.01
loops = 50
step_per_loop = 100

# AVB: the parameters (that I usually use)
a = 0.05
r0 = 2.0*a
kBT = 4.0e-6
vwf_type = 0
monomer_mass = 0.01

# System setup
system = System(box_l = [32.0, 32.0, 32.0])
system.set_random_state_PRNG()
np.random.seed(seed = system.seed)
system.cell_system.skin = 0.4

mpc = 20 # The number of monomers has been set to be 20 as default
         # Change this value for further simulations



# Fene interaction
fene = interactions.FeneBond(k=0.04, d_r_max=0.3)
system.bonded_inter.add(fene)


# Setup polymer of part_id 0 with fene bond
# AVB: Notice the mode, max_tries and shield parameters for pruned self-avoiding random walk algorithm 
polymer.create_polymer(N_P=1, MPC=mpc, bond=fene, bond_length=r0,
                       start_pos = [16.0, 16.0, 16.0], mode=2, max_tries=100, shield=0.6*r0)


# AVB: setting the type of particles and changing mass of each monomer to 0.01
system.part[:].type = vwf_type
system.part[:].mass = monomer_mass

# AVB: I suggest to add Lennard-Jones interaction between the monomers
# AVB: to reproduce hydrophobicity
# AVB: parameters for the potential (amplitude and cut-off redius)
amplVwfVwf = 4.0*kBT    # sometimes we change this to 2.0*kBT  
rcutVwfVwf = 2.0*r0
# AVB: the potential
system.non_bonded_inter[vwf_type,vwf_type].lennard_jones.set_params(
      epsilon = amplVwfVwf, sigma = r0/1.122,
      shift = "auto", cutoff = rcutVwfVwf, min = r0*0.6)

print("Warming up the polymer chain.")
## For longer chains (>100) an extensive 
## warmup is neccessary ...
system.time_step = 0.002
system.thermostat.set_langevin(kT=4.0e-6, gamma=1.0)
# AVB: Here the Langevin thermostat is needed, because we have not yet initialized the LB-fluid.
# AVB: And somehow it is necessary so that the polymer adopts the equilibrium conformation of the globule.
# AVB: you may skip this step

for i in range(100):
    system.force_cap = float(i) + 1
    system.integrator.run(100)

print("Warmup finished.")
system.force_cap = 0
system.integrator.run(100)
system.time_step = time_step
system.integrator.run(500)

# AVB: the following command turns the Langevin thermostat on in line 49
system.thermostat.turn_off()

# AVB: This command sets the velocities of all particles to zero
system.part[:].v = [0,0,0]

# AVB: The density was too small here. I have set 1.0 for now.
# AVB: It would be necessary to recalculate, but the density of the liquid should not affect the movements of the polymer (this is how our physical model works).
lbf = espressomd.lb.LBFluid(agrid=1, dens=1.0, visc=1.0e2, tau=time_step, fric=0.01)
system.actors.add(lbf)
system.thermostat.set_lb(kT= 4.0e-6)

print("Warming up the system with LB fluid.")
system.integrator.run(5000)
print("LB fluid warming finished.")
# AVB: after this you should have a completely collapsed polymer globule
# AVB: If you want to watch the process of globule formation in Paraview, just change 5000 to 0 in line 100


# configure correlators
com_pos = ComPosition(ids=(0,))
c = Correlator(obs1 = com_pos, tau_lin=16, tau_max=loops*step_per_loop, delta_N=1,
        corr_operation="square_distance_componentwise", compress1="discard1")
system.auto_update_accumulators.add(c)

print("Sampling started.")
for i in range(loops):
    system.integrator.run(step_per_loop)
    system.analysis.append()
    lbf.print_vtk_velocity(outdir+"/fluid%04i.vtk" %i)
    system.part.writevtk(outdir+"/vwf_all%04i.vtk" %i) 
    sys.stdout.write("\rSampling: %05i"%i)
    sys.stdout.flush()

c.finalize()
corrdata = c.result()
corr = zeros((corrdata.shape[0],2))
corr[:,0] = corrdata[:,0]
corr[:,1] = (corrdata[:,2] + corrdata[:,3] + corrdata[:,4]) / 3

savetxt(outdir+"/msd_nom"+str(mpc)+".dat", corr)

with open(outdir+"/rh_out.dat","a") as datafile:
    rh = system.analysis.calc_rh(chain_start=0, number_of_chains=1, chain_length=mpc-1)
    datafile.write(str(mpc)+ "    " + str(rh[0])+"\n")
