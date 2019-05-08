from __future__ import print_function, division
from espressomd import System, interactions, lb, polymer, lbboundaries, shapes
from espressomd.observables import ComPosition
from espressomd.accumulators import Correlator

from numpy import savetxt, zeros
import numpy as np
import sys
import espressomd.lb

	
import os


def calc(var, loopsnumber):

    # AVB: Create an output directory for this to store the output files
    outdir = "./Noelle/out-adsorpVWF_ww/SR=" + str(var)
    
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    
    # Setup constant
    time_step = 0.01
    loops = loopsnumber
    step_per_loop = 100
    
    # AVB: the parameters (that I usually use)
    fscale= 1e-9
    tscale= 1e-4
    lscale= 1e-6
 
    # per um^2
    collagen_density = 3

    a = 0.05
    r0 = 2.0*a
    kBT = 4.0e-6
    vwf_type = 0
    collagen_type = 1
    wall_type = 2
    monomer_mass = 0.01
    N_vwf = 1    

    box_l = 8.0
    bot_wall_dist = 1.0
    top_wall_dist = box_l-1.0
    print("Shear velocity:"+ str(var) + "/s")
    shear_rate = var
    vy = box_l * (shear_rate * tscale)
    #vy = var
    v = [0, vy, 0]
    
    # System setup
    
    system=0
    
    system = System(box_l = [box_l, box_l, box_l])
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
    polymer.create_polymer(N_P=1, MPC=mpc, bond=fene, bond_length=r0, start_pos = [bot_wall_dist+2*a, box_l*0.5, box_l*0.5], mode=2, max_tries=100, shield=0.6*r0)
    
    
    next_part_id = mpc * N_vwf  
    print(" -- next particle id after VWF creation: ")
    print(next_part_id)  

    # AVB: setting the type of particles and changing mass of each monomer to 0.01
    system.part[:].type = vwf_type
    system.part[:].mass = monomer_mass
    
    # AVB: Shifting the polymer closer to the wall
    print(system.part[:].pos)
    px=[ ]
    for p in system.part[:].pos:
        px.append(p[0])
    minx = np.min(px)
    print("closest x coord = %f" % minx)
    print("\nmoving polymer to closer to the wall")
    for pid in range(next_part_id):
        p = system.part[pid].pos
        posx = p[0] - minx +bot_wall_dist
        posy = p[1]
        posz = p[2]
        print([posx, posy, posz])
        system.part[pid].pos = [posx, posy, posz]
        

    # AVB: I suggest to add Lennard-Jones interaction between the monomers
    # AVB: to reproduce hydrophobicity
    # AVB: parameters for the potential (amplitude and cut-off redius)
    amplVwfVwf = 4.0*kBT    # sometimes we change this to 2.0*kBT  
    rcutVwfVwf = 1.5*r0
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
    

    # Setup boundaries
    walls = [lbboundaries.LBBoundary() for k in range(2)]
    walls[0].set_params(shape=shapes.Wall(normal=[1,0,0], dist = bot_wall_dist), velocity=[0,0,0])
    walls[1].set_params(shape=shapes.Wall(normal=[-1,0,0], dist = -1.0*top_wall_dist), velocity=v)
    for wall in walls:
        system.lbboundaries.add(wall)
    

    # AVB:  Put some adhesive sites on the bottom wall     
    # AVB:  N is a number of collagen binding sites on a surface
    N = int(collagen_density * box_l* box_l)
    x_coord = np.array([bot_wall_dist+a] * N)
    y_coord = np.array(box_l*np.random.rand(N))
    z_coord = np.array(box_l*np.random.rand(N))   



    # AVB: assign collagen type
    for i in range(N):
        system.part.add(id = next_part_id+i, pos = np.array([x_coord[i], y_coord[i], z_coord[i]]), v = np.array([0, 0, 0]), type = collagen_type)  
        system.part[next_part_id+i].fix = (1, 1, 1)
    next_part_id = N+next_part_id

    print("\n")
    print(" -- next particle id after Collagen creation: ")
    print(next_part_id)  

    # AVB: Setting adsorption potential vwf-coll
    ### AVB: parameters for the potential (amplitude and cut-off redius)
    amplVwfCol = 500.0*kBT    # COLLAGEN ADSORPTION 
    rcutVwfCol = 1.2*r0
    ### And the potential
    system.non_bonded_inter[vwf_type,collagen_type].lennard_jones.set_params(
              epsilon = amplVwfCol, sigma = r0/1.122,
              shift = "auto", cutoff = rcutVwfCol, min = r0*0.6) 
    
    # Setting repulsive potential vwf-wall
    
    system.constraints.add(shape=shapes.Wall(normal=[1,0,0], dist = (bot_wall_dist-2*a)), particle_type=wall_type, penetrable=True)       
    
    ### AVB: parameters for the potential (amplitude and cut-off radius)
    amplVwlCol = 1e-5  
    rcutVwlCol = 2.0*a
    ### And the potential
    system.non_bonded_inter[vwf_type,wall_type].soft_sphere.set_params(
    a = amplVwlCol, n = 1.2,cutoff = rcutVwlCol, offset = -0.5*a)     
    
    
    # configure correlators
    com_pos = ComPosition(ids=(0,))
    c = Correlator(obs1 = com_pos, tau_lin=16, tau_max=loops*step_per_loop, delta_N=1,
            corr_operation="square_distance_componentwise", compress1="discard1")
    system.auto_update_accumulators.add(c)
    
    
    print("Sampling started.")
    print("length after warmup:")
    print(system.analysis.calc_re(chain_start=0, number_of_chains=1, chain_length=mpc-1)[0])
    
    lengths = []
    
    ylengths=[]
    
    for i in range(loops):
        system.integrator.run(step_per_loop)
        system.analysis.append()
        lengths.append(system.analysis.calc_re(chain_start=0, number_of_chains=1, chain_length=mpc-1)[0]) 
        lbf.print_vtk_velocity(outdir+"/"+str(vy)+"fluid%04i.vtk" %i)
        system.part.writevtk(outdir+"/"+str(vy)+"coll%04i.vtk" %i, types=[collagen_type])  
        system.part.writevtk(outdir+"/"+str(vy)+"vwf%04i.vtk" %i, types=[0])  
        cor = list(system.part[:].pos)
        y = []
        for l in cor:
            y.append(l[1])
        ylengths.append(max(y) - min(y))
        
        sys.stdout.write("\rSampling: %05i"%i)
        sys.stdout.flush()
    
    sys.stdout.write("\n")
    # AVB: Removing LBB with shear
    for wall in walls:
        system.lbboundaries.remove(wall) 
    # AVB: Setting new LBB without shear
    walls[0].set_params(shape=shapes.Wall(normal=[1,0,0], dist = bot_wall_dist), velocity=[0,0,0])
    walls[1].set_params(shape=shapes.Wall(normal=[-1,0,0], dist = -1.0*top_wall_dist), velocity=[0,0,0])    
    for wall in walls:
        system.lbboundaries.add(wall)   
        
    for i in range(100):
        system.integrator.run(step_per_loop)
        lengths.append(system.analysis.calc_re(chain_start=0, number_of_chains=1, chain_length=mpc-1)[0])
    
    system.part.writevtk(outdir+"/"+str(vy)+"vwf_all[r0=2,kBT=4]intheEND.vtk")
        
    with open(outdir+"/lengths"+str(vy)+".dat","a") as datafile:
        datafile.write("\n".join(map(str,lengths)))
        
    with open(outdir+"/lengthsY"+str(vy)+".dat","a") as datafile:
                datafile.write("\n".join(map(str, ylengths)))  
                
    mean_vy = [(vy /tscale) / 32, sum(ylengths)/len(ylengths)]
    
    print("\n mean_vy")
    print(mean_vy)
    
    
                
    with open(outdir+"/mean_vy" + ".dat","a") as datafile:
                datafile.write(" ".join(map(str, mean_vy)))     
    
    c.finalize()
    corrdata = c.result()
    corr = zeros((corrdata.shape[0],2))
    corr[:,0] = corrdata[:,0]
    corr[:,1] = (corrdata[:,2] + corrdata[:,3] + corrdata[:,4]) / 3
    savetxt(outdir+"/msd_nom"+str(mpc)+".dat", corr)
    
    with open(outdir+"/rh_out.dat","a") as datafile:
        rh = system.analysis.calc_rh(chain_start=0, number_of_chains=1, chain_length=mpc-1)
        datafile.write(str(mpc)+ "    " + str(rh[0])+"\n")
    

# main course

calc(int(sys.argv[1]), 1000)
   
