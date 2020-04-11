################### 1D Laser Wakefield with envelope
ups = 1 #upscaling of simulation in space 
dx = 1/ups 
#dtrans = 3.
dt = 0.8*dx 
nx = 192*ups #number of cells 
Lx = nx * dx # length of simulation box
npatch_x = 32*ups
laser_fwhm = 2 
center_laser = Lx-2.*laser_fwhm # the temporal center here is the same as waist position, but in principle they can differ
time_start_moving_window =  0.


Main(
    geometry = "1Dcartesian",
    interpolation_order = 2, #defines particle shape function, 2:3-point stencil
    timestep = dt, 
    simulation_time = ups*350.*dt,
    cell_length  = [dx],
    grid_length = [Lx],
    number_of_patches =[npatch_x], #the cells are split inte npatch_x number of patches. 
                                   #Each patch will be calculated and stored in differnt locations
    # clrw = nx/npatch_x, 
    EM_boundary_conditions = [ ["silver-muller"] ], # absorbing boundary conditions
    solve_poisson = False, # Decides if Poisson correction must be applied or not initially
    print_every = 100, # output on screen
    random_seed = smilei_mpi_rank #random seed chooser
)

MovingWindow(
    time_start = time_start_moving_window,
    velocity_x = 1 
)

# to make the simulation run smoother
LoadBalancing(
    initial_balance = False,
        every = 20,
    cell_load = 1.,
    frozen_particle_load = 0.1
)

Species(
    name = "electron",
    position_initialization = "regular",
    momentum_initialization = "cold",
    particles_per_cell = 8,
    c_part_max = 1.0, #???
    ponderomotive_dynamics = True, # = this species interacts with laser envelope
    mass = 1.0, 
    charge = -1.0,
    charge_density = polygonal(xpoints=[center_laser+2.*laser_fwhm,
                                        center_laser+2.1*laser_fwhm,
                                        15000,
                                        20000],
                               xvalues=[0.,0.0045,0.0045,0.]), #initial density
    mean_velocity = [0.0, 0.0, 0.0],
    temperature = [0.0],
    pusher = "ponderomotive_boris", # pusher to interact with envelope / scheme to advance momenta
    time_frozen = 0.0,
    boundary_conditions = [
       ["remove", "remove"],
    ], # particles are removed at boundary
)

LaserEnvelopePlanar1D( # linear regime of LWFA
    a0              = 0.1, #normalized vectorpotential  
    time_envelope   = tgaussian(duration = Lx, # default simulation time
                                center=center_laser, # default center of puls
                                fwhm=laser_fwhm), #full width half max
    envelope_solver = 'explicit',
    Envelope_boundary_conditions = [ ["reflective"],],
)


'''
For restart
Checkpoints(
    dump_step = 0,
    dump_minutes = 0.0,
    exit_after_dump = False,
)
'''

# collect data of fields
list_fields = ['Ex', # Electric field component
               'Rho', #Total density
               'Env_A_abs', #laser vector potential component along the polarization direction
               'Env_Chi', #susceptability
               'Env_E_abs' # Electric field component along the polarization direction
               ]
DiagFields(
   every = 50,
        fields = list_fields
)


DiagScalar(every = 10, vars=['Env_A_absMax','Env_E_absMax'])
