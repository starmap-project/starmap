=================================================================
StaRMAP
=================================================================
A second order staggered grid finite difference solver for linear
hyperbolic moment approximations to the equations of radiative
transfer in 1D, 2D and 3D geometry.

Version 2.0
Copyright (c) 06/28/2022 Benjamin Seibold, Martin Frank, and
                         Rujeko Chinomona
http://www.math.temple.edu/~seibold
https://www.scc.kit.edu/personen/martin.frank.php
Contributers: Edgar Olbrant (v1.0), Kerstin Kuepper (v1.5,v2.0),

StaRMAP project website:
https://github.com/starmap-project

Older versions of StaRMAP:
http://math.temple.edu/~seibold/research/starmap/


=================================================================
Files
=================================================================
starmap_closure_pn.m           Creates P_N moment matrices
starmap_closure_spn.m          Creates SP_N moment matrices
starmap_1d_ex_gauss.m          Example file: Gaussian initial condition
starmap_1d_ex_planesource.m    Example file: plane source test 
starmap_2d_create_beam.m       Creates file starmap_2d_ex_beam_auto.m
 -> starmap_2d_ex_beam_auto.m  Example file: beam in void and medium
starmap_2d_create_hemisphere_filter_convergence.m
                               Creates file starmap_2d_ex_hemisphere_
                               filter_convergence_auto.m
 -> starmap_2d_ex_hemisphere_filter_convergence_auto.m
                               Example file: FPN convergence for
                               hemisphere test
starmap_2d_create_mms.m        Creates file starmap_ex_mms_auto.m
 -> starmap_2d_ex_mms_auto.m   Example file: manufactured solution
starmap_2d_ex_asymptotic_preserving.m
                               Example file: asymptotic preserving 
                                             properties
starmap_2d_ex_boxes.m          Example file: boxes test case
starmap_2d_ex_controlrod.m     Example file: control rod test
starmap_2d_ex_gauss_filter_convergence.m
                               Example file: FPN convergence for
                               Gauss test
starmap_2d_ex_homogeneous.m    Example file: homogeneous material 
                                             coefficients
starmap_2d_ex_l2norm.m         Example file: l2norm of solution
starmap_2d_ex_lattice.m        Example file: lattice/checkerboard
                               test
starmap_2d_ex_lattice_filter.m Example file: lattice/checkerboard
                               test with filtering
starmap_2d_ex_lattice_filter_convergence.m
                               Example file: FPN convergence for
                               lattice/checkerboard test
starmap_2d_ex_linesource.m     Example file: line source test
starmap_2d_ex_linesource_filter.m
                               Example file: line source test with
                               filtering
starmap_2d_ex_linesource_filter_reconstruction.m
                               Example file: line source test with
                               angular reconstruction in a point
starmap_3d_ex_lattice.m        Example file: lattice test in 3D
starmap_3d_ex_pointsource.m    Example file: point source test 
starmap_init.m                 Initializes data structures compatible with 
                               solver
starmap_solver.m               The computational code
ChangeLog.txt                  Changes from previous versions
README.txt                     This file

=================================================================
Operation
=================================================================
To run a simulation with StaRMAP, the file starmap_init.m must first
be called with a struct prob in which the problem parameters are defined.
Within starmap_init.m, variables in user-defined struct prob are transformed 
into data structures compatible with the solver, defaults are set and 
closure matrices are defined. This is stored in a new struct par. The 
starmap_solver.m must then be called with the struct par. 
The meaning of most parameters can be seen in the
"Parameter Defaults" part in starmap_init.m, in the provided
example files starmap_ex_*.m, and in the paper
http://arxiv.org/abs/1211.2205.
Below explained is the syntax for the parameter functions that
define the problem. General philosophy: do not prescribe
parameters when not needed; for special cases (e.g.,
time-independence, isotropy) the code will run faster.
=================================================================

=================================================================
Syntax for functions (which parameters need to be provided)
=================================================================
Initial conditions: par.ic: 
1D: (x) or (x,k)
2D: (x,y) or (x,y,k)
3D: (x,y,z) or (x,y,z,k)
 - If only spatial coordinates specified, then initial conditions 
are only nonzero for k=1 (i.e., zeroth moment)

Absorption coefficient: sigma_a: 
1D: (x) or (x,t)
2D: (x,y) or (x,y,t)
3D: (x,y,z) or (x,y,z,t)
 - If only spatial coordinates specified, the function is only evaluated once
   (in the first time step).

Isotropic scattering coefficient: par.sigma_s0: 
1D: (x) or (x,t)
2D: (x,y) or (x,y,t)
3D: (x,y,z) or (x,y,z,t)
 - Zeroth moment of scattering kernel (i.e., isotropic part).
 - If only spatial coordinates specified, the function is only evaluated once.

Anisotropic scattering: par.sigma_sm: 
1D: (x,m) or (x,m,t)
2D: (x,y,m) or (x,y,m,t)
3D: (x,y,z,m) or (x,y,z,m,t)
 - Higher moments of scattering kernel.
 - Evaluated for m>=1, where m = par.mom_order(k) moment order
   and k = system component.
 - Vector par.mom_order must exist if this function is specified.
 - If only spatial coordinates and moment order specified, 
  the function is only evaluated once (in the first time step).

Source: par.source: 
1D: (x,t) or (x,t,k)
2D: (x,y,t) or (x,y,t,k)
3D: (x,y,z,t) or (x,y,z,t,k)
 - If only spatial coordinates and time specified, then source 
   only active for k=1 (i.e., zeroth moment)
 - Via the vector par.source_ind one can prescribe in which
   moments the source term is active (only those will be
   evaluated).

Filter: par.filterfunction: 
1D: (par,dt,k,x) or (par,dt,k,x,t)
2D: (par,dt,k,x,y) or (par,dt,k,x,y,t)
3D: (par,dt,k,x,y,z) or (par,dt,k,x,y,z,t)
- Multiplication factor for decay coefficients.

- Vector par.mom_order must exist if this function is specified.
- If t not specified, the function is only evaluated once
  (in the first time step).
- Via par.f_position one can define if the filter is applied
  after each 'substep' or 'full' time-step ('' turns off the
  filer).
=================================================================
