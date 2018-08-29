Single Species
===============

File Structure
---------------
**Main file**

    * /diffusion2d/single/single2d.m
        main function for numeric simulation of the dynamics of one single species in 2D

**ALso include**

    * /diffusion2d/single/minmod.m
        a generalized minmod limiter function
    * /diffusion2d/single/drho.m
        gives the time derivative of the density function at given configuration
    * /diffusion2d/single/draw.m
        a help function for 3D visualization
    * /diffusion2d/single/SSPRK3.m
        SSP-RK3 method
    * /diffusion2d/single/Euler.m
        Euler method

**Example**

    * /diffusion2d/single/Example1.m
        Example 1


Functions
------------
.. topic:: single2d

    .. code-block:: matlab

        rho = single2d(rho0, L, W, dt, T)


    SINGLE2D:
        numeric simulation of the dynamics of one single species in 2D

    Input:

        * rho0:        initial density, N by N matrix
        * L:          domain [-L,L] x [-L,L]
        * W:          W(x) interacting potentials
        * dt:         time step
        * T:          simulation time. Total #iterations = T/dt

    Output:

        * rho:        density at t = T, N by N matrix

    Optional parm:

    .. code-block:: matlab
        :emphasize-lines: 3-4,9-10,15, 20, 25-26

        rho = single2d(.. ,H);

        %   H is a symbolic function for internal energy as a
        %   function of the density. Default H(r) = 0
        --------------------------------------------------------------------------

        rho = single2d(.. ,V);

        %   optionally sets the environmental confinement potential V,
        %   which is a NxNmatrix. Default: 0
        --------------------------------------------------------------------------

        rho = single2d(.. ,e)

        %   sets the diffusion coefficient for some e > 0. Default e = 0
        --------------------------------------------------------------------------

        rho = single2d(.. ,'v') or rho = single2d(.. ,'V')

        %   enables visual display during the simulation. Default disabled.
        --------------------------------------------------------------------------

        rho = single2d(.. ,'solver')

        %   where 'solver' sets the numeric method used for ODE.
        %   Possible options: 'Euler', 'SSPRK3'. Default 'SSPRK3'.
        --------------------------------------------------------------------------

.. topic:: drho

    .. code-block:: matlab

        dr = drho(r, K, dH, V, dx, e)


    DRHO:
            gives time derivative of the density function evaluated at r

    Input:

            * r:          density distribution
            * K:          convolution matrix
            * dH:         a symbolic function of the derivative of H
            * V:          matrix of confinement potential
            * dx:         time step
            * e:          diffusion coefficient

    Output:

            * dr:        time derivative of the density function
