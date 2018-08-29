Two Species
===============

File Structure
---------------
**Main file**

    * /diffusion2d/two/two2d.m
        main function for numeric simulation of the dynamics of two species in 2D

**ALso include**

    * /diffusion2d/two/minmod.m
        a generalized minmod limiter function
    * /diffusion2d/two/drho2.m
        gives the time derivatives of the density functions at given configuration
    * /diffusion2d/two/draw.m
        a help function for 3D visualization
    * /diffusion2d/two/SSPRK3.m
        SSP-RK3 method
    * /diffusion2d/two/Euler.m
        Euler method

**Examples**

    * /diffusion2d/two/Example_two_pulses.m
        :ref:`two_pul`
    * /diffusion2d/two/Example_three_pulses.m
        :ref:`three_pul`


Functions
------------
.. topic:: two2d

    .. code-block:: matlab

        [r1,r2] = two2d(r1_0,r2_0, L, W, dt, T)


    TWO2D:
        numeric simulation of the dynamics of two species in 2D

    Input:

        * r1_0,r2_0:  initial density of the two species, N by N matrix
        * L:          domain [-L,L] x [-L,L]
        * W:          W(x) interacting potentials
        * dt:         time step
        * T:          simulation time. Total #iterations = T/dt

    Output:

        * r1, r2:        densities at t = T, N by N matrix

    Optional parm:

    .. code-block:: matlab
        :emphasize-lines: 3-4,9-10,15, 20, 25-26

        [r1,r2] = two2d(.. ,H);

        %   H is a symbolic function for internal energy as a
        %   function of the density. Default H(r) = 0
        --------------------------------------------------------------------------

        [r1,r2] = two2d(.. ,V);

        %   optionally sets the environmental confinement potential V,
        %   which is a NxNmatrix. Default: 0
        --------------------------------------------------------------------------

        [r1,r2] = two2d(.. ,e)

        %   sets the diffusion coefficient for some e > 0. Default e = 0
        --------------------------------------------------------------------------

        [r1,r2] = two2d(.. ,'v') or rho = single2d(.. ,'V')

        %   enables visual display during the simulation. Default disabled.
        --------------------------------------------------------------------------

        [r1,r2] = two2d(.. ,'solver')

        %   where 'solver' sets the numeric method used for ODE.
        %   Possible options: 'Euler', 'SSPRK3'. Default 'SSPRK3'.
        --------------------------------------------------------------------------

.. topic:: drho2

    .. code-block:: matlab

        dr = drho2(r, K11,K12,K21,K22, dH, V, dx, e)


    DRHO2:
            gives time derivatives of the density functions evaluated at r

    Input:

            * r:          density distribution
            * K:          convolution matrix
            * dH:         a symbolic function of the derivative of H
            * V:          matrix of confinement potential
            * dx:         time step
            * e:          diffusion coefficient

    Output:

            * dr:        time derivative of the density function
