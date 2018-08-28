=======================================================
Finite Volume Method for Nonlinear Nonlocal Equations
=======================================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

Overview
---------

**Model Description:** For N groups in a 2D space, let :math:`\rho_i(\mathbf{x}, t)` denotes the mass density distribution of group *i*
at time *t*, :math:`i=1,2,...,N`. The problem we are interested in:

.. math::

    \partial_t\rho_i = \nabla\cdot\Big[\rho_i\nabla\Big(H'(\rho)+V(\mathbf{x})+W_{ii}\ast\rho_i+\sum\limits_{j\neq{i}}{W_{ij}\ast\rho_j}+\epsilon\rho\Big)\Big], \mathbf{x}\in\mathbb{R}^2, t\gt{0}


:math:`\mathbf{\rho_i}` \: mass density of group *i*

:math:`\mathbf{\rho}`\: :math:`\sum{\rho_i}` , sum of densities of all N groups

:math:`\mathbf{H(\rho)}`\: density of internal energy

:math:`\mathbf{V(\mathbf{x})}`\: environmental confinement potential

:math:`\mathbf{W_{ii}}`\: self-interaction potential (intraspecific interaction potential)

:math:`\mathbf{W_{ij}}(j\neq{i})`\: cross-interaction potential (interspecific interaction potential)

:math:`\mathbf{\epsilon}`\: diffusion coefficient

with initial condition :math:`\rho_i(\mathbf{x},0)=\rho_{i0}(\mathbf{x})`.


**Numeric Analysis:** This numeric scheme is developed based on the finite volume method described in `Links`_ [2]. The paper proposes a
finite method scheme for nonlinear nonlocal system for a single group in one and two spacial dimensions. We extend the scheme to include more than one group,
with cross-interaction mechanism mentioned in `Links`_ [3].


.. topic:: System Requirements

    ================================== ==============================
     Software Requires                  MATLAB Release Compatibility
    ================================== ==============================
     MATLAB                             Created with R2018a
    ================================== ==============================



Links
----------

#. Source Code: `Github <github.com/uidmice/finiteVolumeMethod>`_
#. A Finite-Volume Method for Nonlinear Nonlocal Equations with a Gradient Flow Structure: `arXiv:1402.4252v2 <http://de.arxiv.org/abs/1402.4252v2>`_
#. Zoology of a Nonlocal Cross-Diffusion Model for Two Species: `Permalink <https://doi.org/10.1137/17M1128782>`_

Support
-------

If you are having issues, please email `<yydxhhuyi@gmail.com>`_

License
-------

.. The project is licensed under the BSD license.
