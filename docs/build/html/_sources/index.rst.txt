=======================================================
Finite Volume Method for Nonlinear Nonlocal Equations
=======================================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

Overview
---------


This MATLAB project explores the following system:

.. math::

    \partial_t\rho_i = \nabla\cdot[\rho_i\nabla(H'(\rho_i)+V(x)+W\ast\rho_i)]

where :math:`\rho_i` describes a probability measure of one species in the system

implements the finite volume scheme described in `this paper <http://de.arxiv.org/abs/1402.4252v2>`_, and extending the scheme to include more than one species.

.. topic:: System Requirements

    ================================== ==============================
     Software Requires                  MATLAB Release Compatibility
    ================================== ==============================
     MATLAB                             Created with R2018a          
    ================================== ==============================

Links
----------

- Source Code: `Github <github.com/uidmice/finiteVolumeMethod>`_
- A Finite-Volume Method for Nonlinear Nonlocal Equations with a Gradient Flow Structure: `arXiv:1402.4252v2 <http://de.arxiv.org/abs/1402.4252v2>`_


Support
-------

If you are having issues, please email `<yydxhhuyi@gmail.com>`_

License
-------

.. The project is licensed under the BSD license.
