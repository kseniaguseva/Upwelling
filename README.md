This code produces the results of the paper *Guseva et al (2019)*,
"The effect of intermittent upwelling events on plankton blooms."

The code consists of a semi-Lagrangian algorithm to solve the
reaction-advectio-diffusion system. It is written in Python, with
additional functions imported from C++ using Boost.Python. To run the
code make sure to have installed Python 3 and Boost.Python.


# File list:
--------


1. **flow.py**: The hydrodynamic model in agreement with *Jung et al (1993)*
and *Sandulescu et al (2006)*.  The default parameters are taken from
*Sandulescu et al (2006)*.

2. **population.py**: The population model from *Chakraborty et al
   (2014)*, the default parameters are also the ones from this
   manuscript.  The time and length scales of the flow field are used
   to make these parameters dimensionless.  Additional functions are
   included to evolve the state of the system and computes averages
   of concentrations.

3. **space.py**: Defines the space grid. Contains functions that help to
   evolve the species concentrations at grid points: diffusion,
   interpolation.

4. **upwelling.py**: Builds the time series of the upwelling events.
   Defines the upwelling region.

5. **intermit.py**: Generates the intermittent time series using the
   dynamical system from *Platt et al (1993)*.

6. **set_output.py**: Functions used to save the output.

7. **sys1.py**: Parameters used for system (I) of the paper.

8. **sys2.py**: Parameters used for system (II) of the paper.

9. **upw_loops.cc**: Function used for interpolation, diffusion and
to compute space averages.

10. **numpy_bind.hh**, **demangle.hh**, **demangle.cc**: Functions for
moving arrays between Numpy and C++, requires Boost.Python.


# How to run the program.
-----------------------


Run make. Create two folders for the output data:

* Data_space;
* Data_time series.

To run the code use **main_timeseries.py**.

