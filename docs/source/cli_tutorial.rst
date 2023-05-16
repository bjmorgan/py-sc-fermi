CLI tutorial
===============

Inputs
-------

``py-sc-fermi`` features a command line tool for calculating self consistent Fermi
energies. In this case defect system is defined by a ``.yaml`` file structured like so::

    bandgap: 1 # replace with your calculated bandgap
    temperature: 300 # the temperature to be used in the solution for Fermi energy
    nelect: 18 # number of electrons in the density of states calculation carried out on the unit cell
    defect_species:
      - V_Na: # name of first defect species
          nsites: 1 # site degeneracy of the defect species
          charge_states:
            0: # charge of first charge state
              charge: 0
              formation_energy: 2 # formation energy of first charge state
              degeneracy: 1 # degeneracy of first charge state
            -1:
              charge: -1
              formation_energy: 1
              degeneracy: 2
      ... # repeat for each defect in your system

The dimension of the unit cell can be given within a ``.cif`` file (or a number of other 
common structure files), or as a file which is structured as::

    1 # scaling factor for lattice vector components
    4.1 0.0 0.0
    0.0 4.1 0.0
    0.0 0.0 4.1

where the matrix above defines the lattice vector of the unit cell. Alternatively,
we can skip the definition of the full unit cell, and simply include ``volume: x`` (where
``x`` is the unit cell volume) in the ``.yaml`` file described above.

The density of states can either be specified in the ``.yaml`` file as::

    ...
    edos: [energy values of the dos]
    dos: [array of total dos values]
    ...

or the density of states can be read directly from a ``vasprun.xml`` file.

sc_fermi_solve
---------------

the command to solve for the self consistent Fermi energy once you have specifed 
the inputs as above is ``sc_fermi_solve [input .yaml]``. If the ``.yaml`` file specifies
the volume and dos information for your system, this is all you need to do. Otherwise, there
are some additional arguments accepted by this command:

   - ``-s, --structure_file`` path to the structure file which defines the volume of the unit cell
   - ``-d, --dos_file`` path to the file which defines the dos
       - if this argument is specified, you must also specify ``-b, --band_gap`` which gives
         the bulk band-gap of the system.

frozen-concentration defects 
-----------------------------

You are able to specify the concentrations of different defects in the solver, just as with the
API. In this case, you may either fix the concentration of the defect species by specifying e.g::

    defect_species:
        - V_Na:
            fixed_concentration: 1e20
            ...

and specify the concentration of a defect charge state like so::

        defect_species:
        - V_Na:
            nsites: 1
            charge_states:
                -1:
                    charge: -1
                    formation_energy: 1
                    fixed_concentration: 1e20
                    degeneracy: 1
        ...

within the input ``.yaml`` file. If you do so, you must add the flag ``-f`` or 
``--frozen_defects`` when you call ``sc_fermi_solve``.
            

       



