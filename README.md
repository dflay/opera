# opera

Tosca (Opera) analysis code and command input files for tosca modeling. 
This code is focused on the SBS GMn and GEn-II experiments at Jefferson Lab 
in Hall A.

# Prerequisites 

For analysis and visualization of the Tosca/Opera magnetic field maps, the user should 
have the following installed:  
- The CERN [ROOT framework](https://root.cern.ch)
- The utility library [util_df](https://github.com/dflay/util_df)

# Directories 

## C++ Files 

The directories `include`, `src`, and `input` contain code and text files used for 
the C++ files (`*.cxx`).  These scripts work in the CERN ROOT framework and are used 
for visualizing and analyzing data produced by Tosca/Opera. 

## opc 

The `opc` directory has the 3D Opera model files (`*.opc`) that are used to generate 
the particle trajectories and magnetic field maps from the BigBite and SBS magnet 
configurations for the GMn and GEn experiments.  The subdirectory `bogdan` has 
the builds created by Bogdan, while the `flay` directory has those same builds modified 
by David Flay to resolve geometry meshing conflicts and updates to data storage levels 
and magnetic potential definitions of all materials to improve the accuracy of the 
output of the Opera simulation.  
 
## comi 

The `comi` (Command Input Files) directory has files that contain Tosca/Opera code 
that can generate a 3D model of a magnet system when loaded into Opera.  These were 
built by David Flay, but ultimately not as complete as the `*.opc` files in the `opc` 
directory since they focus on the SBS magnet only, and do not have the BigBite magnet nor 
the bracing structures for both the BigBite and SBS magnets.  

