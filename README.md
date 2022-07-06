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

## comi 

The `comi` (Command Input Files) directory has files that contain Tosca/Opera code 
that can generate a 3D model of a magnet system when loaded into Opera. 

## C++ Files 

The directories `include`, `src`, and `input` contain code used for the C++ files 
(`*.cxx`).  These scripts work in the CERN ROOT framework and are used for visualizing 
and analyzing data produced by Tosca/Opera.  
