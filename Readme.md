
# Abstract
The Mathematica package MasterTwo allows the calculation of one* and two*loop Feynman integrals which can be reduced to scalar integrals independent of external momenta and depending on up to two different masses.
It consists of two subpackages, Fermions~ and Integrals.  Fermions covers the standard Dirac Algebra to transform the integrals in scalar integrals. Integrals perform the subsequent Taylor expansion, partial fraction, tensor reduction and the integration of the thus achieved scalar integrals.



# Structure
The main programme  `MasterTwo.m` calls these two subpackages

* `Fermions.m`:  containing functions for standard Dirac algebra, traces etc...

* `Integrals.m`: containing functions for Taylor expansion, tensor reduction, partial fraction

# Help

The package is equipped with an online help to each command.
`MasterTwoInfo[]` produces a list with all the available
commands and `
?command `prints a short information on syntax and effect of command. A detailed description of all functionalities is given in the manual being part of this distribution. 

# Installation under Linux

The installer MasterTwoInstall copies the files into the correct path, updates the init.m file in the .Mathematica/Autoload/directory in the home directory,
To install the package,

 * copy the zip file MasterTwo-1.0.zip to your disk,

 * unpack it with unzip MasterTwo-1.0.zip.

 * Change directory to MasterTwo-1.0

 * Change the permission of the installation script:
  ` chmod +x MasterTwoInstall`

 * Execute it: `./MasterTwoInstall `

Follow the instructions... 

 
# Uninstallation under Linux

Run the program ./MasterTwoUninstall in the installation directory of MasterTwo.  

# Installatiion under  MACOS and Windows

To install the package

* Quit Mathematica if it is running. Uninstall the older version of MasterTwo,
if present.
* Start Mathematica and type $Path in order to find out the exact path of the
   Autoload directories of your Mathematica distribution

* Copy the files
  * `Fermions.m`
  *  `Integrals.m `
  *  `MasterTwo.m `

in one of the Autoload directories.

*  Close Mathematica.
*  In your next Mathematica session you can start the package by typing in
`<<MasterTwo` `


