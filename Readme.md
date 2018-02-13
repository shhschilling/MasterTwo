
# Purpose

The calculation of loop decays is a tedious work, which can rarely be performed in a reasonable time by calculating all arising diagrams by hand.
The Mathematica package `MasterTwo `allows the automated  calculation of  all one- and two-loop Feynman integrals reducable to scalar integrals independent of external momenta and depending on up to two different masses. In contrast to other programmes like `Reduce`and `Form` it works completely inside `Mathematica`. Compared to other programmes like `HIP`, `Tracer` and `FeynArts` , `MasterTwo` is
much smaller.

It consists of two subpackages, `Fermions` and `Integrals`.  `Fermions` covers the standard Dirac Algebra allowing the transformation of the integrals in scalar integrals. `Integrals` perform the subsequent Taylor expansion, partial fraction, tensor reduction and the integration of the thus achieved scalar integrals.


# Structure
The main programme  `MasterTwo.m` calls these two subpackages

* `Fermions.m`:  containing functions for standard Dirac algebra  allowing the transformation of the integrals in scalar integrals.

* `Integrals.m`: containing functions for  Taylor expansion, partial fraction, tensor reduction and the integration of the thus achieved scalar integrals.

# Help

The package is equipped with an online help to each command.
`MasterTwoInfo[]` produces a list with all the available
commands and `
?command `prints a short information on syntax and effect of command. A detailed description of all functionalities is given in the manual being part of this distribution.

# Installation under Linux

The installer `MasterTwoInstall` copies the files into the correct path, updates the `init.m` file in the ` .Mathematica/Autoload/directory` in the home directory,
To install the package,
* copy the zip file `MasterTwo-1.0.zip `to your disk,
* unpack it with `unzip MasterTwo-1.0.zip`.
* Change directory to `MasterTwo-1.0`
* Change the permission of the installation script:
  `chmod +x MasterTwoInstall`
* Execute it: `./MasterTwoInstall `
Follow the instructions.

# Uninstallation under Linux

Run the program `./MasterTwoUninstall` in the installation directory of `MasterTwo`.

# Installation under  MACOS and Windows

To install the package

* Quit `Mathematica` . Eventually uninstall  older version of `MasterTwo`.
* Start Mathematica
* Type $Path: Lists the path(s)  of the Autoload directories of your Mathematica distribution
* Copy the files
  * `Fermions.m`
  *  `Integrals.m `
  *  `MasterTwo.m `
in one of the Autoload directories.
*  Close Mathematica.
*  In your next Mathematica session you can start the package by typing in
<<MasterTwo`


