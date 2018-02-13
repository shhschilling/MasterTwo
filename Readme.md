The Mathematica package MasterTwo allows the calculation of one- and two-loop integrals which can be reduced to scalar integrals independent of external momenta and depending on up to two masses. 

The main programme MasterTwo.m calls the two subpackages

- Fermions.m:  containing functions for standard Dirac algebra, traces etc...

- Integrals.m: containing functions for Taylor expansion, tensor reduction, partial fraction

The package is equipped with an on-line help to each

command. MasterTwoInfo[] produces a list with all the available

commands and ?command prints a short information on syntax and effect

of command.

A detailed description of all functionalities is given in the manual being part of this distribution. 

INSTALLATION UNDER LINUX 

There is an installation program with the Package. It copies the files where they belong, updates the init.m

file in the .Mathematica/Autoload/directory in the home directory, so that you can load the package without having to give the whole path.

To install the package,

 - copy the zip file MasterTwo-1.0.zip to your disk,

 - unpack it with unzip MasterTwo-1.0.zip.

 - Change directory to MasterTwo-1.0

 - Change the permission of the installation script:chmod +x MasterTwoInstal

 - Execute it: ./MasterTwoInstall 

Follow the instructions... 

 
UNINSTALLATION

Run the program ./MasterTwoUninstall in the installation directory of MasterTwo.  

INSTALLATION UNDER MACOS and Windows

To install the package

- Start Mathematica and type $Path in order to find out the exact path of the 

  Autoload directories of your Mathematica distribution

- Copy the files

                 Fermions.m

                 Integrals.m

                 MasterTwo.m  


