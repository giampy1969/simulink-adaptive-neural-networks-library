CONTENTS:

-------------------------------------------------------------------------

ann.mdl 	- adaptive neural network library
anndemo.mdl	- demo file

contents.m 	- info file for the Matlab Help and Ver commands
info.xml 	- info file for the Matlab Start Menu
demos.m 	- info file for the Matlab Demos Menu

emran8.dll	- dll file for grbf network
dcsgl2.dll	- dll file for the gdcs network 
vrmult.dll	- dll file for real matrices multiplication

source	- folder containing source code
training	- folder containing example on grbf training

-------------------------------------------------------------------------

INSTALLATION:

Unzip the file in a folder of choice, usually under the "toolbox" folder
under the Matlab main installation folder.

From MATLAB®, open the "Set Path" interface from the "File" menu,
and use the "Add with Subfolders" feature to add the new folder 
to the matlab path. Then Save and exit the Set Path dialog box.

-------------------------------------------------------------------------

HOW TO USE THE LIBRARY:

From MATLAB, use the command "ann" to open the library, and double click 
on the "ann demo" block to open the main demo. This should open a simulink 
scheme capable of running "out of the box" with any matlab version from
R11 to R2007a, that is from Matlab 5.3 to 7.4.

To use the library blocks, just drag and drop them into any simulink scheme.
The folder "training" contains instructions on how to train the GRBF network.

-------------------------------------------------------------------------

ACKNOWLEDGMENTS:

Partial support for the latter development phases of this library 
was given by the NASA Grant # NCC5-685

-------------------------------------------------------------------------

G.Campa, West Virginia University, October 2007
