The adaptive Neural Network Library (MATLAB® 5.3.1 and later) is a collection of blocks that implement several Adaptive Neural Networks featuring different adaptation algorithms.

It was developed mainly in June-July 2001 by Giampiero Campa (West Virginia University) and Mario Fravolini (Perugia University). Later improvements were partially supported by the NASA Grant NCC5-685.

There are blocks that implement basically these kinds of neural networks:

Adaptive Linear Networks (ADALINE)

Multilayer Layer Perceptron Networks

Generalized Radial Basis Functions Networks

Dynamic Cell Structure (DCS) Networks with gaussian or conical basis functions

Also, a Simulink example regarding the approximation of a scalar nonlinear function is included.

Finally, the folder "training" includes step by step instrucions on how to train the GRBF network and the supporting example.


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
