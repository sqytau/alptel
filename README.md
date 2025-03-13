# Geant4 simulation of ALPIDE telescope.
## Installation
Requires Geant4 v11.

Assuming Geant4 is installed:
```
cmake ./
make
```
## Running
```
./alptelsim telescope_sim.mac 1
```
This will simulate 5k cosmic muons, parameter 1 means 1 thread.
It is build in multithreading mode but not tested...

Runnig it without parameters will open Qt GUI, of course, if Qt is installed and Geat4 link with it.
```
./alptelsim
```
Then in command line of the GUI:
```
/control/execute vis_ev_e_v1.mac
```
will visualize the geometry.

