Geant4 simulation of ALPIDE telescope.
Requires Geant4 v11.
Assuming Geant4 is installed
cmake ./
make

./alptelsim telescope_sim.mac 1
This will simulate 5k cosmic muons, 1 means 1 thread.
It is build in multithreading mode but not tested...

Runnig it without parameter will open Qt GUI, of course, if Qt is installed and Geat4 link with it.

./alptelsim
Then in command linre of the GUI:
/control/execute vis_ev_e_v1.mac
will open the geometry.

