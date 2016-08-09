# NuSim - example scripts for Geant4 simulation for neutrinos

```
git clone git@github.com:hep-skku/NuSim.git
source setup.shrc
mkdir run
cd run
cmake ..
cmake --build .
./gen.exe event.hepmc 100 # Generate 100 events
#./gen.exe event2.hepmc 100 1000 # one can add more events with different initial seed and event number
./sim.exe sim.in
```
