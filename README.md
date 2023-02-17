# A Monte Carlo Simulation package for STIX


# Requirements
- ROOT version > 6.0
- cmake
- Geant4.10.05, other geant4 versions may also work but have not been tested
# Compile and run
```sh
cmake .
make
./g4STIX  -m test.mac -o test.root -Ba133
```
# Materials and positions 
see https://docs.google.com/spreadsheets/d/1_cStJZBpLF8T0uwpfPfBadFIsnr4oUwPPqxzYmvimJ4/edit?usp=sharing



#history
2023-02-17
  added pads to caliste, pins, 
  attenuator stays in the model when att is out
  added PCB, approximation of PCB


