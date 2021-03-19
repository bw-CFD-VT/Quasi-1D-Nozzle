# Quasi-1D-Nozzle
Brendan Walsh
AOE 6145
Quasi 1D Nozzle FVM Code 

------- Version 1.0 (03/14/21) ----------- 
- Even # of cell definition required to ensure interface at throat
- Central Flux Quadrature Approximation F(U_i+1/2)
- Artificial Dissipation (JST damping based on implementation discussed in lecture + Swanson, Radespiel and Turkel (1945) )
- 1 Ghost Cell at Inflow and Outflow
