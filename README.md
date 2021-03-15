# Quasi-1D-Nozzle
Quasi 1D Nozzle FVM Code in support of HW #2

// Solver Notes
//-------------------------------------//
// Even # of cell definition required to ensure interface at throat
// Flux Quadrature -> Centered Approximation F(U_i+1/2)
// JST Artificial Dissipation based on implementation discussed in lecture and Swanson,Radespiel and Turkel (1945)
// 1 Ghost Cell at Inflow and Outflow
