# Drift Wave Toy Model (Vlasov)

This repository contains a simple kinetic toy model used to test the stability of a density gradient in a drift-wave configuration using the Gkeyll Vlasov solver.

## Model Description

The setup is intentionally simplified to study gradient stability.  
The background plasma density follows a smooth spatial profile (cosine/tanh-type gradient), and the electron and ion masses are taken moderate(`mi/me = 25`) to keep the computational cost manageable.

The goal of these tests was to determine whether the imposed density gradient remains stable under the Vlasov evolution.

After successive runs, the gradient appears to be stable for a fair amount of simulation time.

## Running the Simulation

All simulations were executed using the pre-built Gkeyll binary available on the Dartmouth cluster:

/dartfs/rc/lab/E/EPaCO/libs/gkylsoft/gkeyll/bin/gkeyll drift_wave_instability.lua


## Gkeyll

The Gkeyll plasma simulation framework can be downloaded from the official project page:

https://gkeyll.readthedocs.io

## Notes

This repository contains only the **input scripts used for the simulation tests**.  
The solver itself is not included.

The three folders contain Lua input scripts for Gkeyll simulations and Python scripts for visualization; the goal was to start from a 2D2V Vlasov–Maxwell solver benchmark (the Weibel instability) and progressively add the relevant physics to study the fully kinetic drift-wave instability, whose reduced fluid description is captured by the Hasegawa–Wakatani model (also included), with the Python analysis scripts requiring Postgkyl (documentation available at https://gkeyll.readthedocs.io).
P.S. Vlasov Maxwell drift wave instability runs aren't fully consistent with existing studies and is still under development.
