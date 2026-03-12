# Drift Wave Toy Model (Vlasov)

This repository contains a simple kinetic toy model used to test the stability of a density gradient in a drift-wave configuration using the Gkeyll Vlasov solver.

## Model Description

The setup is intentionally simplified to study gradient stability.  
The background plasma density follows a smooth spatial profile (cosine/tanh-type gradient), and the electron and ion masses are taken equal (`mi/me = 1`) to keep the computational cost manageable.

The goal of these tests was to determine whether the imposed density gradient remains stable under the Vlasov evolution.

During testing it was observed that the gradient diffused rapidly at low resolution. After increasing the spatial resolution and strengthening the guide magnetic field, the density profile remains significantly more stable.

## Running the Simulation

All simulations were executed using the pre-built Gkeyll binary available on the Dartmouth cluster:

/dartfs/rc/lab/E/EPaCO/libs/gkylsoft/gkeyll/bin/gkeyll drift_wave_instability.lua


## Gkeyll

The Gkeyll plasma simulation framework can be downloaded from the official project page:

https://gkeyll.readthedocs.io

## Notes

This repository contains only the **input scripts used for the simulation tests**.  
The solver itself is not included.
