
This v51 branch has been ported to foam-extend-v5.1. 
It compiles and runs against a slightly modified version of v5.1.

However, significant issues still remain.

1. The bdf2 ddtScheme was not ported.
2. All partitioned DNA loops seem to seg-fault in transportTemperature::setCoupledEqns()
3. I only got two tutorials to run using v5.1:
   conjugateHeatTransfer/flowOverHeatedPlate for monolithic coupling
   multiphaseFlow/2dRisingBubble for monolithic coupling (I do not try running the 3d bubbles)

greg.burgreen@msstate.edu
March 26, 2025

