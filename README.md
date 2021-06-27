# SSFB_cost_calculation
Calculates system-level capital cost of a SSFB system

**Calculation files**
1. Capital_cost.m
Calculates capital cost for different discharge durations using cost of power calculated from ssfbcostanalysis.m and cost of energy

2. ssfbcostanalysis.m 
Optimizes channel height for lowest cost of power for three cost estimates (upper, mid, lower)

3. nonlinconstr_SSFB.m
Optimization constraint called by ssfbcostanalysis.m file

**Input files**
1. Catholyte_4p2CB 11MnO2.mat (_upper.mat, _lower.mat)
Cost of energy and other inputs to calculate cost of energy using "Derived constants.m"

2. Cost components for power_mid.mat (_upper.mat, _lower.mat)
Cost of power inputs

**Instruction**
1. Download the folder and run Capital_cost.m file to find overnight capital cost of Zn-Mn SSFB with 4.2v% CB + 11 v% MnO2 as positive semi-solid electrode.
2. Costs of other electrode compositions can be calculated by updating the input files accordingly.
