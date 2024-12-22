# python-state-estimation
Power system state estimation examples and power flow in Python using the Gauss Newton algorithm.

# Example program

The example program "stateEstimator.py" estimates bus voltage and phase and power flows for a simple three-bus network with a line and transformer. The state estimation and power flow algorithms are implemented in the stateEstimate() and powerFlow() methods of the PowerNetwork class defined in the ".estimator" package. 

Network parameters are added programmatically to the model defined in the PowerNetwork object as per below.

```
# Import custom packages
from estimator import PowerNetwork

# Create a power network object
net = PowerNetwork()

# Add buses to the network
bus1 = net.addBus("POI Bus", 275, slack = True)
bus2 = net.addBus("MV Tx Bus", 33)
bus3 = net.addBus("MV Load Bus", 33)

# Base values for per-unit system
Vbase = 275
Sbase = 250
Zbase = (Vbase*1000)**2/(Sbase*1000000)

# Base impedances
Zbase1 = (33*1000)**2/(Sbase*1000000)
# ^ actually for transformer in pscad
Zbase2 = (33*1000)**2/(1*1000000)

# Add a transformer between buses
net.addTransformer(
    "MV Tx1",
    bus1, # Near bus 
    bus2, # Far bus
    0, # Series conductance
    0, # Shunt conductance
    1/(-0.1 * Zbase1 / Zbase2), # Series suseptance
    0, # Shunt suseptance
    275, # Near nominal voltage (for tap ratio) 
    33, # Far nominal voltage (for tap ratio)
    30) # Phase shift

# Add a line between buses
net.addLine("MV line", bus2, bus3, 0, 0, 1/(-2*math.pi*1*50/Zbase2), 0)

# Add measurements
net.addMeasurement("POI voltage", "V", 275.123/275, bus1)
net.addMeasurement("HV P", "P", 0.9900, bus1, bus2)
net.addMeasurement("HV Q", "Q", 0.3206, bus1, bus2)
net.addMeasurement("MV Q", "P", 0.9909, bus2, bus3)
net.addMeasurement("MV Q", "Q", 0.3113, bus2, bus3)
```

The state (bus voltage magnitude and phase) and power flow solution can then be calculated using the methods implemented in the PowerNetwork class:

```
# Estimate the state (bus voltage magnitudes and angles) of the network
busVoltages = net.stateEstimate()

# Calculate the power flows from the estimated state
powerFlows = net.powerFlow()
```

The resulting output for the three bus network example is:

```
STATE ESTIMATION RESULTS
--------------------------------------------------------------------------------
 Description                    |V|        Phase
--------------------------------------------------------------------------------
 POI Bus                  275.123 V    0.000 deg
 MV Tx Bus                 33.011 V  -30.023 deg
 MV Load Bus               31.493 V  -47.447 deg

POWER FLOW RESULTS
--------------------------------------------------------------------------------
 Description              Near Bus    Far Bus            P             Q
--------------------------------------------------------------------------------
 MV Tx1                          0          1     0.990 MW    0.321 MVAr
 MV line                         1          2     0.991 MW    0.311 MVAr
```
