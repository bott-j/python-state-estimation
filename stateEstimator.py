# Import built-in packages
import math

# Import custom packages
from estimator import PowerNetwork

# Authorship information
__author__ = "James Bott"
__copyright__ = "Copyright 2024, James Bott"
__credits__ = ["James Bott"]
__license__ = "MIT"
__version__ = "0.0.1"
__maintainer__ = "James Bott"
__email__ = "https://github.com/bott-j"
__status__ = "Development"

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

# Estimate the state (bus voltage magnitudes and angles) of the network
busVoltages = net.stateEstimate()

# Calculate the power flows from the estimated state
powerFlows = net.powerFlow()

# Display the state estimates for the network
print("STATE ESTIMATION RESULTS")
print('-'*80)
print(' Description', ' ' * 18, '|V|', ' ' * 6, 'Phase')
print('-'*80)
for bus in busVoltages.values():
    print(" {:25}{:7.3f} V {:8.3f} deg".format(bus["name"], bus["vnominal"] * bus["magnitude"], 180*bus["phase"]/math.pi))

# Display the power flow estimates for the network
print("")
print("POWER FLOW RESULTS")
print('-'*80)
print(' Description', ' ' * 12, 'Near Bus' , ' ' * 2, 'Far Bus', ' ' * 10, 'P', ' ' * 11, 'Q')
print('-'*80)
for eq in powerFlows.values():
    print(" {:22}{:11}{:11}{:10.3f} MW {:8.3f} MVAr".format(eq["name"], eq["nearbus"], eq["farbus"], eq["p"], eq["q"]))
