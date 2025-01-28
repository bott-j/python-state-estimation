#!/usr/bin/env python
"""test_PowerNetwork.py: unit tests for PowerNetwork class."""

# Import built-in modules
import math
import sys
import pytest

# Import custom modules
sys.path.append('../')
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

def test_StateEstimateLine():
    """Verify filterpoint() method deadband filter functionality."""

    # Create power network object
    net = PowerNetwork()

    # Create buses
    busPoc = net.addBus("PoC Bus", 330, True)
    busHv1 = net.addBus("HV Bus 1", 330)
    
    # Create line
    Sbase = 100
    Zbase = (330000)**2/(Sbase*1000000)
    net.addLine("Line 1", busPoc, busHv1, 0, 0, 1/(-2*math.pi*1*50/Zbase), 0)

    # Add measurements
    net.addMeasurement("PoC", "V", 1, busPoc)
    net.addMeasurement("PoC", "P", 0.9, busPoc, busHv1)
    net.addMeasurement("PoC", "Q", 0.3, busPoc, busHv1)
    
    # Estimate states
    states = net.stateEstimate()

    # Correct type
    assert(type(states) == dict)
    # Correct length
    assert(len(states) == 2)

    # Check for expected states
    Expected0 = {'name': 'PoC Bus', 'vnominal': 330, 'magnitude': 1.0, 'phase': 0.0}
    assert(states[0] == pytest.approx(Expected0, rel=1e-6, abs=1e-12))
    Expected1 = {'name': 'HV Bus 1', 'vnominal': 330, 'magnitude': 0.9496367096114429, 'phase': -0.2769313269844059}
    assert(states[1] == pytest.approx(Expected1, rel=1e-6, abs=1e-12))

    return
