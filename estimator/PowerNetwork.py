#!/usr/bin/env python
"""PowerNetwork.py: Power system estimation class."""

# Import built-in packages
import collections
import math

# Import third-party packages
import numpy as np

# Authorship information
__author__ = "James Bott"
__copyright__ = "Copyright 2024, James Bott"
__credits__ = ["James Bott"]
__license__ = "MIT"
__version__ = "0.0.1"
__maintainer__ = "James Bott"
__email__ = "https://github.com/bott-j"
__status__ = "Development"

# Named tupple holds a bus definition
Bus = collections.namedtuple('Bus',
                            ['id',
                             'name',
                             'vNominal',
                             'slack'])

# Named tupple holds a equipment definition
Equipment = collections.namedtuple('Equipment',
                            ['id',
                             'name',
                             'busNear',
                             'busFar',
                             'gSeries',
                             'gShunt',
                             'bSeries',
                             'bShunt',
                             'voltageNear',
                             'voltageFar',
                             'phaseShift'])
        
# Named tupple holds a equipment definition
Measurement = collections.namedtuple('Measurement',
                            ['name',
                             'type',
                             'value',
                             'busNear',
                             'busFar'])

""" The PowerNetwork class provides a model of the network and
    automates functions for state estimation and power flow calculation.
"""
class PowerNetwork():

    def __init__(
        self,
        maxIterations = 1000,
        tolerance = 0.001
    ):
        """ Constructor for the power network model object. """

        # Maximum number of iterations for Gauss-Newton state estimation
        self._maxIterations = maxIterations
        # Voltage tolerance for convergance of Gauss-Newton state estimation
        self._tolerance = tolerance

        # Index of the slack bus
        self._slackBusIndex = None 

        # Holds list of buses
        self._buses = list()
        # Holds list of equipment (line, transformer) connecting buses
        self._equipment = list()
        # Measurements
        self._measurements = list()
        # Voltage magnitude states as complex quantity 
        self._states = None

        return

    def addMeasurement(
        self, 
        name : str, # Name of this measurement
        type : str, # Type of measurement as 'V', 'P' or 'Q'
        value : float, # Value for this measurement
        busNear : Bus, # The near bus object
        busFar : Bus = None # The far bus object
        ):
        """ Adds a new measurement value. """

        # Create measurement named tuple 
        measurement = Measurement(name, type, value, busNear, busFar)
        
        # Add to internal list
        self._measurements += [measurement]
        
        return

    def addBus(
        self, 
        name : str, # Name of the bus
        vNominal : float, # Nominal base voltage of the bus 
        slack : bool = False # If true, this is the slack bus
        ):
        """ Adds a bus to the network model. """

        # Create a new bus named tuple
        bus = Bus(len(self._buses), name, vNominal, slack)

        # Will need to know later which bus is the slack bus
        if(slack):
            self._slackBusIndex = bus.id

        # Add the new bus to the list
        self._buses += [bus]
        
        return bus
    
    def addLine(
            self, 
            name : str, # Name of the line
            busNear : Bus, # Near bus named tuple
            busFar : Bus, # Far bus named tuple
            gSeries : float, # Series conductance
            gShunt : float, # Shunt conductance
            bSeries : float, # Series susceptance
            bShunt : float # Shunt susceptance
            ):
        """ Adds a line joining buses to the network model. """
                
        # Create a new equipment object
        eq = Equipment(len(self._equipment), name, busNear, busFar, gSeries, gShunt, bSeries, bShunt, busNear.vNominal, busFar.vNominal, 0)
        
        # Add to list
        self._equipment += [eq]

        return eq

    def addTransformer(
            self, 
            name : str, # Name of the transformer
            busNear : Bus,  # Near bus object
            busFar : Bus, # Far bus object
            gSeries : float, # Series conductance
            gShunt : float, # Shunt conductance
            bSeries : float, # Series susceptance
            bShunt : float, # Shunt susceptance 
            vNear : float, # Nominal voltage of near bus for tap ratio
            vFar : float, # Nominal voltage of far bus for tap ratio
            phaseShift : float = 30 # Phase shift of transformer
            ):
        """ Adds a transformer joining buses to the network model. """
        
        # Create a new equipment object
        eq = Equipment(len(self._equipment), name, busNear, busFar, gSeries, gShunt, bSeries, bShunt, vNear, vFar, phaseShift)
        
        # Add to list
        self._equipment += [eq]
        
        return eq

    def _dPijVk(self, vs, ts, gs, bs, gss, tau, i, j, k):
        """ Calculates partial derivative of real power measurement to voltage magnitude state. """ 

        result = 0
        if(k == i):
            result = -vs[i]*(gs[i][j]*math.cos(tau[i][j] + ts[i] - ts[j]) + bs[i][j]*math.sin(tau[i][j] + ts[i] - ts[j])) + 2*(gs[i][j] + gss[i])*vs[i]
        if(k == j):
            result = -vs[i]*(gs[i][j]*math.cos(tau[i][j] + ts[i] - ts[j]) + bs[i][j]*math.sin(tau[i][j] + ts[i] - ts[j])) 
        return result

    def _dPijTk(self, vs, ts, gs, bs, tau, i, j, k):
        """ Calculates partial derivative of real measurement to phase magnitude state. """ 
        result = 0
        if(k == i):
            result = vs[i]*vs[j]*(gs[i][j]*math.sin(tau[i][j] + ts[i] - ts[j]) - bs[i][j]*math.cos(tau[i][j] + ts[i] - ts[j]))
        if(k == j):
            result = -vs[i]*vs[j]*(gs[i][j]*math.sin(tau[i][j] + ts[i] - ts[j]) - bs[i][j]*math.cos(tau[i][j] + ts[i] - ts[j])) 
        return result

    def _dQijVk(self, vs, ts, gs, bs, bss, tau, i, j, k):
        """ Calculates partial derivative of reactive power measurement to voltage magnitude state. """ 

        result = 0
        if(k == i):
            result = -vs[i]*(gs[i][j]*math.sin(tau[i][j] + ts[i] - ts[j]) - bs[i][j]*math.cos(tau[i][j] + ts[i] - ts[j])) - 2*(bs[i][j] + bss[i])*vs[i]
        if(k == j):
            result = -vs[i]*(gs[i][j]*math.sin(tau[i][j] + ts[i] - ts[j]) - bs[i][j]*math.cos(tau[i][j] + ts[i] - ts[j])) 
        return result

    def _dQijTk(self, vs, ts, gs, bs, tau, i, j, k):
        """ Calculates partial derivative of reactive power measurement to voltage magnitude phase. """ 

        result = 0
        if(k == i):
            result = -vs[i]*vs[j]*(gs[i][j]*math.cos(tau[i][j] + ts[i] - ts[j]) + bs[i][j]*math.sin(tau[i][j] + ts[i] - ts[j]))
        if(k == j):
            result = vs[i]*vs[j]*(gs[i][j]*math.cos(tau[i][j] + ts[i] - ts[j]) + bs[i][j]*math.sin(tau[i][j] + ts[i] - ts[j])) 
        return result

    def _Pij(self, vs, ts, gs, bs, gss, tau, i, j):
        """ Calculates real power flow from state. """ 
        return vs[i]**2*(gss[i]+gs[i][j])-vs[i]*vs[j]*(gs[i][j]*math.cos(tau[i][j] + ts[i]-ts[j])+bs[i][j]*math.sin(tau[i][j] + ts[i]-ts[j]))

    def _Qij(self, vs, ts, gs, bs, bss, tau, i, j):
        """ Calculates reactive power flow from state. """ 

        return -vs[i]**2*(bss[i]+bs[i][j])-vs[i]*vs[j]*(gs[i][j]*math.sin(tau[i][j] + ts[i]-ts[j])-bs[i][j]*math.cos(tau[i][j] + ts[i]-ts[j]))

    def _buildJacobian(self, vs, ts, gs, gss, bs, bss, tau):
        """ Build the Jacobian matrix of partial derivatives used in state estimation. """

        # Number of measurements
        rows = len(self._measurements)

        # States include |V| and theta_V for each bus
        nBuses = len(self._buses)
        cols = 2*nBuses
        
        # Create an empty Jacobian
        J = np.zeros([rows, cols])

        # For each measurement
        for m in range(0, rows):
            # For each pair (magnitude and angle) of state variables
            for k in range(0, nBuses):
                
                # If this is a real power measurement
                if(self._measurements[m].type == 'V'
                  and self._measurements[m].busNear.id == k):
                    J[m][k] = 1
                if(self._measurements[m].type == 'P'):                    
                    # The voltage magnitude partial derivative
                    J[m][k] = self._dPijVk(
                        vs,
                        ts,
                        gs,
                        bs,
                        gss,
                        tau,
                        self._measurements[m].busNear.id, 
                        self._measurements[m].busFar.id,
                        k)
                    # The voltage angle partial derivative
                    J[m][k+nBuses] = self._dPijTk(
                        vs,
                        ts,
                        gs,
                        bs,
                        tau,
                        self._measurements[m].busNear.id, 
                        self._measurements[m].busFar.id,
                        k)
                # If this is a reactive power measurement
                elif(self._measurements[m].type == 'Q'):
                    # The voltage magnitude partial derivative
                    J[m][k] = self._dQijVk(
                        vs, 
                        ts, 
                        gs, 
                        bs, 
                        bss, 
                        tau,
                        self._measurements[m].busNear.id, 
                        self._measurements[m].busFar.id,
                        k)
                    # The voltage angle partial derivative
                    J[m][k+nBuses] = self._dQijTk(
                        vs,
                        ts,
                        gs,
                        bs,
                        tau,
                        self._measurements[m].busNear.id, 
                        self._measurements[m].busFar.id,
                        k)
        return J

    def _estimateMeasurements(self, x, vs, ts, gs, bs, gss, bss, tau):
        """ The measurement function h(.) estimates the measurements from the state. """
        h_vk = list()

        # For each measurement
        for measurement in self._measurements:
            if(measurement.type == 'V'):
                h_vk += [x[measurement.busNear.id]]
            elif(measurement.type == 'P'):
                h_vk += [self._Pij(vs, ts, gs, bs, gss, tau, measurement.busNear.id, measurement.busFar.id)]
            elif(measurement.type == 'Q'):
                h_vk += [self._Qij(vs, ts, gs, bs, bss, tau, measurement.busNear.id, measurement.busFar.id)] 
        
        return np.asarray(h_vk)

    def stateEstimate(self):
        """ Perform state estimation of bus voltage magnitudes and phases. """
        x_next = None
        nBuses = len(self._buses)

        # Initialize states as flat start, with voltage magnitudes have a value of 1
        # and phases having a value of 0
        x = np.append(np.ones([nBuses,1]), np.zeros([nBuses,1]), axis=0) 
        states = x.T
        
        # Calculate network parameters
        gs, bs, gss, bss, tau = self._calculateNetParameters()

        # Setup measurement matrix
        z = np.transpose(np.asarray([[measurement.value for measurement in self._measurements]]))

        # Solve nonlinear system of equations using Gauss-Newton method 
        for i in range (0, self._maxIterations):
            vs = x[:nBuses] # voltage slice of the state
            ts = x[nBuses:] # phase slice of the state

            J = self._buildJacobian(vs, ts, gs, gss, bs, bss, tau)

            # Estimate measurement from state
            h_vk = self._estimateMeasurements(x, vs, ts, gs, bs, gss, bss, tau) 

            # Singular value decomposition of the Jacobian
            U,S,Vt = np.linalg.svd(J, full_matrices=False)
            S = np.diag(S)
            V = np.transpose(Vt)
            
            # Apply the Gauss-Newton update equation 
            # Pure GN, no SVD
            #x_next = np.array(x) + np.linalg.inv(np.transpose(J)@J)@np.transpose(J)@(np.array(z) - np.array(h_vk)) 
            # With SVD
            x_next = np.array(x) + V@np.linalg.inv(np.transpose(S)@S)@np.transpose(S)@np.transpose(U)@(np.array(z) - np.array(h_vk))
            
            # Calculate the maximum change in state for convergence evaluation
            if(x_next is not None):
                max_change = max(abs((x_next - x)))

            x = x_next

            # Save state iterations for plotting
            states = np.append(states, x.T, axis=0)

            # Error threshold for convergence
            if(max_change <= self._tolerance):
                break

        # Offset the slack bus to zero phase 
        phaseOffset = x[self._slackBusIndex + nBuses][0]
        
        # Return results as dictionary
        results = {
            bus.id : 
                {
                    "name" : bus.name, 
                    "vnominal" : bus.vNominal,
                    "magnitude" : x[i][0], 
                    "phase" : x[i + nBuses][0] - phaseOffset
                } 
                for i, bus in zip(range(0, nBuses), self._buses)
            }

        self._states = x

        return results

    def _calculateNetParameters(self):
        """ Calculates network parameters. """
        nBuses = len(self._buses)

        # Setup impedance related matrices
        gs = np.zeros([nBuses, nBuses])
        bs = np.zeros([nBuses, nBuses])
        gss = [0]*nBuses
        bss = [0]*nBuses
        # Setup phase offset matrix
        tau = np.zeros([nBuses, nBuses])

        for eq in self._equipment:
            tau[eq.busNear.id][eq.busFar.id] = -math.radians(eq.phaseShift)
            tau[eq.busFar.id][eq.busNear.id] = math.radians(eq.phaseShift)
            gs[eq.busNear.id][eq.busFar.id] = eq.gSeries
            gs[eq.busFar.id][eq.busNear.id] = eq.gSeries
            bs[eq.busNear.id][eq.busFar.id] = eq.bSeries
            bs[eq.busFar.id][eq.busNear.id] = eq.bSeries

        return gs, bs, gss, bss, tau 

    def powerFlow(self):
        """ Calculates power flow for network model. """
        nBuses = len(self._buses) 

        # Check states have been estimated
        if(self._states is None):
            raise ValueError("State estimation must be run before power flow.")

        # Take slices of state
        vs = self._states[:nBuses] # voltage slice of the state
        ts = self._states[nBuses:] # phase slice of the state
            
        # Calculate network parameters
        gs, bs, gss, bss, tau = self._calculateNetParameters()

        # Return results as dictionary
        results = {
            eq.id : 
                {
                    "name" : eq.name, 
                    "nearbus" : eq.busNear.id,
                    "farbus" : eq.busFar.id,
                    "nearbusname" : eq.busNear.name,
                    "farbusname" : eq.busFar.name,
                    "p" : self._Pij(vs, ts, gs, bs, gss, tau, eq.busNear.id, eq.busFar.id)[0], 
                    "q" : self._Qij(vs, ts, gs, bs, bss, tau, eq.busNear.id, eq.busFar.id)[0]
                } 
                for eq in self._equipment
            }

        return results


