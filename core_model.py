import numpy as np
from dataclasses import dataclass
from typing import Dict, NamedTuple

@dataclass # Holds all the basic properties of our neuron

#Parameters:
class HHParameters:
    C_m: float = 1.0      # Membrane capacitance (μF/cm²)
    g_Na: float = 120.0   # Maximum sodium conductance (mS/cm²)
    g_K: float = 36.0     # Maximum potassium conductance (mS/cm²)
    g_L: float = 0.3      # Leak conductance (mS/cm²)
    E_Na: float = 55.0    # Sodium reversal potential (mV)
    E_K: float = -72.0    # Potassium reversal potential (mV)
    E_L: float = -49.0    # Leak reversal potential (mV)
    dt: float = 0.01      # Time step (ms)

# Dynamic variables:
class StateVector(NamedTuple):
    V: float      # Membrane potential
    m: float      # Na+ activation gate
    h: float      # Na+ inactivation gate
    n: float      # K+ activation gate
    I_Na: float   # Sodium current
    I_K: float    # Potassium current
    I_L: float    # Leak current

# Core model:
class HodgkinHuxleyModel:
    def __init__(self):
        # Default parameters
        self.params = HHParameters()
        
        # Initial conditions
        self.V = -65.0    # Initial membrane potential (mV) ; neuron at rest
        self.m = 0.05     # Initial m gate value
        self.h = 0.6      # Initial h gate value
        self.n = 0.32     # Initial n gate value

# The following functions control how fast the Na+ activation gates open and close
# This is essentially how the gates respond to voltage changes
    @staticmethod
    def alpha_m(V: float) -> float:
        return 0.1 * (V + 40.0) / (1.0 - np.exp(-(V + 40.0) / 10.0))

    @staticmethod
    def beta_m(V: float) -> float:
        return 4.0 * np.exp(-(V + 65.0) / 18.0)

# Similar functions for h (Na+ inactivation gate)
    @staticmethod
    def alpha_h(V: float) -> float:
        return 0.07 * np.exp(-(V + 65.0) / 20.0)

    @staticmethod
    def beta_h(V: float) -> float:
        return 1.0 / (1.0 + np.exp(-(V + 35.0) / 10.0))
    
# And n (K+ activation gate)...
    @staticmethod
    def alpha_n(V: float) -> float:
        return 0.01 * (V + 55.0) / (1.0 - np.exp(-(V + 55.0) / 10.0))

    @staticmethod
    def beta_n(V: float) -> float:
        return 0.125 * np.exp(-(V + 65.0) / 80.0)


# Calculates how much current is flowing through each type of ion channel
    def ionic_currents(self, V: float, m: float, h: float, n: float) -> tuple:
       
        # Sodium current (Na+) - flows into the cell
        I_Na = self.params.g_Na * m**3 * h * (V - self.params.E_Na)

        # Potassium current (K+) - flows out of the cell
        I_K = self.params.g_K * n**4 * (V - self.params.E_K)

        # Leak current - a general trickle of ions
        I_L = self.params.g_L * (V - self.params.E_L)

        return I_Na, I_K, I_L

# Updates the position of all gates based on the current voltage
    def update_states(self, V: float, m: float, h: float, n: float) -> tuple:
        
        # Calculate rate constants
        am, bm = self.alpha_m(V), self.beta_m(V)
        ah, bh = self.alpha_h(V), self.beta_h(V)
        an, bn = self.alpha_n(V), self.beta_n(V)

        # Update gating variables
        m_new = m + self.params.dt * (am * (1 - m) - bm * m)
        h_new = h + self.params.dt * (ah * (1 - h) - bh * h)
        n_new = n + self.params.dt * (an * (1 - n) - bn * n)

        return m_new, h_new, n_new

# Takes one step forward in time
    def step(self, I_ext: float = 0) -> StateVector:
        # Calculate ionic currents
        I_Na, I_K, I_L = self.ionic_currents(self.V, self.m, self.h, self.n)

        # Update membrane potential
        dV = (I_ext - I_Na - I_K - I_L) / self.params.C_m
        V_new = self.V + self.params.dt * dV

        # Update gating variables
        m_new, h_new, n_new = self.update_states(self.V, self.m, self.h, self.n)

        # Update state
        self.V, self.m, self.h, self.n = V_new, m_new, h_new, n_new

        return StateVector(V_new, m_new, h_new, n_new, I_Na, I_K, I_L)

# Runs the full simulation for a specified amount of time
    #duration: Simulation duration (in ms)
    #I_ext: External current being applied in μA/cm²
            

    def simulate(self, duration: float, I_ext: float = 0) -> Dict[str, np.ndarray]:
       
        steps = int(duration / self.params.dt)
        
        # Start w/ empty arrays for better performance
        time = np.arange(steps) * self.params.dt
        V = np.zeros(steps)
        m = np.zeros(steps)
        h = np.zeros(steps)
        n = np.zeros(steps)
        I_Na = np.zeros(steps)
        I_K = np.zeros(steps)
        I_L = np.zeros(steps)

        # Run the simulation
        for i in range(steps):
            state = self.step(I_ext)
            V[i] = state.V
            m[i] = state.m
            h[i] = state.h
            n[i] = state.n
            I_Na[i] = state.I_Na
            I_K[i] = state.I_K
            I_L[i] = state.I_L

# Return recordings
        return {
            't': time,
            'V': V,
            'm': m,
            'h': h,
            'n': n,
            'I_Na': I_Na,
            'I_K': I_K,
            'I_L': I_L
        }
