import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

class HodgkinHuxleyModel:
    def __init__(self):
        # Default parameters (at 6.3°C)
        self.params = {
            'C_m': 1.0,      # Membrane capacitance (μF/cm²)
            'g_Na': 120.0,   # Maximum sodium conductance (mS/cm²)
            'g_K': 36.0,     # Maximum potassium conductance (mS/cm²)
            'g_L': 0.3,      # Leak conductance (mS/cm²)
            'E_Na': 55.0,    # Sodium reversal potential (mV)
            'E_K': -72.0,    # Potassium reversal potential (mV)
            'E_L': -49.0,    # Leak reversal potential (mV)
            'dt': 0.01       # Time step (ms)
        }
        
        # Initial conditions
        self.V = -65.0    # Initial membrane potential (mV)
        self.m = 0.05     # Initial m gate value
        self.h = 0.6      # Initial h gate value
        self.n = 0.32     # Initial n gate value

    def alpha_m(self, V):
        return 0.1 * (V + 40.0) / (1.0 - np.exp(-(V + 40.0) / 10.0))

    def beta_m(self, V):
        return 4.0 * np.exp(-(V + 65.0) / 18.0)

    def alpha_h(self, V):
        return 0.07 * np.exp(-(V + 65.0) / 20.0)

    def beta_h(self, V):
        return 1.0 / (1.0 + np.exp(-(V + 35.0) / 10.0))

    def alpha_n(self, V):
        return 0.01 * (V + 55.0) / (1.0 - np.exp(-(V + 55.0) / 10.0))

    def beta_n(self, V):
        return 0.125 * np.exp(-(V + 65.0) / 80.0)

    def simulate(self, duration, I_ext=0):
        # Calculate number of time steps
        steps = int(duration / self.params['dt'])
        
        # Initialize arrays to store results
        time = np.arange(steps) * self.params['dt']
        V_record = np.zeros(steps)
        m_record = np.zeros(steps)
        h_record = np.zeros(steps)
        n_record = np.zeros(steps)
        
        # Set initial values
        V = self.V
        m = self.m
        h = self.h
        n = self.n
        
        # Simulation loop
        for i in range(steps):
            # Calculate current gate parameters
            alpha_m = self.alpha_m(V)
            beta_m = self.beta_m(V)
            alpha_h = self.alpha_h(V)
            beta_h = self.beta_h(V)
            alpha_n = self.alpha_n(V)
            beta_n = self.beta_n(V)
            
            # Update gate variables
            m = m + self.params['dt'] * (alpha_m * (1 - m) - beta_m * m)
            h = h + self.params['dt'] * (alpha_h * (1 - h) - beta_h * h)
            n = n + self.params['dt'] * (alpha_n * (1 - n) - beta_n * n)
            
            # Calculate ionic currents
            I_Na = self.params['g_Na'] * m**3 * h * (V - self.params['E_Na'])
            I_K = self.params['g_K'] * n**4 * (V - self.params['E_K'])
            I_L = self.params['g_L'] * (V - self.params['E_L'])
            
            # Update membrane potential
            dV = (I_ext - I_Na - I_K - I_L) / self.params['C_m']
            V = V + self.params['dt'] * dV
            
            # Store results
            V_record[i] = V
            m_record[i] = m
            h_record[i] = h
            n_record[i] = n
        
        return {
            't': time,
            'V': V_record,
            'm': m_record,
            'h': h_record,
            'n': n_record
        }

def plot_simulation(duration=100, I_ext=10):
    # Create model and run simulation
    model = HodgkinHuxleyModel()
    results = model.simulate(duration, I_ext)
    
    # Create figure with secondary y-axis
    fig = make_subplots(
        rows=2, cols=1,
        subplot_titles=('Membrane Potential', 'Gating Variables'),
        vertical_spacing=0.15
    )
    
    # Add traces for membrane potential
    fig.add_trace(
        go.Scatter(
            x=results['t'],
            y=results['V'],
            name='Membrane Potential',
            line=dict(color='royalblue')
        ),
        row=1, col=1
    )
    
    # Add traces for gating variables
    fig.add_trace(
        go.Scatter(
            x=results['t'],
            y=results['m'],
            name='m (Na⁺ activation)',
            line=dict(color='red')
        ),
        row=2, col=1
    )
    
    fig.add_trace(
        go.Scatter(
            x=results['t'],
            y=results['h'],
            name='h (Na⁺ inactivation)',
            line=dict(color='green')
        ),
        row=2, col=1
    )
    
    fig.add_trace(
        go.Scatter(
            x=results['t'],
            y=results['n'],
            name='n (K⁺ activation)',
            line=dict(color='purple')
        ),
        row=2, col=1
    )
    
    # Update layout
    fig.update_layout(
        height=800,
        showlegend=True,
        title_text="Hodgkin-Huxley Model Simulation",
    )
    
    # Update axes labels
    fig.update_xaxes(title_text="Time (ms)", row=2, col=1)
    fig.update_yaxes(title_text="Membrane Potential (mV)", row=1, col=1)
    fig.update_yaxes(title_text="Gating Value", row=2, col=1)
    
    return fig

# Run simulation ; display plot
if __name__ == "__main__":
    fig = plot_simulation(duration=100, I_ext=10)
    fig.show()