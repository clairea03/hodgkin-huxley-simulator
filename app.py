import streamlit as st
import numpy as np
import time
from core_model import HodgkinHuxleyModel
import plotly.graph_objects as go
from plotly.subplots import make_subplots

st.set_page_config(page_title="Hodgkin-Huxley Model", layout="wide")

st.title("Interactive Hodgkin-Huxley Neuron Model")
st.markdown("""
This app simulates the Hodgkin-Huxley model of neuronal action potentials.
Adjust the parameters on the sidebar and see how they affect the neuron's behavior.
""")

st.sidebar.header("Model Parameters")

# Simulation parameters
st.sidebar.subheader("Simulation Settings")
st.sidebar.info("Controls the basic simulation parameters.")
duration = st.sidebar.slider("Duration (ms)", 10, 500, 100, 
                           help="Total time to simulate in milliseconds")
I_ext = st.sidebar.slider("Stimulus Current (μA/cm²)", 0.0, 50.0, 10.0, 
                         help="External current applied to the neuron (higher values trigger action potentials)")
animation_speed = st.sidebar.slider("Animation Speed", 1, 10, 5, 
                                  help="Controls how fast the simulation animation plays")

# Membrane properties
st.sidebar.subheader("Membrane Properties")
st.sidebar.info("The fundamental electrical properties of the neuron's cell membrane.")
C_m = st.sidebar.slider("Membrane Capacitance (μF/cm²)", 0.1, 5.0, 1.0, 
                       help="Represents the cell membrane's ability to store charge (higher values slow down voltage changes)")

# Channel conductances
st.sidebar.subheader("Maximum Conductances")
st.sidebar.info("Controls how readily ions can flow through each type of channel when fully open.")
g_Na = st.sidebar.slider("Sodium (mS/cm²)", 0.0, 200.0, 120.0, 
                        help="Maximum conductance of sodium channels (higher values increase spike amplitude)")
g_K = st.sidebar.slider("Potassium (mS/cm²)", 0.0, 100.0, 36.0, 
                       help="Maximum conductance of potassium channels (higher values speed up repolarization)")
g_L = st.sidebar.slider("Leak (mS/cm²)", 0.0, 1.0, 0.3, 0.01, 
                       help="Conductance of non-gated leak channels (affects resting potential)")

# Reversal potentials
st.sidebar.subheader("Reversal Potentials")
st.sidebar.info("The equilibrium voltage for each ion type, determined by concentration gradients.")
E_Na = st.sidebar.slider("Sodium (mV)", 20.0, 80.0, 55.0, 
                        help="Sodium reversal potential (more positive values increase spike height)")
E_K = st.sidebar.slider("Potassium (mV)", -100.0, -40.0, -72.0, 
                       help="Potassium reversal potential (more negative values lower resting potential)")
E_L = st.sidebar.slider("Leak (mV)", -80.0, -20.0, -49.0, 
                       help="Leak reversal potential (affects resting membrane potential)")

# Initial conditions
st.sidebar.subheader("Initial Conditions")
st.sidebar.info("Starting values for all variables at the beginning of the simulation.")
V_init = st.sidebar.slider("Membrane Potential (mV)", -80.0, 0.0, -65.0, 
                          help="Initial voltage across the membrane (typically around -65mV when at rest)")
m_init = st.sidebar.slider("m gate", 0.0, 1.0, 0.05, 
                          help="Initial value of m (Na+ activation gate), normally low at rest")
h_init = st.sidebar.slider("h gate", 0.0, 1.0, 0.6, 
                          help="Initial value of h (Na+ inactivation gate), normally high at rest")
n_init = st.sidebar.slider("n gate", 0.0, 1.0, 0.32, 
                          help="Initial value of n (K+ activation gate), intermediate value at rest")

# Initialize model with custom parameters
model = HodgkinHuxleyModel()
model.params.C_m = C_m
model.params.g_Na = g_Na
model.params.g_K = g_K
model.params.g_L = g_L
model.params.E_Na = E_Na
model.params.E_K = E_K
model.params.E_L = E_L

# Set initial conditions
model.V = V_init
model.m = m_init
model.h = h_init
model.n = n_init

# Initialize model with custom parameters
model = HodgkinHuxleyModel()
model.params.C_m = C_m
model.params.g_Na = g_Na
model.params.g_K = g_K
model.params.g_L = g_L
model.params.E_Na = E_Na
model.params.E_K = E_K
model.params.E_L = E_L

# Set initial conditions
model.V = V_init
model.m = m_init
model.h = h_init
model.n = n_init

# Run simulation
results = model.simulate(duration, I_ext)

# Create placeholders for figures
plot_placeholder = st.empty()
status_text = st.empty()

# Set up plot
def create_figure(results, time_window=None):
    fig = make_subplots(
        rows=3, cols=1,
        subplot_titles=('Membrane Potential', 'Gating Variables', 'Ionic Currents'),
        vertical_spacing=0.1,
        row_heights=[0.4, 0.3, 0.3]
    )
    
    if time_window is not None:
        start_idx, end_idx = time_window
        t = results['t'][start_idx:end_idx]
        V = results['V'][start_idx:end_idx]
        m = results['m'][start_idx:end_idx]
        h = results['h'][start_idx:end_idx]
        n = results['n'][start_idx:end_idx]
        I_Na = results['I_Na'][start_idx:end_idx]
        I_K = results['I_K'][start_idx:end_idx]
        I_L = results['I_L'][start_idx:end_idx]
        x_range = [t[0], t[-1]]
    else:
        t = results['t']
        V = results['V']
        m = results['m']
        h = results['h']
        n = results['n']
        I_Na = results['I_Na']
        I_K = results['I_K']
        I_L = results['I_L']
        x_range = [0, duration]
    
    # Add traces for membrane potential
    fig.add_trace(
        go.Scatter(
            x=t,
            y=V,
            name='Membrane Potential',
            line=dict(color='royalblue', width=2)
        ),
        row=1, col=1
    )
    
    # Add threshold line at 0 mV
    fig.add_trace(
        go.Scatter(
            x=x_range,
            y=[0, 0],
            name='Threshold',
            line=dict(color='gray', width=1, dash='dash')
        ),
        row=1, col=1
    )
    
    # Add traces for gating variables
    fig.add_trace(
        go.Scatter(
            x=t,
            y=m,
            name='m (Na⁺ activation)',
            line=dict(color='red')
        ),
        row=2, col=1
    )
    
    fig.add_trace(
        go.Scatter(
            x=t,
            y=h,
            name='h (Na⁺ inactivation)',
            line=dict(color='green')
        ),
        row=2, col=1
    )
    
    fig.add_trace(
        go.Scatter(
            x=t,
            y=n,
            name='n (K⁺ activation)',
            line=dict(color='purple')
        ),
        row=2, col=1
    )
    
    # Add traces for ionic currents
    fig.add_trace(
        go.Scatter(
            x=t,
            y=I_Na,
            name='I_Na',
            line=dict(color='orange')
        ),
        row=3, col=1
    )
    
    fig.add_trace(
        go.Scatter(
            x=t,
            y=I_K,
            name='I_K',
            line=dict(color='cyan')
        ),
        row=3, col=1
    )
    
    fig.add_trace(
        go.Scatter(
            x=t,
            y=I_L,
            name='I_L',
            line=dict(color='brown')
        ),
        row=3, col=1
    )
    
    # Update layout
    fig.update_layout(
        height=800,
        showlegend=True,
        title_text="Hodgkin-Huxley Model Simulation Results"
    )
    
    # Update axes labels
    fig.update_xaxes(title_text="Time (ms)", row=3, col=1)
    fig.update_yaxes(title_text="Membrane Potential (mV)", row=1, col=1)
    fig.update_yaxes(title_text="Gating Value", row=2, col=1)
    fig.update_yaxes(title_text="Current (μA/cm²)", row=3, col=1)
    
    return fig

fig = create_figure(results)
plot_placeholder.plotly_chart(fig, use_container_width=True)

# Basic explanation of the model
with st.expander("About the Hodgkin-Huxley Model"):
    st.markdown("""
    ### The Hodgkin-Huxley Model
    
    The Hodgkin-Huxley model describes how action potentials in neurons are initiated and propagated. It was first introduced by Alan Hodgkin and Andrew Huxley in 1952 to explain the electrical properties of squid giant axons.
    
    #### Key Components:
    
    1. **Membrane Potential (V)**: The voltage difference across the cell membrane.
    
    2. **Ion Channels**:
       - **Sodium (Na⁺)**: Controlled by activation (m) and inactivation (h) gates
       - **Potassium (K⁺)**: Controlled by activation (n) gates
       - **Leak**: Represents passive ion flow
    
    3. **Gating Variables**:
       - **m**: Na⁺ activation gate (opens quickly)
       - **h**: Na⁺ inactivation gate (closes slowly)
       - **n**: K⁺ activation gate (opens slowly)
    
    #### Action Potential Phases:
    
    1. **Resting**: Membrane potential is around -65 mV
    2. **Depolarization**: Na⁺ channels open, causing rapid rise in voltage
    3. **Repolarization**: K⁺ channels open and Na⁺ channels inactivate, bringing voltage back down
    4. **Hyperpolarization**: K⁺ channels remain open, driving potential below resting
    5. **Recovery**: Return to resting potential
    """)

st.markdown("---")
st.caption("Created by Claire Alverson | View source code at https://github.com/clairea03/hodgkin-huxley-simulator.git")