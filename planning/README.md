Interactive Hodgkin-Huxley Neural Simulator


## What is the Hodgkin-Huxley Model?
The Hodgkin-Huxley model, published in 1952, revolutionized neuroscience by mathematically describing how neurons generate action potentials through voltage-gated ion channels. 

Using precise measurements from the giant axon of a squid, Alan Hodgkin and Andrew Huxley created a series of equations that accurately predicted how neurons generate and propagate electrical signals. 

Their work was so fundamental to our understanding of neural function that it earned them a Nobel Prize in 1963, and their equations remain the gold standard for detailed neural modeling today. This project brings this model to life, letting the user:

    - Visualize membrane potential dynamics and action potential propagation
    - Interact with Na⁺ and K⁺ voltage-gated ion channels
    - Explore the interplay between ionic currents and membrane voltage
    - Experiment with different neuronal parameters and stimulation patterns

# Core Features

    ## Real-time Neural Dynamics

        - Full implementation of HH differential equations
        - Voltage-dependent Na⁺ activation (m), Na⁺ inactivation (h), and K⁺ activation (n) gates
        - Accurate ionic current calculations (Na⁺, K⁺, leak currents)

    ## Interactive Interface

        - Dynamic parameter manipulation during simulation
        - Real-time visualization of:
            - Membrane potential traces
            - Gating variable dynamics
            - Individual ionic currents
            - Phase plane trajectories


    ## Analysis Tools

        - Phase plane analysis for understanding dynamics
        - Spike timing and frequency measurements
        - Parameter sensitivity analysis
        - Data export for external analysis

# Technical Implementation
The project is built in four phases:

    ## Phase 1: Core Development

    [HH Model Implementation] --> [Basic Simulation]


    ## Phase 2: Interface

    [UI Development] --> [Interactive Controls]

    [Real-time Visualization]
    

    ## Phase 3: Features

    [Model Variations] --> [Analysis Tools]

    ## Phase 4: Polish

    [Optimization] --> [Documentation]


# Mathematical Core
The model centers around four key differential equations:

    Membrane Potential (V):
        C_m * dV/dt = I_ext - (I_Na + I_K + I_L)
            where:
                I_Na = g_Na * m³ * h * (V - E_Na)
                I_K = g_K * n⁴ * (V - E_K)
                I_L = g_L * (V - E_L)

    Gating Variables (m, h, n):
        dx/dt = αx(V)(1-x) - βx(V)x, where x ∈ {m,h,n}

