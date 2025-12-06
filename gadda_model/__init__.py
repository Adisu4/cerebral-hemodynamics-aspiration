"""
Intracranial pressure model package

Computational model of intracranial pressure dynamics and cerebral venous outflow.

Core modules:
    parameters - Model parameter definitions
    equations  - ODE system for the 16-state model
    solver     - High-level simulation interface
    aspiration - Venous aspiration protocol definitions

Basic usage:
    >>> from gadda_model import ModelParameters, run_simulation
    >>> params = ModelParameters()
    >>> result = run_simulation(params, duration_s=900)
    >>> print(f"Final ICP: {result.ICP[-1]:.1f} mmHg")
"""

from gadda_model.parameters import ModelParameters
from gadda_model.equations import (
    gadda_ode_system_FIXED as gadda_ode_system,
    jugular_conductance,
    cerebral_arterioles_resistance,
    cerebral_veins_resistance,
    cerebral_veins_capacitance
)
from gadda_model.solver import (
    SimulationResult,
    run_simulation,
    run_full_protocol
)
from gadda_model.aspiration import AspirationProtocol

__version__ = '1.0.0'

__all__ = [
    # Parameters
    'ModelParameters',
    
    # Equations
    'gadda_ode_system',
    'jugular_conductance',
    'cerebral_arterioles_resistance',
    'cerebral_veins_resistance',
    'cerebral_veins_capacitance',
    
    # Solver
    'SimulationResult',
    'run_simulation',
    'run_full_protocol',
    
    # Aspiration
    'AspirationProtocol'
]
