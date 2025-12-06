"""
Model parameter definitions

Core parameter set for the intracranial and venous system model.
"""

import numpy as np

class ModelParameters:
    def __init__(self, posture='supine'):
        self.posture = posture
        
        # --- INTRACRANIAL (Table 1) ---
        self.Cpan = 0.205
        self.Cpa1 = 2.87
        self.Cpa2 = 0.164
        self.Pa = 100.0
        self.Picn = 9.5
        self.Ppa = 58.9
        self.Pv = 14.1
        self.Pv1 = -2.5
        self.Pvs = 6.0
        self.Qn = 12.5
        self.R0 = 526.3
        self.Rf = 2.38e3
        self.Rla = 0.6
        self.Rpv = 0.880
        self.Rvs1 = 0.366
        self.Gaut = 3.0
        self.tau_aut = 20.0
        self.xaut = 2.16e-4
        self.kE = 0.077
        self.kR = 13.1e3
        self.kven = 0.155
        
        # --- EXTRACRANIAL CAPACITIES (Table 4) ---
        self.Cvs = 0.5
        self.Cjr3 = 1.0
        self.Cjl3 = 1.0
        self.Cjr2 = 2.5
        self.Cjl2 = 2.5
        # Table 4 does not include Cjr1/Cjl1; J1 is modeled as resistive, not capacitive.
        
        self.Cc3 = 0.7
        self.Cc2 = 1.4
        self.Csvc = 20.0
        self.Cazy = 0.5
        self.Cvv = 0.5
        
        # --- CONDUCTANCE PARAMETERS (published baseline values) ---
        # Baseline conductances from the reference publication
        self.kjr3 = 11.0   # Right J3 baseline conductance
        self.kjl3 = 11.0   # Left J3 baseline conductance
        self.kjr2 = 13.0   # Right J2 conductance
        self.kjl2 = 13.0   # Left J2 conductance
        self.kjr1 = 6.9    # Right J1 conductance
        self.kjl1 = 6.9    # Left J1 conductance
        
        # Starling resistor collapse sensitivity (published)
        self.A = 0.8       # Published value for Starling collapse
        
        self.G0 = 1.0 / self.R0
        self.Gvs = 1.0 / self.Rvs1
        self.Gex = 5.0 / (100.0 - 6.0)
        
        # Collaterals (Table 3 derived)
        if posture == 'supine':
            self.Gc3 = 0.0
        else:
            self.Gc3 = 3.5
        
        self.Gc2 = 6.0  / (6.0 - 5.85)
        self.Gc1 = 5.0 / (5.85 - 5.0)
        
        # Collateral conductances (REDUCED for acute stenosis modeling)
        # Original values represent fully developed collaterals
        # For IIH stenosis studies, use 10% of original values
        self.Gcjr3 = 0.1 * (1.0 / (6.0 - 5.85))  # ~0.67 instead of 6.67
        self.Gcjl3 = 0.1 * (1.0 / (6.0 - 5.85))
        self.Gcjr2 = 0.1 * (1.0 / (5.85 - 5.70))
        self.Gcjl2 = 0.1 * (1.0 / (5.85 - 5.70))
        
        self.Gvvr = 0.4 / (6.0 - 5.8)
        self.Gvvl = 0.4 / (6.0 - 5.8)
        self.Gvv2 = 0.4 / (5.8 - 5.0)
        
        self.Gazy1 = 0.4 / (5.8 - 5.5)
        self.Glv = 0.13 / (5.5 - 5.0)
        self.Grv = 0.27 / (5.5 - 5.0)
        self.Gazy2 = 1e-6
        
        self.Gsvc1 = 15.74 / (5.2 - 5.0)
        self.Pcv = 5.0
        
        # External Pressures
        if posture == 'supine':
            self.Pj3ext = 0.0
            self.Pj2ext = 0.0
            self.Pj1ext = -6.5
            self.h_j3 = 0.0
            self.h_j2 = 0.0
            self.h_j1 = 0.0
            self.h_vs = 0.0
        else:
            self.Pj3ext = 0.0
            self.Pj2ext = 0.0
            self.Pj1ext = -6.5
            self.h_j3 = 35.0
            self.h_j2 = 25.0
            self.h_j1 = 10.0
            self.h_vs = 42.0
        
        self.blood_density = 1.05
        self.g = 980.0
        
        # Catheter (Reinfusion is direct, no separate state)
        self.blood_viscosity = 3.5e-3
        self.catheter_length = 20.0
        self.catheter_radius = 0.10
        self.R_catheter = self.calculate_catheter_resistance(self.catheter_radius, self.catheter_length)
        
        self.aspiration_protocol = None
    
    def calculate_catheter_resistance(self, radius_cm, length_cm=None):
        L_cm = length_cm if length_cm is not None else self.catheter_length
        L_m = L_cm / 100.0
        r_m = radius_cm / 100.0
        R_si = (8.0 * self.blood_viscosity * L_m) / (np.pi * r_m**4)
        return R_si * 0.00750062 / 1e6
    
    def get_initial_conditions(self):
        pic_val = getattr(self, 'Pic0', None)
        if pic_val is None: pic_val = 9.5
        
        # 15 State Variables (No P_reinfuse)
        return np.array([
            pic_val,          # 0: Pic
            58.9,             # 1: Ppa
            14.1,             # 2: Pv
            6.0,              # 3: Pvs
            5.85,             # 4: Pjr3
            5.85,             # 5: Pjl3
            5.70,             # 6: Pjr2
            5.70,             # 7: Pjl2
            6.00,             # 8: Pc3
            5.85,             # 9: Pc2
            5.80,             # 10: Pvv
            5.50,             # 11: Pazy
            5.20,             # 12: Psvc
            self.xaut,        # 13: xaut
            self.Cpan         # 14: Cpa
        ])

__all__ = ['ModelParameters']