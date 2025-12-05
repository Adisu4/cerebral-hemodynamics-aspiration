"""
Cerebral hemodynamics ODE system

Implements cerebral autoregulation, venous collapse, jugular collapsibility,
and intracranial pressure-volume dynamics.
"""
import numpy as np

def jugular_conductance(P_internal, P_external, k, A):
    P_transmural = P_internal - P_external
    G = k * (1 + (2/np.pi) * np.arctan(P_transmural / A))**2
    return max(G, 1e-9)

def cerebral_arterioles_resistance(Cpa, Ppa, Pic, params):
    driving_pressure = max(Ppa - Pic, 0.1)
    Cpa_safe = max(Cpa, 1e-6)
    Rpa = params.kR * (params.Cpan**2) / (driving_pressure**2 * Cpa_safe**2)
    return Rpa

def cerebral_veins_resistance(Pv, Pic, Pvs, params):
    """
    Gadda Eq 14 - Starling resistor (ROBUST IMPLEMENTATION).
    """
    # Calculate transmural pressure (distending pressure)
    # If Pic > Pv, the vein collapses.
    transmural = Pv - Pic
    
    # 1. Normal state: Pv > Pvs (Flow downstream)
    # 2. Backpressure state: Pv <= Pvs (Flow stops or reverses - handled by driving pressure in flow eq)

    
    if Pv > Pvs:
        if transmural > 0.1:
            # Open or partially open
            Rvs = params.Rvs1 * (Pv - Pvs) / transmural
        else:
            # Collapsed / Collapsing state
            # As transmural pressure approaches zero, resistance approaches infinity,
            # causing numerical instability. Physiologically, veins maintain residual
            # channels preventing complete occlusion. Limit Rvs to 100x baseline to
            # represent this physiological constraint and maintain numerical stability.
            R_max = params.Rvs1 * 100.0
            Rvs = R_max
    else:
        # If Pvs > Pv, the "Starling" effect acts differently or is moot because flow reverses.
        # Standard model assumes baseline resistance for backflow or low pressure gradient not driven by collapse.
        Rvs = params.Rvs1
        
    return Rvs

def cerebral_veins_capacitance(Pv, Pic, Pv1):
    kven = 0.155
    denominator = kven * (Pv - Pic - Pv1)
    # Clamp capacitance to avoid singularity or negative capacitance
    return 1.0 / denominator if denominator > 1e-4 else 1e9

def gadda_ode_system_FIXED(t, y, params, flow_rate_callback=None):
    # Unpack state
    Pic, Ppa, Pv, Pvs, Pjr3, Pjl3, Pjr2, Pjl2, Pc3, Pc2, Pvv, Pazy, Psvc, xaut, Cpa = y
    
    # Hydrostatic pressure corrections (supine configuration)
    h_hydro_vs = h_hydro_j3 = h_hydro_j2 = h_hydro_j1 = 0.0
    if params.posture == 'upright':
        pass  # Upright posture hydrostatic corrections applied in parameters module
    
    Pvs_int = Pvs - h_hydro_vs
    Pjr3_int = Pjr3 - h_hydro_j3
    Pjl3_int = Pjl3 - h_hydro_j3
    Pjr2_int = Pjr2 - h_hydro_j2
    Pjl2_int = Pjl2 - h_hydro_j2
    
    # Conductances
    Gjr3 = jugular_conductance(Pjr3_int, params.Pj3ext, params.kjr3, params.A)
    Gjl3 = jugular_conductance(Pjl3_int, params.Pj3ext, params.kjl3, params.A)
    Gjr2 = jugular_conductance(Pjr2_int, params.Pj2ext, params.kjr2, params.A)
    Gjl2 = jugular_conductance(Pjl2_int, params.Pj2ext, params.kjl2, params.A)
    Gjr1 = jugular_conductance(Pjr2, params.Pj1ext, params.kjr1, params.A)
    Gjl1 = jugular_conductance(Pjl2, params.Pj1ext, params.kjl1, params.A)
    
    Rpa = cerebral_arterioles_resistance(Cpa, Ppa, Pic, params)
    Rvs = cerebral_veins_resistance(Pv, Pic, Pvs, params)
    Cvi = cerebral_veins_capacitance(Pv, Pic, params.Pv1)
    Cic = 1.0 / (params.kE * max(Pic, 0.1))
    
    # Flows
    Pc_num = (Pv / params.Rpv) + (Ppa / (Rpa / 2.0)) + (Pic / params.Rf)
    Pc_den = (1.0 / params.Rpv) + (1.0 / (Rpa / 2.0)) + (1.0 / params.Rf)
    Pc = Pc_num / Pc_den
    
    Qf = max(0.0, (Pc - Pic) / params.Rf)
    Q0 = (Pic - Pvs) * params.G0 if Pic > Pvs else 0.0
    Q_in_pa = (params.Pa - Ppa) / (params.Rla + Rpa/2.0)
    Q_out_pa = (Ppa - Pc) / (Rpa / 2.0)
    Q_in_v = (Pc - Pv) / params.Rpv
    Q_out_v = (Pv - Pvs) / Rvs
    
    # Autoregulation
    Q = Q_in_pa
    dxaut_dt = (1/params.tau_aut) * (-xaut + params.Gaut * (Q - params.Qn) / params.Qn)
    if xaut < 0:
        kCpa = params.Cpa1 / 4.0
        Cpa_target = params.Cpan + (params.Cpa1 - params.Cpan) / (1.0 + np.exp(-xaut / kCpa))
    else:
        kCpa = params.Cpa2 / 4.0
        Cpa_target = params.Cpa2 + (params.Cpan - params.Cpa2) / (1.0 + np.exp(xaut / kCpa))
    dCpa_dt = (Cpa_target - Cpa) / params.tau_aut
    
    # ASPIRATION
    Q_asp_Pv = Q_asp_Pvs = Q_asp_Jr3 = Q_asp_Jl3 = Q_asp_Jr2 = Q_asp_Jl2 = 0.0
    
    if flow_rate_callback:
        fr = flow_rate_callback(t, Pic)
        if fr > 0:
            prot = getattr(params, "aspiration_protocol", {})
            site = prot.get('site', 'Pvs')
            
            if site == "Pvs": Q_asp_Pvs = fr
            elif site == "Pv": Q_asp_Pv = fr
            elif site == "J3": 
                Q_asp_Jr3 = fr / 2.0
                Q_asp_Jl3 = fr / 2.0
            elif site == "J2":
                Q_asp_Jr2 = fr / 2.0
                Q_asp_Jl2 = fr / 2.0
            elif site == "J1":
                Q_asp_Jr2 = fr / 2.0
                Q_asp_Jl2 = fr / 2.0

    # Derivatives
    d_trans_pa = (1/Cpa) * (Q_in_pa - Q_out_pa)
    d_trans_v = (1/Cvi) * (Q_in_v - Q_out_v - Q_asp_Pv)
    dPic = (1/Cic) * (Cpa*d_trans_pa + Cvi*d_trans_v + dCpa_dt*(Ppa-Pic) + Qf - Q0)
    dPpa = d_trans_pa + dPic
    dPv = d_trans_v + dPic
    
    dPvs = (1/params.Cvs) * (
        Q_out_v + Q0 - Q_asp_Pvs -
        (Pvs-Pjr3)*Gjr3 - (Pvs-Pjl3)*Gjl3 - (Pvs-Pc3)*params.Gc3 - 
        (Pvs-Pvv)*params.Gvvl - (Pvs-Pvv)*params.Gvvr
    )
    dPjr3 = (1/params.Cjr3)*((Pvs-Pjr3)*Gjr3 + (Pc3-Pjr3)*params.Gcjr3 - (Pjr3-Pjr2)*Gjr2 - Q_asp_Jr3)
    dPjl3 = (1/params.Cjl3)*((Pvs-Pjl3)*Gjl3 + (Pc3-Pjl3)*params.Gcjl3 - (Pjl3-Pjl2)*Gjl2 - Q_asp_Jl3)
    dPjr2 = (1/params.Cjr2)*((Pjr3-Pjr2)*Gjr2 + (Pc2-Pjr2)*params.Gcjr2 - (Pjr2-Psvc)*Gjr1 - Q_asp_Jr2)
    dPjl2 = (1/params.Cjl2)*((Pjl3-Pjl2)*Gjl2 + (Pc2-Pjl2)*params.Gcjl2 - (Pjl2-Psvc)*Gjl1 - Q_asp_Jl2)
    dPc3 = (1/params.Cc3)*((params.Pa-Pc3)*params.Gex - (Pc3-Pvs)*params.Gc3 - (Pc3-Pjr3)*params.Gcjr3 - (Pc3-Pjl3)*params.Gcjl3 - (Pc3-Pc2)*params.Gc2)
    dPc2 = (1/params.Cc2)*((Pc3-Pc2)*params.Gc2 - (Pc2-Pjr2)*params.Gcjr2 - (Pc2-Pjl2)*params.Gcjl2 - (Pc2-params.Pcv)*params.Gc1)
    dPvv = (1/params.Cvv)*((Pvs-Pvv)*params.Gvvl + (Pvs-Pvv)*params.Gvvr - (Pvv-Pazy)*params.Gazy1 - (Pvv-params.Pcv)*params.Gvv2)
    dPazy = (1/params.Cazy)*((Pvv-Pazy)*params.Gazy1 - (Pazy-Psvc)*params.Gazy2 - (Pazy-params.Pcv)*params.Glv - (Pazy-params.Pcv)*params.Grv)
    dPsvc = (1/params.Csvc)*((Pjr2-Psvc)*Gjr1 + (Pjl2-Psvc)*Gjl1 + (Pazy-Psvc)*params.Gazy2 - (Psvc-params.Pcv)*params.Gsvc1)
    
    return [dPic, dPpa, dPv, dPvs, dPjr3, dPjl3, dPjr2, dPjl2, dPc3, dPc2, dPvv, dPazy, dPsvc, dxaut_dt, dCpa_dt]