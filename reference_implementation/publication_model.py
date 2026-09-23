"""Source-traceable supine cerebral venous aspiration model.

The state equations follow Gadda et al. (2015), Eqs. 3-32.  The model is
deliberately limited to supine simulations.  Aspiration adds an external
withdrawal term at one of four explicitly represented sites.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Mapping

import numpy as np

from publication_parameters import PublicationParameters


STATE_NAMES = (
    "Pic", "Ppa", "Pv", "Pvs", "Pjr3", "Pjl3", "Pjr2", "Pjl2",
    "Pc3", "Pc2", "Pvv", "Pazy", "Psvc", "xaut",
)
PRESSURE_STATE_COUNT = 13
VALID_ASPIRATION_SITES = ("Pv", "Pvs", "J3", "J2")


class ModelDomainError(RuntimeError):
    """Raised when the published constitutive law leaves its valid domain."""


@dataclass(frozen=True)
class Aspiration:
    site: str = "Pv"
    flow_ml_s: float = 0.0

    def __post_init__(self) -> None:
        if self.site not in VALID_ASPIRATION_SITES:
            raise ValueError(f"site must be one of {VALID_ASPIRATION_SITES}")
        if self.flow_ml_s < 0:
            raise ValueError("flow_ml_s must be nonnegative")


def autoregulation_capacity(
    xaut: float, p: PublicationParameters
) -> tuple[float, float]:
    """Return Cpa and dCpa/dxaut from Gadda Eqs. 8-10."""
    amplitude = p.Cpa1 if xaut < 0.0 else p.Cpa2
    h = np.tanh(2.0 * xaut / amplitude)
    return p.Cpan - 0.5 * amplitude * h, -(1.0 - h * h)


def jugular_conductance(
    upstream_pressure: float, external_pressure: float, k: float, A: float
) -> float:
    """Jugular tube law (Gadda Eq. 2 and Eqs. 27-32)."""
    if A <= 0.0:
        raise ValueError("A must be positive")
    value = k * (1.0 + (2.0 / np.pi) * np.arctan(
        (upstream_pressure - external_pressure) / A
    )) ** 2
    return float(value)


def jugular_conductances(y: np.ndarray, p: PublicationParameters) -> dict[str, float]:
    """Return all six conductances using the source-defined upstream node."""
    Pvs, Pjr3, Pjl3, Pjr2, Pjl2 = y[3], y[4], y[5], y[6], y[7]
    return {
        "Gjr3": jugular_conductance(Pvs, p.Pj3ext, p.kjr3, p.A),
        "Gjl3": jugular_conductance(Pvs, p.Pj3ext, p.kjl3, p.A),
        "Gjr2": jugular_conductance(Pjr3, p.Pj2ext, p.kjr2, p.A),
        "Gjl2": jugular_conductance(Pjl3, p.Pj2ext, p.kjl2, p.A),
        "Gjr1": jugular_conductance(Pjr2, p.Pj1ext, p.kjr1, p.A),
        "Gjl1": jugular_conductance(Pjl2, p.Pj1ext, p.kjl1, p.A),
    }


def cerebral_arterioles_resistance(
    Cpa: float, Ppa: float, Pic: float, p: PublicationParameters
) -> float:
    """Pial arteriolar resistance from Gadda Eq. 11."""
    if Cpa <= 0.0 or Ppa <= Pic:
        raise ModelDomainError("Eq. 11 requires Cpa > 0 and Ppa > Pic")
    return p.kR * p.Cpan**2 / ((Ppa - Pic) * Cpa) ** 2


def cerebral_veins_resistance(
    Pv: float, Pic: float, Pvs: float, p: PublicationParameters
) -> tuple[float, str]:
    """Terminal cerebral-vein resistance from Gadda Eq. 14.

    No unpublished cap is applied.  The open branch requires positive
    transmural pressure; a violation is reported instead of silently changing
    the constitutive law.
    """
    if p.terminal_resistance_mode == "fixed_Rvs1":
        return p.Rvs1, "fixed_comparator"
    if p.terminal_resistance_mode != "source_eq14":
        raise ValueError("terminal_resistance_mode must be 'source_eq14' or 'fixed_Rvs1'")
    if Pv > Pvs:
        transmural = Pv - Pic
        if transmural <= 0.0:
            raise ModelDomainError(
                "Eq. 14 open branch requires Pv-Pic > 0; "
                f"received {transmural:.6g} mmHg"
            )
        return p.Rvs1 * (Pv - Pvs) / transmural, "open_starling"
    return p.Rvs1, "reverse_or_equal"


def cerebral_veins_capacitance(
    Pv: float, Pic: float, p: PublicationParameters
) -> float:
    """Intracranial venous capacitance from Gadda Eq. 6."""
    denominator = p.kven * (Pv - Pic - p.Pv1)
    if denominator <= 0.0:
        raise ModelDomainError(
            "Eq. 6 requires kven*(Pv-Pic-Pv1) > 0; "
            f"received {denominator:.6g}"
        )
    return 1.0 / denominator


def superior_vena_cava_junction(
    y: np.ndarray, conductances: Mapping[str, float], p: PublicationParameters
) -> float:
    """Algebraic Psvc1 satisfying mass balance at the jugular confluence."""
    Pjr2, Pjl2, Psvc = y[6], y[7], y[12]
    Gjr1 = conductances["Gjr1"]
    Gjl1 = conductances["Gjl1"]
    return (
        Gjr1 * Pjr2 + Gjl1 * Pjl2 + p.Gsvc1 * Psvc
    ) / (Gjr1 + Gjl1 + p.Gsvc1)


def _aspiration_sinks(site: str, flow_ml_s: float) -> dict[str, float]:
    sinks = {name: 0.0 for name in ("Pv", "Pvs", "Jr3", "Jl3", "Jr2", "Jl2")}
    if flow_ml_s <= 0.0:
        return sinks
    if site == "Pv":
        sinks["Pv"] = flow_ml_s
    elif site == "Pvs":
        sinks["Pvs"] = flow_ml_s
    elif site == "J3":
        sinks["Jr3"] = sinks["Jl3"] = flow_ml_s / 2.0
    elif site == "J2":
        sinks["Jr2"] = sinks["Jl2"] = flow_ml_s / 2.0
    else:
        raise ValueError(f"site must be one of {VALID_ASPIRATION_SITES}")
    return sinks


def evaluate(
    y: np.ndarray,
    p: PublicationParameters,
    aspiration: Aspiration | None = None,
    return_to_svc_ml_s: float = 0.0,
) -> tuple[np.ndarray, dict[str, float | str]]:
    """Evaluate the ODE and return derivatives plus auditable auxiliaries."""
    if p.posture != "supine":
        raise ValueError("The publication model is supine-only")
    if len(y) != len(STATE_NAMES):
        raise ValueError(f"Expected {len(STATE_NAMES)} states, received {len(y)}")
    if return_to_svc_ml_s < 0.0:
        raise ValueError("return_to_svc_ml_s must be nonnegative")

    (
        Pic, Ppa, Pv, Pvs, Pjr3, Pjl3, Pjr2, Pjl2,
        Pc3, Pc2, Pvv, Pazy, Psvc, xaut,
    ) = np.asarray(y, dtype=float)

    Cpa, dCpa_dx = autoregulation_capacity(xaut, p)
    Rpa = cerebral_arterioles_resistance(Cpa, Ppa, Pic, p)
    Rvs, rvs_branch = cerebral_veins_resistance(Pv, Pic, Pvs, p)
    Cvi = cerebral_veins_capacitance(Pv, Pic, p)
    if Pic <= 0.0:
        raise ModelDomainError("Eq. 16 requires Pic > 0")
    Cic = 1.0 / (p.kE * Pic)

    Pc = (
        Pv / p.Rpv + Ppa / (Rpa / 2.0) + Pic / p.Rf
    ) / (1.0 / p.Rpv + 1.0 / (Rpa / 2.0) + 1.0 / p.Rf)
    Qf = max((Pc - Pic) / p.Rf, 0.0)
    Q0 = max((Pic - Pvs) / p.R0, 0.0)
    Qin_pa = (p.Pa - Ppa) / (p.Rla + Rpa / 2.0)
    Qout_pa = (Ppa - Pc) / (Rpa / 2.0)
    Qin_v = (Pc - Pv) / p.Rpv
    Qout_v = (Pv - Pvs) / Rvs

    dxaut = (-xaut + p.Gaut * (Qin_pa - p.Qn) / p.Qn) / p.tau_aut
    dCpa = dCpa_dx * dxaut

    aspiration = aspiration or Aspiration()
    sinks = _aspiration_sinks(aspiration.site, aspiration.flow_ml_s)

    dtrans_pa = (
        Qin_pa - Qout_pa - dCpa * (Ppa - Pic)
    ) / Cpa
    dtrans_v = (Qin_v - Qout_v - sinks["Pv"]) / Cvi
    dPic = (
        Cpa * dtrans_pa
        + Cvi * dtrans_v
        + dCpa * (Ppa - Pic)
        + Qf
        - Q0
    ) / Cic
    dPpa = dtrans_pa + dPic
    dPv = dtrans_v + dPic

    G = jugular_conductances(y, p)
    Psvc1 = superior_vena_cava_junction(y, G, p)

    dPvs = (
        Qout_v + Q0 - sinks["Pvs"]
        - (Pvs - Pjr3) * G["Gjr3"]
        - (Pvs - Pjl3) * G["Gjl3"]
        - (Pvs - Pc3) * p.Gc3
        - (Pvs - Pvv) * p.Gvvl
        - (Pvs - Pvv) * p.Gvvr
    ) / p.Cvs
    dPjr3 = (
        (Pvs - Pjr3) * G["Gjr3"]
        - (Pjr3 - Pc3) * p.Gcjr3
        - (Pjr3 - Pjr2) * G["Gjr2"]
        - sinks["Jr3"]
    ) / p.Cjr3
    dPjl3 = (
        (Pvs - Pjl3) * G["Gjl3"]
        - (Pjl3 - Pc3) * p.Gcjl3
        - (Pjl3 - Pjl2) * G["Gjl2"]
        - sinks["Jl3"]
    ) / p.Cjl3
    dPjr2 = (
        (Pjr3 - Pjr2) * G["Gjr2"]
        - (Pjr2 - Pc2) * p.Gcjr2
        - (Pjr2 - Psvc1) * G["Gjr1"]
        - sinks["Jr2"]
    ) / p.Cjr2
    dPjl2 = (
        (Pjl3 - Pjl2) * G["Gjl2"]
        - (Pjl2 - Pc2) * p.Gcjl2
        - (Pjl2 - Psvc1) * G["Gjl1"]
        - sinks["Jl2"]
    ) / p.Cjl2
    dPc3 = (
        (Pvs - Pc3) * p.Gc3
        + (Pjr3 - Pc3) * p.Gcjr3
        + (Pjl3 - Pc3) * p.Gcjl3
        + (p.Pa - Pc3) * p.Gex
        - (Pc3 - Pc2) * p.Gc2
    ) / p.Cc3
    dPc2 = (
        (Pc3 - Pc2) * p.Gc2
        + (Pjr2 - Pc2) * p.Gcjr2
        + (Pjl2 - Pc2) * p.Gcjl2
        - (Pc2 - p.Pcv) * p.Gc1
    ) / p.Cc2
    dPsvc = (
        (Psvc1 - Psvc) * p.Gsvc1
        + (Pazy - Psvc) * p.Gazy2
        - (Psvc - p.Pcv) * p.Gsvc2
        + return_to_svc_ml_s
    ) / p.Csvc
    dPvv = (
        (Pvs - Pvv) * p.Gvvl
        + (Pvs - Pvv) * p.Gvvr
        - (Pvv - Pazy) * p.Gazy1
        - (Pvv - p.Piv) * p.Gvv2
    ) / p.Cvv
    dPazy = (
        (Pvv - Pazy) * p.Gazy1
        + (p.Plv - Pazy) * p.Glv
        - (Pazy - Psvc) * p.Gazy2
    ) / p.Cazy

    dy = np.array(
        [
            dPic, dPpa, dPv, dPvs, dPjr3, dPjl3, dPjr2, dPjl2,
            dPc3, dPc2, dPvv, dPazy, dPsvc, dxaut,
        ],
        dtype=float,
    )
    aux: dict[str, float | str] = {
        "Cpa": Cpa,
        "Cvi": Cvi,
        "Cic": Cic,
        "Rpa": Rpa,
        "Rvs": Rvs,
        "Rvs_branch": rvs_branch,
        "Pc": Pc,
        "Psvc1": Psvc1,
        "Q": Qin_pa,
        "Qout_pa": Qout_pa,
        "Qin_v": Qin_v,
        "Qout_v": Qout_v,
        "Qf": Qf,
        "Q0": Q0,
        "min_intracranial_transmural": min(Ppa - Pic, Pv - Pic - p.Pv1, Pv - Pic),
    }
    aux.update(G)
    return dy, aux


def rhs(
    p: PublicationParameters,
    aspiration_callback: Callable[[float, np.ndarray], Aspiration] | None = None,
    return_callback: Callable[[float, np.ndarray], float] | None = None,
) -> Callable[[float, np.ndarray], np.ndarray]:
    """Build a scipy-compatible right-hand-side function."""

    def function(t: float, y: np.ndarray) -> np.ndarray:
        aspiration = aspiration_callback(t, y) if aspiration_callback else Aspiration()
        returned = return_callback(t, y) if return_callback else 0.0
        return evaluate(y, p, aspiration, returned)[0]

    return function


def flow_snapshot(
    y: np.ndarray,
    p: PublicationParameters,
    aspiration: Aspiration | None = None,
    return_to_svc_ml_s: float = 0.0,
) -> dict[str, float | str]:
    """Return named flows and constitutive diagnostics at one state."""
    _, aux = evaluate(y, p, aspiration, return_to_svc_ml_s)
    P = dict(zip(STATE_NAMES, np.asarray(y, dtype=float)))
    G = {name: float(aux[name]) for name in ("Gjr3", "Gjl3", "Gjr2", "Gjl2", "Gjr1", "Gjl1")}
    Psvc1 = float(aux["Psvc1"])
    flows: dict[str, float | str] = {
        "Q_cerebral": float(aux["Q"]),
        "Q_terminal_cerebral": float(aux["Qout_v"]),
        "Qjr3": (P["Pvs"] - P["Pjr3"]) * G["Gjr3"],
        "Qjl3": (P["Pvs"] - P["Pjl3"]) * G["Gjl3"],
        "Qjr2": (P["Pjr3"] - P["Pjr2"]) * G["Gjr2"],
        "Qjl2": (P["Pjl3"] - P["Pjl2"]) * G["Gjl2"],
        "Qjr1": (P["Pjr2"] - Psvc1) * G["Gjr1"],
        "Qjl1": (P["Pjl2"] - Psvc1) * G["Gjl1"],
        "Qc3": (P["Pvs"] - P["Pc3"]) * p.Gc3,
        "Qc2": (P["Pc3"] - P["Pc2"]) * p.Gc2,
        "Qc1": (P["Pc2"] - p.Pcv) * p.Gc1,
        "Qvvr": (P["Pvs"] - P["Pvv"]) * p.Gvvr,
        "Qvvl": (P["Pvs"] - P["Pvv"]) * p.Gvvl,
        "Qazy1": (P["Pvv"] - P["Pazy"]) * p.Gazy1,
        "Qvv2": (P["Pvv"] - p.Piv) * p.Gvv2,
        "Qlv": (p.Plv - P["Pazy"]) * p.Glv,
        "Qazy2": (P["Pazy"] - P["Psvc"]) * p.Gazy2,
        "Qsvc1": (Psvc1 - P["Psvc"]) * p.Gsvc1,
        "Qsvc2": (P["Psvc"] - p.Pcv) * p.Gsvc2,
    }
    flows.update(aux)
    return flows


__all__ = [
    "Aspiration", "ModelDomainError", "PRESSURE_STATE_COUNT", "STATE_NAMES",
    "VALID_ASPIRATION_SITES", "autoregulation_capacity",
    "cerebral_arterioles_resistance", "cerebral_veins_capacitance",
    "cerebral_veins_resistance", "evaluate", "flow_snapshot",
    "jugular_conductance", "jugular_conductances", "rhs",
    "superior_vena_cava_junction",
]
