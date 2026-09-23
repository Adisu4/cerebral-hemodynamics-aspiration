"""Physiological and numerical parameters for the hemodynamic model.

The intracranial and basal extracranial values follow Gadda et al. (2015),
Tables 1-4 and Appendix equations 17-32.  Three conductances omitted from
the 2015 tables (Gc3, Gvv2 and Glv) are taken from the authors' subsequent
validation paper (Gadda et al., AJNR 2016, Table 1).  The fixed boundary
pressures Piv and Plv are then calculated from the 2015 basal flows so the
printed equations close at the stated 0.40/0.13/0.27 mL/s split.

Only the supine configuration is supported.  TBI-like changes are isolated
in :func:`make_tbi_parameters` and are never mixed into the reference model.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import numpy as np


@dataclass
class ModelParameters:
    posture: str = "supine"
    terminal_resistance_mode: str = "source_eq14"

    # Intracranial model (Gadda Table 1 / cited Ursino model).
    Cpan: float = 0.205
    Cpa1: float = 2.87
    Cpa2: float = 0.164
    Pa: float = 100.0
    Picn: float = 9.5
    Ppa0: float = 58.9
    Pv0: float = 14.1
    Pv1: float = -2.5
    Pvs0: float = 6.0
    Qn: float = 12.5
    R0: float = 526.3
    Rf: float = 2380.0
    Rla: float = 0.6
    Rpv: float = 0.880
    Rvs1: float = 0.366
    Gaut: float = 3.0
    tau_aut: float = 20.0
    xaut0: float = 2.16e-4
    kE: float = 0.077
    kR: float = 13.1e3
    kven: float = 0.155

    # Capacities (Gadda Table 4), mL/mmHg.
    Cvs: float = 0.5
    Cjr3: float = 1.0
    Cjl3: float = 1.0
    Cjr2: float = 2.5
    Cjl2: float = 2.5
    Cc3: float = 0.7
    Cc2: float = 1.4
    Csvc: float = 20.0
    Cazy: float = 0.5
    Cvv: float = 0.5

    # Jugular sigmoid parameters (Gadda H221 and Eqs. 27-32).
    kjr3: float = 11.0
    kjl3: float = 11.0
    kjr2: float = 13.0
    kjl2: float = 13.0
    kjr1: float = 6.9
    kjl1: float = 6.9
    A: float = 0.8

    # Constant extracranial conductances, mL/(mmHg s).
    Gex: float = 5.0 / (100.0 - 6.0)
    # c3 is latent at the basal state because Pvs == Pc3, not deleted.
    # Gadda et al. 2016, Table 1B (the 2015 paper did not tabulate it).
    Gc3: float = 21.43
    Gc2: float = 3.0 / (6.0 - 5.85)
    # Table 6 gives Qc1=1.00 mL/s.
    Gc1: float = 1.0 / (5.85 - 5.0)
    Gcjr3: float = 1.0 / (6.0 - 5.85)
    Gcjl3: float = 1.0 / (6.0 - 5.85)
    Gcjr2: float = 1.0 / (5.85 - 5.70)
    Gcjl2: float = 1.0 / (5.85 - 5.70)
    Gvvr: float = 0.4 / (6.0 - 5.8)
    Gvvl: float = 0.4 / (6.0 - 5.8)
    # Gadda et al. 2016, Table 1B.  These two conductances do not appear in
    # the numerical tables of the 2015 paper.
    Gvv2: float = 0.83
    Gazy1: float = 0.4 / (5.8 - 5.5)
    Glv: float = 0.89
    # At equilibrium Qazy2 = Qazy1 + Qlv = 0.53 mL/s.
    Gazy2: float = (0.4 + 0.13) / (5.5 - 5.2)
    # Source basal jugular-confluence flow is approximately 15.7 mL/s.
    Gsvc1: float = 15.7 / (5.4 - 5.2)
    # Lower SVC receives Qsvc1 plus Qazy2.
    Gsvc2: float = (15.7 + 0.4 + 0.13) / (5.2 - 5.0)

    Pcv: float = 5.0
    # Fixed boundary pressures implied by the reported basal flows and the
    # later published conductances: Qvv2=0.40 and Qlv=0.13 mL/s.
    Piv: float = 5.8 - 0.4 / 0.83
    Plv: float = 5.5 + 0.13 / 0.89
    Pj3ext: float = 0.0
    Pj2ext: float = 0.0
    Pj1ext: float = -6.5

    aspiration_site: str = "Pv"

    @property
    def G0(self) -> float:
        return 1.0 / self.R0

    def initial_state(self) -> np.ndarray:
        if self.posture != "supine":
            raise ValueError("Only the supine configuration is implemented.")
        return np.array(
            [
                self.Picn,
                self.Ppa0,
                self.Pv0,
                self.Pvs0,
                5.85,
                5.85,
                5.70,
                5.70,
                6.00,
                5.85,
                5.80,
                5.50,
                5.20,
                self.xaut0,
            ],
            dtype=float,
        )

    def serializable(self) -> dict:
        result = asdict(self)
        result["G0"] = self.G0
        result["model_scope"] = "supine_only"
        result["lower_body_boundaries"] = (
            "Piv and Plv fixed from 2015 basal flows using conductances in "
            "Gadda et al. 2016 Table 1B"
        )
        result["Gc3_provenance"] = "Gadda et al. 2016 Table 1B; sensitivity tested"
        return result


def make_reference_parameters() -> ModelParameters:
    return ModelParameters()


def make_tbi_parameters() -> ModelParameters:
    """Construct the prespecified TBI-like state without altering base topology."""
    p = ModelParameters()
    p.R0 = 1800.0
    p.Gaut = 0.3
    p.Rvs1 *= 3.0
    p.kjr3 *= 0.15
    p.kjl3 *= 0.15
    return p


__all__ = ["ModelParameters", "make_reference_parameters", "make_tbi_parameters"]
