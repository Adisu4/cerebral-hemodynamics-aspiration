"""Equation-level verification tests for the hemodynamic model."""

from __future__ import annotations

import copy
import unittest

import numpy as np
from scipy.integrate import solve_ivp

from cerebral_hemodynamics_aspiration import model
from cerebral_hemodynamics_aspiration.parameters import (
    make_reference_parameters,
    make_tbi_parameters,
)


class PublicationModelChecks(unittest.TestCase):
    def test_capacity_matches_published_logistic_and_derivative(self) -> None:
        p = make_reference_parameters()
        for x in (-0.5, -0.01, 0.0, 0.01, 0.1):
            amplitude = p.Cpa1 if x < 0.0 else p.Cpa2
            expected = p.Cpan - 0.5 * amplitude * np.tanh(2.0 * x / amplitude)
            capacity, derivative = model.autoregulation_capacity(x, p)
            self.assertAlmostEqual(capacity, expected, places=13)
            h = 1e-7
            numeric = (
                model.autoregulation_capacity(x + h, p)[0]
                - model.autoregulation_capacity(x - h, p)[0]
            ) / (2.0 * h)
            self.assertAlmostEqual(derivative, numeric, places=7)

    def test_jugular_conductances_use_source_upstream_nodes(self) -> None:
        p = make_reference_parameters()
        y = p.initial_state()
        original = model.jugular_conductances(y, p)

        downstream_changed = y.copy()
        downstream_changed[4] += 1.0  # Pjr3 is downstream of the J3 conductance.
        changed = model.jugular_conductances(downstream_changed, p)
        self.assertEqual(original["Gjr3"], changed["Gjr3"])
        self.assertNotEqual(original["Gjr2"], changed["Gjr2"])

        upstream_changed = y.copy()
        upstream_changed[3] += 1.0  # Pvs is upstream of both J3 segments.
        changed = model.jugular_conductances(upstream_changed, p)
        self.assertNotEqual(original["Gjr3"], changed["Gjr3"])
        self.assertNotEqual(original["Gjl3"], changed["Gjl3"])

    def test_svc_junction_has_zero_mass_residual(self) -> None:
        p = make_reference_parameters()
        y = p.initial_state()
        G = model.jugular_conductances(y, p)
        junction = model.superior_vena_cava_junction(y, G, p)
        residual = (
            (y[6] - junction) * G["Gjr1"]
            + (y[7] - junction) * G["Gjl1"]
            - (junction - y[12]) * p.Gsvc1
        )
        self.assertAlmostEqual(residual, 0.0, places=12)

    def test_reported_lower_body_flow_split_closes(self) -> None:
        p = make_reference_parameters()
        flows = model.flow_snapshot(p.initial_state(), p)
        self.assertAlmostEqual(float(flows["Qazy1"]), 0.40, places=12)
        self.assertAlmostEqual(float(flows["Qvv2"]), 0.40, places=12)
        self.assertAlmostEqual(float(flows["Qlv"]), 0.13, places=12)
        self.assertAlmostEqual(float(flows["Qazy2"]), 0.53, places=12)
        self.assertAlmostEqual(
            float(flows["Qazy1"]) + float(flows["Qlv"]),
            float(flows["Qazy2"]),
            places=12,
        )

    def test_whole_network_mass_balance_for_every_site(self) -> None:
        p = make_reference_parameters()
        y = p.initial_state()
        capacities = np.array(
            [p.Cvs, p.Cjr3, p.Cjl3, p.Cjr2, p.Cjl2,
             p.Cc3, p.Cc2, p.Cvv, p.Cazy, p.Csvc]
        )
        for site in model.VALID_ASPIRATION_SITES:
            aspiration = model.Aspiration(site, 4.0)
            dy, aux = model.evaluate(y, p, aspiration, return_to_svc_ml_s=1.25)
            storage = dy[0] / (p.kE * y[0]) + capacities @ dy[3:13]
            external = (
                float(aux["Q"])
                + (p.Pa - y[8]) * p.Gex
                + (p.Plv - y[11]) * p.Glv
                + 1.25
                - (y[9] - p.Pcv) * p.Gc1
                - (y[10] - p.Piv) * p.Gvv2
                - (y[12] - p.Pcv) * p.Gsvc2
                - 4.0
            )
            self.assertAlmostEqual(storage, external, places=10)

    def test_starling_law_is_exact_and_domain_is_explicit(self) -> None:
        p = make_reference_parameters()
        value, branch = model.cerebral_veins_resistance(14.0, 9.0, 6.0, p)
        self.assertEqual(branch, "open_starling")
        self.assertAlmostEqual(value, p.Rvs1 * (14.0 - 6.0) / (14.0 - 9.0))
        value, branch = model.cerebral_veins_resistance(5.0, 4.0, 6.0, p)
        self.assertEqual(value, p.Rvs1)
        self.assertEqual(branch, "reverse_or_equal")
        with self.assertRaises(model.ModelDomainError):
            model.cerebral_veins_resistance(10.0, 10.0, 6.0, p)
        p.terminal_resistance_mode = "fixed_Rvs1"
        value, branch = model.cerebral_veins_resistance(14.0, 9.0, 6.0, p)
        self.assertEqual(value, p.Rvs1)
        self.assertEqual(branch, "fixed_comparator")

    def test_reference_equilibrium_matches_gadda_targets(self) -> None:
        p = make_reference_parameters()
        sol = solve_ivp(
            model.rhs(p),
            (0.0, 8000.0),
            p.initial_state(),
            method="BDF",
            rtol=1e-8,
            atol=1e-10,
            max_step=4.0,
        )
        self.assertTrue(sol.success)
        y = sol.y[:, -1]
        flows = model.flow_snapshot(y, p)
        self.assertAlmostEqual(y[0], 9.44, delta=0.03)
        self.assertAlmostEqual(y[3], 6.00, delta=0.03)
        self.assertAlmostEqual(
            float(flows["Qvvr"]) + float(flows["Qvvl"]), 0.79, delta=0.03
        )
        self.assertLess(np.max(np.abs(model.rhs(p)(8000.0, y)[:13])) * 60.0, 1e-4)

    def test_tbi_changes_are_isolated(self) -> None:
        reference = make_reference_parameters()
        tbi = make_tbi_parameters()
        changed = {
            name for name in reference.__dataclass_fields__
            if getattr(reference, name) != getattr(tbi, name)
        }
        self.assertEqual(changed, {"R0", "Gaut", "Rvs1", "kjr3", "kjl3"})

    def test_scope_and_site_guards(self) -> None:
        p = make_reference_parameters()
        upright = copy.deepcopy(p)
        upright.posture = "upright"
        with self.assertRaises(ValueError):
            model.evaluate(upright.initial_state(), upright)
        with self.assertRaises(ValueError):
            model.Aspiration("J1", 1.0)


if __name__ == "__main__":
    unittest.main()
