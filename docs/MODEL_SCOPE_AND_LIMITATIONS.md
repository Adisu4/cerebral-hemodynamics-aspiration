# Model scope and limitations

## What the implementation supports

The implementation is a deterministic, supine-only lumped-parameter model of
intracranial hemodynamics, cerebrospinal-fluid dynamics, cerebral venous
outflow, and an extracranial venous network. Withdrawal can be placed at the
aggregate cerebral-vein node, venous-sinus node, or bilateral J3 and J2 nodes.

The TBI-like state changes five model parameters simultaneously: CSF outflow
resistance, autoregulatory gain, terminal venous resistance, and the right and
left J3 conductance coefficients. It is a constructed scenario rather than a
validated representation of a TBI population.

## Interpretation limits

The primary intervention removes flow from a model with fixed arterial and
central venous pressure boundaries. The system has no finite whole-body blood
volume and no systemic pressure response. At the accepted endpoints, the
240 mL/min and 480 mL/min cerebral-vein simulations integrate to approximately
31.08 L and 86.16 L of withdrawal, respectively. Those values expose the
open-reservoir assumption; they are not treatment volumes.

The large upstream/downstream difference is partly structural. In the active
terminal-vein branch,

```text
Rvs = Rvs1 * (Pv - Pvs) / (Pv - Pic)
```

so terminal outflow simplifies to `(Pv - Pic) / Rvs1` and no longer depends
directly on sinus pressure. Withdrawal at the cerebral-vein node also enters
the intracranial volume balance directly, whereas downstream sinks act through
other pressure states. The resulting site ordering is a property of the
implemented topology and constitutive law and needs physiological testing.

The model does not represent catheter geometry, local vessel collapse or wall
contact, regional cerebral anatomy, oxygen delivery, hemolysis, thrombosis,
coagulation, embolism, blood-volume depletion, baroreflexes, or patient
heterogeneity. It must not be used to select an aspiration rate, duration, or
catheter position for a patient.

## Evidence level

The reference-state comparison reproduces selected resting values from the
source model. This is an implementation benchmark, not independent validation
of the aspiration extension or the TBI-like parameter state. The saved output
and solver comparisons establish numerical reproducibility within the same
equations. They do not establish biological accuracy, efficacy, or safety.
