# Declared cofactor pools and net production

6 September 2026. `gembench/cofactor.py` packages two small checks used in the Putida pathway investigation. These are applications of established conservation and metabolite-production methods.

`pool_balance_certificate(model, weights)` computes a weighted sum of explicitly named metabolite rows. It records every incident reaction, its pool coefficients and bounds, and whether its allowed direction could supply the pool. Rational arithmetic on the stored floating-point coefficients determines whether the sum is exactly zero; small discrepancies are retained. The check also verifies that the operative solver rows match the stored reaction stoichiometry and impose zero steady-state balances. Missing metabolites, invalid weights or altered solver rows cause an error.

For example, equal weights on `q8_c` and `q8h2_c` test the represented quinone pool. A zero row sum certifies conservation within the specified model. It does not prove which chemical pool the organism uses, or that a growth requirement for it is appropriate. Atom-mapped conservation analysis can establish chemical moiety identity more directly; that is beyond this helper's scope. [Haraldsdóttir and Fleming, 2016](https://doi.org/10.1371/journal.pcbi.1004999). Later work further characterizes conserved and reacting molecular substructures. [Rahou et al., 2026](https://doi.org/10.1016/j.jtbi.2025.112348).

`probe_metabolite_production(model, identifiers)` separately maximizes an outward demand for each requested metabolite. It retains the current medium, maintenance, growth bounds and additional constraints, then restores the original model. The report distinguishes the production threshold from the artificial probe's capacity. A nonoptimal or nonfinite solution raises an error. An optimal zero means that net production is blocked under the supplied constraints. A direction-compatible reaction alone does not establish feasible production: its substrates may be unavailable.

```python
from gembench.cofactor import pool_balance_certificate, probe_metabolite_production

# Configure the intended medium and bounds before probing.
certificate = pool_balance_certificate(model, {"q8_c": 1.0, "q8h2_c": 1.0})
production = probe_metabolite_production(model, ["q8_c", "q8h2_c"])
```

MEMOTE already includes metabolite producibility and stoichiometric-consistency tests. This repository-specific helper adds an explicit small-pool certificate and retained-constraint reporting to the current experimental workflow; it makes no novelty claim for conservation analysis or demand maximization. [MEMOTE consistency API](https://memote.readthedocs.io/en/latest/autoapi/memote/support/consistency/index.html).

Twenty focused tests cover catalytic cycling without net supply, a real precursor route, blocked medium, reversed and disabled supply directions, tiny stoichiometric leaks, missing metabolites, custom solver supplies, failed or nonfinite solves, and restoration of existing constraints and objectives. The prepared Putida baseline reproduces the original ten-reaction conservation certificate.

The biological implication still requires care. A cofactor can be recycled in steady-state FBA while its de novo synthesis remains blocked. Adding a growth demand exposes that gap, but an unconditional demand may misrepresent condition-dependent cofactor use. The dilution mechanism and this limitation are established prior work. [Benyamini et al., 2010](https://pmc.ncbi.nlm.nih.gov/articles/PMC2884546/).
