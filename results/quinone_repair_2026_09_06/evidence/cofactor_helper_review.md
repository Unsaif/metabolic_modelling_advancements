# Independent review of the cofactor helper

The review identified and corrected two issues before freezing the repair study:

1. A zero-bounded solver balance did not guarantee that its expression matched the stored stoichiometric row. A custom supply variable added directly to a quinone balance produced flux 1 while the original helper reported a conserved pool from reaction metadata. The certificate now rejects altered solver-row coefficients or additional balance terms. Additional constraints elsewhere remain permitted.
2. The saturation flag used the production classification threshold as its distance from capacity. When threshold equaled capacity, a blocked zero was reported as reaching capacity. Saturation now uses a separate relative comparison with no positive absolute tolerance.

The 20 tests in `tests/test_cofactor.py` pass. Added regressions cover custom solver supply, modified balance coefficients, preservation of additional constraints and a custom objective, nonfinite solution cleanup, and the saturation edge case. The prepared Putida baseline still returns an exact conserved-pool certificate for `q8_c + q8h2_c`, involving 10 reactions.

The production probe keeps the current medium, maintenance, growth bounds, and additional constraints. Each demand is removed before the next probe, including on failure. Nonoptimal or nonfinite results raise an error rather than becoming zero production. The helper distinguishes a direction-compatible source reaction from feasible production: a nonconserved pool can still be blocked by precursor or medium constraints. Its conclusions remain properties of the declared model and constraints, without claims of physiological cofactor abundance or gene essentiality.
