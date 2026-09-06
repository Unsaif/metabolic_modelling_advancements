# One-condition numerical diagnosis after the v2 failure

The first curated case, L-Arginine on MOPS, was reconstructed from the frozen
source loader, recipe and v2 preparation helper. Its stoichiometry, bounds,
objective, GPRs and chemical metadata match the legacy formulation exactly.
Four methods and their settings were recorded before these diagnostic solves.
No numeric fitness data were read. The primary failed run remains unchanged.

| Method | Status | Maximum absolute S·v residual | Original 1e−8 gate |
|---|---|---:|---|
| Native GLPK, already auto-scaled | Optimal | 7.13e−14 | Pass |
| GLPK with advanced basis | Optimal | 7.14e−14 | Pass |
| HiGHS simplex | Optimal | 7.18e−8 | Fail |
| HiGHS IPM | 60-second time limit | No valid primal vector | Fail |

The reproduced HiGHS simplex vector violates the original gate in `h_c`
(7.1820e−8) and `adp_c` (1.9215e−8). Independent compensated summation gives
7.3009e−8 and 1.9111e−8 respectively, so changing the accumulation method alone
does not make the vector acceptable. HiGHS reports an internal maximum primal
infeasibility of 7.7221e−11; that reported scalar does not replace checking the
exported, unscaled vector against the original equations.

Its vector contains many opposing fluxes at the inherited ±999999 bounds.
The proton row's absolute term sum is approximately 1.6004e7, making small
absolute residuals sensitive to cancellation and numerical accuracy. The native
GLPK vector's maximum absolute flux is only 46.2567. Both methods' objective
values differ by 2.9974e−9, within the objective-agreement tolerance; objective
agreement is insufficient to accept an infeasible primal vector.

Each valid returned vector contains all 2,935 reaction fluxes. The records retain
every metabolite residual, the largest contributing terms, bound violations and
large fluxes. IPM returned no valid primal vector, which is recorded explicitly.
All 164 input records and eight pre-solve artifacts remained unchanged, including
the 159 frozen study inputs. The copied diagnostic source and runtime versions
are included. No successful alternative HiGHS method is established by this
first diagnostic, and success on one condition would not establish reliability
across the remaining panel.
