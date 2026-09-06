# Sources and notices for the quinone repair models

This attribution covers `runs/**/model.xml.gz` below this study directory, including completed and preserved attempts. It was added on 6 September 2026 after the numerical inputs were frozen. No SBML, compressed-model, source-code, recipe, or manifest bytes were changed to add this documentation.

The generated files are research variants of the prepared Putida draft. They are not unmodified iJN1463 releases. Their precise interventions, gene aliases, and metadata normalizations are recorded in each arm's `intervention.json`; the corresponding attempt's `manifest.json` records frozen input hashes and runtime versions.

| Component | Source and attribution |
|---|---|
| Prepared baseline | Local `models/gapfilled/Putida.xml.gz`, SHA256 `378f67c3c8d4642c135942917cd46f3b6109c53c16f7987c9751fb643ade5739`, followed by the frozen universe v0.1, model v0.4, and GPR v0.2/v0.3/v0.4 patches. |
| Earlier draft and gap filling | `models/gapfilled/Putida_gapfill.json` records the historical EMBL/CarveMe draft `models/embl/Pseudomonas_putida_KT2440.xml.gz` (SHA256 `5b69bf07cf15fdecb7cf5a6c87e35271cd94e97c7d16747aee62f4c7594fd7c1`) and CarveMe universe `external/carveme/universe_bacteria.xml.gz` (SHA256 `b983aec83c4f6c30ec58ded6684d15e90355e7d6c9d573747dde9d0bbae35d19`). These are distinct provenance layers; the BiGG notice below is not asserted to replace their applicable terms. |
| Transferred quinone content | Local `models/bigg/iJN1463.xml`, SHA256 `d573833328ffae0dfa752a1fa3262ed939ed5862288beab287fca30d0fefb4a1`, downloaded through the [BiGG iJN1463 model page](https://bigg.ucsd.edu/models/iJN1463). The exact five source reactions are `OHPHM`, `OMPHHX`, `OMBZLM`, `OMMBLHX`, and `DMQMT`; transfer arms retain their source stoichiometry and bounds. Source gene aliases and the declared sequence-informed overrides are separately recorded. |
| Biomass coefficients | The template coefficient comes from the frozen gap-filled Putida biomass. The curated coefficient, used only in its named arm, comes from `BIOMASS_KT2440_WT3` in the frozen BiGG source. |

The source model's [official license page](https://bigg.ucsd.edu/license), linked directly from its model page, grants educational, research and non-profit uses subject to retaining the specified notice. It does **not** state CC-BY-SA 4.0. The complete official notice is included below and in [the source model's NOTICE](../../models/bigg/NOTICE). When distributing an individual generated model separately from this repository, include this document and the complete notice with the copy. No broader license or commercial permission for the composite model is claimed by this attribution.

Please cite the BiGG resource: King ZA et al. (2016), *BiGG Models: A platform for integrating, standardizing, and sharing genome-scale models*, Nucleic Acids Research 44(D1):D515-D522, [doi:10.1093/nar/gkv1049](https://doi.org/10.1093/nar/gkv1049).

The related reconstruction paper is Nogales J et al. (2020), *High-quality genome-scale metabolic modelling of Pseudomonas putida highlights its broad metabolic capabilities*, Environmental Microbiology 22(1):255-269, [doi:10.1111/1462-2920.14843](https://doi.org/10.1111/1462-2920.14843). The publisher calls its supplementary reconstruction **iJN1462**. Its article is marked open access; an alternative license for the exact BiGG iJN1463 bytes was not established in this check. No article-access label is substituted for the BiGG download terms.

The derivatives preserve existing target metabolite metadata and normalize newly copied cytosol labels from `c` to `C_c`. They therefore are not source-identical metadata copies. The studies and their results remain provisional development analyses; attribution does not imply endorsement by BiGG or the reconstruction authors.

## Required BiGG notice

The following notice was extracted directly from [the official license page](https://bigg.ucsd.edu/license) on 6 September 2026. HTML entities and paragraph whitespace were normalized; its wording was preserved.

Copyright © 2019 The Regents of the University of California

All Rights Reserved

Permission to use, copy, modify and distribute any part of BiGG Models for educational, research and non-profit purposes, without fee, and without a written agreement is hereby granted, provided that the above copyright notice, this paragraph and the following three paragraphs appear in all copies.

Those desiring to incorporate BiGG Models into commercial products or use for commercial purposes should contact the Technology Transfer & Intellectual Property Services, University of California, San Diego, 9500 Gilman Drive, Mail Code 0910, La Jolla, CA 92093-0910, Ph: (858) 534-5815, FAX: (858) 534-7345, e-mail: invent@ucsd.edu.

In no event shall the University of California be liable to any party for direct, indirect, special, incidental, or consequential damages, including lost profits, arising out of the use of this bigg database, even if the University of California has been advised of the possibility of such damage.

The BiGG Models provided herein is on an "as is" basis, and the University of California has no obligation to provide maintenance, support, updates, enhancements, or modifications. The University of California makes no representations and extends no warranties of any kind, either implied or express, including, but not limited to, the implied warranties of merchantability or fitness for a particular purpose, or that the use of the BiGG Models will not infringe any patent, trademark or other rights.
