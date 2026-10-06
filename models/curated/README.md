# Published hand-curated models (Paper 1 reference comparison)

Downloaded 6 October 2026 through the browser on Tim's computer (the cloud workspace cannot reach these hosts) and
checked by sha256. Each `*_bigg_view.xml.gz` is the same model with BiGG identifiers
(`scripts/translate_curated_model.py`, relabelling only; `*_translation.json` lists every change).

| Directory | Model | Source | Licence |
|---|---|---|---|
| `iSO783/` | iSO783, *S. oneidensis* MR-1 (Pinchuk et al. 2010, PLoS Comput Biol 6:e1000822) | BioModels MODEL1507180036 (`MODEL1507180036_url.xml`) | BioModels: CC0 |
| `iGD1575/` | iGD1575, *S. meliloti* 1021 (diCenzo et al. 2016, Nat Commun 7:12219) | Supplementary Data 6 (`41467_2016_BFncomms12219_MOESM594_ESM.zip`, contains `ncomms12219-s7.xml`) | article and supplement CC BY 4.0 |
| `iAH991/` | iAH991, *B. thetaiotaomicron* VPI-5482 (Heinken et al. 2013, Gut Microbes 4:28-40) | No public SBML. Rebuilt from Supplementary Tables S10a/S10b of the supplement `2012GUTMICROBES0043R-Sup.pdf` (Taylor & Francis; supplied by Tim on 6 October 2026; sha256 49502d765c7d049208aeb03a140aa7973d4569b7cddf22b9ca2ea4fa019975b9; not redistributed here) with `scripts/rebuild_iAH991_from_pdf.py`; see `iAH991/REBUILD.md` for the checks against the paper's own predictions | rebuilt model shared for research comparison; credit the original authors |
