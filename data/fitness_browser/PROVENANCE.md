# Fitness Browser downloads

Source: Fitness Browser, Arkin lab / LBNL (Price et al. 2018, Nature 557:503-509), https://fit.genomics.lbl.gov
Downloaded: 5 September 2026, through the Claude desktop browser pane, one file per link on each organism page:

- `fit_logratios.tsv`  — https://fit.genomics.lbl.gov/cgi-bin/createFitData.cgi?orgId=<orgId>          (gene fitness, log2 ratios; genes x experiments)
- `fit_t.tsv`          — https://fit.genomics.lbl.gov/cgi-bin/createFitData.cgi?orgId=<orgId>&t=1      (t-like statistics, same layout)
- `experiments.tsv`    — https://fit.genomics.lbl.gov/cgi-bin/createExpData.cgi?orgId=<orgId>          (experiment metadata: group, media, condition, concentration, quality metrics)
- `genes.tsv`          — https://fit.genomics.lbl.gov/cgi-bin/orgGenes.cgi?orgId=<orgId>               (locus table)
- `specific_phenotypes.tsv` — https://fit.genomics.lbl.gov/cgi-bin/spec.cgi?orgId=<orgId>&download=1 (genes with specific phenotypes)

| orgId | organism (Fitness Browser name, Sept 2026) | why it is here |
|---|---|---|
| Keio | Escherichia coli BW25113 | reference organism; Bernstein et al. 2023 regression test (their copy is from 2022; this is current) |
| Btheta | Bacteroides thetaiotaomicron VPI-5482 | gut commensal; in AGORA2 and with a curated model (iAH991); 519 experiments |
| Bvulgatus_CL09T03C04 | Phocaeicola vulgatus CL09T03C04 | gut commensal; AGORA2 has P. vulgatus strains; 22 experiments |
| Koxy | Klebsiella grimontii M5al (formerly K. michiganensis) | gut-associated Enterobacteriaceae; 208 experiments |
| Putida | Aquipseudomonas alloputida KT2440 (Pseudomonas putida KT2440) | curated BiGG model iJN1463; 314 experiments; ChatGEM's organism |
| MR1 | Shewanella oneidensis MR-1 | curated model iSO783; 176 experiments |
| SynE | Synechococcus elongatus PCC 7942 | curated BiGG model iJB785 (photoautotroph); 129 experiments |
| Smeli | Sinorhizobium meliloti 1021 | curated model iGD1575; 90 experiments |
| DvH | Nitratidesulfovibrio vulgaris Hildenborough JW710 | curated model iJF744; 757 experiments |

Terms: the Fitness Browser data are publicly available; cite Price et al. 2018 and check https://fit.genomics.lbl.gov/cgi-bin/help.cgi before redistribution.
Not downloaded for every organism: t scores and specific phenotypes (fetched for Keio, Btheta, Bvulgatus, Koxy, Putida only).
