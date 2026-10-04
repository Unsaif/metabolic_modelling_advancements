"""Physiological and diet constraints for whole-body models, ported from the COBRA Toolbox.

runIEM_HH.m (Harvey/Harvetta branch, toolbox commit 67c790d) prepares the model with
    sex = model.sex; standardPhysiolDefaultParameters;
    model = physiologicalConstraintsHMDBbased(model, IndividualParameters);
    EUAverageDietNew; model = setDietConstraints(model, Diet);
before the global reaction constraints and the IEM loop. This module ports those three steps for
the default call (Type 'HMDB', Biofluid 'all', no exclusion list, setDefault 1, ExclMet 0, factor 1).
The data inputs (HMDB concentration tables, blood-flow percentages, organ weights, the diet and
the AGORA essential-metabolite list) are read from data/iem/wbm_constraint_inputs_v0.3.json,
written by scripts/extract_wbm_constraint_inputs.py with the toolbox file hashes.

Arithmetic follows the MATLAB expressions term by term (same operand order), so that re-applying
the port to a model built by the Toolbox reproduces its bounds exactly where nothing changed.
Each function returns new bound arrays and a change log; inputs are never modified in place.
"""
from __future__ import annotations

import json
import os
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import scipy.sparse as sp

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_INPUTS = os.path.join(ROOT, "data", "iem", "wbm_constraint_inputs_v0.3.json")

# OrganLists.m, OrgansListExt
ORGANS_EXT = {
    "male": ["Heart", "Muscle", "Lung", "Skin", "Stomach", "sIEC", "Colon", "Urinarybladder", "Retina", "Scord",
             "Brain", "Adipocytes", "Liver", "Gall", "Kidney", "Pancreas", "Spleen", "Agland", "Thyroidgland",
             "Pthyroidgland", "Testis", "Prostate", "Bcells", "CD4Tcells", "Nkcells", "Monocyte", "Platelet", "RBC",
             "BBB", "Diet", "SI", "GI", "LI", "BileDuct", "Excretion"],
    "female": ["Heart", "Muscle", "Lung", "Skin", "Stomach", "sIEC", "Colon", "Urinarybladder", "Retina", "Scord",
               "Brain", "Adipocytes", "Liver", "Gall", "Kidney", "Pancreas", "Spleen", "Agland", "Thyroidgland",
               "Pthyroidgland", "Ovary", "Uterus", "Breast", "Cervix", "Bcells", "CD4Tcells", "Nkcells", "Monocyte",
               "Platelet", "RBC", "BBB", "Diet", "SI", "GI", "LI", "BileDuct", "Excretion"],
}

# physiologicalConstraintsHMDBbased.m: concentration thresholds (uM) for setting lower bounds
MIN_CONC_CONSTRAINT = 5
MAX_CONC_CONSTRAINT = 50
# kidney exchanges whose filtration (upper bound) is set to 0
KIDNEY_NO_FILTRATION = ["Kidney_EX_na1(e)_[bc]", "Kidney_EX_hco3(e)_[bc]", "Kidney_EX_urea(e)_[bc]",
                        "Kidney_EX_k(e)_[bc]", "Kidney_EX_cl(e)_[bc]", "Kidney_EX_ca2(e)_[bc]",
                        "Kidney_EX_HC02172(e)_[bc]", "Kidney_EX_avite1(e)_[bc]"]
# BBB uptakes left unconstrained ("otw infeasible - Nov 2017")
BBB_UNCONSTRAINED = ["BBB_NH4[CSF]upt", "BBB_CHOL[CSF]upt", "BBB_PI[CSF]upt", "BBB_STRDNC[CSF]upt",
                     "BBB_HC00250[CSF]upt", "BBB_PYDXN[CSF]upt", "BBB_5MTHF[CSF]upt", "BBB_SO3[CSF]upt"]
# CSF metabolites without an enforced export (lower bound 0)
CSF_NO_LOWER_BOUND = ["na1[csf]", "cl[csf]", "k[csf]", "h2o[csf]", "sucsal[csf]", "ca2[csf]", "ser_D[csf]"]
# urine exchanges without an enforced excretion (lower bound 0); the list in the source, duplicates included
URINE_NO_LOWER_BOUND = [
    "EX_na1[u]", "EX_cl[u]", "EX_k[u]", "EX_ca2[u]", "EX_C05767[u]", "EX_C05770[u]", "EX_C05302[u]", "EX_trypta[u]",
    "EX_ppbng[u]", "EX_13dampp[u]", "EX_mhista[u]", "EX_tym[u]", "EX_2hyoxplac[u]", "EX_pmtcrn[u]", "EX_dheas[u]",
    "EX_34dhphe[u]", "EX_srtn[u]", "EX_gthrd[u]", "EX_pcholhep_hs[u]", "EX_pcholste_hs[u]", "EX_pcholn204_hs[u]",
    "EX_3moxtyr[u]", "EX_aldstrn[u]", "EX_tststerone[u]", "EX_pydxn[u]", "EX_sphgn[u]", "EX_sphings[u]", "EX_csn[u]",
    "EX_arab_L[u]", "EX_tststerone[u]", "EX_tststerone[u]", "EX_pydxn[u]", "EX_mma[u]", "EX_tsul[u]", "EX_5htrp[u]",
    "EX_7dhchsterol", "EX_etoh[u]", "EX_gsn[u]", "EX_5aop[u]", "EX_uri[u]", "EX_dad_2[u]", "EX_ocdca[u]", "EX_gua[u]",
    "EX_dcyt[u]", "EX_glyleu[u]", "EX_acald[u]", "EX_HC02191[u]"]
CREATININE_MW = 113.1179  # g/mol

# setDietConstraints.m
DIET_MISSING_COMPOUNDS = ["Diet_EX_asn_L[d]", "Diet_EX_gln_L[d]", "Diet_EX_chol[d]", "Diet_EX_crn[d]",
                          "Diet_EX_elaid[d]", "Diet_EX_hdcea[d]", "Diet_EX_dlnlcg[d]", "Diet_EX_adrn[d]",
                          "Diet_EX_hco3[d]", "Diet_EX_sprm[d]", "Diet_EX_carn[d]", "Diet_EX_7thf[d]",
                          "Diet_EX_Lcystin[d]", "Diet_EX_hista[d]", "Diet_EX_orn[d]", "Diet_EX_ptrc[d]",
                          "Diet_EX_creat[d]"]
DIET_MISSING_COMPOUNDS_2 = ["Diet_EX_cytd[d]", "Diet_EX_so4[d]"]
DIET_MICRONUTRIENTS = ["Diet_EX_adocbl[d]", "Diet_EX_vitd2[d]", "Diet_EX_vitd3[d]", "Diet_EX_psyl[d]",
                       "Diet_EX_gum[d]", "Diet_EX_bglc[d]", "Diet_EX_phyQ[d]", "Diet_EX_fol[d]", "Diet_EX_5mthf[d]",
                       "Diet_EX_q10[d]", "Diet_EX_retinol_9_cis[d]", "Diet_EX_pydxn[d]", "Diet_EX_pydam[d]",
                       "Diet_EX_pydx[d]", "Diet_EX_pheme[d]", "Diet_EX_ribflv[d]", "Diet_EX_thm[d]",
                       "Diet_EX_avite1[d]", "Diet_EX_pnto_R[d]"]
DIET_IONS = ["Diet_EX_na1[d]", "Diet_EX_cl[d]", "Diet_EX_k[d]", "Diet_EX_pi[d]", "Diet_EX_zn2[d]", "Diet_EX_cu2[d]"]
DIET_SO4 = ["Diet_EX_so4[d]"]


def load_inputs(path: str = DEFAULT_INPUTS) -> dict:
    with open(path) as fh:
        return json.load(fh)


@dataclass
class IndividualParameters:
    """standardPhysiolDefaultParameters.m (reference man or woman)."""
    sex: str
    body_weight_kg: float
    organ_weights_g: Dict[str, Optional[float]]
    blood_flow_fraction: Dict[str, Optional[float]]   # organ -> fraction of cardiac output (None = empty cell)
    heart_rate: float = 67
    stroke_volume: float = 80
    hematocrit: float = 0.4
    creatinine_urine_max_mg_dl: float = 1.2
    creatinine_urine_min_mg_dl: float = 0.5
    mcon_default_bc: float = 20
    mcon_default_csf: float = 20
    mcon_default_ur_max: float = 20
    mcon_default_ur_min: float = 0
    csf_flow_rate: float = 0.35
    csf_blood_flow_rate: float = 0.52
    urine_flow_rate: float = 2000
    gfr: float = 90
    # Not a Toolbox parameter: the CSF export upper bound uses CSFBloodFlowRate in the current code and
    # CSFFlowRate in physiologicalConstraintsHMDBbased_old.m. None means the current code.
    csf_export_ub_flow_rate: Optional[float] = None

    @property
    def cardiac_output(self) -> float:
        return self.heart_rate * self.stroke_volume


def default_parameters(inputs: dict, sex: str) -> IndividualParameters:
    sex = sex.lower()
    if sex not in ("male", "female"):
        raise ValueError("sex must be 'male' or 'female'")
    ow = inputs["organ_weights"][sex]
    flow = {organ: rec[sex] for organ, rec in inputs["blood_flow"]["organs"].items()}
    return IndividualParameters(sex=sex, body_weight_kg=ow["body_weight_g"] / 1000,
                                organ_weights_g=_first_by_name(ow["organs"]),
                                blood_flow_fraction=flow)


def _first_by_name(organs: List[dict]) -> Dict[str, Optional[float]]:
    out: Dict[str, Optional[float]] = {}
    for o in organs:
        out.setdefault(o["name"], o["weight_g"])   # find(ismember(...)) then {1}: the first row
    return out


def legacy_gfr(params: IndividualParameters) -> float:
    """GFR as computed by physiologicalConstraintsHMDBbased_old.m: 20 percent of the renal plasma flow."""
    bk = params.blood_flow_fraction["Kidney"]
    renal_flow_rate = bk * params.cardiac_output * (1 - params.hematocrit)
    return renal_flow_rate * 0.2


def plasma_flow_rates(params: IndividualParameters) -> Dict[str, float]:
    """Plasma flow (ml/min) for each organ of OrgansListExt, as in the blood-flow loop of the Toolbox."""
    co, hct = params.cardiac_output, params.hematocrit
    out = {}
    for organ in ORGANS_EXT[params.sex]:
        if organ in params.blood_flow_fraction:
            b = params.blood_flow_fraction[organ]
            out[organ] = b * co * (1 - hct) if b is not None else 0.01 * co * (1 - hct)
        elif organ == "BBB":
            brain, scord = params.blood_flow_fraction.get("Brain"), params.blood_flow_fraction.get("Scord")
            if brain is None or scord is None:
                raise ValueError("BBB plasma flow needs Brain and Scord blood-flow fractions")
            out[organ] = (brain + scord) * co * (1 - hct)
        else:
            out[organ] = 0.01 * co * (1 - hct)
    return out


@dataclass
class ConstraintLog:
    changes: List[Tuple[str, str, str, float, float]] = field(default_factory=list)  # section, rxn, bound, old, new
    warnings: List[str] = field(default_factory=list)

    def set(self, arr: np.ndarray, i: int, value: float, section: str, rxn: str, bound: str) -> None:
        old = float(arr[i])
        arr[i] = value
        if not (old == value or (np.isnan(old) and np.isnan(value))):
            self.changes.append((section, rxn, bound, old, float(value)))

    def summary(self) -> Dict[str, int]:
        out: Dict[str, int] = {}
        for section, *_ in self.changes:
            out[section] = out.get(section, 0) + 1
        return out


class _Model:
    """Reaction ids, the stoichiometric matrix and working bound arrays."""

    def __init__(self, rxns: Sequence[str], mets: Sequence[str], S: sp.spmatrix, lb: np.ndarray, ub: np.ndarray):
        self.rxns = [str(r) for r in rxns]
        self.mets = [str(m) for m in mets]
        self.S = sp.csc_matrix(S)
        self.lb = np.array(lb, dtype=float).copy()
        self.ub = np.array(ub, dtype=float).copy()
        self.pos: Dict[str, int] = {}
        for i, r in enumerate(self.rxns):
            self.pos.setdefault(r, i)

    def mets_with_sign(self, j: int, positive: bool) -> List[str]:
        col = self.S[:, j]
        rows, vals = col.indices, col.data
        return [self.mets[r] for r, v in zip(rows, vals) if (v > 0 if positive else v < 0)]


def _first_index(table: List[dict], prefix: str, suffix: str) -> Dict[str, int]:
    """Key -> first row index (MATLAB uses the first matching row if an id occurs twice)."""
    out: Dict[str, int] = {}
    for k, r in enumerate(table):
        out.setdefault(prefix + r["vmh"] + suffix, k)
    return out


def _first_row(table: List[dict], keys: Dict[str, int], mets: Sequence[str]) -> Optional[dict]:
    """MATLAB: X = find(ismember(strcat(ids, suffix), ExM)); value = data{X, col} takes the first X."""
    hits = sorted(keys[m] for m in mets if m in keys)
    return table[hits[0]] if hits else None


def physiological_constraints(rxns, mets, S, lb, ub, params: IndividualParameters, inputs: dict
                              ) -> Tuple[np.ndarray, np.ndarray, ConstraintLog]:
    """physiologicalConstraintsHMDBbased(model, IndividualParameters) with Type 'HMDB'."""
    m = _Model(rxns, mets, S, lb, ub)
    log = ConstraintLog()
    blood, csf, urine = inputs["hmdb"]["blood"], inputs["hmdb"]["csf"], inputs["hmdb"]["urine"]
    blood_keys = _first_index(blood, "", "[bc]")
    blood_rbc_keys = _first_index(blood, "RBC_", "[bc]")
    csf_keys = _first_index(csf, "", "[csf]")
    urine_keys = _first_index(urine, "", "[u]")
    pf = plasma_flow_rates(params)
    gfr = params.gfr
    organs = ORGANS_EXT[params.sex]

    # --- blood: uptake limits from plasma flow; kidney filtration from GFR ------------------------------------
    csf_upt = [j for j, r in enumerate(m.rxns) if "[CSF]upt" in r]
    for organ in organs:
        if organ.startswith("BBB"):
            exr = csf_upt
        else:
            prefix = organ + "_EX_"
            exr = [j for j, r in enumerate(m.rxns) if r.startswith(prefix)]
        for j in exr:
            r = m.rxns[j]
            if not (("[bc]" in r or "[CSF]upt" in r) and "_o2s(e)" not in r and "_h2o(e)" not in r
                    and "_H2O[CSF]upt" not in r):
                continue
            exm = m.mets_with_sign(j, positive=True)
            if not exm:
                continue
            if organ == "Kidney" and "[bcK]" not in r:
                row = _first_row(blood, blood_keys, exm)
                mmin, mmax = (0, params.mcon_default_bc) if row is None else (row["min"], row["max"])
                rate = 1 * (mmin / 1000) * gfr * 60 * 24 / 1000
                if r not in KIDNEY_NO_FILTRATION:
                    log.set(m.ub, j, -rate if mmax >= MAX_CONC_CONSTRAINT else 0, "blood_kidney", r, "ub")
                else:
                    log.set(m.ub, j, 0, "blood_kidney", r, "ub")
                rate = 1 * (mmax / 1000) * gfr * 60 * 24 / 1000
                log.set(m.lb, j, -rate, "blood_kidney", r, "lb")
            else:
                row = _first_row(blood, blood_keys, exm)
                if row is None:
                    row = _first_row(blood, blood_rbc_keys, exm)
                mcon = params.mcon_default_bc if row is None else row["max"]
                if r in BBB_UNCONSTRAINED:
                    continue
                if m.lb[j] < 0:
                    if mcon > 1e-3:
                        rate = ((mcon) / 1000) * pf[organ] * 60 * 24 / 1000
                    else:
                        rate = ((1e-3) / 1000) * pf[organ] * 60 * 24 / 1000
                    log.set(m.lb, j, -1 * rate, "blood_uptake", r, "lb")

    # --- CSF export across the blood-brain barrier (the Toolbox repeats this loop once per organ; it is
    # idempotent, so one pass gives the same bounds) -------------------------------------------------------------
    for j, r in enumerate(m.rxns):
        if not (r.startswith("BBB_") and "[CSF]" in r and "exp" in r and "_o2(e)" not in r and "_o2s(e)" not in r
                and "_co2(e)" not in r):
            continue
        exm = m.mets_with_sign(j, positive=False)
        if not exm:
            continue
        row = _first_row(csf, csf_keys, exm)
        mmin, mmax = (0, params.mcon_default_csf) if row is None else (row["min"], row["max"])
        lb_rate = (mmin / 1000) * params.csf_blood_flow_rate * 60 * 24 / 1000
        if mmin >= MIN_CONC_CONSTRAINT and mmax >= MAX_CONC_CONSTRAINT and not any(x in exm for x in CSF_NO_LOWER_BOUND):
            log.set(m.lb, j, lb_rate, "csf_export", r, "lb")
        else:
            log.set(m.lb, j, 0, "csf_export", r, "lb")
        ub_flow = params.csf_blood_flow_rate if params.csf_export_ub_flow_rate is None else params.csf_export_ub_flow_rate
        log.set(m.ub, j, (mmax / 1000) * ub_flow * 60 * 24 / 1000, "csf_export", r, "ub")

    # --- urine excretion ------------------------------------------------------------------------------------
    cr_max = params.creatinine_urine_max_mg_dl * 10 / CREATININE_MW
    cr_min = params.creatinine_urine_min_mg_dl * 10 / CREATININE_MW
    for j, r in enumerate(m.rxns):
        if not (r.startswith("EX_") and "[u]" in r and "_o2(e)" not in r and "_o2s(e)" not in r and "_co2(e)" not in r):
            continue
        exm = m.mets_with_sign(j, positive=False)
        if not exm:
            continue
        row = _first_row(urine, urine_keys, exm)
        mmin, mmax = (params.mcon_default_ur_min, params.mcon_default_ur_max) if row is None else (row["min"], row["max"])
        lb_rate = (mmin / 1000) * cr_min * params.urine_flow_rate / 1000
        ub_rate = (mmax / 1000) * cr_max * params.urine_flow_rate / 1000
        if m.ub[j] > 0:
            if mmax >= MAX_CONC_CONSTRAINT and mmin >= MIN_CONC_CONSTRAINT and r not in URINE_NO_LOWER_BOUND:
                log.set(m.lb, j, lb_rate, "urine", r, "lb")
            else:
                log.set(m.lb, j, 0, "urine", r, "lb")
            log.set(m.ub, j, ub_rate, "urine", r, "ub")

    # --- no milk production -------------------------------------------------------------------------------
    for j, r in enumerate(m.rxns):
        if "(miB)_[mi]" in r:
            log.set(m.lb, j, 0, "milk", r, "lb")
            log.set(m.ub, j, 0, "milk", r, "ub")

    # --- literature constraints (changeRxnBounds warns and skips a missing reaction) --------------------
    def change(rxn: str, value: float, bound: str) -> None:
        j = m.pos.get(rxn)
        if j is None:
            log.warnings.append(f"Reaction {rxn} not in model")
            return
        if bound in ("l", "b"):
            log.set(m.lb, j, value, "literature", rxn, "lb")
        if bound in ("u", "b"):
            log.set(m.ub, j, value, "literature", rxn, "ub")

    def idx(rxn: str) -> int:
        if rxn not in m.pos:   # the MATLAB comparisons below would fail on an empty index
            raise KeyError(f"{rxn} is required by physiologicalConstraintsHMDBbased")
        return m.pos[rxn]

    change("EX_o2[a]", -15000, "u"); change("EX_o2[a]", -25000, "l")
    change("EX_co2[a]", 15000 * 0.8, "l"); change("EX_co2[a]", 25000, "u")
    change("EX_h2o[a]", 47182 * 0.8, "l"); change("EX_h2o[a]", 47182 * 1.2, "u")
    change("EX_h2o[sw]", 36080 * 0.8, "l"); change("EX_h2o[sw]", 36080 * 1.2, "u")
    change("EX_h2o[u]", 77711 * 0.8, "l"); change("EX_h2o[u]", 77711 * 1.2, "u")
    change("Excretion_EX_h2o[fe]", 5550 * 0.8, "l"); change("Excretion_EX_h2o[fe]", 5550 * 1.2, "u")
    change("Gall_H2Ot[bdG]", 1000, "u"); change("Liver_H2Ot[bdL]", 1000, "u")

    bw = params.body_weight_kg

    def daily_mmol(mg_per_min: float, mw: float) -> float:
        met = (mg_per_min * 60 * 24 * bw / 65) / 1000
        return met * 1000 / mw

    j = idx("Muscle_EX_glc_D(e)_[bc]")
    if m.ub[j] >= -0.01 * 1000 and m.lb[j] <= -0.01 * 1000:
        change("Muscle_EX_glc_D(e)_[bc]", -0.01 * 1000, "u")
    elif m.lb[j] > m.ub[j]:
        log.set(m.ub, j, 0, "literature", m.rxns[j], "ub")

    met = daily_mmol(12.5, 89.09)
    j = idx("Muscle_EX_ala_L(e)_[bc]")
    if m.lb[j] < met * 0.80 and m.ub[j] >= met * 0.8:
        change("Muscle_EX_ala_L(e)_[bc]", met * 0.8, "l"); change("Muscle_EX_ala_L(e)_[bc]", met * 1.2, "u")
    elif m.lb[j] > m.ub[j]:
        log.set(m.lb, j, 0, "literature", m.rxns[j], "lb")

    met = daily_mmol(25, 180.16)
    j = idx("RBC_EX_glc_D(e)_[bc]")
    if m.lb[j] < -met * 1.2 and m.ub[j] >= -met * 1.2:
        change("RBC_EX_glc_D(e)_[bc]", -met * 0.8, "u"); change("RBC_EX_glc_D(e)_[bc]", -met * 1.2, "l")
    elif m.lb[j] > m.ub[j]:
        log.set(m.ub, j, 0, "literature", m.rxns[j], "ub")

    met = daily_mmol(80, 180.16)
    j = idx("Brain_EX_glc_D(e)_[csf]")
    if m.lb[j] < -met * 1.20 and m.ub[j] >= -met * 1.2:
        change("Brain_EX_glc_D(e)_[csf]", -met * 0.8, "u"); change("Brain_EX_glc_D(e)_[csf]", -met * 1.2, "l")
    elif m.lb[j] > m.ub[j]:
        log.set(m.ub, j, 0, "literature", m.rxns[j], "ub")

    brain_weight = params.organ_weights_g["Brain"]
    brain_o2 = (156 * 60 * 24 * brain_weight / 100) / 1000
    j = idx("Brain_EX_o2(e)_[csf]")
    if m.lb[j] < -brain_o2 * 1.2 and m.ub[j] >= -brain_o2 * 1.2:
        change("Brain_EX_o2(e)_[csf]", -brain_o2 * 1.2, "l"); change("Brain_EX_o2(e)_[csf]", -brain_o2 * 0.7, "u")
    elif m.lb[j] > m.ub[j]:
        log.set(m.ub, j, 0, "literature", m.rxns[j], "ub")

    met = daily_mmol(12.5, 89.09)
    j = idx("Liver_EX_ala_L(e)_[bc]")
    if m.lb[j] < -met * 1.20 and m.ub[j] >= -met * 1.2:
        change("Liver_EX_ala_L(e)_[bc]", -met * 0.8, "u"); change("Liver_EX_ala_L(e)_[bc]", -met * 1.2, "l")
    elif m.lb[j] > m.ub[j]:
        log.set(m.ub, j, 0, "literature", m.rxns[j], "ub")

    met = daily_mmol(130, 180.16)
    j = idx("Liver_EX_glc_D(e)_[bc]")
    if m.lb[j] < met * 0.80 and m.ub[j] >= met * 0.8:
        change("Liver_EX_glc_D(e)_[bc]", met * 1.2, "u"); change("Liver_EX_glc_D(e)_[bc]", met * 0.8, "l")
    elif m.lb[j] > m.ub[j]:
        log.set(m.lb, j, 0, "literature", m.rxns[j], "lb")

    met = daily_mmol(12, 92.09)
    j = idx("Adipocytes_EX_glyc(e)_[bc]")
    if m.lb[j] < met * 0.80 and m.ub[j] >= met * 0.8:
        change("Adipocytes_EX_glyc(e)_[bc]", met * 1.2, "u")
    elif m.lb[j] > m.ub[j]:
        log.set(m.lb, j, 0, "literature", m.rxns[j], "lb")

    # no carbon fixation, except brain, liver, lung and kidney
    for j, r in enumerate(m.rxns):
        if "EX_co2(e)_[bc]" in r:
            log.set(m.lb, j, 0, "literature", r, "lb")
    change("Brain_EX_co2(e)_[csf]", -10000, "l")
    change("Liver_EX_co2(e)_[bc]", -10000, "l")
    change("Lung_EX_co2(e)_[bc]", -10000, "l")
    change("Kidney_EX_co2(e)_[bc]", -10000, "l")
    change("Brain_DM_atp_c_", 0, "l")
    change("Heart_DM_atp_c_", 6000, "l")
    return m.lb, m.ub, log


def diet_constraints(rxns, lb, ub, diet: Sequence[Sequence[str]], agora_essential: Sequence[str], factor: float = 1
                     ) -> Tuple[np.ndarray, np.ndarray, ConstraintLog]:
    """setDietConstraints(model, Diet) with microbiotaEnabling = 1; Diet as [rxn, value text] pairs."""
    rxns = [str(r) for r in rxns]
    lb = np.array(lb, dtype=float).copy(); ub = np.array(ub, dtype=float).copy()
    log = ConstraintLog()
    where: Dict[str, List[int]] = {}
    for i, r in enumerate(rxns):
        where.setdefault(r, []).append(i)

    def set_bound(arr, names, value, section, bound):
        for name in names:
            for i in where.get(name, []):
                log.set(arr, i, value, section, name, bound)

    diet_rxns = [i for i, r in enumerate(rxns) if r.startswith("Diet_EX_")]
    for i in diet_rxns:
        log.set(lb, i, 0, "diet_reset", rxns[i], "lb")
        log.set(ub, i, 0, "diet_reset", rxns[i], "ub")
    in_diet = {row[0] for row in diet}
    missing_uptakes = sorted(set(agora_essential) - in_diet)
    set_bound(lb, missing_uptakes, -0.1, "diet_agora_essential", "lb")
    set_bound(lb, DIET_MISSING_COMPOUNDS, -50, "diet_missing_compounds", "lb")
    set_bound(lb, DIET_MISSING_COMPOUNDS_2, -50, "diet_missing_compounds", "lb")
    set_bound(lb, ["Diet_EX_chol[d]"], -41.251, "diet_missing_compounds", "lb")
    for rxn, text in diet:
        v = float(text)
        if rxn in DIET_MICRONUTRIENTS and v <= 10:
            value = -10 * factor
        elif rxn in DIET_MICRONUTRIENTS and v > 0.1:
            value = -1.2 * v * 100 * factor
        elif rxn in DIET_IONS:
            value = -1.2 * v * 100 * factor
        elif rxn in DIET_SO4:
            value = -1000 * factor
        else:
            value = -1.2 * v * factor
        set_bound(lb, [rxn], value, "diet", "lb")
    for rxn, text in diet:
        set_bound(ub, [rxn], -0.8 * float(text) * factor, "diet", "ub")
    absent = [rxn for rxn, _ in diet if rxn not in where]
    if absent:
        log.warnings.append(f"diet reactions absent from the model: {absent}")
    return lb, ub, log


def runiem_model_setup(rxns, mets, S, lb, ub, sex: str, inputs: Optional[dict] = None
                       ) -> Tuple[np.ndarray, np.ndarray, Dict[str, object]]:
    """The Harvey/Harvetta preparation in runIEM_HH.m: physiological constraints, then the EU average diet."""
    inputs = inputs or load_inputs()
    params = default_parameters(inputs, sex)
    lb1, ub1, log1 = physiological_constraints(rxns, mets, S, lb, ub, params, inputs)
    lb2, ub2, log2 = diet_constraints(rxns, lb1, ub1, inputs["diet_eu_average"], inputs["agora_essential"])
    report = {"physiology_changes": log1.summary(), "diet_changes": log2.summary(),
              "warnings": log1.warnings + log2.warnings, "changes": log1.changes + log2.changes,
              "parameters": {"sex": params.sex, "body_weight_kg": params.body_weight_kg,
                             "brain_weight_g": params.organ_weights_g["Brain"], "cardiac_output_ml_min": params.cardiac_output,
                             "toolbox_commit": inputs["provenance"]["toolbox_commit"]}}
    return lb2, ub2, report
