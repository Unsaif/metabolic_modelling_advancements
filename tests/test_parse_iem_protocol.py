import json
from pathlib import Path

from scripts.parse_iem_protocol import extract_protocol, strip_matlab_comments


def test_matlab_comments_preserve_quoted_percent_and_escaped_quote():
    text = "'100% deficiency' % ignored\n'name''s marker' % comment\n%{\nblock comment\n%}\n"
    assert strip_matlab_comments(text).splitlines()[:2] == ["'100% deficiency' ", "'name''s marker' "]
    assert "block comment" not in strip_matlab_comments(text)


def test_disabled_iems_commented_biomarkers_and_recon_branch_are_excluded():
    source = """
if 1
  if 0 % disease excluded by source author
    model = modelO;
    R = '_DISABLED';
    for i = 1:3
      if 1
        x = i;
      end
    end
    BiomarkerRxns = {'EX_no' 'Increased'};
    [IEMSol_DISABLED] = checkIEM_WBM(model,IEMRxns,BiomarkerRxns,minRxnsFluxHealthy);
  end
  model = modelO;
  R = '_ACTIVE';
  if ~strcmp(modelName,'Recon3D')
    model = addDemandReaction(model, 'b[bc]');
    model = addDemandReaction(model, 'extra[bc]');
    BiomarkerRxns = {
      'DM_b[bc]' 'Increased (blood)'
      % 'EX_comment' 'Decreased'
    };
  else
    BiomarkerRxns = {'EX_recon' 'Decreased'};
  end
  % [IEMSol_COMMENT] = checkIEM_WBM(model,IEMRxns,BiomarkerRxns,minRxnsFluxHealthy);
  [IEMSol_ACTIVE] = checkIEM_WBM(model,IEMRxns,BiomarkerRxns,minRxnsFluxHealthy);
end
"""
    protocol = extract_protocol(source)
    assert len(protocol) == 1
    block = protocol[0]
    assert block["iem"] == "ACTIVE"
    assert block["call_index"] == 2  # retains identity of source call, not filtered ordinal
    assert block["include_patterns"] == ["_ACTIVE"]
    assert block["biomarkers"] == [("DM_b[bc]", "Increased (blood)")]
    assert block["demand_metabolites"] == ["b[bc]", "extra[bc]"]


def test_disabled_branch_else_is_active():
    source = """
if 0
  R = '_BAD';
else
  R = '_GOOD';
  BiomarkerRxns = {'EX_good' 'Increased'};
  [IEMSol_GOOD] = checkIEM_WBM(model,IEMRxns,BiomarkerRxns,minRxnsFluxHealthy);
end
"""
    assert extract_protocol(source)[0]["include_patterns"] == ["_GOOD"]


def test_corrected_protocol_inventory_and_provenance():
    root = Path(__file__).resolve().parents[1]
    protocol = json.loads((root / "data/iem/iem_protocol_v0.2.json").read_text())
    provenance = json.loads((root / "data/iem/iem_protocol_v0.2.provenance.json").read_text())
    assert len(protocol) == 57
    assert sum(len(p["biomarkers"]) for p in protocol) == 252
    assert len({p["iem"] for p in protocol}) == 57
    assert not {"CRFD", "GSD6", "HFI", "OROA", "MNGIE"} & {p["iem"] for p in protocol}
    assert next(p for p in protocol if p["iem"] == "3MGA")["biomarkers"] == [["EX_3ivcrn[u]", "Increased (urine)"]]
    assert provenance["source_git_blob_sha1"] == "31d3d27cc297e4e59feefe2c9b2018d5e3688e58"
    assert next(p for p in protocol if p["iem"] == "EP")["bound_tweaks"] == [
        {"pattern": "_XYLUR", "bound": "lb", "value": 0.0},
        {"pattern": "_r0784", "bound": "lb", "value": 0.0},
        {"pattern": "_r0784", "bound": "ub", "value": 0.0},
    ]
