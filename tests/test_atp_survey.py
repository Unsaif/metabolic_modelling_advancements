"""Offline checks for survey evidence, failure denominators, and provenance."""
import csv
import gzip
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout
from unittest.mock import patch

from gembench.survey import screen_sbml, summarize_rows

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("atp_survey_cli", ROOT / "scripts/survey_embl_atp_synthase.py")
cli = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cli)


def model(reaction_id="renamed", name="", reverse=False, extra=""):
    substrate = '<speciesReference species="M_adp_c"/><speciesReference species="M_pi_c"/><speciesReference species="M_h_p" stoichiometry="4"/>'
    product = '<speciesReference species="M_atp_c"/><speciesReference species="M_h2o_c"/><speciesReference species="M_h_c" stoichiometry="3"/>'
    if reverse:
        substrate, product = product, substrate
    return f'''<sbml xmlns="http://www.sbml.org/sbml/level3/version1/core"><model>
    <notes><p id="R_ATPS_NOTE">irrelevant</p></notes>
    <listOfReactions><reaction id="{reaction_id}" name="{name}">
    <listOfReactants>{substrate}{extra}</listOfReactants><listOfProducts>{product}</listOfProducts>
    </reaction></listOfReactions></model></sbml>'''.encode()


class SurveyTests(unittest.TestCase):
    def test_renamed_reaction_is_an_equation_candidate_not_identifier_hit(self):
        result = screen_sbml(model())
        self.assertFalse(result["has_atps"])
        self.assertEqual(result["atps_ids"], "")
        self.assertEqual(result["ion_coupled_atp_candidates"], "renamed")
        self.assertEqual(result["n_reactions"], 1)

    def test_names_and_ids_are_separate_evidence_and_reverse_equation_is_detected(self):
        result = screen_sbml(model("R_ATPS4rpp", "ATP synthase", reverse=True))
        self.assertEqual(result["atps_ids"], "R_ATPS4rpp")
        self.assertEqual(result["atp_synthase_name_candidates"], "R_ATPS4rpp")
        self.assertEqual(result["ion_coupled_atp_candidates"], "R_ATPS4rpp")

    def test_other_solute_transport_does_not_count_as_atp_synthase_equation(self):
        result = screen_sbml(model(extra='<speciesReference species="M_glc__D_p"/>'))
        self.assertEqual(result["ion_coupled_atp_candidates"], "")

    def test_namespace_prefix_and_single_quotes_do_not_break_identifier_detection(self):
        xml = b"<s:sbml xmlns:s='urn:sbml'><s:model><s:listOfReactions><s:reaction id='ATPS4rpp'/></s:listOfReactions></s:model></s:sbml>"
        self.assertEqual(screen_sbml(xml)["atps_ids"], "ATPS4rpp")

    def test_non_sbml_or_empty_document_is_not_an_absence(self):
        for xml in (b"<html/>", b"<sbml><model/></sbml>"):
            with self.subTest(xml=xml), self.assertRaises(ValueError):
                screen_sbml(xml)

    def test_failures_are_excluded_from_denominator(self):
        rows = [
            dict(atps_ids="R_ATPS4rpp;", bytes="123", n_reactions="1"),
            dict(atps_ids="", bytes="123", n_reactions="1"),
            dict(atps_ids="DOWNLOAD_FAILED: timeout", bytes="0", n_reactions="0"),
            dict(atps_ids="", bytes="123", n_reactions="0", status="failed"),
        ]
        result = summarize_rows(rows)
        self.assertEqual(result["successful_models"], 2)
        self.assertEqual(result["failed_models"], 2)
        self.assertEqual(result["fraction_without_atps_identifier"], 0.5)
        self.assertIsNone(summarize_rows(rows[2:])["fraction_without_atps_identifier"])

    def test_archived_surveys_reproduce_52_present_448_absent_and_4_present_29_absent(self):
        for filename, present, absent in [("embl_atp_synthase_survey_sample500.tsv", 52, 448),
                                          ("embl_atp_synthase_survey.tsv", 4, 29)]:
            with self.subTest(filename=filename), (ROOT / "results" / filename).open() as stream:
                summary = summarize_rows(list(csv.DictReader(stream, delimiter="\t")))
                self.assertEqual(summary["with_atps_identifier"], present)
                self.assertEqual(summary["without_atps_identifier"], absent)
                self.assertEqual(summary["failed_models"], 0)

    def test_offline_summary_never_downloads(self):
        filename = ROOT / "results/embl_atp_synthase_survey_sample500.tsv"
        with patch.object(cli.sys, "argv", ["survey", "--summarize", str(filename)]), \
             patch.object(cli.urllib.request, "urlopen", side_effect=AssertionError("network")), \
             redirect_stdout(io.StringIO()) as output:
            cli.main()
        self.assertEqual(json.loads(output.getvalue())["without_atps_identifier"], 448)

    def test_new_survey_pins_sources_hashes_files_and_records_parse_failure(self):
        manifest = b"assembly_accession\ttaxid\torganism_name\tfile_path\nA\t1\tSpecies one\tmodels/a.xml.gz\nB\t2\tSpecies two\tmodels/b.xml.gz\nC\t3\tSpecies three\tmodels/c.xml.gz\n"
        files = {"model_list.tsv": manifest, "models/a.xml.gz": gzip.compress(model("R_ATPS4rpp")),
                 "models/b.xml.gz": gzip.compress(model("other", extra='<speciesReference species="M_glc__D_p"/>')),
                 "models/c.xml.gz": gzip.compress(b"<html/>")}
        urls = []

        def download(url, timeout):
            urls.append(url)
            prefix = f"https://raw.githubusercontent.com/cdanielmachado/embl_gems/{cli.DEFAULT_REF}/"
            self.assertTrue(url.startswith(prefix))
            return io.BytesIO(files[url.removeprefix(prefix)])

        with tempfile.TemporaryDirectory() as directory:
            out = Path(directory) / "screen.tsv"
            with patch.object(cli.sys, "argv", ["survey", "--n", "3", "--output", str(out)]), \
                 patch.object(cli.urllib.request, "urlopen", side_effect=download), redirect_stdout(io.StringIO()):
                cli.main()
            meta = json.loads(out.with_suffix(".metadata.json").read_text())
            self.assertTrue(meta["complete"])
            self.assertEqual(meta["upstream_commit"], cli.DEFAULT_REF)
            self.assertEqual(meta["model_list_sha256"], hashlib.sha256(manifest).hexdigest())
            self.assertEqual(meta["output_sha256"], hashlib.sha256(out.read_bytes()).hexdigest())
            self.assertEqual(meta["summary"]["failed_models"], 1)
            self.assertEqual(meta["summary"]["successful_models"], 2)
            self.assertEqual(meta["summary"]["fraction_without_atps_identifier"], 0.5)
            with out.open() as stream:
                rows = list(csv.DictReader(stream, delimiter="\t"))
            for row in rows:
                self.assertEqual(row["sha256"], hashlib.sha256(files[row["source_path"]]).hexdigest())
            self.assertEqual(len(urls), 4)

    def test_existing_results_and_mutable_references_are_rejected_before_network(self):
        with tempfile.TemporaryDirectory() as directory:
            out = Path(directory) / "screen.tsv"
            out.write_text("preserve me")
            for arguments in (["--output", str(out)], ["--ref", "master"]):
                with self.subTest(arguments=arguments), \
                     patch.object(cli.sys, "argv", ["survey", *arguments]), \
                     patch.object(cli.urllib.request, "urlopen", side_effect=AssertionError("network")), \
                     redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                    cli.main()
            self.assertEqual(out.read_text(), "preserve me")


if __name__ == "__main__":
    unittest.main()
