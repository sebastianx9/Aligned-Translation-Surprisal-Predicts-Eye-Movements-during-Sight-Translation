import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]


class AnalysisSubmissionTests(unittest.TestCase):
    def read(self, relative_path):
        return (REPO / relative_path).read_text(encoding="utf-8")

    def test_rq2_interaction_cv_has_exactly_two_models(self):
        script = self.read("RQ2/rq2_interaction_kfold.R")
        self.assertIn('common = f("c_mono + cond:c_mono + c_nmt")', script)
        self.assertIn("stage_specific =", script)
        self.assertIn("ptw$stage_specific - ptw$common", script)
        self.assertNotIn("M1 = f(", script)
        self.assertNotIn("N1 = f(", script)

    def test_rq1_joint_script_reports_all_three_comparisons(self):
        script = self.read("RQ1/rq1_joint_surprisal_kfold.R")
        for contrast in (
            '"M_nmt - M_mono"',
            '"M_both - M_mono"',
            '"M_both - M_nmt"',
        ):
            self.assertIn(contrast, script)
        for cache_name in ("c_mono", "c_nmt", "c_mono_nmt"):
            self.assertIn(f'fit_kfold("{cache_name}"', script)

    def test_rq1_uses_all_six_lim_source_attention_features(self):
        script = self.read("RQ1/rq1_kfold_elpd.R")
        extractor = self.read(
            "data-extraction/extract_attention_features_norm.py"
        )
        self.assertIn("f_self=attn_self", script)
        self.assertIn('"c_fself"', script)
        self.assertIn('"attn_self"', extractor)

    def test_submission_graph_uses_refactored_jobs(self):
        submit = self.read("hpc/submit_core_jobs.sh")
        runner = self.read("hpc/csf3_analysis.sbatch")
        self.assertIn("rq1_joint_predictive", submit)
        self.assertIn("rq2_interaction_cv", submit)
        self.assertNotIn("submit rq2_beyond", submit)
        self.assertNotIn("submit rq2_cv", submit)
        self.assertIn("RQ1/rq1_joint_surprisal_kfold.R", runner)
        self.assertIn("RQ2/rq2_interaction_kfold.R", runner)
        self.assertIn("EXPECTED_GIT_COMMIT", submit)
        self.assertIn("EXPECTED_MANIFEST_SHA256", submit)
        self.assertIn("rq2_joint_diagnostics", submit)
        self.assertIn("rq3_rrt_diagnostics", submit)
        self.assertIn("rq3coef_v3_total_RRT.rds", submit)

    def test_csf_preflight_verifies_the_input_manifest(self):
        check = self.read("hpc/csf3_check.sbatch")
        analysis = self.read("hpc/csf3_analysis.sbatch")
        self.assertIn("sha256sum -c MANIFEST.sha256", check)
        self.assertIn("EXPECTED_MANIFEST_SHA256", check)
        self.assertIn("EXPECTED_MANIFEST_SHA256", analysis)
        self.assertIn("#SBATCH --ntasks=1", check)
        self.assertIn("#SBATCH --cpus-per-task=4", analysis)

    def test_rq3_primary_contrast_locates_total_cnmt_association(self):
        script = self.read("RQ3/rq3_kfold_elpd.R")
        self.assertIn(
            'compare("c_nmt_total", "c_nmt", "base", "focal")',
            script,
        )
        self.assertNotIn('fitkf("c_mono"', script)
        self.assertNotIn('fitkf("c_mono_nmt"', script)
        self.assertNotIn("c_nmt_unique", script)
        self.assertIn("(1 + c_nmt | participant)", script)
        self.assertIn("rq3_coefficient_results.csv", script)
        self.assertIn('Go_past=run_outcome(', script)

    def test_rq3_figure_uses_matching_total_tfd_contrast(self):
        script = self.read("RQ3/rq3_gd_rrt.R")
        self.assertIn('"rq1_kfold_elpd.rds"', script)
        self.assertIn('predictor == "c_nmt"', script)
        self.assertIn('contrast == "c_nmt_total"', script)


if __name__ == "__main__":
    unittest.main()
