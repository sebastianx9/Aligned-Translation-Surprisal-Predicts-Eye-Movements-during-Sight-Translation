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

    def test_submission_graph_uses_refactored_jobs(self):
        submit = self.read("hpc/submit_core_jobs.sh")
        runner = self.read("hpc/csf3_analysis.sbatch")
        self.assertIn("rq1_joint_predictive", submit)
        self.assertIn("rq2_interaction_cv", submit)
        self.assertNotIn("submit rq2_beyond", submit)
        self.assertNotIn("submit rq2_cv", submit)
        self.assertIn("RQ1/rq1_joint_surprisal_kfold.R", runner)
        self.assertIn("RQ2/rq2_interaction_kfold.R", runner)


if __name__ == "__main__":
    unittest.main()
