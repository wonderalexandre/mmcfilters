"""Check the workflow exception and the required aggregate CI result."""

import copy
import json
import os
from pathlib import Path
import subprocess
import unittest

import yaml


WORKFLOWS = Path(__file__).resolve().parents[1] / "workflows"


def lookup(context, path):
    value = context
    for key in path.split("."):
        if not isinstance(value, dict):
            return None
        value = value.get(key)
    return value


def should_validate(expression, context):
    """Evaluate the policy's limited OR/inequality syntax, including toJSON."""
    def operand(text):
        text = text.strip()
        if text.startswith("'") and text.endswith("'"):
            return text[1:-1]
        if text.startswith("toJSON(") and text.endswith(")"):
            return json.dumps(lookup(context, text[7:-1]))
        if not text.startswith("github."):
            raise AssertionError(f"Unsupported policy operand: {text}")
        return lookup(context, text)

    for clause in expression.split("||"):
        left, right = map(operand, clause.strip().split(" != "))
        # Actions compares strings without regard to case.
        if isinstance(left, str) and isinstance(right, str):
            left, right = left.casefold(), right.casefold()
        if left != right:
            return True
    return False


class CiPolicyTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.ci = yaml.safe_load((WORKFLOWS / "ci.yml").read_text())["jobs"]
        cls.condition = cls.ci["repository-promotion"]["if"]
        cls.sync = {
            "github": {
                "event_name": "pull_request",
                "repository": "wonderalexandre/mmcfilters-dev",
                "event": {
                    "pull_request": {
                        "base": {"ref": "main"},
                        "head": {
                            "ref": "sync/public-main",
                            "repo": {"full_name": "wonderalexandre/mmcfilters-dev"},
                        },
                        "user": {"login": "mmcfilters-repository-promotion[bot]"},
                        "changed_files": 0,
                    }
                },
            }
        }

    def test_only_trusted_content_free_sync_skips_validation(self):
        self.assertFalse(should_validate(self.condition, self.sync))
        cases = [
            ("github.repository", "wonderalexandre/mmcfilters"),
            ("github.event.pull_request.base.ref", "other-base"),
            ("github.event.pull_request.head.ref", "feature"),
            ("github.event.pull_request.head.repo.full_name", "other/fork"),
            ("github.event.pull_request.user.login", "wonderalexandre"),
            ("github.event.pull_request.changed_files", 1),
            ("github.event.pull_request.changed_files", None),
            ("github.event.pull_request.changed_files", "0"),
            ("github.event.pull_request.changed_files", False),
        ]
        for path, value in cases:
            with self.subTest(path=path, value=value):
                context = copy.deepcopy(self.sync)
                parent, key = path.rsplit(".", 1)
                lookup(context, parent)[key] = value
                self.assertTrue(should_validate(self.condition, context))

    def test_missing_changed_files_runs_validation(self):
        context = copy.deepcopy(self.sync)
        del context["github"]["event"]["pull_request"]["changed_files"]
        self.assertTrue(should_validate(self.condition, context))

    def test_non_pr_events_run_validation_without_pr_payload(self):
        for event in ("push", "schedule", "workflow_dispatch"):
            with self.subTest(event=event):
                self.assertTrue(should_validate(
                    self.condition, {"github": {"event_name": event}}
                ))

    def test_all_validation_jobs_share_policy_and_feed_result(self):
        gate = self.ci["ci-result"]
        jobs = set(self.ci) - {"ci-result"}
        self.assertEqual(set(gate["needs"]), jobs)
        self.assertEqual(
            " ".join(gate["if"].split()), f"always() && ( {self.condition} )"
        )
        for job in jobs:
            with self.subTest(job=job):
                self.assertEqual(self.ci[job]["if"], self.condition)
                self.assertFalse(self.ci[job].get("continue-on-error", False))
        for filename, job in (("docs.yml", "doxygen"),
                              ("package.yml", "build-distributions")):
            workflow = yaml.safe_load((WORKFLOWS / filename).read_text())
            self.assertEqual(workflow["jobs"][job]["if"], self.condition)

    def test_aggregate_result_requires_every_job_to_succeed(self):
        gate = self.ci["ci-result"]
        script = gate["steps"][0]["run"]
        success = {job: {"result": "success"} for job in gate["needs"]}
        cases = [(success, True), ({}, False)]
        for result in ("failure", "cancelled", "skipped", "unknown", None):
            for job in gate["needs"]:
                jobs = copy.deepcopy(success)
                jobs[job]["result"] = result
                cases.append((jobs, False))
        malformed = copy.deepcopy(success)
        malformed[gate["needs"][0]] = {}
        cases.append((malformed, False))
        for jobs, expected_success in cases:
            with self.subTest(jobs=jobs):
                result = subprocess.run(
                    ["bash", "-e", "-o", "pipefail", "-c", script],
                    env={**os.environ, "CI_NEEDS": json.dumps(jobs)},
                    capture_output=True, text=True,
                )
                self.assertEqual(result.returncode == 0, expected_success,
                                 result.stdout + result.stderr)


if __name__ == "__main__":
    unittest.main()
