"""Exercise publication with real local Git repositories and an in-memory GitHub API."""

import copy
from contextlib import redirect_stdout
import io
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from urllib.parse import parse_qs, urlsplit

from repository_promotion import (
    Git, GitHub, PRIVATE, PUBLIC, PROMOTION_PREFIX, SYNC_BRANCH, PromotionError,
    command, main, metadata, promote, promotion_body, synchronize,
)


class LocalGit(Git):
    def __init__(self, root, remotes):
        super().__init__(root)
        self.remotes = remotes
        self.pushes = []

    def url(self, repo):
        return str(self.remotes[repo])

    def push(self, repo, branch, commit, previous=None):
        super().push(repo, branch, commit, previous)
        self.pushes.append((repo, branch, commit))


class FakeGitHub(GitHub):
    def __init__(self, git):
        self.git = git
        self.pulls = {PRIVATE: [], PUBLIC: []}
        self.sync_requests = 0
        self.queue_resumes = 0
        self.auto_merges = []

    def branch_tip(self, repo, branch):
        return self.git.text("--git-dir=" + self.git.url(repo), "rev-parse", branch)

    def api(self, path, method="GET", data=None):
        url = urlsplit(path)
        query = {k: v[0] for k, v in parse_qs(url.query).items()}
        parts = url.path.split("/")
        repo = "/".join(parts[1:3])
        resource = parts[3:]
        if resource[0] == "commits":
            rows = [p for p in self.pulls[repo] if p.get("merge_commit_sha") == resource[1]]
        elif len(resource) == 2:
            row = next(p for p in self.pulls[repo] if p["number"] == int(resource[1]))
            if method == "PATCH":
                row.update(data)
            if row["state"] == "open":
                row["head"]["sha"] = self.branch_tip(repo, row["head"]["ref"])
            return copy.deepcopy(row)
        elif method == "POST":
            number = len(self.pulls[repo]) + 1
            row = dict(
                number=number, title=data["title"], body=data["body"], state="open",
                merged_at=None, merge_commit_sha=None, auto_merge=None,
                head=dict(ref=data["head"], sha=self.branch_tip(repo, data["head"])),
                base=dict(ref=data["base"]), html_url=f"https://github.com/{repo}/pull/{number}",
            )
            self.pulls[repo].append(row)
            return copy.deepcopy(row)
        else:
            rows = self.pulls[repo]
            if query.get("state") not in (None, "all"):
                rows = [p for p in rows if p["state"] == query["state"]]
            if query.get("head"):
                rows = [p for p in rows if p["head"]["ref"] == query["head"].split(":", 1)[1]]
            if query.get("base"):
                rows = [p for p in rows if p["base"]["ref"] == query["base"]]
        start = (int(query.get("page", 1)) - 1) * int(query.get("per_page", 100))
        return copy.deepcopy(rows[start:start + int(query.get("per_page", 100))])

    def auto_merge(self, repo, pull, method):
        if not pull.get("auto_merge"):
            self.auto_merges.append((repo, pull["number"], method))
            next(p for p in self.pulls[repo] if p["number"] == pull["number"])["auto_merge"] = {
                "merge_method": method.upper(),
            }

    def request_sync(self):
        self.sync_requests += 1

    def resume_queue(self):
        self.queue_resumes += 1


class PromotionTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        remotes = {PUBLIC: self.root / "public.git", PRIVATE: self.root / "private.git"}
        for remote in remotes.values():
            command(["git", "init", "--bare", "--initial-branch=main", str(remote)])
        work = self.root / "work"
        command(["git", "init", "--initial-branch=main", str(work)])
        self.git = LocalGit(work, remotes)
        self.git.run("config", "user.name", "Promotion test")
        self.git.run("config", "user.email", "test@example.com")
        # Ignore machine-specific signing settings in real-commit tests.
        self.git.run("config", "commit.gpgsign", "false")
        self.initial = self.commit("value.txt", "zero\n", "initial")
        for repo in (PUBLIC, PRIVATE):
            self.git.push(repo, "main", self.initial)
        self.git.pushes.clear()
        self.github = FakeGitHub(self.git)

    def commit(self, filename, content, message):
        path = Path(self.git.root) / filename
        path.write_text(content)
        self.git.run("add", filename)
        self.git.run("commit", "-m", message)
        return self.git.text("rev-parse", "HEAD")

    def tip(self, repo):
        return self.github.branch_tip(repo, "main")

    def private_change(self, content, filename="value.txt"):
        base = self.tip(PRIVATE)
        self.git.checkout("private-work", base)
        source = self.commit(filename, content, "fix: private change")
        self.git.push(PRIVATE, "main", source, base)
        number = len(self.github.pulls[PRIVATE]) + 1
        row = dict(
            number=number, title=f"fix: private change {number}", body="Private review notes",
            state="closed", merged_at="2026-10-27T00:00:00Z", merge_commit_sha=source,
            head=dict(ref=f"feature/{number}", sha=source), base=dict(ref="main"),
        )
        self.github.pulls[PRIVATE].append(row)
        return row

    def external_change(self, filename="external.txt", content="public\n"):
        base = self.tip(PUBLIC)
        self.git.checkout("public-work", base)
        head = self.commit(filename, content, "fix: independent public change")
        self.git.push(PUBLIC, "main", head, base)
        return head

    def publish(self):
        return promote(self.git, self.github, self.tip(PRIVATE), self.tip(PUBLIC))

    def sync(self):
        return synchronize(self.git, self.github, self.tip(PRIVATE), self.tip(PUBLIC))

    def merge_pull(self, repo, number):
        pull = next(p for p in self.github.pulls[repo] if p["number"] == number)
        base = self.tip(repo)
        self.git.checkout("merge-work", base)
        head = self.github.branch_tip(repo, pull["head"]["ref"])
        if repo == PUBLIC:
            self.git.run("merge", "--squash", head)
            self.git.run("commit", "-m", f"Public squash #{number}")
        else:
            self.git.run("merge", "--no-ff", head, "-m", f"Private return #{number}")
        commit = self.git.text("rev-parse", "HEAD")
        self.git.push(repo, "main", commit, base)
        pull.update(state="closed", merged_at="2026-10-27T01:00:00Z", merge_commit_sha=commit)
        return commit

    def contents(self, commit, filename="value.txt"):
        return self.git.text("show", commit + ":" + filename)

    def test_round_trip_and_reruns_do_not_publish_sync_pr(self):
        source = self.private_change("one\n")["merge_commit_sha"]
        self.publish()
        public_pr = self.github.pulls[PUBLIC][0]
        self.assertNotIn("Private review notes", public_pr["body"])
        self.assertEqual(metadata(public_pr)["source"], source)
        count = len(self.git.pushes)
        self.publish()
        self.assertEqual(len(self.git.pushes), count)
        self.assertEqual(len(self.github.auto_merges), 1)
        public_merge = self.merge_pull(PUBLIC, 1)
        self.sync()
        sync_pr = self.github.pulls[PRIVATE][-1]
        self.assertEqual(sync_pr["head"]["ref"], SYNC_BRANCH)
        self.assertEqual(self.github.auto_merges[-1][2], "merge")
        private_merge = self.merge_pull(PRIVATE, sync_pr["number"])
        self.assertTrue(self.git.ancestor(public_merge, private_merge))
        self.assertEqual(self.contents(private_merge), "one")
        self.assertIn("All merged", self.publish())
        self.sync()
        self.assertEqual(self.github.queue_resumes, 1)
        self.assertEqual(len(self.github.pulls[PUBLIC]), 1)

    def test_later_private_work_survives_return_and_is_published_next(self):
        self.private_change("one\n")
        self.publish()
        later = self.private_change("two\n")["merge_commit_sha"]
        self.publish()
        self.assertEqual(len(self.github.pulls[PUBLIC]), 1)
        self.merge_pull(PUBLIC, 1)
        self.sync()
        returned = self.github.branch_tip(PRIVATE, SYNC_BRANCH)
        self.assertEqual(self.contents(returned), "two")
        self.assertTrue(self.git.ancestor(later, returned))
        self.merge_pull(PRIVATE, self.github.pulls[PRIVATE][-1]["number"])
        self.publish()
        second = self.github.pulls[PUBLIC][-1]
        self.assertEqual(metadata(second)["source"], later)
        self.assertEqual(self.contents(second["head"]["sha"]), "two")

    def test_external_public_change_is_kept_when_pending_public_pr_is_rebased(self):
        self.private_change("one\n")
        self.publish()
        old = self.github.pulls[PUBLIC][0]["head"]["sha"]
        public = self.external_change()
        self.assertIn("Waiting for public main", self.publish())
        self.assertEqual(self.github.sync_requests, 1)
        self.sync()
        self.merge_pull(PRIVATE, self.github.pulls[PRIVATE][-1]["number"])
        self.publish()
        pull = self.github.pulls[PUBLIC][0]
        self.assertNotEqual(pull["head"]["sha"], old)
        self.assertEqual(metadata(pull)["public_base"], public)
        self.assertEqual(self.contents(pull["head"]["sha"], "external.txt"), "public")
        self.merge_pull(PUBLIC, 1)
        self.sync()
        self.assertEqual(self.contents(self.github.branch_tip(PRIVATE, SYNC_BRANCH)), "one")

    def test_successive_returns_reuse_branch_but_create_distinct_prs(self):
        for index in range(3):
            self.external_change(filename=f"external-{index}.txt")
            self.sync()
            pull = self.github.pulls[PRIVATE][-1]
            self.merge_pull(PRIVATE, pull["number"])
        self.assertEqual(len(self.github.pulls[PRIVATE]), 3)
        self.assertTrue(self.git.ancestor(self.tip(PUBLIC), self.tip(PRIVATE)))

    def test_repeated_sync_keeps_the_open_pr_and_its_commit(self):
        self.external_change()
        self.sync()
        first = self.github.branch_tip(PRIVATE, SYNC_BRANCH)
        pushes = len(self.git.pushes)
        self.sync()
        self.assertEqual(self.github.branch_tip(PRIVATE, SYNC_BRANCH), first)
        self.assertEqual(len(self.github.pulls[PRIVATE]), 1)
        self.assertEqual(len(self.git.pushes), pushes)

    def test_delayed_return_integrates_public_followup_without_losing_private_work(self):
        self.private_change("one\n")
        self.publish()
        self.private_change("two\n")
        self.merge_pull(PUBLIC, 1)
        public = self.external_change()
        self.sync()
        head = self.github.branch_tip(PRIVATE, SYNC_BRANCH)
        self.assertEqual(self.contents(head), "two")
        self.assertEqual(self.contents(head, "external.txt"), "public")
        self.assertTrue(self.git.ancestor(public, head))
        self.merge_pull(PRIVATE, self.github.pulls[PRIVATE][-1]["number"])
        self.publish()
        self.merge_pull(PUBLIC, 2)
        self.sync()
        head = self.github.branch_tip(PRIVATE, SYNC_BRANCH)
        self.assertEqual(self.contents(head), "two")
        self.assertEqual(self.contents(head, "external.txt"), "public")

    def test_conflicting_external_change_does_not_push_partial_sync(self):
        self.private_change("private\n")
        self.external_change(filename="value.txt", content="external\n")
        before = len(self.git.pushes)
        with self.assertRaisesRegex(PromotionError, "conflict"):
            self.sync()
        self.assertEqual(len(self.git.pushes), before)
        self.assertEqual(self.contents(self.tip(PRIVATE)), "private")

    def test_manual_changes_on_public_branch_are_not_merged_or_overwritten(self):
        self.private_change("one\n")
        self.publish()
        pull = self.github.pulls[PUBLIC][0]
        old = pull["head"]["sha"]
        self.git.checkout("manual-work", old)
        changed = self.commit("unexpected.txt", "manual\n", "Manual public edit")
        self.git.push(PUBLIC, pull["head"]["ref"], changed, old)
        pull["head"]["sha"] = changed
        before = len(self.git.pushes)
        with self.assertRaisesRegex(PromotionError, "manual changes"):
            self.publish()
        self.assertEqual(len(self.git.pushes), before)
        self.assertEqual(self.github.branch_tip(PUBLIC, pull["head"]["ref"]), changed)

    def test_tampered_public_merge_is_never_marked_as_already_present(self):
        self.private_change("one\n")
        self.publish()
        pull = self.github.pulls[PUBLIC][0]
        old = pull["head"]["sha"]
        self.git.checkout("tamper-work", old)
        changed = self.commit("unexpected.txt", "extra\n", "Extra change")
        self.git.push(PUBLIC, pull["head"]["ref"], changed, old)
        pull["head"]["sha"] = changed
        self.merge_pull(PUBLIC, 1)
        before = len(self.git.pushes)
        with self.assertRaisesRegex(PromotionError, "beyond the approved"):
            self.sync()
        self.assertEqual(len(self.git.pushes), before)

    def test_closed_public_pr_blocks_later_publications(self):
        self.private_change("one\n")
        self.publish()
        self.github.pulls[PUBLIC][0]["state"] = "closed"
        self.private_change("two\n")
        with self.assertRaisesRegex(PromotionError, "closed without merging"):
            self.publish()
        self.assertEqual(len(self.github.pulls[PUBLIC]), 1)

    def test_run_interrupted_after_push_reuses_the_existing_branch(self):
        self.private_change("one\n")
        self.publish()
        head = self.github.pulls[PUBLIC][0]["head"]["sha"]
        self.github.pulls[PUBLIC].clear()
        count = len(self.git.pushes)
        self.publish()
        self.assertEqual(len(self.git.pushes), count)
        self.assertEqual(self.github.pulls[PUBLIC][0]["head"]["sha"], head)

    def test_unmerged_private_branch_is_not_exposed(self):
        base = self.tip(PRIVATE)
        self.git.checkout("unmerged", base)
        private = self.commit("draft.txt", "internal\n", "Unmerged draft")
        self.git.push(PRIVATE, "feature/draft", private)
        self.assertIn("All merged", self.publish())
        self.assertEqual(self.github.pulls[PUBLIC], [])
        self.assertEqual(self.tip(PUBLIC), self.initial)

    def test_non_squash_private_merge_requires_review(self):
        base = self.tip(PRIVATE)
        self.git.checkout("feature", base)
        feature = self.commit("value.txt", "one\n", "Feature")
        self.git.checkout("merge-private", base)
        self.git.run("merge", "--no-ff", feature, "-m", "Private merge")
        merge = self.git.text("rev-parse", "HEAD")
        self.git.push(PRIVATE, "main", merge, base)
        self.github.pulls[PRIVATE].append(dict(
            number=1, title="feature", state="closed", merged_at="now",
            merge_commit_sha=merge, head=dict(ref="feature"), base=dict(ref="main"),
        ))
        with self.assertRaisesRegex(PromotionError, "must use squash"):
            self.publish()

    def test_lease_preserves_a_concurrent_remote_update(self):
        base = self.tip(PRIVATE)
        later = self.private_change("one\n")["merge_commit_sha"]
        with self.assertRaises(PromotionError):
            self.git.push(PRIVATE, "main", base, base)
        self.assertEqual(self.tip(PRIVATE), later)

    def test_queue_order_follows_history_instead_of_pr_creation_number(self):
        first = self.private_change("one\n")
        first["number"] = 9
        second = self.private_change("two\n")
        second["number"] = 1
        self.github.pulls[PRIVATE].reverse()
        self.publish()
        self.assertEqual(self.github.pulls[PUBLIC][0]["head"]["ref"], PROMOTION_PREFIX + "9")


class MetadataTests(unittest.TestCase):
    def test_rejects_missing_invalid_or_option_like_object_ids(self):
        for body in ("", "<!-- mmcfilters-promotion: {invalid} -->",
                     promotion_body(dict(private_pr=1, source="--help",
                                         public_base="a" * 40, tree="b" * 40))):
            with self.subTest(body=body), self.assertRaises(PromotionError):
                metadata(dict(body=body))

    def test_lists_more_than_one_page(self):
        class PagedGitHub(GitHub):
            def api(self, path, method="GET", data=None):
                page = int(parse_qs(urlsplit(path).query)["page"][0])
                return list(range((page - 1) * 100, min(page * 100, 207)))
        self.assertEqual(PagedGitHub().pages("repos/test/example/pulls"), list(range(207)))


class PublicControllerTests(unittest.TestCase):
    def test_public_merge_waits_for_checks(self):
        pull = dict(number=1, head=dict(ref=PROMOTION_PREFIX + "4", sha="a" * 40))
        with patch("repository_promotion.command") as run:
            GitHub().auto_merge(PUBLIC, pull, "squash")
            self.assertEqual(run.call_args.args[0], [
                "gh", "pr", "merge", "1", "--repo", PUBLIC, "--squash",
                "--match-head-commit", "a" * 40, "--auto",
            ])
            run.reset_mock()
            pull["auto_merge"] = {"merge_method": "SQUASH"}
            GitHub().auto_merge(PUBLIC, pull, "squash")
            run.assert_not_called()

    def test_private_sync_merges_without_enabling_auto_merge(self):
        pull = dict(number=5, head=dict(ref=SYNC_BRANCH, sha="b" * 40))
        with patch("repository_promotion.command") as run:
            GitHub().auto_merge(PRIVATE, pull, "merge")
            self.assertEqual(run.call_args.args[0], [
                "gh", "pr", "merge", "5", "--repo", PRIVATE, "--merge",
                "--match-head-commit", "b" * 40,
            ])

    def test_direct_merge_is_restricted_to_private_sync(self):
        for repo, branch, method in ((PRIVATE, "feature", "merge"),
                                     (PRIVATE, SYNC_BRANCH, "squash"),
                                     ("other/repository", SYNC_BRANCH, "merge")):
            with self.subTest(repo=repo, branch=branch, method=method), \
                    patch("repository_promotion.command") as run:
                pull = dict(number=1, head=dict(ref=branch, sha="c" * 40))
                with self.assertRaisesRegex(PromotionError, "Only the reserved"):
                    GitHub().auto_merge(repo, pull, method)
                run.assert_not_called()

    def test_all_dispatches_target_public_actions(self):
        class RecordingGitHub(GitHub):
            def __init__(self):
                self.calls = []

            def api(self, path, method="GET", data=None):
                self.calls.append((path, method, data))

        github = RecordingGitHub()
        github.resume_queue()
        github.request_sync()
        self.assertEqual(github.calls, [
            (f"repos/{PUBLIC}/dispatches", "POST", {"event_type": "resume-promotion"}),
            (f"repos/{PUBLIC}/dispatches", "POST", {"event_type": "sync-private"}),
        ])

    def test_private_execution_is_rejected_before_accessing_git(self):
        for mode in ("promote", "sync"):
            with self.subTest(mode=mode), \
                    patch.dict(os.environ, {"GITHUB_REPOSITORY": PRIVATE}), \
                    patch("sys.argv", ["repository_promotion.py", mode]), \
                    patch("repository_promotion.Git") as git:
                with self.assertRaisesRegex(PromotionError, "must run in the public"):
                    main()
                git.assert_not_called()

    def test_both_operations_execute_from_public(self):
        for mode, operation in (("promote", "promote"), ("sync", "synchronize")):
            with self.subTest(mode=mode), \
                    patch.dict(os.environ, {"GITHUB_REPOSITORY": PUBLIC,
                                            "GITHUB_STEP_SUMMARY": ""}), \
                    patch("sys.argv", ["repository_promotion.py", mode]), \
                    patch("repository_promotion.Git") as git_type, \
                    patch("repository_promotion.GitHub") as github_type, \
                    patch("repository_promotion." + operation, return_value="finished") as run, \
                    redirect_stdout(io.StringIO()):
                git = git_type.return_value
                git.text.return_value = ""
                git.fetch_main.side_effect = ["private-tip", "public-tip"]
                main()
                run.assert_called_once_with(git, github_type.return_value,
                                            "private-tip", "public-tip")


if __name__ == "__main__":
    unittest.main()
