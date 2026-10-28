#!/usr/bin/env python3
"""Publish approved private PRs and reconcile the two repository histories."""

import argparse
import json
import os
from pathlib import Path
import re
import subprocess
from urllib.parse import urlencode


PUBLIC = "wonderalexandre/mmcfilters"
PRIVATE = "wonderalexandre/mmcfilters-dev"
SYNC_BRANCH = "sync/public-main"
PROMOTION_PREFIX = "promotion/private-pr-"
MARKER = re.compile(r"<!-- mmcfilters-promotion: (\{[^\n]+\}) -->")


class PromotionError(RuntimeError):
    pass


def command(args, *, cwd=None, data=None, check=True):
    result = subprocess.run(args, cwd=cwd, input=data, text=True, capture_output=True)
    if check and result.returncode:
        raise PromotionError(result.stderr.strip() or result.stdout.strip())
    return result


class GitHub:
    def api(self, path, method="GET", data=None):
        args = ["gh", "api", "--method", method, path]
        if data is not None:
            args += ["--input", "-"]
        result = command(args, data=json.dumps(data) if data is not None else None)
        return json.loads(result.stdout) if result.stdout.strip() else None

    def pages(self, path, **query):
        rows = []
        for page in range(1, 10001):
            batch = self.api(path + "?" + urlencode(dict(query, per_page=100, page=page)))
            rows.extend(batch)
            if len(batch) < 100:
                return rows
        raise PromotionError("Pagination limit exceeded.")

    def pull(self, repo, number):
        return self.api(f"repos/{repo}/pulls/{number}")

    def find_pull(self, repo, branch, state="all"):
        rows = self.pages(
            f"repos/{repo}/pulls", state=state, base="main",
            head=f"{repo.split('/')[0]}:{branch}",
        )
        if len(rows) > 1:
            raise PromotionError(f"Multiple PRs use {repo}:{branch}; resolve them first.")
        return rows[0] if rows else None

    def auto_merge(self, repo, pull, method):
        if not pull.get("auto_merge"):
            command([
                "gh", "pr", "merge", str(pull["number"]), "--repo", repo,
                "--auto", "--" + method, "--match-head-commit", pull["head"]["sha"],
            ])

    def request_sync(self):
        self.api(f"repos/{PUBLIC}/dispatches", "POST", {"event_type": "sync-private"})

    def resume_queue(self):
        self.api(f"repos/{PUBLIC}/dispatches", "POST", {"event_type": "resume-promotion"})


class Git:
    def __init__(self, root="."):
        self.root = root

    def run(self, *args, check=True):
        return command(["git", *args], cwd=self.root, check=check)

    def text(self, *args):
        return self.run(*args).stdout.strip()

    def url(self, repo):
        return f"https://github.com/{repo}.git"

    def fetch_main(self, repo, name):
        self.run("fetch", "--no-tags", self.url(repo),
                 f"+refs/heads/main:refs/remotes/{name}/main")
        return self.text("rev-parse", f"refs/remotes/{name}/main")

    def ancestor(self, older, newer):
        result = self.run("merge-base", "--is-ancestor", older, newer, check=False)
        if result.returncode not in (0, 1):
            raise PromotionError(result.stderr.strip())
        return result.returncode == 0

    def tree(self, commit):
        return self.text("rev-parse", commit + "^{tree}")

    def parents(self, commit):
        return self.text("show", "-s", "--format=%P", commit).split()

    def remote_branch(self, repo, branch):
        result = self.run("ls-remote", "--exit-code", "--heads",
                          self.url(repo), "refs/heads/" + branch,
                          check=False)
        if result.returncode == 2:
            return None
        if result.returncode:
            raise PromotionError(result.stderr.strip())
        sha = result.stdout.split()[0]
        self.run("fetch", "--no-tags", self.url(repo),
                 "refs/heads/" + branch)
        return sha

    def checkout(self, branch, base):
        self.run("checkout", "-B", branch, base)

    def replay(self, source, base):
        self.checkout("automation-replay", base)
        result = self.run("cherry-pick", "-x", source, check=False)
        if result.returncode:
            self.run("cherry-pick", "--abort", check=False)
            raise PromotionError(
                "The approved change cannot be applied to public main. "
                "Resolve the conflict or an already-applied change before retrying.\n"
                + result.stderr.strip()
            )
        return self.text("rev-parse", "HEAD")

    def push(self, repo, branch, commit, previous=None):
        args = ["-c", "push.followTags=false", "push"]
        # A lease also prevents overwriting a branch created after the initial read.
        args += [f"--force-with-lease=refs/heads/{branch}:{previous or ''}"]
        args += [self.url(repo), f"{commit}:refs/heads/{branch}"]
        self.run(*args)


def metadata(pull):
    match = MARKER.search(pull.get("body") or "")
    if not match:
        raise PromotionError("The public PR has no promotion metadata.")
    try:
        result = json.loads(match[1])
    except json.JSONDecodeError as error:
        raise PromotionError("Invalid promotion metadata.") from error
    if set(result) != {"private_pr", "source", "public_base", "tree"}:
        raise PromotionError("Invalid promotion metadata.")
    if type(result["private_pr"]) is not int or result["private_pr"] <= 0:
        raise PromotionError("Invalid private PR number.")
    for key in ("source", "public_base", "tree"):
        if not re.fullmatch(r"[0-9a-f]{40}", result[key]):
            raise PromotionError("Invalid promotion object ID.")
    return result


def promotion_body(record):
    return (
        "Publish changes approved in the private development repository.\n\n"
        "The public checks must pass before this PR can merge.\n\n"
        "<!-- mmcfilters-promotion: " + json.dumps(record, sort_keys=True) + " -->\n"
    )


def merged(pull):
    return bool(pull.get("merged_at"))


def promote(git, github, private_tip, public_tip):
    pulls = github.pages(f"repos/{PRIVATE}/pulls", state="closed", base="main")
    pulls = [p for p in pulls if merged(p) and p["head"]["ref"] != SYNC_BRANCH]
    order = git.text("rev-list", "--first-parent", "--reverse", private_tip).splitlines()
    positions = {sha: i for i, sha in enumerate(order)}
    if any(p["merge_commit_sha"] not in positions for p in pulls):
        raise PromotionError("A merged private PR is outside the main history.")
    pulls.sort(key=lambda p: positions[p["merge_commit_sha"]])

    for source_pr in pulls:
        number = source_pr["number"]
        source = source_pr["merge_commit_sha"]
        branch = PROMOTION_PREFIX + str(number)
        public_pr = github.find_pull(PUBLIC, branch)
        old_record = None
        if public_pr:
            old_record = metadata(public_pr)
            if old_record["source"] != source or old_record["private_pr"] != number:
                raise PromotionError("The public PR refers to another private change.")
            if merged(public_pr):
                continue
            if public_pr["state"] != "open":
                raise PromotionError(f"Public PR #{public_pr['number']} was closed without merging.")
            previous = git.remote_branch(PUBLIC, branch)
            if (previous != public_pr["head"]["sha"]
                    or git.parents(previous) != [old_record["public_base"]]
                    or git.tree(previous) != old_record["tree"]):
                raise PromotionError("The publication branch has manual changes; preserve and review them.")

        if not git.ancestor(public_tip, private_tip):
            github.request_sync()
            return "Waiting for public main to be reconciled into private main."
        if len(git.parents(source)) != 1:
            raise PromotionError(f"Private PR #{number} must use squash merge.")
        if public_pr and old_record["public_base"] == public_tip:
            github.auto_merge(PUBLIC, public_pr, "squash")
            return f"Waiting for public PR #{public_pr['number']} to merge."

        candidate = git.replay(source, public_tip)
        record = dict(private_pr=number, source=source,
                      public_base=public_tip, tree=git.tree(candidate))
        previous = git.remote_branch(PUBLIC, branch)
        if previous:
            expected = old_record or record
            if (git.parents(previous) != [expected["public_base"]]
                    or git.tree(previous) != expected["tree"]):
                raise PromotionError("The publication branch has manual changes; preserve and review them.")
            if not public_pr:
                candidate = previous  # Recover a run interrupted between push and PR creation.
        if candidate != previous:
            git.push(PUBLIC, branch, candidate, previous)
        if public_pr:
            body = MARKER.sub(lambda _: MARKER.search(promotion_body(record))[0],
                              public_pr["body"])
            github.api(f"repos/{PUBLIC}/pulls/{public_pr['number']}", "PATCH", {"body": body})
            public_pr = github.pull(PUBLIC, public_pr["number"])
        else:
            public_pr = github.api(f"repos/{PUBLIC}/pulls", "POST", {
                "head": branch, "base": "main", "title": source_pr["title"],
                "body": promotion_body(record),
            })
        github.auto_merge(PUBLIC, public_pr, "squash")
        return f"Publication pending: {public_pr['html_url']}"
    return "All merged private PRs have been published."


def verify_published_change(git, github, commit, private_tip, pull):
    """An equivalent squash is already present in private history, even with later work."""
    record = metadata(pull)
    source_pr = github.pull(PRIVATE, record["private_pr"])
    source = record["source"]
    if (not merged(source_pr) or source_pr["base"]["ref"] != "main"
            or source_pr["head"]["ref"] == SYNC_BRANCH
            or source_pr["merge_commit_sha"] != source
            or pull["head"]["ref"] != PROMOTION_PREFIX + str(record["private_pr"])
            or git.parents(commit) != [record["public_base"]]
            or len(git.parents(source)) != 1
            or not git.ancestor(source, private_tip)
            or not git.ancestor(record["public_base"], private_tip)):
        raise PromotionError("The public squash does not match approved private history.")
    expected = git.text("merge-tree", "--write-tree",
                        "--merge-base=" + git.parents(source)[0],
                        record["public_base"], source).splitlines()[0]
    if expected != record["tree"] or git.tree(commit) != expected:
        raise PromotionError("The public merge contains changes beyond the approved private PR.")


def synchronize(git, github, private_tip, public_tip):
    if git.ancestor(public_tip, private_tip):
        github.resume_queue()
        return "Private main already includes public main."
    git.checkout("automation-sync", private_tip)
    commits = git.text("rev-list", "--first-parent", "--reverse",
                       private_tip + ".." + public_tip).splitlines()
    for commit in commits:
        associated = github.pages(f"repos/{PUBLIC}/commits/{commit}/pulls")
        promotions = [p for p in associated if merged(p)
                      and p["merge_commit_sha"] == commit and p["base"]["ref"] == "main"
                      and p["head"]["ref"].startswith(PROMOTION_PREFIX)]
        if len(promotions) > 1:
            raise PromotionError("Ambiguous public promotion history.")
        args = ["merge", "--no-ff", "-m", "Sync public main: " + commit, commit]
        if promotions:
            verify_published_change(git, github, commit, git.text("rev-parse", "HEAD"), promotions[0])
            # The content was verified above. Record ancestry without reverting later private work.
            args[1:1] = ["--strategy=ours"]
        result = git.run(*args, check=False)
        if result.returncode:
            git.run("merge", "--abort", check=False)
            raise PromotionError("Public changes conflict with private work; resolve them manually.\n"
                                 + result.stderr.strip())
    candidate = git.text("rev-parse", "HEAD")
    previous = git.remote_branch(PRIVATE, SYNC_BRANCH)
    if (previous and git.ancestor(private_tip, previous) and git.ancestor(public_tip, previous)
            and git.tree(previous) == git.tree(candidate)):
        candidate = previous
    if candidate != previous:
        git.push(PRIVATE, SYNC_BRANCH, candidate, previous)
    sync_pr = github.find_pull(PRIVATE, SYNC_BRANCH, state="open")
    if not sync_pr:
        sync_pr = github.api(f"repos/{PRIVATE}/pulls", "POST", {
            "head": SYNC_BRANCH, "base": "main", "title": "chore: sync public main",
            "body": "Reconcile public main into private main.\n\n"
                    "Use a merge commit to retain the public ancestry. "
                    "This reserved branch is excluded from publication.\n",
        })
    else:
        sync_pr = github.pull(PRIVATE, sync_pr["number"])
    github.auto_merge(PRIVATE, sync_pr, "merge")
    return f"Synchronization pending: {sync_pr['html_url']}"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("promote", "sync"))
    args = parser.parse_args()
    if os.environ.get("GITHUB_REPOSITORY") != PUBLIC:
        raise PromotionError("Repository automation must run in the public repository.")
    git, github = Git(), GitHub()
    if git.text("status", "--porcelain", "--untracked-files=no"):
        raise PromotionError("The automation checkout must be clean.")
    private_tip = git.fetch_main(PRIVATE, "private")
    public_tip = git.fetch_main(PUBLIC, "public")
    message = (promote if args.mode == "promote" else synchronize)(
        git, github, private_tip, public_tip)
    print(message)
    if os.environ.get("GITHUB_STEP_SUMMARY"):
        with Path(os.environ["GITHUB_STEP_SUMMARY"]).open("a", encoding="utf-8") as stream:
            stream.write(message + "\n")


if __name__ == "__main__":
    try:
        main()
    except PromotionError as error:
        raise SystemExit(str(error))
