# Reviewing and Contributing to Someone Else's Pull Request

Pull requests (PRs) to `oomph-lib`'s main branch `oomph-lib/oomph-lib` have to be reviewed by at least one maintainer. Usually there are issues that have to be dealt with before the PR can be merged. This is generally done by annotating code on GitHub and/or communicating via email. However, sometimes it is easier for the maintainer to get involved, just because it's often quicker to "just do it yourself". A similar problem arises if the maintainer wants to "try out" the newly developed functionality before approving it. This can obviously not be done on GitHub but requires the maintainer to download (and possibly modify) the PR in its entirety. This document (the original draft of which was created by Copilot) explains how to do this. 

Here's the scenario:

- A contributor (JoeContributor, say) has submitted a PR to `oomph-lib/oomph-lib`.
- You, the maintainer (JosephineMaintainer, say), want to test the code locally.
- You may want to make changes yourself.
- If so, you want to pass the changes back to the contributor (for approval or otherwise) via a "PR on the PR" rather than pushing them directly to the contributor's branch (which they may not want to share with you directly).
- The contributor may continue working on their branch while you are reviewing it.

To give specific (almost-)copy-and-paste instructions, the example below uses the following assumptions:

- The original repository is `oomph-lib/oomph-lib`. This represents the GitHub repository at https://github.com/oomph-lib/oomph-lib.
- Original PR: `#231`. 

   **[Replace '231' with whatever number GitHub has assigned to the PR that the contributor has submitted to the `oomph-lib` repo.]**
- The contributor repository is `JoeContributor/oomph-lib`. This represents the GitHub repository (a fork of `oomph-lib`'s main GitHub repo) at https://github.com/JoeContributor/oomph-lib. 

   **[Replace 'JoeContributor' with the GitHub username of the actual contributor.]**
- The contributor wishes to merge their branch `bug_fix` into the `oomph-lib` main branch. 

   **[Replace `bug_fix` with the name of the actual branch.]**
- Your own fork of  `oomph-lib`'s GitHub repo is `JosephineMaintainer/oomph-lib`. This represents the GitHub repo at https://github.com/JosephineMaintainer/oomph-lib. 

   **[Replace JosephineMaintainer with your GitHub username.]**

---

## One-time Setup

Check out a clean version of the `oomph-lib` repository
```bash
git clone git@github.com:oomph-lib/oomph-lib.git
```

Add your fork as an additional remote in your local clone:

```bash
cd oomph-lib
git remote add myfork git@github.com:JosephineMaintainer/oomph-lib.git
```

Check your configuration:

```bash
git remote -v
```

Expected output:

```text
myfork	git@github.com:JosephineMaintainer/oomph-lib.git (fetch)
myfork	git@github.com:JosephineMaintainer/oomph-lib.git (push)
origin	git@github.com:oomph-lib/oomph-lib.git (fetch)
origin	git@github.com:oomph-lib/oomph-lib.git (push)
```

This setup only needs to be done once.

---

## Step 1: Check Out the Original PR Locally

Start from the `oomph-lib` main branch:

```bash
git switch main
git pull origin main
```

Fetch PR #231:

```bash
git fetch origin pull/231/head:pr231
```

Switch to the PR branch:

```bash
git switch pr231
```

You are now looking at exactly the code contained in PR #231.

---

## Step 2: Create Your Own Working Branch

Do **not** edit `pr231` directly.

Create a branch for your own work:

```bash
git switch -c pr231-josephine
```

Current branch structure:

```text
main
  \
   pr231
      \
       pr231-josephine
```

At this point `pr231` and `pr231-josephine` point to the same commit.

---

## Step 3: Make and Commit Your Changes

Edit the files as required.

Check status:

```bash
git status
```

Commit your changes:

```bash
git add <files>
git commit -m "Fix xyz in PR #231"
```

The branch structure is now:

```text
main
  \
   pr231
      \
       pr231-josephine
            \
             your commits
```

---

## Step 4: Push Your Branch to Your Fork

Push the branch to GitHub:

```bash
git push -u myfork pr231-josephine
```

This creates the branch

```text
JosephineMaintainer/oomph-lib
    └── pr231-josephine
```

on your fork.

Subsequent pushes require only:

```bash
git push
```

---

## Step 5: Create a Pull Request Against the Contributor's Branch

Create a new PR on GitHub.

Use:

```text
Base repository: JoeContributor/oomph-lib
Base branch:     bug_fix

Head repository: JosephineMaintainer/oomph-lib
Compare branch:  pr231-josephine
```

In words:

> Please merge my extra commits into the branch that underlies PR #231.

This creates a "PR on the PR".

---

## Step 6: Contributor Adds More Commits

Suppose the contributor pushes additional commits to their `bug_fix` branch.

The original PR #231 has now changed.

To update your local copy of the original PR:

```bash
git switch pr231

git fetch origin pull/231/head
git reset --hard FETCH_HEAD
```

Your local `pr231` branch now matches the latest state of PR #231.

---

## Step 7: Pull the Contributor's Updates into Your Branch

Switch back to your branch:

```bash
git switch pr231-josephine
```

Merge the updated PR branch:

```bash
git merge pr231
```

This produces:

```text
main
  \
   pr231 (updated)
      \
       pr231-josephine
            \
             your commits
              +
             contributor updates
```

Resolve any merge conflicts if necessary.

---

## Step 8: Update Your PR

Push your updated branch:

```bash
git push
```

GitHub automatically updates your PR against the contributor's branch.

No further action is required.

---

## Useful Commands

Show the branch structure:

```bash
git log --graph --oneline --decorate --all
```

Compare your changes against the original PR:

```bash
git diff pr231
```

Update local copy of PR #231:

```bash
git switch pr231
git fetch origin pull/231/head
git reset --hard FETCH_HEAD
```

Update your branch with contributor changes:

```bash
git switch pr231-josephine
git merge pr231
```

Push updates to GitHub:

```bash
git push
```

---

## Recommended Mental Model

Treat the original PR branch as read-only:

```text
pr231
```

Treat your branch as:

```text
pr231
   +
your commits
```

Whenever the contributor adds more commits:

```text
update pr231
       ↓
merge pr231 into pr231-josephine
       ↓
git push
```

Your PR automatically tracks both:

- the contributor's ongoing work; and
- your additional commits.

This allows maintainer and contributor to collaborate without either party pushing directly to the other's branch. When you're happy with the version you've jointly created in the contributor's branch, merge the original PR and thus pull it into the `oomph-lib` main repository. Done!