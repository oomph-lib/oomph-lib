# oomph-lib's GitHub workflow
**A guide to working with and contributing to the (now-GitHub-hosted) `oomph-lib` repository.**

_Notation:_ We prefix any command line input with "`>>>`" and generally show the resulting output from Git underneath. Lengthy output is sometimes truncated and omitted parts are then indicated by "`[...]`". Comments for specific commands are prefixed with "`#`".
## Basic setup (only to be done once)

We assume that you have created a GitHub account, and for the purpose of this document assume that your GitHub home page is https://github.com/JoeCoolDummy. Unless your name is Joe Cool Dummy, you'll have to change the name to your own.

Contributing to `oomph-lib` involves three separate repositories:
- The official repository, https://github.com/oomph-lib/oomph-lib.

     This is a remote repository, hosted on GitHub. In Git terminology this repository is known as `upstream`.

- Your remote forked version of the official repository, https://github.com/JoeCoolDummy/oomph-lib.

     This is a remote repository, hosted on GitHub. In Git terminology this repository is known as `origin` because it is the repository that you clone (as described in the next step) onto your local computer. You create the remote forked repository by going to the GitHub page for the official repository, https://github.com/oomph-lib/oomph-lib, and clicking on the fork button in the top right corner of that page.

    ![](doc/README/fork_button.png)

     Pressing that button will create a deep copy of the repository in your own account.

     Make sure you enable the Actions. These ensure that, once you have committed any changes to your remote forked repository, the self-tests are run automatically on a variety of operating systems. These tests must pass before you can issue a pull request to merge your changes into the official repository. By default the Actions are disabled, so click on the Actions button

     ![](doc/README/actions_button.png)

     and enable the workflows:

     ![](doc/README/enable_workflows.png)


- The cloned repository on your computer (obtained by cloning your forked repository).

    This repository is local to your computer and is cloned from your remote forked repository. It is where you do all your work before ultimately commiting it, via the procedure described below, to the GitHub-hosted remote repositories.

    You create this repository from the command line on your computer, using
   ```bash
   >>> git clone git@github.com:JoeCoolDummy/oomph-lib.git
   ```
   This command will create a new directory, `oomph-lib` which contains all the code and the relevant Git information. Have a look around:
   ```bash
   cd oomph-lib
   ls -l
   ```

The official (`upstream`) repository has two key branches: `main` and `development`. These have now made it onto your computer. You can list the branches as follows:
```bash
>>> git branch
  development
  main
```
You should only ever interact with the `development` branch. The `main` branch is only updated (and never by you!) when a new stable version with significant new features is available and needs to be shared with the whole wide world.

Your local forked repository is automatically aware of the GitHub-hosted remote counterpart it was cloned from (i.e. your forked one). The `git remote` command lists the URLs associated with the `origin` repository:
```bash
>>> git remote -v
origin	git@github.com:JoeCoolDummy/oomph-lib.git (fetch)
origin	git@github.com:JoeCoolDummy/oomph-lib.git (push)
```
You can therefore use `origin` as a shorthand for the full URL. It is useful to add a similar shorthand for the official `oomph-lib` repository (the `upstream` one):
```bash
>>> git remote add upstream git@github.com:oomph-lib/oomph-lib.git
```
You can check that this is now known within your local forked repository:
```bash
>>> git remote -v
origin	git@github.com:JoeCoolDummy/oomph-lib.git (fetch)
origin	git@github.com:JoeCoolDummy/oomph-lib.git (push)
upstream	git@github.com:oomph-lib/oomph-lib.git (fetch)
upstream	git@github.com:oomph-lib/oomph-lib.git (push)
```

## The workflow
We will assume that you have updated your remote forked repository from the official one, and have updated the `development` branch on your computer before doing any new
work. (We show below how this is done as part of the overall workflow but want
to make sure we know where to start; if you have just completed the steps described in 'Basic setup' above, you are now in this blissful state).

The end-goal is to emulate the Git branching model described
in https://nvie.com/posts/a-successful-git-branching-model/.

The idea is as follows: you add your changes to your *local forked forked*
repository, and then push these changes to your *remote forked* repository. Once
you are happy to share your additions with the *upstream* repository, you can
create a contribution via a mechanism called a "pull request", from your *remote forked repository* on GitHub.

This involves the following steps:

1. Create a new branch on your computer.
2. Do some work and add/commit it on your computer (i.e. locally).
3. Check it.
4. Check it again.
5. Push the changes from your local forked repository (the one on your computer) to your remote forked repository (on GitHub).
6. Wait for the automated self-tests to pass (check the Actions tab of your remote forked repository on GitHub for progress updates and log files). If they don't, resolve the issues and repeat from step 2.
7. Create a pull request to merge the changes created in that branch into the `development` branch of the official (`oomph-lib/oomph-lib`) repository (`upstream`).
8. Once the pull request has been accepted (and your changes have thus been merged into the official repository), update your remote forked repository (`origin`) on GitHub.
9. Update the `development` branch on your computer from your remote forked repository (`origin`).

## The steps in detail:


1. Start from the branch you want to work on, e.g. `development`,
   ```bash
   >>> git checkout development
   Switched to branch 'development'
   Your branch is up-to-date with 'origin/development'.
   ```
   and check that everything is clean:
   ```bash
   >>> git status
   On branch development
   Your branch is up to date with 'origin/development'.

   nothing to commit, working tree clean
   ```


2. Make a new branch to do some work (here `feature/add-new-important-headers`) and switch to it at the same time:
   ```bash
   >>> git checkout -b feature/add-new-important-headers
   Switched to a new branch 'feature/add-new-important-headers'
   ```
   (This combines the two separate commands `git branch feature/add-new-important-headers` and `git checkout feature/add-new-important-headers`).

   **Note:** the "/" is just part of the branch name and doesn't play a special role; it could in principle be replaced by any other character like another hyphen. There's also no need to prefix the branch with "feature", though this is common practice. Using the "/" has the advantage that it can be perceived as a creating a hierarchical subset of branches (so `upstream/feature/*` refers to all the feature branches at `upstream`).

   Now check (if you don't believe it) that you're really on the new branch:
   ```bash
   >>> git status
   On branch feature/add-new-important-headers
   nothing to commit, working tree clean
   ```

3. Now do some work and add the files you wish to keep (i.e. stage them for the next commit):
   ```bash
   # ...do some work...

   >>> git status
   On branch feature/add-new-important-headers
   Untracked files:
   (use "git add <file>..." to include in what will be committed)
         new-file1.h
         new-file2.h
         new-file3.h
         new-file4.h

   nothing added to commit but untracked files present (use "git add" to track)

   # Don't want to add new-file3.h and new-file4.h, so only add the others
   >>> git add new-file1.h new-file2.h

   >>> git status
   On branch feature/add-new-important-headers
   Changes to be committed:
   (use "git restore --staged <file>..." to unstage)
         new file:   new-file1.h
         new file:   new-file2.h

   Untracked files:
   (use "git add <file>..." to include in what will be committed)
         new-file3.h
         new-file4.h
   ```
   (Note: If you only want to add files that were changed and that are
   already under Git version control, you can use `git add -u`; this command will not add any newly created files.)

4. Now commit. A one-line commit message can be provided from the command line.
   The first 50 characters should be self-contained and written in "imperative"
   voice (i.e. "Update .clang-format" rather than "Updated .clang-format"). If
   no message is specified on the command line (with the `-m` flag), an editor
   will open a file for a possibly longer commit message (to which the same rules
   apply; specifically the first line should be no longer than 50 characters and
   be self-contained). You can (and in fact are encouraged to) format any longer
   messages nicely using Markdown; see
   https://github.com/adam-p/markdown-here/wiki/Markdown-Cheatsheet.
   ```bash
   >>> git commit -m "Add new-file1.h and new-file2.h to repository."
   [feature/add-new-important-headers 688ef2cfe8] Add new-file1.h and new-file2.h to repository.
   2 files changed, 0 insertions(+), 0 deletions(-)
   create mode 100644 new-file1.h
   create mode 100644 new-file2.h
   ```
   **IMPORTANT:**
   You can switch between branches at any point (e.g. `git checkout main` will get you onto the `main` branch in your local repository) but if you have not committed your work, the changes will automatically be moved across to the branch you are switching to. This is unlikely to be the desired outcome as we are specifically working on a separate branch to keep our new work isolated from the rest of the code. If you don't want to commit your changes before switching to another branch, you can use
   ```bash
   git stash
   ```
   to stash away your changes; see `git help stash` for details on how to reapply the stashed changes.

5. Now push the commit to your (GitHub-hosted) remote forked repository (`origin`).

   The first time you do this you have to this explicitly state that you want the branch to be added to the `origin` repository, using the `--set-upstream` flag,
   ```bash
   >>> git push --set-upstream origin feature/add-new-important-headers
   ```
   (where `--set-upstream` can be abbreviated to `-u`)

   Subsequent pushes will no longer require the `--set-upstream` flag, so you can just do
   ```bash
   git push origin feature/add-new-important-headers
   ```

6. Now go to the GitHub webpage for your remote forked repository
   (https://github.com/JoeCoolDummy/oomph-lib/) and click on the button with a branch symbol and the text "`main`":

    ![](doc/README/main_button.png)

    then click on the `feature/add-new-important-headers` branch in the drop-down menu that appears. When you reach the new page, click the green "Compare & pull request" button. Carefully check all changes and select reviewers, using the button on the right hand side of the screen.

    ![](doc/README/reviewers_button.png)

   By default, the pull request will attempt to merge into the `main` branch of
   the base repository (`oomph-lib/oomph-lib`) as indicated in this box:

   ![](doc/README/base_repository_button.png)

   Since commits must only be made to the `development` branch, not the `main` branch, click the "`base: main`" button and click "`development`" from the dropdown menu:

   ![](doc/README/base_repository_button_development.png)


   Provide a meaningful description of your changes in the textbox, press the button "Create pull request" and wait for somebody to merge your changes in (or get back to you with comments/requests). Note that subsequent changes (in response to discussions/requests, say) can simply be submitted to the same branch (repeating everything from step 3 above); they will automatically be added/included to the same pull request.

7. Once the pull request has been accepted and the changes made have been merged
   into the official repository, update the `development` branch on your remote forked branch. This is done most easily via the webpage: go to the `development` branch for the remote forked repository, i.e. go to https://github.com/JoeCoolDummy/oomph-lib and click on the button with a branch symbol and the text "`main`":

    ![](doc/README/main_button.png)

    Select the `development` branch from the dropdown menu.
     GitHub will show you that "This branch is n commits behind `oomph-lib:development`."
     The "Fetch upstream" button allows a comparison ("Compare") and "Fetch and merge". Go for the latter when happy. GitHub will now announce that "This branch is even with `oomph-lib:development`".

8. Now the `development` branch on your computer (cloned from the forked
   remote repository) is out of sync with the updated version at GitHub. (And if you do `git difftool feature/add-new-important-headers..development` you'll see that the `development` branch obviously(!) does not yet contain the changes we've made in the local forked `feature/add-new-important-headers` branch.) So pull the updated
   `development` branch from the `origin` repository, as follows:

   First go back to `development` branch (assuming we're still on the
   `feature/add-new-important-headers` one)

   ```bash
   >>> git checkout development
   ```

   Now `git pull` (= `git fetch` + `git merge`) the changes from the
   `development` branch in the remote forked repository (`origin`):

   ```bash
   >>> git pull origin development
   ```

   Assuming this pull doesn't conflict with any other changes made to your
   `development` branch, the local forked `development` and
   `feature/add-new-important-headers` branches should now agree. Check by
   running

   ```bash
   >>> git diff feature/add-new-important-headers..development
   ```

   The above command should not produce an output. We can then get rid of the branch on which we did the work, both locally,

   ```bash
   >>> git branch -d feature/add-new-important-headers
   ```

   and remotely on the remote forked repository (`origin`),

   ```bash
   >>> git push origin --delete feature/add-new-important-headers
   ```


## Advanced approach to pulling in upstream changes

Described below is an alternative way to pull changes from the official repository (`upstream`) into your local forked and remote forked repositories using the command-line.

1. Switch to the `development` branch
   ```bash
   >>> git checkout development
   ```

2. If you haven't done this already, add the upstream repository as a new remote repository named `upstream`
   ```bash
   >>> git remote add upstream https://github.com/oomph-lib/oomph-lib
   ```

3. Confirm that a remote named `upstream` has been added to your list of remotes
   ```bash
   >>> git remote
   ```

4. Fetch the upstream changes
   ```bash
   >>> git fetch upstream
   ```

5. Merge the changes in the upstream `development` branch into your local forked `development` branch
   ```bash
   >>> git merge upstream/development
   ```

6. Resolve any merge conflicts you may have. Strictly speaking, you
   should only be using the `development` branch to stay in sync with the upstream repository i.e. you shouldn't be adding changes to it directly yourself, so there should never be any merge conflicts.

7. Now push the changes in your local forked repository up to your remote forked
   repository (`origin`)
   ```bash
   >>> git push origin development
   ```
   Now both your local forked repository and your remote forked repository are in
   sync with the upstream (i.e. official) repository. Hurray!
