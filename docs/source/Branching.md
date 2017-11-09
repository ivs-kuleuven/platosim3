# Branching Strategy for PlatoSim3 {#branching}

We have adopted the <a href="http://nvie.com/posts/a-successful-git-branching-model/">branching strategy of Vincent Driessen</a>, which means that the following two branches will be used permanently:

* the <b>master</b> branch, with the stable (i.e. tested) code, which is considered safe to use
* and the <b>develop</b> branch, with the latest development, which might be (highly) experimental at times (use at your own risk).


---


## Branches

To switch to a specific branch, use the command:

\code
git checkout BRANCH_NAME
\endcode

To grab all branches and get an overview, use the following commands:

\code
git fetch
git branch
\endcode

The first command pulls all remote branches (not only the one you are currently working on) and the second command gives you an overview of all branches you have on your system.

---

## Release Candidates and Releases

Release candidates and releases correspond to tagged versions of the master branch.  To start using a specific release or release candidate, you have to check out the version of the master branch with a specific tag to a new branch, like this:

\code
git checkout -b NEW_BRANCH_NAME TAG_NAME
\endcode

We will send around the tag name for new release candidates and releases once they become available.  To get an overview of the available tags, use

\code
git tag -l
\endcode

---

## Switching between Branches

If you switch to another branch and want to run simulations with the current branch, you will have to re-build the software, as described @ref building "here".

It is of course always possible that you come across problems when using either of these branches!  Please, tell use about them via the GitHub @ref #issueTracking "issue tracking" system. 