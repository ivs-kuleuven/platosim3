# Issue Tracking {#IssueTracking}

You are encouraged to use the issue tracking system of GitHub to report any problem you would come across, rather than sending one of the PlatoSim developers an email.  This helps you and us to better keep track of which problems arise and what their status is.

Browse to the PlatoSim3 repository on GitHub and select the <code>Issues</code> tab (see Fig. 1).

@image html /images/issueTracking.png "Figure 1: Issue tracking on GitHub."

The issue were are currently working on are listed on this page.

You can raise a new issue by clicking the green <code>New issue</code> button.  This will direct you to a page as shown in Fig. 2.

@image html /images/newIssue.png "Figure 2: Reporting a new issue on GitHub."





## <a name=whichInformationToProvide>Which Information to Provide?

In the upper part (where it says "Title") you can enter a concise description of the problem.  A more elaborate description can be entered in the bottom part (where it says "Leave a comment").  In order to save some time, it is suggested to provide us with the following information in the issue:

- Which version of the PlatoSim software are you using?  This can be easily retrieved by typing <code>git describe</code>
on the command line in your installation directory.
- How can the problem be reproduced?  You can, e.g., attach an example script to the issue.
- In case the results are different than you expected, please explain us which numbers you expected and why.
- Which system are you using?
	- operating system
	- gcc/clang (<code>gcc --version</code> or <code>clang --version</code>)
	- CMake (<code>cmake --version</code>)
- What have you tried to solve the problem?
- Error and log messages.