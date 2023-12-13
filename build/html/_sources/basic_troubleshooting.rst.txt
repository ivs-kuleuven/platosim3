Troubleshooting
===============

In case of problems and/or questions related to PlatoSim software we encourage you to use the `GitHub issue tracking system <https://github.com/IvS-KULeuven/PlatoSim3/issues>`_, rather than sending one of the PlatoSim developers an email. This helps you and us to better keep track of *active* problems and what their status are.

Browse to the PlatoSim repository on GitHub and select the Issues tab (see figure below). The issues we are currently working on are listed on this page, hence before opening a new issue, first check if your issue is not reported already.

.. figure:: ../figures/github_issueTracking.png
   :align: center
   :width: 700
   
You can raise a new issue by clicking the green ``New Issue`` button. This will direct you to a page (shown in the figure below).
 
.. figure:: ../figures/github_newIssue.png
   :align: center
   :width: 700

To make the issue tracking faster and smoother please provide the following information:

* PlatoSim branch and version (use command: ``platosim --version``)
* Appropriate issue label (see right-hand menu bottom ``Labels``)
* Concise explanation of the problem
* What you have tried to solve the problem
* Please provide the ``inputfile.yaml`` and ``log.txt`` file (in debug mode: ``<logLevel>=3``)  
* In case of a installation or execution problem, please further provide:

  * Your operating system
  * Version of gcc (``gcc --version``) or clang (``clang --version``)
  * Version of CMake (``cmake --version``)
