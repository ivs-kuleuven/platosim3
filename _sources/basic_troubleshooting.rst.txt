Troubleshooting
===============

In case of any problems and questions related to PlatoSim we encourage you to use the `GitHub issue tracking system <https://github.com/IvS-KULeuven/PlatoSim3/issues>`_, rather than sending one of the PlatoSim developers an email. This helps you and us to better keep track of which problems arise and what their status are.

Browse to the PlatoSim repository on GitHub and select the Issues tab (see figure below). The issue were are currently working on are listed on this page.

.. figure:: ../figures/github_issueTracking.png
   :align: center
   :width: 700
   
You can raise a new issue by clicking the green New issue button. This will direct you to a page (shown in the figure below).
 
.. figure:: ../figures/github_newIssue.png
   :align: center
   :width: 700

To make the issue tracking faster and smoother please provide the following information:

* PlatoSim version (bash command: ``platosim --version``)
* Appropiate issue **label** (see right-hand menu)
* Concise explanation of the problem
* What you tried to investigate the problem
* Please provide the ``inputfile.yaml`` and ``log.txt`` files (in debug mode: ``<logLevel>=3``)  
* In case of installation or execution problem please provide:

  * Your operating system
  * Version of gcc (``gcc --version``) or clang (``clang --version``)
  * Version of CMake (``cmake --version``)
