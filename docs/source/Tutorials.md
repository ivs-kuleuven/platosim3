# Tutorials {#Tutorials}

Under <code>/docs/tutorials</code> you can find a number of <a href="http://jupyter.org/">Jupyter Notebooks</a> (previously iPython Notebooks), each of which serves as a (documented) tutorial for a specific use case.  To view such a tutorial, <code>cd</code> to the <code>/docs/tutorial</code> directory (or one of its sub-directories) and type

\code
jupyter notebook
\endcode


Open a browser and go to <a href="http://localhost:8888/tree">http://localhost:8888/tree</a> (this may happen automatically).  This will show the directory structure of the <code>/docs/tutorials</code>.  Each of the sub-directories (e.g. <code>FineGuidanceStars</code>) contains (at least) one IPYNB file (and a set of images that are included in such notebooks).  By clicking on the name of such an IPYNB file, a notebook will be shown in your browser.

The individual cells of the notebook can be run one-by-one by pressing the <code>Play</code> button in the toolbar or by selecting <code>Cell > Run Cells</code>.



## Jupyter

If you have installed <a href="https://docs.continuum.io/anaconda/install">Anaconda</a>, you can install <code>Jupyter</code> by typing the following command:

\code
conda install jupyter
\endcode

Experience has learnt us that this command causes problems in using Spyder (the GUI that comes with Anaconda).  This can be solved by executing the following series of commands:

\code
spyder --reset
conda update spyder
conda update matplotlib
conda update pytables
conda update h5py
\endcode


Launching <code>Jupyter</code> while having <code>Spyder</code> open in the same environment crashed the kernel in <code>Spyder</code>.  This can be solved by cloning the <code>platosim</code> environment and launching <code>Jupyter</code> and <code>Spyder</code> in two different environments:

\code
conda create -n platosimjupyter --clone platosim
\endcode