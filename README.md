# PlatoSim's Sphinx documentation


## Setting up a Sphinx documentation

The documentation was set up following:

* Install Spinx: `conda install sphinx`
* [The Sphinx documentation](https://www.sphinx-doc.org/en/master/index.html). 
* [Installing a theme](https://sphinx-themes.org/sample-sites/sphinx-rtd-theme/).

Note that sphinx 3.1 is required for the [markdown integration](https://www.sphinx-doc.org/en/master/usage/markdown.html) to work with the theme 'sphinx-rtd-theme'. However, note that there is a bug with sphinx-rtd-theme for the required sphinx version to use markdown.

To keep the nice layout of the rtd-theme you need to downgrade sphinx to `pip install sphinx==1.8`

## How to update the documentaion page:

Changes needs to be made to the files within `source/`. Note all figures are placed within the `figures/` directory. 

If you just want to see the changes in a locally you first need to:

	make html
	firefox build/html/index.rst

We here use firefox as a browser to visualise the html pages.

When you're happy with the changes simply run:

	make github

This will both create the html files and copy them to root. Now commit your changes:

	git add *.html *.js *.inv build/ source/ _sources/
	git commit -m "<message-of-update>

Note that the above command add all files that have changed. Now push it to your fork:

	git push origin sphinx
	
Compare your fork `sphinx` with the upstream branch `sphinx` and merge the branches.

Wait 1-5 min for the webpage to change. 

## How to change GitHub documentation branch

You can always update the branch for which the GitHub looks for the documentation:

* Go to `Settings` and `Pages`
* Under Branch, change to your new branch.
