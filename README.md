# Setting up a SPHINX documentation

See the [documentation](https://www.sphinx-doc.org/en/master/index.html)

Install nice theme [here](https://sphinx-themes.org/sample-sites/sphinx-rtd-theme/).

Note that sphinx 3.1 is required for the [markdown integration](https://www.sphinx-doc.org/en/master/usage/markdown.html) to work with the theme 'sphinx-rtd-theme'. However, note that there is a bug with sphinx-rtd-theme for the required sphinx version to use markdown.

To keep the nice layout of the rtd-theme you need to downgrade sphinx to `pip install sphinx==1.8`

Finally to update your changes run:

make html

and see the changes in your browser (e.g. firefox):

firefox build/html/index.html
