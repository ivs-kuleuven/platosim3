@mainpage





<!-- ************************************ -->
<!-- Welcome to the PLATO Simulator Pages -->
<!-- ************************************ -->

## <a name="welcome"></a>Welcome to the PLATO Simulator Pages

The PLATO Simulator is an end-to-end software tool, designed to perform realistic simulations of the expected observations of the PLATO mission. It can, however, easily be adapted to similar types of missions.

Our simulator models and simulates time series of CCD images by including models of the CCD and its electronics, the telescope optics, the stellar field, the jitter movements of the spacecraft, and all important natural sources of noise.

Many aspects concerning the design trade-off of a space-based instrument and its performance can best be tackled through realistic simulations of the expected observations. The complex interplay of various noise sources in the course of the observations make such simulations an indispensable part of the assessment and design study of any space-based mission.


---

## <a name="gettingStarted"></a>Getting Started












<!-- Running PlatoSim3 -->
<!-- ***************** -->

## <a name="runningPlatoSim3"></a>Run PlatoSim3

 For the simulation itself, only one input file with configuration parameters (e.g. <code>/inputfiles/inputfile.yaml</code>) is required as input. 

To initiate a simulation, <code>cd</code> to the <code>/build</code> directory and type:

\code
./platosim <input file> <non-existing output file> [<log file>]
\endcode

The structure of the input file and the meaning of the individual configuration parameters are described @ref InputFileDescription "here".

Note that - before running PlatoSim3 - you must have an environment variable <code>PLATO_PROJECT_HOME</code>, set to the base folder of PlatoSim3.

If you want to use realistic PSF models instead of a Gaussian, you can download these from <a href="ftp://plato:miSotalP@ftp.ster.kuleuven.be/psf.hdf5
">our FTP server</a>.  A convenient place to store this file is in the <code>/inputfiles</code> directory.

<!-- Accessing the Output -->
<!-- ********************* -->

### <a name="output"></a>Output

The output of the simulations with PlatoSim3 is stored as an HDF5 file.  The structure of such files and how to access (and visualise) the content is described @ref OutputFileDescription "here".

Note that the x-axis is defined along the serial readout register and corresponds to the columns of the detector (and therefore of the pixel and sub-pixel map).  The y-axis corresponds to the rows of the detector (and therefore of the pixel and sub-pixel map).  In Python, for example, the <code>imshow()</code> method of <code>matplotlib</code> transposes the image as it tries to mimic Matlab.

The <code>Armadillo</code> arrays (that are used internally to store the maps) are column-major rather than row-major.



---

<!-- ******************* -->
<!-- In Case of Problems -->
<!-- ******************* -->

## <a name=inCaseOfProblems>In Case of Problems

In case you would come across problems you cannot solve yourself, please, let us know!  We would like you to use the issue tracking on GitHub rather than sending an email to the developers, as this helps to better keep track of the issues and their status.  How raise issues (and which information you must provide us with) is desribed @ref IssueTracking "here".

---

<!-- ********* -->
<!-- Reference -->
<!-- ********* -->

## <a name="reference"></a>Reference

We kindly ask you to refer to <a href="http://arxiv.org/abs/1404.1886">this work</a> in any publication the PLATO Simulator software contributes to.




<style>
ul {
    list-style-type: none;
    margin: 0;
    padding: 0;
    overflow: hidden;
    background-color: #333;
}

li {
    float: left;
}

li a {
    display: block;
    color: white;
    text-align: center;
    padding: 14px 16px;
    text-decoration: none;
}

li a:hover {
    background-color: #111;
}
</style>

<body>

<ul>
  <li><a class="active" href="#home">Home</a></li>
  <li><a href="#news">News</a></li>
  <li><a href="#contact">Contact</a></li>
  <li><a href="#about">About</a></li>
</ul>

</body>
