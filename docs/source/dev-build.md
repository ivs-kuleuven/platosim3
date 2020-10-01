# Building the Code {#dev-build}

## Software Changes

In case of code changes (after retrieving the latest version from GitHub or after introducing changes to the code yourself)

In case you have updated the PlatoSim3 code but the dependencies remain unchanged, you only have to re-build the software but not resolve the dependencies again.  This saves you a tremendous amount of time.  The way to do this, is:

    $ cmake ..
    $ (make clean)
    $ make -j 4

---

## Updated Dependencies

At some stage, we will want to update (some of) the dependencies.  You will be notified by the developer team in case this happens.  You will then have to:

* clear the <code>/dependencies/Installs</code> directory
* and run the install script again (as described above for the first-time installation).


---

## Running the Test Harnesses

In order to be able to run the test harnesses, you must first build the code (see above) and export the required environment variables, as explained @ref ReqsRun "here".

The actual command to run the tests must be executed in the <code>/build</code> directory:

    $ ./testplatosim

---

## Troubleshooting

### Not Using the System Default Compiler

If you want to use a different compiler than the system default to execute the steps described above, you have to export two additional environment variables, <code>CCX</code> and <code>CC</code>, as follows:

    $ export CXX=g++-5
    $ export CC=gcc-5

where you may want to adapt the right-hand side of these two lines to the compiler (version) of your choice.

### Still Experiencing Problems?

If you are still experiencing problems following the instructions above, please, @ref IssueTracking "tell us about it"!  It's convenient if you can send us the full error log, which you can get hold of as follows:

    $ ./install.sh > output.txt 2> errors.txt


---

## Switching between Branches

As you can read @ref dev-branching "here", we don't just use the master branch.  If you switch to another branch and want to run simulations with the current branch, you will have to re-build the software.
