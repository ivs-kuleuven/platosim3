Introduction (VSC)
==================

.. warning::

   Under construction!

**Introduction to High Performance Computing (HPC) on the VSC**

* Infrastructure
* Parallel computing frameworks
* Data storage on cluster
* Setup project and create job script
* Transfer data


.. _hpc_installation_software:

Infrastructure
--------------

The good thing about Tier-2 infrastructure is that you will start with 2000 introduction credits, you pay as you go, there is no computaion limits nor application deadlines, and lastly the available run time is up to 1 week. For bigger simulations you will need to turn to Tier-1 and adhere to application deadlines (February 3rd, June 5th, and October 2nd) and a lower/upper computation limit. The aforementioned introductionary session (HPC-Intro) is highly recommended for a broader review on the infrastructures and the different clusters, however, since Tier-2 will be your first go-to cluster infrastructure, we present here a small list of available **CPU types**, *clusters*, and their associated node types (i.e. Haswells, Skylake, etc.)

- **Compute (thin) nodes:** 
  - *ThinKing:* 
    - Haswells: 48n/86n x 24c with 64/128 GB
  - *Genius:*
    - Skylake: 92n x 36c with 192 GB
	- CascadeLake: 120n x 36c with 192 GB
- **Large memory nodes:**
  - *Genius:*
    - Skylake: 10n x 36c with 768 GB
	- Superdome: 8n x 14c with 6 TB
- **GPU nodes:**
  - *Genius:*
    - Skylake: 20n x 36c with 192 GB + nodes have 4xP100 devices with 16 GB
	- CascadeLake: 2s x 36c with 192 GB + nodes have 8xV100 devices with 32 GB

Here to keep the format short, "n" and "c" refers to nodes and cores, respectively. Keep in mind that depending on which cluster you want to run your simulations on, you might need to login to a different login-node or import special modules. Check these requirements at the [VSC Hardware specifications](https://vlaams-supercomputing-centrum-vscdocumentation.readthedocs-hosted.com/en/latest/hardware.html#hardware).  




Data storage
------------

When logged in, you can quick-access your internal directories using the enviroment variables

- **home:** `$VSC_HOME`
- **data**: `$VSC_DATA`
- **scratch:** `$VSC_SCRATCH`
- **node scratch:** `$VSC_SCRATCH_NODE` (Only while running jobs)
- **staging:** `/staging/leuven/stg_XY` (KU Leuven access only)
- **archive:** `/archive/leuven/arc_XY` (KU Leuven access only)

Since the **home** folder (3 GB) is only for small files (such as configuration files), and since the **scratch** (100 GB) are not backed-up, best practice is to clone and build PlatoSim3 in the **data** directory (75 GB). A word concerning data storage is, **data** is well suited for initial/final results, but ill suited for intensive or parallel read-and-write (Input/Output) jobs, however, **scratch** are optimized toward I/O situations and storage can be extended for free by request. *At runtime only* a **node scratch** (200 GB) is an availble scratch on the compute node that provide the optimal I/O workflow, since no information is inter-changed over the network. The **staiging** and **archive** foldes are on demand and can only be requested by KU Leuven associates. As the name suggest, **staging** (from 1 TB) is well suited for staging your jobs/simulations as it is accessable from your login and the compute nodes, as well as for sharing data with collaborations (since you need to grand access to collaborators for the deafult folders). Like **scratch**, **staging** share the same capabilitiesis of fast file system optimized for parallel I/O simulations. Lastly, **archive** (from 1 TB) is a more slow backup and storage folder that likewise are easy platform for sharing data with collaborators, but only available from login nodes. **Remember to keep track of your storage using the command `myquota`.**

The above mentioned path-variables becomes very useful when running simulations, however, notice that these variables are not back-linked to your local machine, meaning that from your local machine, you'll need to specify the full path, e.g. `/<folder>/leuven/<X>/vsc<XY>`). For more information have a look at [KU Leuven Storage](https://vlaams-supercomputing-centrum-vscdocumentation.readthedocs-hosted.com/en/latest/leuven/tier2_hardware/kuleuven_storage.html) and [Where to store what kind of data](https://vlaams-supercomputing-centrum-vscdocumentation.readthedocs-hosted.com/en/latest/access/where_can_i_store_what_kind_of_data.html). If you would need more storage you can [Request extra storage](https://admin.kuleuven.be/icts/onderzoek/hpc/hpc-storage).




.. _hpc_tutorial_frameworks:

Frameworks
----------

: worker and atools

Starting from the simple job script we would immediately like to take advantage of multiple cores, but increasing the number of cores within the default PBS/Torque cluster environment will still only make our job run once on the first node of our allocation. To start processes on other nodes, we either choose to use tools like ``pbsdsh`` or change framework to e.g., ``worker`` or ``atools`` (notice that the MPI tool ``mpirun`` are specifically designed to start a distributed memory which is not yet applicable for PlatoSim3). 

Altough in theory starting a single-core program on each assigned core should be straight forward, the scheduler do not always serve you justice regarding the job CPU capacity, resulting in a inefficiency use of job resources. The ``pbsdsh`` command assures that all programs will execute under the full control of the resource manager on the cores allocated to your job. We refer the reader to the page `Starting programs in a job <https://vlaams-supercomputing-centrum-vscdocumentation.readthedocs-hosted.com/en/latest/jobs/starting_programs_in_a_job.html?highlight=pbsdsh#starting-a-single-core-program-on-each-assigned-core>`_ for the use of ``pbsdsh``. 

`How to conveniently run many similar computations <https://vlaams-supercomputing-centrum-vscdocumentation.readthedocs-hosted.com/en/latest/jobs/worker_or_atools.html#worker-or-atools>`_ are optimal using the ``worker`` and ``atools`` frameworks, since they are specifically designed to automatically handle the delegation of cores from a large bundle of sigle-core jobs. Advantageous compared to ``pbsdsh``, this also means that one can specify the numbers of cores explicit (e.g., ``ppn=10``) and the scheduler will automatically try to give you those resources. A computation for an individual input file, or, equivalently, an array id, is called a work item in both frameworks. For clarification, note that the ``worker`` and ``atools`` frameworks has been developed for *parameter variation* and *array-like patterns* jobs, which we will see in the next section are the main use of PlatoSim3.

Choosing between the ``worker`` and ``atools`` framework depends on the following items

- **Single Work Item Walltime (SWIW):** The time efficiency for ``worker`` and ``atools`` in respect to the walltime per work item, is by anology the same problem observational astronomers face regarding exposure time vs. overheads (from CCD readout): When work items take only a short time (< 1 min) to complete, the overhead for starting new work items will be a considerable time fraction compared to the total computational run. Since ``worker`` do not rely on the scheduler to start individual work items, and ``atools`` does, ``worker`` can efficiently be used for work items (``1s < SWIW < 1min``), where ``atools`` should only be considered for a sigle work item walltime of > 1 min. However, both ``worker`` and ``atools`` peak in performance using (``1 min < SWIW < 1 day``).

- **Number Of Work Items (NOWI):** While computing a very high number (in the hundres) of SWIW, ``worker`` is the better choice as it run these as a single job, while ``atools`` runs very efficiently for jobs in the bold park of tens to a few hundred. Hence, again the efficiency comes down to how oppressed the scheduler will be upon your job's NOWI.

- **Cluster job policy:** Lastly, the choice of framework depends heavilly on if the desireable cluster has accounting enabled, a single user policy, or a shared user policy. Generally, it is more favorable to use ``atools`` on clusters with a sigle or shared policy (as ``worker`` more agressively block nodes and allows less flexibility for the scheduler), whereas the single job workflow of ``worker`` are more favorable from an accounting (pay as you go) perspective, as charges will be made for each individual work item.

From the items presented above we might conclude that initially ``worker`` is the better framework for PlatoSim3 since our main aim is to simulate a large amount of images/imagttes as observed from all of the 6 x 4 normal cameras. Thus, a choice of framework will specifically depend on the number of cameras included, the baseline of the timeseries, and if parameter variantion is desired along the timeseries (e.g. to investigate instrumental effects, mission efficiency degration, etc.). For optimal use have a look at the commands ``tail``, ``watch``, ``wload``, ``timedrun``, and ``wresume``.

There might be special cases where simulations favorably can be submitted as ``atools`` jobs, however to keep the price to a minimum, if you run simulations on an accounting cluster (like *Genius*) try always to maximumize the use of node-cores before including more nodes. For more information have a look at `worker quick start guide <https://vlaams-supercomputing-centrum-vscdocumentation.readthedocs-hosted.com/en/latest/jobs/worker_framework.html#worker-framework>`_,the `worker documentation <https://worker.readthedocs.io/en/latest/>`_, and the `atools documentation <https://atools.readthedocs.io/en/latest/>`_.


.. _hpc_intro_job:

Job scripts
-----------

To grasp the information to come, let's quickly take a look on how to declare nodes, cores, memory, etc. in a so-called "job-script" (``<JobScriptName>.pbs``) that can be used for the HPC facility to read and execute our simulations. The ingredients of a simple job script using the default *PBS/Torque* framework are

.. code-block:: shell

   #!/bin/bash
   #PBS -l nodes=1:ppn=1
   #PBS -l pmem=1gb
   #PBS -l walltime=00:00:01
   #PBS -N <JobName>
   #PBS -A <ProjectCredits>
   #PBS -m abe -M <Email>
   cd $VSC_DATA/PlatoSim3

- **She-bang:** The first line simply indicates that we are writing a bash script and shall use syntax of such, a so-called *she-bang*.

- **Resources:** The following lines of code starting with ``#PBS`` are the *job resources* for the scheduler and for this simple example we show:
  - *Nodes and Cores:* The ``nodes`` and ``ppn`` parameters refers nodes and cores (or processors per node; CPUs), and explicit we choose here the number of nodes and number of cores per node, respectively.
  - *Memory:* is the memory RAM needed for the job. Units can be given in (``kb, mb, gb, ..``). Notice that the memory for most HPC clusters nodes comes in a few hundred of GBs. 
  - *Walltime:* is your estimated time for the computations to finish. If your job exceeds the specified walltime, it will be killed, so the walltime should not be underestimated. 
  - *Job name:* For a greater overview it is very handy to give your job a name using ``-N``.
  - *Project credits:* First, use the identifier from your *introductionary credits* and later the *project credits* belonging to a VSC group (created by your supervisor/project leader). Upon submission this resource is mandatory and your script will not run without specifying it.
  - *Notification:* ``#PBS -m abe -M <Email>`` tells the scheduler to send a notification to the specified email (``-M``) when the job begins (``-m b``), ends normally (``-m e``), or aborts due to errors (``-m a``). Seen here you can specify everything altogether by ``-m abe``.
  
- **Working directory:** should be specified by ``cd $PBS_O_WORKDIR`` if you submit your job within the folder where you script(s) will be running.

A more ellaborate example on how to create a simple job script, can be found in the VSC section [Running Jobs](https://vlaams-supercomputing-centrum-vscdocumentation.readthedocs-hosted.com/en/latest/jobs/running_jobs.html). Also, more information can be found on how to [Specify Job Resources](https://vlaams-supercomputing-centrum-vscdocumentation.readthedocs-hosted.com/en/latest/jobs/specifying_resources.html#resource-specification). 

The job script "skeleton" presented here do not take advantage of the great capabilities of running computations in parallel, hence, do take a look at the availble resources within the Tier-2 infrastructure regarding nodes, cores, memory, etc. The concept of parallelism is equal to employing multiple processors (cores/CPUs) for a single problem, and good instruction are to be found on how to run [Applications in parallel](https://vlaams-supercomputing-centrum-vscdocumentation.readthedocs-hosted.com/en/latest/software/parallel_software.html). 

Another solution of parallization is to use *GPU partitions*. On *Genius* the choice is between P100 (*skylake*) or V100 (*cascadelake*) GPUs, and the following needs to be added to your job script

.. code-block:: shell

   #PBS -l partition=gpu
   #PBS -l nodes=1:ppn=1:gpus=1:<GPU-node-name>

Remember that there is a lower limit/upper limit on used GPUs for each GPU partition. 

In cases where a huge amount simulations needs to be performed, you can use the *large memory nodes* for which there are 3 ways of requesting larger memory than normal

- Use thin nodes but lower ``ppn`` and increase ``pmem``
- Use large memory node ``skylake`` with: ``#PBS -l partition=bigmem``
- Use stand-alone large memory node ``superdome``: First ``module load superdome`` and add the following lines to your job script: ``#PBS -l partition=superdome -L tasks=1:place=numanode=1:lproces=1``

For more information take a look at the [Genius quick start guide](https://vlaams-supercomputing-centrum-vscdocumentation.readthedocs-hosted.com/en/latest/leuven/genius_quick_start.html?highlight=submit%20to%20a%20GPU%20node).

Submitting and managing jobs
----------------------------

Next, one can sumit the job to the cluster scheduler by `qsub <JobName>.pbs`. The job-ID that is returned after a submission is your unique identifier to follow a job status (`qstat [-f <JobID>]`), delete a job from the que (`qdel <JobID>`), etc. Further information can be found on [Submitting and Managing Jobs](https://vlaams-supercomputing-centrum-vscdocumentation.readthedocs-hosted.com/en/latest/jobs/submitting_and_managing_jobs_with_torque_and_moab.html) but let us sumarize a small list of the most-used commands

.. code-block:: shell

   qsub <JobName>.pbs
   qstat <JobID>
   qdel <JobID>
   showq
   checkjob -v <JobID>
   showstart <JobID>
   pbstop
   monitor <JobScript/Software>

The the job output is by default saved in two files called `<JobName>.o<JobID>` and `<JobName>.e<JobID>`, which follow default Linux `stdout` (standard output) and `stderr` (standard error) output, respectively. Notice besides these default output files, the PlatoSim3 HDF5 output file will likewise be created upon a successful run. If you receive an email that your job terminated abnormally, have a look at the errors written to your `<JobName>.e<JobID>`. For trouble shooting have a look at the following pages: [Time line](https://vlaams-supercomputing-centrum-vscdocumentation.readthedocs-hosted.com/en/latest/jobs/why_doesn_t_my_job_start.html) and [Missing output](https://vlaams-supercomputing-centrum-vscdocumentation.readthedocs-hosted.com/en/latest/jobs/what_if_jobs_fail_after_starting_successfully.html). For a more ellaborate description on the output file format have a look at [Specifying job name, output files and notifications](https://vlaams-supercomputing-centrum-vscdocumentation.readthedocs-hosted.com/en/latest/jobs/specifying_output_files_and_notifications.html#specifying-output-files-and-notifications).


Transfering files
-----------------

Since the access from your local machine and your VSC account is protected by a firewall, tranfering files relies on your ssh authentication. FileZilla is an easy installable GUI for mvoing data across local to server platforms and is available for both Linux, Mac, and Windows users. However, if you prefer to use bash, files can be transfered by the `rsync`, `scp`, `sftp`, or `sshfs` commands. For more details see the [rsync documentation](https://rsync.samba.org/documentation.html) or the [scp/sftp VSC tutorial](https://vlaams-supercomputing-centrum-vscdocumentation.readthedocs-hosted.com/en/latest/access/data_transfer_with_scp_sftp.html).

For your convenience we here present a bash script called `sshfs-vsc.sh` that allows to *mount* and *dismount* your VSC account to your local machine through a folder of your choice. This script is made available to you in the `$PLATO_PROJECT_HOME/python/examples/highPerformanceComputing/` folder. Within this script, change the *User parameters* being the directory name of your choise, `VSC=<FolderName>`, and the `X=<X>` and `Y=<Y>` user ID


.. code-block:: shell

   #!/usr/bin/env bash

   # User parameters
   VSC="<folder name>"
   X=<X>
   Y=<Y>

   # Provide directory structure
   if [ ! -d "$VSC" ]; then
       echo "Creating VSC file structure locally"
       mkdir $VSC
       mkdir $VSC/home
       mkdir $VSC/data
       mkdir $VSC/scratch
   fi

   # If no arguments are given write usage message
   if [ -z "$1" ]; then
       echo "Usage: server-sshfs <option>"
       echo "       mount    :   <option> = m"
       echo "       dismount :   <option> = d"
       exit 1

   else
   
       # Mount VSC account
       if [ $1 = "m" ]; then
           echo "Mounting VSC account"
	   sshfs vsc$X$Y@login.hpc.kuleuven.be:/user/leuven/$X/vsc$X$Y $VSC/home
	   sshfs vsc$X$Y@login.hpc.kuleuven.be:/data/leuven/$X/vsc$X$Y $VSC/data
	   sshfs vsc$X$Y@login.hpc.kuleuven.be:/scratch/leuven/$X/vsc$X$Y $VSC/scratch
       fi

       # Dismount VSC account
       if [ $1 = "d" ]; then
           echo "Dismounting VSC account"
	   fusermount -u $VSC/home
	   fusermount -u $VSC/data
	   fusermount -u $VSC/scratch
       fi
   fi

By running the script without any input argument `./server-sshfs.sh`, a usage message will appear and remind you to use either the mount (`m`) or dismount (`d`) option, for which the latter obviously does not do anything (except of reporting an error) unless you actually did mount your VSC account first. Upon first execution, the script will automatically create a folder `<FolderName>` and within setup a file system of folders (**home**, **data**, and **scratch**) alike your VSC account. Hopefully seen quite intuatively from the code, you can also add your potential **staging** or **archive** folders. 

Files can now simply be moved across your local-server mount, but remember not to dismount before making sure that files are fully moved/copied, which can take some time. Advancing to a rare problem, since the script here works on the assumption that your distribution independently can select free ports locally, the script will not work if the number of folders on your VSC account are greater than local free ports.

---

Words on scheduling
-------------------

Before showing real simulation examples, we need to add a few words on scheduling. Creating a proper job script for your simulations depends on your preparation. The following can be used as a good-practice check list before creating and submitting a job

- **Check available storage** on your VSC account using the command `myquota`. Make sure that you have enough RAM to run your job(s).
  
- **Check the availability** of nodes, cores, and GPUs for the different clusters (`clusterview`); check the overview of the different partitions (`queueview`); check the que from the schedular's perspective (`showq`). Some clusters are faster than others, but if these are occupied, it might be better to use one that is available already. You can specify the cluster by name on the same line as `#PBS -l nodes=1:ppn=1:<cluster>`.
  
- **Number of nodes and cores/GPUs** that is available depends on which cluster you want to run your simulations on. Notice that jobs requesting a large number of nodes typically spend quite some time in the queue, compared to equal computational heavy jobs, that use less nodes but an increased number of cores. In short for normal applications; try to use all cores on a single node before adding more nodes to your computation.
  
- **Specify the walltime** by measuring the PlatoSim3 simulation time of a single/subset of images on your local machine, and multiply with the appropiate factor in correspondance to the total amount of image output. While using parallelization remember to take this into account (i.e. divided the run time with the number of cores), however, please do not underestimate the walltime as your job will terminate.
  
- **How much RAM memory** per core is needed for your job. It is important to check that job memory never exceeds the maximum RAM memory of each cluster node (`nmax`). Use this inequality as a check: `ppn * pmem < nmax - 8 GB`, where the 8 GB is to insure the operating system/other services can function properly.



Debugging and test jobs
-----------------------

Debugging and testing job scripts can nicely be handled by two dedicated nodes on *Genius* (one GPU and one CPU node) and *ThinKing* (two CPU nodes). Thus, remember to take advantage of these nodes if your job fits this description, as the walltime are reduced to 30 minutes. Simply specify *quality of service* to debugging (`-l qos=debugging`) and a walltime of maximum 30 minutes in your job script, e.g. `-l walltime=30:00`. 


Resuming jobs
-------------

Since calculations can get rather large quickly (involving a lot of nodes and cores) something will most likely go wrong at some point. If not all work items are finished during your specified walltime, your job can be resumed with ``wresume -jobid <jobID>``, where you most likely also will add a longer walltime or maybe more nodes or cores to keep the walltime. If some work items failed, after correcting the mistake, it is possible to resume and complete these failed items using ``wresume -jobid <jobID> - retry``.
