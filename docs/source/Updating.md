# Updating PlatoSim3 {#updating}

If you have [downloaded the software via a <code>git clone</code>](#downloading), you can update the software by executing the command

\code git pull origin master \endcode 

in the directory in which you installed PlatoSim3.

However, this will only work smoothly if you did not change any of the PlatoSim3 files or added files to the PlatoSim3 folders. The only exceptions are the <code>/inputfiles</code> and the <code>/build</code> folder, where you can add files.  Please, do not modify the original files in the <code>/inputfiles</code> folder, as this might cause problems when updating the software.  We recommend that you copy the <code>inputfile.yaml</code> file and modify the copy rather than the original file.

Please note that you have to re-build the code each time you fetch software changes. How to do this is explained @ref building "here".