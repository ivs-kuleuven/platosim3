
This variant of PlatoSim3 should enable the user of PlatoSim to:

- set the size and position of the subfield
- set the next jitter step
- receive a finished imagette
- end platosim at a non determined point of time

during the runtime of the program via the respective socket connection.

---------------------------------------------------------------------------------------------------------

Changes:

- introduction of the zmq library during the installation of platosim

- addition of the ControlTcpConnection group within the inputfiles, where it is specified whether the new functions and which socket addresses are to be used

- a new jitterFromNetwork class which inherits from jitterGenerator and governs the receiving of new jitter steps from a server

- a running simulation will wait for a new jitter step indefinitly (Maybe there should be a Timeout)

- before a new exposure generation is started, it is checked whether there is a request to change the position of the subfield. this will not stop the simulation, except for the first exposure. In that case, PlatoSim will wait indefinitly for a message from a server

- since the window position and size is no longer static the psf and some additional maps (flatfield etc.) have to be changed during a running simulation. This causes problems when the maps are written twice to the hdf5 file. This is why I disabled that for the moment  

- whenever a imagette is finished it will be send to client

- the simulation can now be stopped via a jitter step message. When the jitter source is "FromNetwork" this is the only way to stop it

---------------------------------------------------------------------------------------------------------

messages:


The used package for each message is a zmq::message_t object which is easily converted from and to a std::string object. A message string consists of the following values:


jitterMessage:

<bool endOfSimulation> <double timeStamp> <double yaw> <double pitch> <double row>

winPositionMessage:

<uint subFieldRows> <uint subFieldCols> <uint subFieldZeroPointRow> <uint subFieldZeroPointCol> <uint orientationAngle>

imagetteMessage:

<uint imagetteNumber> <uint imagetteRows> <uint imagetteCols> <int each imagette element in row major order>


