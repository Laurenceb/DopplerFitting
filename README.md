# DopplerFitting

Doppler Trajectory fitting tools

Doppler trajectory fitting related tools. Currently this is designed to use the airspy SDR with modified host code to add timestamping. For acceptable performance it is necessary to use a GPSDO, but in the future it might be possible to get away with NTPD only. 
Ground tracking is designed to run independently at remote ground station sites. Pseudorange files can then be sent to a central location via your favourite file transfer method for the trajectory to be produced (the inversion step). Currently the inversion is in need of a better atmospheric model, but has been tested successfully with only 2 simulated station, giving >95% reliability. With three stations very reliable performance is obtained, which high accuracy in the range of meters maximum error if the pseudrange data is trimmed to the point around apogee ~+-50s.
Once files are retrieved from the remote stations, the config.txt file should be edited with the position of each receiver antenna, the best known location of the balloon atlaunch, and the carrier frequency of the transmitter. The station files should be named "station<n>.dat"for n=1,2 etc, and stored in GNU-Octave binary format at the remote stations. A helper script to simplify remote station file processing is in progress...
