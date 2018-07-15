# CSIworks
my CSI(Channel State Infomation) works platform in python2.7

Modified from DanielHaimanot's share, add:
* MUSIC estimate function
* member functions for:

  socket: get packet from AP by socket

  csi smooth: for better eig, more info in the paper SpotFi

  csi extend: Intel5300 NIC gathers 30 subcarriers in a packet, actually 56 with the subcarriers' index is [-28..-1] and [1...28], this fun is to restore the subcarriers info by LP

  multipath remove: zero the tailed FFT components

  MUSIC aoa estimate

  pyqtgraph realtime plot
