#!/usr/local/bin/python2.7
# _*_ coding: utf-8 _*_
# file: test.py
# time: 2018/7/1 下午5:05
# version: 1.0
# __author__: ChengChen
# contact: saicc4869@163.com
from Functions.WiFiArrMUSICAnalysis import *
from matplotlib import cm
import matplotlib.pyplot as plt
from sklearn.neighbors import kde
from pyqtgraph.Qt import QtWidgets, QtCore
import pyqtgraph as pg


"""
read bytes stream from server by socket
"""
csi_trace = IntelDeviceService(Intel_IWL5300_API())
csi_trace.live_stream(8090)

""" 
read from a existed binary file
"""
# csi_trace.open('data/20m-25-1.dat')
# MUSIC = WiFi_MUSIC_API(csi_trace)
#
# degList = MUSIC.csi_facto_aoa()
#
# fig = plt.figure('Kernel Density Estimate')
# X_plot = np.arange(-90, 91, 1)[:, np.newaxis]
# kde = kde.KernelDensity(kernel='gaussian', bandwidth=1).fit(degList)
# log_dens = kde.score_samples(X_plot)
# plt.plot(X_plot[:, 0], np.exp(log_dens), '-')
# plt.show()


# fig = plt.figure()
# plt.plot(abs(b.T))
# plt.show()

# print Pmu
# fig = plt.figure()
# ax = Axes3D(fig)
# rads = np.arange(-90, 90.5, 0.5)  # the resolution of degree
# tofs = np.arange(1, 301, 1)  # the resolution of ToF
# rads, tofs = np.meshgrid(rads, tofs)
# ax.plot_surface(rads, tofs, Pmu.T, rstride=1, cstride=1, cmap=cm.jet, linewidth=0, antialiased=False)
# plt.show()

# fig = plt.figure()
# plt.imshow(Pmu.T, cmap=cm.jet)
# fig.waitforbuttonpress()

# if __name__ == '__main__':
#     import sys
#
#     csi_trace = IntelDeviceService(Intel_IWL5300_API())
#     csi_trace.live_stream(8090)
#
#     if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
#         QtWidgets.QApplication.instance().exec_()
