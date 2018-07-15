#!/usr/local/bin/python2.7
# _*_ coding: utf-8 _*_
# file: pyqttest.py
# time: 2018/6/30 下午1:40
# version: 1.0
# __author__: ChengChen
# contact: 18717387276@163.com
# import pyqtgraph.examples
# pyqtgraph.examples.run()

# -*- coding: utf-8 -*-
"""
This example demonstrates many of the 2D plotting capabilities
in pyqtgraph. All of the plots may be panned/scaled by dragging with
the left/right mouse buttons. Right click on any plot to show a context menu.
"""
from pyqtgraph.Qt import QtGui, QtCore
import numpy as np
import pyqtgraph as pg

win = pg.GraphicsWindow(title="Basic plotting examples")
win.resize(1000,600)
win.setWindowTitle('pyqtgraph example: Plotting')

p6 = win.addPlot(title="Updating plot")
curve = p6.plot(pen='y')
def draw():
    global curve, data, ptr
    data = np.random.normal(size=(1,1000))
    curve.setData(data[0])

# def update():
#     global curve, data, ptr, p6
#     curve.setData(data[ptr%10])
#     if ptr == 0:
#         p6.enableAutoRange('xy', False)  ## stop auto-scaling after the first data set is plotted
#     ptr += 1
#     print ptr
timer = QtCore.QTimer()
timer.timeout.connect(draw)
timer.start(0)
QtGui.QApplication.instance().exec_()

## Start Qt event loop unless running in interactive mode or using pyside.
# if __name__ == '__main__':
#     import sys
#     if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
#         QtGui.QApplication.instance().exec_()
