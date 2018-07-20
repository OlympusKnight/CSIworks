#!/usr/local/bin/python2.7
# _*_ coding: utf-8 _*_
# file: WiFiArrMUSICAnalysis.py
# time: 2018/7/1 下午4:51
# version: 1.0
# __author__: ChengChen
# contact: saicc4869@163.com
import copy
import cmath
from Tools.IntelWiFiDevice import *
from scipy.fftpack import fft, ifft
from math import pi
import matplotlib.pyplot as plt


class WiFi_MUSIC_API():

    def __init__(self, symbol):
        self._csi_trace = symbol
        self._c = 3e8
        self._f2 = 2.412e9
        self._f5 = 5.825e9
        self._fs = 312.5e3
        self._d = 0.03

    def csi_facto(self, countdown=50):
        degList = []
        for idx in range(0, self._csi_trace.get_symbol_count()):
            if self._csi_trace.get_receiver_count(idx) == 3 and countdown > 0:
                countdown -= 1
                csi_value = np.array(self._csi_trace.get_csi_symbols()[idx])
                csi_matrix = np.zeros([3, self._csi_trace.device_subcarriers()], dtype='complex64')
                for subcarriers_idx in range(0, self._csi_trace.device_subcarriers()):
                    csi_matrix[:, subcarriers_idx] = csi_value[subcarriers_idx][0]  #only one transmitter
                smoothed_csi = np.matrix(self.__smooth_csi(self.__csi_extend_57(csi_matrix)))
                deg, tof = self.__csi_find_aoa_spotfi(smoothed_csi)
                degList.append(deg)
        return degList

    def csi_facto_aoa(self, countdown=100):
        print self._csi_trace.get_symbol_count()
        degList = []
        for idx in range(0, self._csi_trace.get_symbol_count()):
            if self._csi_trace.get_receiver_count(idx) == 3 and countdown > 0:
                # countdown -= 1
                csi_value = np.array(self._csi_trace.get_csi_symbols()[idx])
                csi_matrix = np.zeros([3, self._csi_trace.device_subcarriers()], dtype='complex64')
                for subcarriers_idx in range(0, self._csi_trace.device_subcarriers()):
                    csi_matrix[:, subcarriers_idx] = csi_value[subcarriers_idx][0]  #only one transmitter
                # csi = np.matrix(self.__csi_extend_57(csi_matrix))
                csi = np.matrix(csi_matrix)
                # plt.plot(np.angle(csi).T)
                # plt.show()
                removedcsi = np.matrix(self.__csi_remove_multipath(csi))  #remove multipath
                deg = [self.__csi_find_aoa(removedcsi)]
                degList.append(deg)
        return degList

    @staticmethod
    def __csi_extend_57(csi_matrix):
        amp30 = abs(csi_matrix)
        fi30 = np.angle(csi_matrix)
        fii30 = copy.copy(fi30)
        for rx in range(0, 3):
            offset = 0
            for i in range(1, len(fi30[rx, :])):
                if abs(fi30[rx, i] - fi30[rx, i-1]) > pi-1:
                    if fi30[rx, i-1] > fi30[rx, i]:
                        offset += 2*pi
                    else:
                        offset -= 2*pi
                fii30[rx, i] = fi30[rx, i] + offset
        amp57 = np.zeros([3, 57], dtype='float64')
        fii57 = np.zeros([3, 57], dtype='float64')
        csi57 = np.zeros([3, 57], dtype='complex64')
        scidx20M = range(0, 27, 2) + range(27, 56, 2) + [56]
        c_scidx20M = range(1, 27, 2) + range(28, 55, 2)
        amp57[:, scidx20M] = amp30
        fii57[:, scidx20M] = fii30
        fii57[:, c_scidx20M] = (fii57[:, [i-1 for i in c_scidx20M]] + fii57[:, [i+1 for i in c_scidx20M]])/2
        # return fii57
        amp57[:, c_scidx20M] = 10*np.log10((10**(amp57[:, [i-1 for i in c_scidx20M]]/10) + 10**(amp57[:, [i+1 for i in c_scidx20M]]/10))/2)
        # return amp57
        for rx in range(0, 3):
            for i in range(0, 57):
                csi57[rx, i] = cmath.rect(amp57[rx, i], fii57[rx, i])
        return csi57

    @staticmethod
    def __smooth_csi(e_csi):
        smoothed_csi = np.zeros([30, 86], dtype='complex64')
        for i in range(0, 15):
            smoothed_csi[i, :] = np.concatenate((e_csi[0, i:(i+43)], e_csi[1, i:(i+43)]), axis=0)
        for i in range(15, 30):
            smoothed_csi[i, :] = np.concatenate((e_csi[1, (i - 15):(i + 28)], e_csi[2, (i - 15):(i + 28)]), axis=0)
        return smoothed_csi

    def __csi_find_aoa_spotfi(self, smoothed_csi):
        MUSIC_S = np.dot(smoothed_csi, smoothed_csi.H)
        value, vector = np.linalg.eigh(MUSIC_S)
        noiseVectorIdx = [i for i in range(0, len(value)) if value[i] < np.max(value) * 1e-4]

        En = np.matrix(vector[:, noiseVectorIdx])
        degrees = np.arange(-90, 90.5, 0.5)  # the resolution of degree
        tofs = np.arange(1e-9, 301e-9, 1e-9)  # the resolution of ToF
        SP = np.zeros([len(degrees), len(tofs)], dtype='float64')
        for deg_idx in range(0, len(degrees)):
            phi = np.exp(-1j*2*pi*self._d*np.sin(degrees[deg_idx]*pi/180)*self._f5/self._c)
            for tof_idx in range(0, len(tofs)):
                omega = np.exp(-1j*2*pi*self._fs*tofs[tof_idx])
                half = [omega**i for i in range(0, 15)]
                a = np.matrix(half + [tmp*phi for tmp in half]).T
                SP[deg_idx, tof_idx] = 1/(np.abs(a.H * En * En.H * a))
        rol, col = divmod(np.argmax(SP), len(tofs))
        return degrees[rol], tofs[col]

    def __csi_find_aoa(self, csi):
        MUSIC_S = csi * csi.H
        value, vector = np.linalg.eigh(MUSIC_S)
        # noiseVectorIdx = [i for i in range(0, len(value)) if value[i] < np.max(value) * 1e-4]
        noiseVectorIdx = [0]
        En = np.matrix(vector[:, noiseVectorIdx])
        degrees = np.arange(-90, 91, 1)  # the resolution of degree
        tmp = np.matrix([0, 1, 2]).T
        # SP = np.zeros([len(degrees), 1], dtype='float64')
        A = np.exp(-1j * 2 * pi * self._d * tmp * np.sin(degrees * pi / 180) * self._f5 / self._c)
        Ena = En.H * A
        den = np.diag(Ena.H * Ena)
        p = np.abs(1 / den).tolist()
        return degrees[p.index(max(p))]

    def __csi_remove_multipath(self, csi):
        cir = ifft(csi)
        # plt.plot(abs(cir.T))
        # plt.show()
        n=1
        numZeros = np.zeros([3, (30-n)], dtype='complex64')
        removedCSI = fft(np.hstack([cir[:, 0:n], numZeros]))
        # plt.plot(abs(np.hstack([cir[:, 0:n], numZeros]).T))
        # plt.show()
        return removedCSI




