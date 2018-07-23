import functools as func
from Tools.ExtractorService import *
import copy
import cmath
import struct
import time
from io import StringIO
from socket import *
from scipy.fftpack import fft, ifft
from math import pi
import pyqtgraph as pg
from pyqtgraph.Qt import QtWidgets
from sklearn.neighbors import kde
from Functions.fastDTW import fastDTW
import matplotlib.pyplot as plt


class Intel_IWL5300_API():
    def __init__(self):
        self.deviceDriverName = 'INTEL_IWL5300'
        self.deviceSubCarriers = 30
        self.deviceSymbolData = []
        self.streamFile = None
        self._c = 3e8
        self._f2 = 2.412e9
        self._f5 = 5.825e9
        self._fs = 312.5e3
        self._d = 0.03
        self._windowLen = 100

    @staticmethod
    def to_db(aValue):
        return 20*np.log10(aValue)

    @staticmethod
    def inverse_db(aValue):
        return 10**float(aValue/10.0)

    @staticmethod
    def _twos_comp(val, bitwidth = 8):
        mask = 2**(bitwidth-1)
        return -(val & mask) + (val & ~mask)

    def open(self, filePath):
        self.__parse_symbol_file(filePath)

    def open_stream(self, source):
        if not self.streamFile:
            try:
                self.streamFile = open(source, "rb")
                st_results = os.stat(source)
                st_size = st_results[6]
                self.streamFile.seek(st_size)
            except IOError as e:
                return False

        self.where = self.streamFile.tell()
        line = self.streamFile.readline()

        while not line:
            self.streamFile.seek(self.where)
            time.sleep(.01)

            self.where = self.streamFile.tell()
            line = self.streamFile.readline()
        return self.__parse_symbol_file(StringIO(line))

    def live_stream(self, port):
        app = QtWidgets.QApplication([])

        mainWindow = QtWidgets.QMainWindow()
        mainWindow.resize(800, 600)
        mainWindow.setWindowTitle('Real-Time packet flow')

        centralWidget = QtWidgets.QWidget()
        mainWindow.setCentralWidget(centralWidget)

        lay = QtWidgets.QVBoxLayout()
        centralWidget.setLayout(lay)

        specWid1 = pg.PlotWidget(title='phase diff')
        specItem1 = specWid1.getPlotItem()

        specWid2 = pg.PlotWidget(title='rssi')
        # specWid2 = pg.PlotWidget(name='spectrum')
        specItem2 = specWid2.getPlotItem()

        lay.addWidget(specWid1)
        lay.addWidget(specWid2)

        specWid1.setYRange(-4, 4)
        specWid2.setYRange(0, 50)

        mainWindow.show()
        while True:
            host = ''
            bufferSize = 1024
            ADDR = (host, port)

            t = socket(AF_INET, SOCK_STREAM)
            t.bind(ADDR)
            t.listen(3)

            broken_perm = 0
            triangle = [1, 3, 6]
            count = 0

            phasediff21 = []
            phasediff32 = []
            phasediff31 = []
            rssi_alist = []
            rssi_blist = []
            rssi_clist = []
            degList1 = []
            degList2 = [[1]]
            X_plot = np.arange(-90, 91, 1)[:, np.newaxis]

            print 'Waiting for connection...'
            tClient, addr = t.accept()
            print 'Connected from :', addr

            while True:
                try:
                    tClient.settimeout(15)
                    size = struct.unpack('>H', tClient.recv(2))[0]
                except timeout:
                    print 'Timeout, please restart the client and connect again.'
                    tClient, addr = t.accept()
                    print 'Connected from :', addr
                    continue
                code = struct.unpack('B', tClient.recv(1))[0]

                if code == 187:
                    bytes = tClient.recv(size-1)

                    if len(bytes) != (size - 1):
                        print len(bytes)
                        break
                        # tClient.close()
                        # return False
                elif size <= bufferSize:
                    tClient.recv(size - 1)
                    continue
                else:
                    continue
                if code == 187:
                    count = count + 1

                    ret = self.__unpack_symbol(bytes)
                    csi = ret['csi']
                    # rssi_a = ret['rssi_a']
                    # rssi_b = ret['rssi_b']
                    # rssi_c = ret['rssi_c']
                    perm = ret['perm']
                    nrx = ret['Nrx']
                    ntx = ret['Ntx']
                    if nrx == 3:
                        if sum(perm) != triangle[nrx - 1]:
                            if broken_perm == 0:
                                broken_perm = 1
                                print("WARN ONCE: Found CSI with NRX=", nrx, " and invalid perm=", perm)
                        else:
                            colByTransmitter = lambda t, csi: [tuple(csi[0:t])] + colByTransmitter(t, csi[t:]) if len(
                                csi) != 0 else []
                            csi = colByTransmitter(ntx, csi)

                            srcList = lambda r, t, csi, cnt: [csi[cnt::r]] + srcList(r, t, csi,
                                                                                     cnt + 1) if r != cnt else []
                            csi = srcList(nrx, ntx, csi, 0)

                            deTuple = lambda t: list(map(deTuple, t)) if isinstance(t, (list, tuple)) else t
                            csi = deTuple(csi)

                            csiMatrix = np.zeros((30 * ntx, nrx), complex)

                            for n in range(len(csi)):
                                csi[n] = list(func.reduce(lambda x, y: x + y, csi[n]))

                            for n in range(nrx):
                                csiMatrix[:, perm[n] - 1] = csi[n]


                            """this is the plot of ifft of csi to get the cir of a packet"""
                            # cir = ifft(csi_matrix.T, 128)
                            # print np.abs(cir)
                            # specItem.plot(np.abs(cir)[0, :], clear=True)
                            # app.processEvents()
                            # if not mainWindow.isActiveWindow():
                            #     break
                            """this is the plot of raw csi phase"""
                            # specItem2.clearPlots()
                            # specItem2.plot(np.unwrap(np.angle(csiMatrix[:, 0])), pen=(255, 0, 0))
                            # specItem2.plot(np.unwrap(np.angle(csiMatrix[:, 1])), pen=(0, 255, 0))
                            # specItem2.plot(np.unwrap(np.angle(csiMatrix[:, 2])), pen=(0, 0, 255))
                            # app.processEvents()
                            # if not mainWindow.isActiveWindow():
                            #     break
                            """this is the plot of raw csi phase difference"""
                            self.__csi_phase_diff(csiMatrix, phasediff21, phasediff32, phasediff31)
                            specItem1.clearPlots()
                            specItem1.plot(phasediff21, pen=(255, 0, 0))
                            specItem1.plot(phasediff32, pen=(0, 255, 0))
                            # specItem1.plot(list(map(lambda x: x[0]-x[1], zip(phasediff31, phasediff21))), pen=(0, 0, 255))
                            specItem1.plot(phasediff31, pen=(0, 0, 255))
                            # app.processEvents()
                            # if not mainWindow.isActiveWindow():
                            #     break

                            """this is the plot of rssi"""
                            self.__get_rssi(csiMatrix, rssi_alist, rssi_blist, rssi_clist)
                            specItem2.clearPlots()
                            specItem2.plot(rssi_alist, pen=(255, 0, 0))
                            specItem2.plot(rssi_blist, pen=(0, 255, 0))
                            # specItem1.plot(list(map(lambda x: x[0]-x[1], zip(phasediff31, phasediff21))), pen=(0, 0, 255))
                            specItem2.plot(rssi_clist, pen=(0, 0, 255))

                            # csi_matrix = np.matrix(csiMatrix)
                            # self.csi_facto_aoa(csi_matrix, degList1, self._windowLen)
                            # """this is the plot of kde of the aoa of the past 100 packet"""
                            # k = kde.KernelDensity(kernel='gaussian', bandwidth=1).fit(degList1)
                            # log_dens = k.score_samples(X_plot)
                            # specItem2.plot(X_plot[:, 0], np.exp(log_dens), clear=True)
                            app.processEvents()
                            if not mainWindow.isActiveWindow():
                                break

                            """foreground dtw print"""
                            # if len(degList1) == self._windowLen:
                            #     dtw = fastDTW(degList1, degList2)
                            #     val, path = dtw.dtw()
                            #     print val
                            #     degList2 = copy.copy(degList1)

            tClient.close()
            t.close()
            break

    @staticmethod
    def __csi_extend_57(csi_matrix):
        amp30 = abs(csi_matrix)
        fi30 = np.angle(csi_matrix)
        fii30 = copy.copy(fi30)
        for rx in range(0, 3):
            offset = 0
            for i in range(1, len(fi30[rx, :])):
                if abs(fi30[rx, i] - fi30[rx, i - 1]) > pi - 1:
                    if fi30[rx, i - 1] > fi30[rx, i]:
                        offset += 2 * pi
                    else:
                        offset -= 2 * pi
                fii30[rx, i] = fi30[rx, i] + offset
        amp57 = np.zeros([3, 57], dtype='float64')
        fii57 = np.zeros([3, 57], dtype='float64')
        csi57 = np.zeros([3, 57], dtype='complex64')
        scidx20M = range(0, 27, 2) + range(27, 56, 2) + [56]
        c_scidx20M = range(1, 27, 2) + range(28, 55, 2)
        amp57[:, scidx20M] = amp30
        fii57[:, scidx20M] = fii30
        fii57[:, c_scidx20M] = (fii57[:, [i - 1 for i in c_scidx20M]] + fii57[:, [i + 1 for i in c_scidx20M]]) / 2
        # return fii57
        amp57[:, c_scidx20M] = 10 * np.log10((10 ** (amp57[:, [i - 1 for i in c_scidx20M]] / 10) + 10 ** (
        amp57[:, [i + 1 for i in c_scidx20M]] / 10)) / 2)
        # return amp57
        for rx in range(0, 3):
            for i in range(0, 57):
                csi57[rx, i] = cmath.rect(amp57[rx, i], fii57[rx, i])
        return csi57

    @staticmethod
    def __smooth_csi(e_csi):
        smoothed_csi = np.zeros([30, 86], dtype='complex64')
        for i in range(0, 15):
            smoothed_csi[i, :] = np.concatenate((e_csi[0, i:(i + 43)], e_csi[1, i:(i + 43)]), axis=0)
        for i in range(15, 30):
            smoothed_csi[i, :] = np.concatenate((e_csi[1, (i - 15):(i + 28)], e_csi[2, (i - 15):(i + 28)]), axis=0)
        return smoothed_csi

    def __get_rssi(self, csiMatrix, rssi_alist, rssi_blist, rssi_clist):
        if len(rssi_alist) == self._windowLen:
            rssi_alist.pop(0)
            rssi_blist.pop(0)
            rssi_clist.pop(0)
        rssi_alist.append(self.to_db(np.abs(csiMatrix[0, 0])))
        rssi_blist.append(self.to_db(np.abs(csiMatrix[0, 1])))
        rssi_clist.append(self.to_db(np.abs(csiMatrix[0, 2])))

    def __csi_phase_diff(self, csiMatrix, phasediff21, phasediff32, phasediff31):
        if len(phasediff21) == self._windowLen:
            phasediff21.pop(0)
            phasediff32.pop(0)
            phasediff31.pop(0)
        phasediff21.append(self.__unwrap(np.angle(csiMatrix[0, 1]) - np.angle(csiMatrix[0, 0])))
        phasediff32.append(self.__unwrap(np.angle(csiMatrix[0, 2]) - np.angle(csiMatrix[0, 1])))
        phasediff31.append(self.__unwrap(np.angle(csiMatrix[0, 2]) - np.angle(csiMatrix[0, 0])))
        # print np.angle(csiMatrix[0, 1]) - np.angle(csiMatrix[0, 0])
        # print self.__unwrap(np.angle(csiMatrix[0, 1]) - np.angle(csiMatrix[0, 0]))

    @staticmethod
    def __unwrap(phase):
        if phase > pi:
            phase -= 2*pi
        elif phase < -pi:
            phase += 2*pi
        return phase

    def csi_facto_aoa(self, csi_matrix, degList, windowLen):
        if len(degList) == windowLen:
            del degList[:(windowLen/5)]
        csi = csi_matrix.T
        removedcsi = np.matrix(self.__csi_remove_multipath(csi))  #remove multipath
        # smoothed_csi = np.matrix(self.__smooth_csi(self.__csi_extend_57(removedcsi)))
        deg = [self.__csi_find_aoa(csi)]
        degList.append(deg)

    def __csi_find_aoa(self, csi):
        MUSIC_S = csi * csi.H
        value, vector = np.linalg.eigh(MUSIC_S)
        # noiseVectorIdx = [i for i in range(0, len(value)) if value[i] < np.max(value) * 1e-4]
        noiseVectorIdx = [0]
        En = np.matrix(vector[:, noiseVectorIdx])
        degrees = np.arange(-90, 91, 1)  # the resolution of degree
        tmp = np.matrix([0, 1, 2]).T  # only use 3 antennas
        # tmp = np.matrix([0, 1, 2]).T  # use 3an
        # SP = np.zeros([len(degrees), 1], dtype='float64')
        A = np.exp(-1j * 2 * pi * self._d * tmp * np.sin(degrees * pi / 180) * self._f5 / self._c)
        Ena = En.H * A
        den = np.diag(Ena.H * Ena)
        p = np.abs(1 / den).tolist()
        return degrees[p.index(max(p))]

    @staticmethod
    def __csi_remove_multipath(csi):
        cir = ifft(csi)
        # print np.abs(cir)
        n=5
        # numZeros1 = np.zeros([3, n-1], dtype='complex64')
        # numZeros2 = np.zeros([3, (30-n)], dtype='complex64')
        # removedCSI = fft(np.hstack([numZeros1, cir[:, (n-1):n], numZeros2]))
        numZeros2 = np.zeros([3, (30 - n)], dtype='complex64')
        removedCSI = fft(np.hstack([cir[:, 0:n], numZeros2]))
        return removedCSI

    def __parse_symbol_file(self, aFilePath):

        # make sure we only open when not streaming
        if isinstance(aFilePath, str):
            try:
                f = open(aFilePath, "rb")
            except IOError as e:
                print(e.errno)
                print(e)
                return False

        # seek end of file
        f.seek(0, os.SEEK_END)
        fileLength = f.tell()
        f.seek(0, os.SEEK_SET)

        cur = 0
        count = 0
        broken_perm = 0
        triangle = [1, 3, 6]
        self.deviceSymbolData = []

        while cur < (fileLength - 1):
            cur = cur + 3
            size = struct.unpack('>H', f.read(2))[0]  # Big Endian, uint16, 2 bytes for the payload size field
            code = struct.unpack('B', f.read(1))[0]  # uint8, only one byte for the code field

            if code == 187:
                bytes = f.read(size - 1)
                cur = cur + size - 1

                if len(bytes) != (size - 1):
                    f.close()
                    return False
            else:  # skip the rest
                f.seek(size - 1)
                cur = cur + size - 1
                pass

            if code == 187:
                count = count + 1

                ret = self.__unpack_symbol(bytes)
                # print ret
                perm = ret['perm']
                nrx = ret['Nrx']
                ntx = ret['Ntx']
                self.deviceSymbolData.append(ret)

                if sum(perm) != triangle[nrx - 1]:
                    if broken_perm == 0:
                        broken_perm = 1
                        print("WARN ONCE: Found CSI with NRX=", nrx, " and invalid perm=", perm)
                else:
                    csi = self.deviceSymbolData[count - 1]['csi']

                    colByTransmitter = lambda t, csi: [tuple(csi[0:t])] + colByTransmitter(t, csi[t:]) if len(
                        csi) != 0 else []
                    csi = colByTransmitter(ntx, csi)

                    srcList = lambda r, t, csi, cnt: [csi[cnt::r]] + srcList(r, t, csi, cnt + 1) if r != cnt else []
                    csi = srcList(nrx, ntx, csi, 0)

                    deTuple = lambda t: list(map(deTuple, t)) if isinstance(t, (list, tuple)) else t
                    csi = deTuple(csi)

                    # create numpy matrix
                    mtrx = np.zeros((30 * ntx, nrx), complex)

                    for n in range(len(csi)):
                        csi[n] = list(func.reduce(lambda x, y: x + y, csi[n]))

                    for n in range(nrx):
                        mtrx[:, perm[n] - 1] = csi[n]

                    mtrx = mtrx.reshape((30, ntx, nrx))
                    self.deviceSymbolData[count - 1]['csi'] = mtrx
                    self.deviceSymbolData[count - 1]['csi'] = self.__scale_csi_to_ref(self.deviceSymbolData[count - 1])

    def __unpack_symbol(self, bytes):
        # print bytes[0], ' ', bytes[1], ' ', bytes[2], ' ', bytes[3]
        timestamp_low = int(bytes[0].encode('hex'), 16) + (int(bytes[1].encode('hex'), 16) << 8) + (int(bytes[2].encode('hex'), 16) << 16) + (int(bytes[3].encode('hex'), 16) << 24)
        nrx = int(bytes[8].encode('hex'), 16)
        ntx = int(bytes[9].encode('hex'), 16)
        rssi_a = int(bytes[10].encode('hex'), 16)
        rssi_b = int(bytes[11].encode('hex'), 16)
        rssi_c = int(bytes[12].encode('hex'), 16)
        # Unsigned, two's complement
        noise = -(int(bytes[13].encode('hex'), 16) & 2 ** (8 - 1)) + (int(bytes[13].encode('hex'), 16) & ~2 ** (8 - 1))
        agc = int(bytes[14].encode('hex'), 16)
        antenna_sel = int(bytes[15].encode('hex'), 16)
        length = int(bytes[16].encode('hex'), 16) + (int(bytes[17].encode('hex'), 16) << 8)
        fake_rate_n_flags = int(bytes[18].encode('hex'), 16) + (int(bytes[19].encode('hex'), 16) << 8)
        calc_len = int((30 * (nrx * ntx * 8 * 2 + 3) + 7) / 8)

        csi = []
        index = 0
        # Starting offset for payload
        ptr = 20
        # Premutation array
        perm = []

        if calc_len != length:
            print("Lengths don't match!")
            return False

        # Compute CSI values
        for i in range(30):
            index += 3
            remainder = index % 8
            for j in range(nrx * ntx):
                # Only care about 1 byte
                tmp1 = (int(bytes[ptr + int(index / 8)].encode('hex'), 16) >> remainder) & 255
                # Python bit width isn't fixed!
                tmp2 = (int(bytes[ptr + int(index / 8) + 1].encode('hex'), 16) << (8 - remainder)) & 255
                re = self._twos_comp(tmp1 | tmp2)

                tmp1 = (int(bytes[ptr + int(index / 8) + 1].encode('hex'), 16) >> remainder) & 255
                tmp2 = (int(bytes[ptr + int(index / 8) + 2].encode('hex'), 16) << (8 - remainder)) & 255
                img = self._twos_comp(tmp1 | tmp2)

                csi.append(complex(re, img))
                index += 16

        # Compute antenna premutation  array
        perm.append(((antenna_sel) & 0x3) + 1)
        perm.append(((antenna_sel >> 2) & 0x3) + 1)
        perm.append(((antenna_sel >> 4) & 0x3) + 1)

        csi_struct = {"timestamp_low": timestamp_low, "bfee_count": 0, "Nrx": nrx, "Ntx": ntx, "rssi_a": rssi_a,
                      "rssi_b": rssi_b, "rssi_c": rssi_c, "noise": noise, "agc": agc, "perm": perm,
                      "rate": fake_rate_n_flags, "csi": csi}

        return csi_struct

    def __convert_to_total_rss(self, symbol):
        rssi_mag = 0

        if symbol['rssi_a'] != 0:
            rssi_mag += self.inverse_db(symbol['rssi_a'])
        if symbol['rssi_b'] != 0:
            rssi_mag += self.inverse_db(symbol['rssi_b'])
        if symbol['rssi_c'] != 0:
            rssi_mag += self.inverse_db(symbol['rssi_c'])

        return 10*np.log10(rssi_mag) - 44 - symbol['agc']

    def __scale_csi_to_ref(self, parsedData):
        csi = parsedData['csi']
        csi_sqr = csi * np.conjugate(csi)
        csi_pwr = np.sum(csi_sqr)
        rssi_pwr = self.inverse_db(self.__convert_to_total_rss(parsedData))
        scale = rssi_pwr / (csi_pwr / 30)

        if parsedData['noise'] == -127:
            noise_db = -92
        else:
            noise_db = parsedData['noise']

        thermal_noise_pwr = self.inverse_db(noise_db)
        quant_error_pwr = scale * (parsedData['Nrx'] * parsedData['Ntx'])
        total_noise_pwr = thermal_noise_pwr + quant_error_pwr

        scaled = csi * np.sqrt(scale / total_noise_pwr)

        if parsedData['Ntx'] == 2:
            parsedData['csi'] = scaled * np.sqrt(2)
        elif parsedData['Ntx'] == 3:
            parsedData['csi'] = scaled * np.sqrt(self.inverse_db(4.5))
        else:
            parsedData['csi'] = scaled

        return parsedData['csi']


# The API implementation for Intel
class IntelDeviceService(ExtractorService):
    def __init__(self, aIntelDevice):
        super(IntelDeviceService, self).__init__()
        self._mDriver = aIntelDevice
        self.streamPtr = 0
        self.streamMtrx = None

    def open(self, filePath):
        self._mDriver.open(filePath)

    def reset_stream(self):
        self.streamPtr = 0
        self.streamMtrx = None

    def open_stream(self, aStreamLocation, mode = 0, sampleRate = 1, bufferSize = 1):

        # mode = 0 --> REPLAY, mode = 1 --> LIVE streaming
        if mode == 1:
            return self._mDriver.open_stream(aStreamLocation)
        elif mode == 0:
            if self.streamMtrx is None:
                self._mDriver.open(aStreamLocation)
                self.streamMtrx = self.convert_to_csi_matrix()
            if self.streamPtr >= self.streamMtrx.shape[3]:
                return None

            start = self.streamPtr
            self.streamPtr = self.streamPtr + bufferSize
            time.sleep(1/sampleRate)
            return self.streamMtrx[:, :, :,  start:self.streamPtr:]
        else:
            raise NotImplementedError

    def live_stream(self, port):
        self._mDriver.live_stream(port)

    def device_driver_type(self):
        return self._mDriver.deviceDriverName

    def device_subcarriers(self):
        return self._mDriver.deviceSubCarriers

    def get_receiver_count(self, idx):
        return self._mDriver.deviceSymbolData[idx]['Nrx']
        # Count must be the same for all symbols

    def get_transmitter_count(self, idx):
        return self._mDriver.deviceSymbolData[idx]['Ntx']

    def get_full_RSSI(self):
        rssi_a = []
        rssi_b = []
        rssi_c = []

        for i in self._mDriver.deviceSymbolData:
            rssi_a.append(i['rssi_a'])
            rssi_b.append(i['rssi_b'])
            rssi_c.append(i['rssi_c'])

        return [['rssi_a',rssi_a], ['rssi_b', rssi_b], ['rssi_c', rssi_c]]

    def get_scaled_RSSI(self):
        temp = []
        for i in self._mDriver.deviceSymbolData:
            scaled = self.__convert_to_total_rss(i)
            temp.append(scaled)
        return temp

    def get_symbol_count(self):
        return len(self._mDriver.deviceSymbolData)

    def get_symbol_agc(self, symbolIndex):
        return self._mDriver.deviceSymbolData[symbolIndex]['agc']

    def get_symbol_noise(self, symbolIndex):
        return self._mDriver.deviceSymbolData[symbolIndex]['noise']

    def get_symbol_rate(self, symbolIndex):
        return self._mDriver.deviceSymbolData[symbolIndex]['rate']

    # def convert_to_csi_matrix(self):
    #     mtrx = csiMatrix(self.get_transmitter_count(), self.get_receiver_count(), self.device_subcarriers(), self.get_symbol_count())
    #
    #     for syms in range(self.get_symbol_count()):
    #         tmp = self._mDriver.deviceSymbolData[syms]
    #         print tmp[0]
    #         for t in range(tmp['Ntx']):
    #             for r in range(tmp['Nrx']):
    #                 mtrx[t, r, :, syms] = tmp['csi'][:, t, r]
    #     return mtrx

    def get_csi_symbols(self):
        tempList = []
        for sym in self._mDriver.deviceSymbolData:
            tempList.append(sym['csi'])
        return tempList

    def __convert_to_total_rss(self, symbol):
        rssi_mag = 0

        if symbol['rssi_a'] != 0:
            rssi_mag += self.InvDb(symbol['rssi_a'])
        if symbol['rssi_b'] != 0:
            rssi_mag += self.InvDb(symbol['rssi_b'])
        if symbol['rssi_c'] != 0:
            rssi_mag += self.InvDb(symbol['rssi_c'])

        return 10*np.log10(rssi_mag) - 44 - symbol['agc']

    def InvDb(self, symbol):
        return 10**(symbol/10.0)
