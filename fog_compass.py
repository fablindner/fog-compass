#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib
# Set the backend for matplotlib.
matplotlib.use("TkAgg")
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import Tkinter

from obspy.core.utcdatetime import UTCDateTime
from obspy import *
from obspy.seedlink.seedlinkexception import SeedLinkException
from obspy.seedlink.slclient import SLClient
from obspy.seedlink.slpacket import SLPacket
import numpy as np
import sys
import traceback
import logging
import math as m
import matplotlib.pyplot as plt
from copy import deepcopy
import threading


# default logger
logger = logging.getLogger('obspy.seedlink')


class MySLClient(SLClient):
    """
    A custom SeedLink client.
    """
    def __init__(self, lock, shared, *args, **kwargs):
        """
        Creates a new instance of SLClient accepting a realtime trace handler.
        """
        self.lock = lock
        self.history = shared['history']
        self.pitch = shared['pitch']
        self.roll = shared['roll']
        self.stream = Stream()
        super(self.__class__, self).__init__(*args, **kwargs)
        self.starttimes_ok = False

    def packetHandler(self, count, slpack):
        """
        Processes each packet received from the SeedLinkConnection.

        This method should be overridden when sub-classing SLClient.

        :type count: int
        :param count:  Packet counter.
        :type slpack: :class:`~obspy.seedlink.SLPacket`
        :param slpack: packet to process.
        :return: Boolean true if connection to SeedLink server should be \
            closed and session terminated, false otherwise.
        """
        # check if not a complete packet
        if slpack is None or (slpack == SLPacket.SLNOPACKET) or \
                (slpack == SLPacket.SLERROR):
            return False

        # get basic packet info
        seqnum = slpack.getSequenceNumber()
        type = slpack.getType()

        # process INFO packets here
        if (type == SLPacket.TYPE_SLINF):
            return False
        if (type == SLPacket.TYPE_SLINFT):
            print "-" * 40
            print "Complete INFO:\n" + self.slconn.getInfoString()
            if self.infolevel is not None:
                return True
            else:
                return False

        # can send an in-line INFO request here
        if (count % 100 == 0):
            infostr = "ID"
            self.slconn.requestInfo(infostr)

        # if here, must be a data blockette
        print "-" * 40
        print self.__class__.__name__ + ": packet seqnum:",
        print str(seqnum) + ": blockette type: " + str(type)
# import ipdb;ipdb.set_trace()
        # process packet data
        
        trace = slpack.getTrace()
        #print trace
        if trace is not None:
            self.stream += trace
            self.stream.merge()
        
        if len(self.stream) != 6:
            return False
        
        if not self.starttimes_ok:
            self.stream.trim(starttime=max([tr.stats.starttime for tr in self.stream]))
            self.starttimes_ok = True
            
           
            #print self.__class__.__name__ + ": blockette contains a trace: ",
            #print trace.id, trace.stats['starttime'],
            #print " dt:" + str(1.0 / trace.stats['sampling_rate']),
            #print " npts:" + str(trace.stats['npts']),
            #print " sampletype:" + str(trace.stats['sampletype']),
            #print " dataquality:" + str(trace.stats['dataquality'])
            
            # Custom: append packet data to RtTrace
            #g_o_check = True    # raises Error on gap or overlap
            #g_o_check = False   # clears RTTrace memory on gap or overlap
            
        t1 = self.stream[0].stats.starttime
        t2 = t1 + 2

        st_work = self.stream.slice(starttime=t1, endtime=t2)
        #print st_work
        for tr in st_work:
            if abs(tr.stats.endtime - t2) > tr.stats.delta:
                return False
                   
        # average acceleration
        ave_accx = np.mean((st_work.select(id="BW.FOG1..EN1"))) * 4.e-6
        ave_accy = np.mean((st_work.select(id="BW.FOG1..EN2"))) * 4.e-6
        ave_accz = np.mean((st_work.select(id="BW.FOG1..EN3"))) * 4.e-6
        
        # average rotation            
        ave_rotx = np.mean((st_work.select(id="BW.FOG1..EJ1"))) * 2.e-7
        ave_roty = np.mean((st_work.select(id="BW.FOG1..EJ2"))) * 2.e-7
        ave_rotz = np.mean((st_work.select(id="BW.FOG1..EJ3"))) * 2.e-7
                    
        # tilt
        # theta
        theta = m.asin(ave_accx/9.81)
        
        # phi
        phi = m.asin(ave_accy/(9.81*m.cos(theta)))
        
        # rotation vector
        rot_vec = np.array([ave_rotx, ave_roty, ave_rotz])
        
        # rotation matrix
        R=np.array([[m.cos(theta)               ,       0       ,                 -m.sin(theta)],
                    [m.sin(-phi)*m.sin(theta)   , m.cos(-phi)   ,      m.sin(-phi)*m.cos(theta)],
                    [m.cos(-phi)*m.sin(theta)   , -m.sin(-phi)  ,      m.cos(-phi)*m.cos(theta)]])
        
        # rotation vector corrected for tilt            
        rot_vec_cor = np.dot(R.transpose() , rot_vec)
        
        # azimuth vs north
        azimuth = -m.atan2(rot_vec_cor[1],rot_vec_cor[0]) *  180/m.pi
        
        # calculate pitch and roll
        pitch = m.degrees(theta)
        roll = m.degrees(phi)

        with self.lock:
            self.history.append(azimuth)
            self.pitch.append(pitch)
            self.roll.append(roll)

        self.stream.trim(starttime=t2)

        return False


class AzimuthPlotter(Tkinter.Tk):
    """
    """
    def __init__(self, lock, shared, *args, **kwargs):
        Tkinter.Tk.__init__(self, *args, **kwargs)

        self.lock = lock
        self.history = shared['history']
        self.pitch = shared['pitch']
        self.roll = shared['roll']
        self.update_interval = 500  # in milliseconds

        self.focus_set()
        self._bind_keys()
        ### size and position
        self.geometry("400"+'x'+"400"+'+'+"100"+'+'+"100")
        w, h, pad = self.winfo_screenwidth(), self.winfo_screenheight(), 3
        self._geometry = ("%ix%i+0+0" % (w - pad, h - pad))

        # main figure
        self.fig = Figure()
        self.ax = self.fig.add_subplot(1, 1, 1, polar=True)

        canvas = FigureCanvasTkAgg(self.fig, master=self)
        canvas.show()
        canvas.get_tk_widget().pack(fill=Tkinter.BOTH, expand=1)
        self.canvas = canvas

        self.plot()

    def _quit(self, event):
        event.widget.quit()

    def _bind_keys(self):
        self.bind('<Escape>', self._quit)
        self.bind('q', self._quit)
        self.bind('f', self._toggle_fullscreen)

    def _toggle_fullscreen(self, event):
        g = self.geometry()
        self.geometry(self._geometry)
        self._geometry = g

    def plot(self):
        with self.lock:
            short_history = deepcopy(self.history[-10:])
            pitch = deepcopy(self.pitch[-10:])
            roll = deepcopy(self.roll[-10:])

        try:
            azimuth = short_history[-1]
            pitch = pitch[-1]
            roll = roll[-1]
        # no data to plot yet
        except IndexError:
            self.after(self.update_interval, self.plot)
            return

        std_history = np.std(short_history)

        #print "Mean History %f\n"%(mean_history)
        print "Standart Deviation History %f\n"%(std_history)

        print "Azimuth[deg]: %f\n"%(azimuth)
        print  "Pitch: %f, Roll: %f\n"%(pitch, roll)

        # Plot
        azimuth_rad = np.deg2rad(azimuth)
        std_rad = np.deg2rad(std_history)
        #std = np.deg2rad(1)
        radii = 1
        width = 2 * std_rad
        
        ax = self.ax
        ax.clear()
        ax.set_theta_direction(-1)
        ax.set_theta_zero_location("N")
        ax.annotate("Azimuth: %.1f\n"%(azimuth),
                    xy=(azimuth_rad,radii),
                    xytext=(0.02, 0.11),
                    textcoords='figure fraction')
        ax.annotate("Std Azimuth: %.1f\n"%(std_history),
                    xy=(azimuth_rad,radii),
                    xytext=(0.02, 0.08),
                    textcoords='figure fraction')
        ax.annotate("pitch: %.1f\n"%(pitch),
                    xy=(azimuth_rad,radii),
                    xytext=(0.02, 0.03),
                    textcoords='figure fraction')
        ax.annotate("roll: %f\n"%(roll),
                    xy=(azimuth_rad,radii),
                    xytext=(0.02, 0.0),
                    textcoords='figure fraction')
        bars = ax.bar(azimuth_rad - 0.5 * width, radii, width=width, alpha=0.2)
        ax.plot([azimuth_rad] * 2, [0, 1], color="r", lw=2)
        self.fig.canvas.draw()

        self.after(self.update_interval, self.plot)


def main():
    # initialize objects shared between threads
    lock = threading.Lock()
    shared = {'history': [], 'pitch': [], 'roll': []}

    # create SeedLink client
    try:
        slClient = MySLClient(lock, shared)
        #
        #slClient.slconn.setSLAddress("erde:18000")
        #slClient.multiselect = ("BW_FOG1:EJ1 EJ2 EJ3 EN1 EN2 EN3")
        slClient.slconn.setSLAddress("localhost:18000")
        slClient.multiselect = ("BW_FOG1:EJ1 EJ2 EJ3 EN1 EN2 EN3")
        #
        #slClient.slconn.setSLAddress("discovery.rm.ingv.it:39962")
        #slClient.multiselect = ("IV_MGAB:BHZ")
        #
        #slClient.slconn.setSLAddress("rtserve.iris.washington.edu:18000")
        #slClient.multiselect = ("AT_TTA:BHZ")
        #
        # set a time window from 2 min in the past to 5 sec in the future
        print "SeedLink date-time range:", slClient.begin_time, " -> ",
        print slClient.end_time
        slClient.verbose = 3
        slClient.initialize()
        # start slClient in a thread
        thread = threading.Thread(target=slClient.run)
        thread.setDaemon(True)
        thread.start()
    except SeedLinkException, sle:
        logger.critical(sle)
        traceback.print_exc()
        raise sle
    except Exception, e:
        sys.stderr.write("Error:" + str(e))
        traceback.print_exc()
        raise e

    # start GUI window
    gui = AzimuthPlotter(lock, shared)
    gui.mainloop()

if __name__ == '__main__':
    main()
