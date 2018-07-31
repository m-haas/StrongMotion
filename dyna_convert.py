#Michael 17.10.17 rewritten to function
#!/usr/bin/env python

"""
dyna-convert.py - v0.3
"""

"""
>>>
>>> usage:
>>> dyna-convert.py /path/to/dyna-file /path/to/out-file FORMAT
>>>
>>> example MSEED:
>>> dyna-convert.py /home/dyna/data/IT.LRG..HNE.D.19820321.094400.C.DIS.ASC out.mseed mseed
>>>
>>> example SAC:
>>> dyna-convert.py /home/dyna/data/IT.LRG..HNE.D.19820321.094400.C.DIS.ASC out.sac sac
>>>
>>> allowed FORMATS:
>>> mseed, sac, ecc. (see obspy.write module for further details)
"""

"""
LICENSE:

dyna-convert.py is released under the following BSD-style license:

Copyright (c) 2014, INGV Milano-Pavia
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies,
either expressed or implied, of the FreeBSD Project.
"""

from sys import argv as sys_argv
from StringIO import StringIO
from obspy.core import Stream, Trace, UTCDateTime, Stats
import numpy as np
import re
import string

#filename_in = sys_argv[1]
#filename_out = sys_argv[2]
#format_out = sys_argv[3]

def toUTCDateTime(value):
    try:
        date, time = value.split('_')
    except ValueError:
        date = value

    year = int(date[0:4])
    month = int(date[4:6])
    day = int(date[6:8])

    hour = int(time[0:2])
    mins = int(time[2:4])
    secs = float(time[4:])

    return UTCDateTime(year, month, day, hour, mins) + secs

def strtofloat(sf):
    try:
        x = float(sf)
    except:
        return None
    return x

def strtoint(sf):
    try:
        x = int(sf)
    except:
        return None
    return x

def read_ems_to_stream(filename_in):
    """
    Reads an ESM file and returns an obspy stream
    """
    headers = {}
    data = StringIO()

    # read file
    with open(filename_in, 'rt') as fh:
        for i in xrange(64):
            key, value = fh.readline().strip().split(':', 1)
            headers[key.strip()] = value.strip()

        # create ObsPy stream object
        stream = Stream()
        header = Stats()
        header['dyna'] = {}

        header['network'] = headers['NETWORK']
        header['station'] = headers['STATION_CODE']
        header['location'] = headers['LOCATION']
        header['channel'] = headers['STREAM']
        try:
            header['starttime'] = toUTCDateTime(headers['DATE_TIME_FIRST_SAMPLE_YYYYMMDD_HHMMSS']) # use toUTCDateTime to convert from DYNA format
        except:
            #Michael: fixed
            header['starttime'] = toUTCDateTime(headers['EVENT_DATE_YYYYMMDD']+'_'+headers['EVENT_TIME_HHMMSS'])
            #header['starttime'] = toUTCDateTime('19700101_000000')
        header['sampling_rate'] = 1/float(headers['SAMPLING_INTERVAL_S'])
        header['delta'] = float(headers['SAMPLING_INTERVAL_S'])
        header['npts'] = int(headers['NDATA'])
        header['calib'] = 1 # not in file header

        ##DYNA dict float data
        header['dyna']['EVENT_LATITUDE_DEGREE'] = strtofloat(headers['EVENT_LATITUDE_DEGREE'])
        header['dyna']['EVENT_LONGITUDE_DEGREE'] = strtofloat(headers['EVENT_LONGITUDE_DEGREE'])
        header['dyna']['EVENT_DEPTH_KM'] = strtofloat(headers['EVENT_DEPTH_KM'])
        header['dyna']['HYPOCENTER_REFERENCE'] = headers['HYPOCENTER_REFERENCE']
        header['dyna']['MAGNITUDE_W'] = strtofloat(headers['MAGNITUDE_W'])
        header['dyna']['MAGNITUDE_L'] = strtofloat(headers['MAGNITUDE_L'])
        header['dyna']['STATION_LATITUDE_DEGREE'] = strtofloat(headers['STATION_LATITUDE_DEGREE'])
        header['dyna']['STATION_LONGITUDE_DEGREE'] = strtofloat(headers['STATION_LONGITUDE_DEGREE'])
        header['dyna']['VS30_M_S'] = strtofloat(headers['VS30_M/S'])
        header['dyna']['EPICENTRAL_DISTANCE_KM'] = strtofloat(headers['EPICENTRAL_DISTANCE_KM'])
        header['dyna']['EARTHQUAKE_BACKAZIMUTH_DEGREE'] = strtofloat(headers['EARTHQUAKE_BACKAZIMUTH_DEGREE'])
        header['dyna']['DURATION_S'] = strtofloat(headers['DURATION_S'])
        header['dyna']['INSTRUMENTAL_FREQUENCY_HZ'] = strtofloat(headers['INSTRUMENTAL_FREQUENCY_HZ'])
        header['dyna']['INSTRUMENTAL_DAMPING'] = strtofloat(headers['INSTRUMENTAL_DAMPING'])
        header['dyna']['FULL_SCALE_G'] = strtofloat(headers['FULL_SCALE_G'])

        # data type is acceleration
        if headers['DATA_TYPE'] == "ACCELERATION" \
        or headers['DATA_TYPE'] == "ACCELERATION RESPONSE SPECTRUM":
            header['dyna']['PGA_CM_S_2'] = strtofloat(headers['PGA_CM/S^2'])
            header['dyna']['TIME_PGA_S'] = strtofloat(headers['TIME_PGA_S'])
        # data type is velocity
        if headers['DATA_TYPE'] == "VELOCITY" \
        or headers['DATA_TYPE'] == "PSEUDO-VELOCITY RESPONSE SPECTRUM":
            header['dyna']['PGV_CM_S'] = strtofloat(headers['PGV_CM/S'])
            header['dyna']['TIME_PGV_S'] = strtofloat(headers['TIME_PGV_S'])
        # data type is displacement
        if headers['DATA_TYPE'] == "DISPLACEMENT" \
        or headers['DATA_TYPE'] == "DISPLACEMENT RESPONSE SPECTRUM":
            header['dyna']['PGD_CM'] = strtofloat(headers['PGD_CM'])
            header['dyna']['TIME_PGD_S'] = strtofloat(headers['TIME_PGD_S'])

        header['dyna']['LOW_CUT_FREQUENCY_HZ'] = strtofloat(headers['LOW_CUT_FREQUENCY_HZ'])
        header['dyna']['HIGH_CUT_FREQUENCY_HZ'] = strtofloat(headers['HIGH_CUT_FREQUENCY_HZ'])

        ##DYNA dict int data
        header['dyna']['STATION_ELEVATION_M'] = strtoint(headers['STATION_ELEVATION_M'])
        header['dyna']['SENSOR_DEPTH_M'] = strtoint(headers['SENSOR_DEPTH_M'])
        header['dyna']['N_BIT_DIGITAL_CONVERTER'] =  strtoint(headers['N_BIT_DIGITAL_CONVERTER'])
        header['dyna']['FILTER_ORDER'] = strtoint(headers['FILTER_ORDER'])

        ##DYNA dict string data
        header['dyna']['EVENT_NAME'] = headers['EVENT_NAME']
        header['dyna']['EVENT_ID'] = headers['EVENT_ID']
        header['dyna']['EVENT_DATE_YYYYMMDD'] = headers['EVENT_DATE_YYYYMMDD']
        header['dyna']['EVENT_TIME_HHMMSS'] = headers['EVENT_TIME_HHMMSS']
        header['dyna']['MAGNITUDE_W_REFERENCE'] = headers['MAGNITUDE_W_REFERENCE']
        header['dyna']['MAGNITUDE_L_REFERENCE'] = headers['MAGNITUDE_L_REFERENCE']
        header['dyna']['FOCAL_MECHANISM'] = headers['FOCAL_MECHANISM']
        header['dyna']['STATION_NAME'] = headers['STATION_NAME']
        header['dyna']['SITE_CLASSIFICATION_EC8'] = headers['SITE_CLASSIFICATION_EC8']
        header['dyna']['MORPHOLOGIC_CLASSIFICATION'] = headers['MORPHOLOGIC_CLASSIFICATION']
        header['dyna']['DATE_TIME_FIRST_SAMPLE_PRECISION'] = headers['DATE_TIME_FIRST_SAMPLE_PRECISION']
        header['dyna']['UNITS'] = headers['UNITS']
        header['dyna']['INSTRUMENT'] = headers['INSTRUMENT']
        header['dyna']['INSTRUMENT_ANALOG_DIGITAL'] = headers['INSTRUMENT_ANALOG/DIGITAL']
        header['dyna']['BASELINE_CORRECTION'] = headers['BASELINE_CORRECTION']
        header['dyna']['FILTER_TYPE'] = headers['FILTER_TYPE']
        header['dyna']['LATE_NORMAL_TRIGGERED'] = headers['LATE/NORMAL_TRIGGERED']
        header['dyna']['HEADER_FORMAT'] = headers['HEADER_FORMAT']
        header['dyna']['DATABASE_VERSION'] = headers['DATABASE_VERSION']
        header['dyna']['DATA_TYPE'] = headers['DATA_TYPE']
        header['dyna']['PROCESSING'] = headers['PROCESSING']
        header['dyna']['DATA_LICENSE'] = headers['DATA_LICENSE']
        header['dyna']['DATA_TIMESTAMP_YYYYMMDD_HHMMSS'] = headers['DATA_TIMESTAMP_YYYYMMDD_HHMMSS']
        header['dyna']['DATA_CITATION'] = headers['DATA_CITATION']
        header['dyna']['DATA_CREATOR'] = headers['DATA_CREATOR']
        header['dyna']['ORIGINAL_DATA_MEDIATOR_CITATION'] = headers['ORIGINAL_DATA_MEDIATOR_CITATION']
        header['dyna']['ORIGINAL_DATA_MEDIATOR'] = headers['ORIGINAL_DATA_MEDIATOR']
        header['dyna']['ORIGINAL_DATA_CREATOR_CITATION'] = headers['ORIGINAL_DATA_CREATOR_CITATION']
        header['dyna']['ORIGINAL_DATA_CREATOR'] = headers['ORIGINAL_DATA_CREATOR']
        header['dyna']['USER1'] = headers['USER1']
        header['dyna']['USER2'] = headers['USER2']
        header['dyna']['USER3'] = headers['USER3']
        header['dyna']['USER4'] = headers['USER4']
        header['dyna']['USER5'] = headers['USER5']

        # read data
        data = np.loadtxt(fh, dtype='float32')
        if headers['DATA_TYPE'][-8:] == "SPECTRUM":
            data_1 = np.array([], dtype=np.float32)
            data_2 = np.array([], dtype=np.float32)
            for j in xrange(len(data)):
                for i in xrange(2):
                    if i == 0:
                        data_1 = np.append(data_1,data[j][i])
                    elif i == 1:
                        data_2 = np.append(data_2,data[j][i])
            stream.append(Trace(data=data_1, header=header))
            stream.append(Trace(data=data_2, header=header))
        else:
            stream.append(Trace(data=data, header=header))

    return stream

#stream.write(filename_out,format=format_out)
