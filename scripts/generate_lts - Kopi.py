# Tool for reading DFS0 or KM2 files and creating LTS files from it
# Created by Emil Nielsen
# Contact:
# E-mail: enielsen93@hotmail.com

# Import of modules
import os  # path module
import sys
import numpy as np  # matrix manipulation
import bisect  # math bisection method
import re  # regex for reading KM2
from collections import OrderedDict  # Dictionairy with order for writing .CSV
from distutils.util import strtobool  # Converting string to boolean

import datetime  # time series management
from datetime import datetime as dtnow  # get time of code
import matplotlib.dates as dates  # time series management
import dateutil  # parse string to date

from jinja2 import Environment  # Templating language
from jinja2.loaders import FileSystemLoader  # Templating language

import configparser  # .ini read/write
import inspect  # global variable
local_vars = {}
import rainreader

# Function that reads KM2


# def readKM2(filename):
    # global local_vars
    # # Read KM2 file as string
    # with open(filename, 'r') as km2:
        # km2Str = km2.readlines()

    # # Pre-compile regex search patterns
    # eventstartlineRE = re.compile(r"^1 \d{8}")
    # eventinfoRE = re.compile(
        # r"^1 ?(\d{8}) {0,}(\d{4}) {1,}\d+ {1,}\d+ {1,}(\d+) {1,}([\d\.]+) {1,}(\d+)")
    # gaugeintRE = re.compile("([\d\.]+)")

    # # Definining vectors for event information
    # eventstarttime = []  # The start time of each event
    # gaugetime = []  # The time vector of the rain gauge
    # gaugeint = []  # The intensity vector of the rain gauge in [mu m/s]
    # timedelay = 0
    # eventrejected = False

    # # Read the KM2 line by line
    # for i, line in enumerate(km2Str):
        # # If the line contains information about the event:
        # if eventstartlineRE.search(line):
            # # Split the information into segments
            # eventinfo = eventinfoRE.match(line)
            # # If it's not rejected ( == 2 ), include it
            # # THIS IS NOW DISABLED: It doesn't appear like this feature works
            # # like it's supposed to in the KM2 files
            # if not eventinfo.group(5) == "4":
                # # Get the start time of the event
                # eventstarttime.append(
                    # dates.date2num(
                        # datetime.datetime.strptime(
                            # eventinfo.group(1) +
                            # " " +
                            # eventinfo.group(2),
                            # "%Y%m%d %H%M")))
                # # Remember that the next line will be the first registrered intensity for the event, so the first measurement can be excluded
                # # It's not rejected, so don't reject the following measurements
                # eventrejected = False
                # if timedelay > 0:
                    # gaugeint.extend([0])
                    # gaugetime.extend([gaugetime[-1] + 1. / 60 / 24])
                    # timedelay = 0
            # # If the event is rejected, remember this
            # else:
                # eventrejected = True
        # # If the line does not contain information about the event, it must contain intensities.
        # # If it's not rejected, read the intensities
        # elif not eventrejected:
            # ints = map(float, gaugeintRE.findall(line))
            # # Exclude the first measurement
            # gaugeint.extend(ints)
            # gaugetime.extend((np.arange(0, len(ints), dtype=float) +
                              # timedelay) / 60 / 24 + eventstarttime[-1])
            # timedelay += len(ints)
    # return np.asarray(gaugetime, dtype=float), np.asarray(gaugeint)

# Function that reads DFS0 file


def readDFS0(filename):
    import mikeio
    
    dfs0 = mikeio.dfs0.Dfs0(filename)
    dfs0_read = dfs0.read()
    gaugetime, gaugeint = dates.date2num(dfs0_read.time), dfs0_read.data[0]
    return gaugetime, gaugeint

# Function that writes DFS0 file
def writeDFS0(gaugetime, gaugeint, outfile):
    global local_vars
    gaugetime = gaugetime.tolist()
    gaugeint = gaugeint.tolist()
    
    import mikeio

    dfs0 = mikeio.dfs0.Dfs0()
    dfs0.write(outfile, data = [gaugeint], start_time = dates.num2date(gaugetime[0]),
           items = [mikeio.eum.ItemInfo("Rainfall Intensity", mikeio.eum.EUMType.Rainfall_Intensity, unit = mikeio.eum.EUMUnit.mu_m_per_sec)], datetimes = dates.num2date(gaugetime))
           
    return

# Function that returns the data period for a rain series


def getDataPeriod(parameters, scriptFolder):
    if ".dfs0" in parameters[0].valueAsText:
        gaugetime, _ = readDFS0(parameters[0].valueAsText)
    else:
        km2 = rainreader.KM2(parameters[0].valueAsText)
        gaugetime = km2.gaugetime
    dataperiod = (gaugetime[-1] - gaugetime[0]) / 365
    daterange = dates.num2date(gaugetime[0]).strftime(
        "%d/%m/%Y") + " - " + dates.num2date(gaugetime[-1]+1).strftime("%d/%m/%Y")
    return dataperiod, daterange

# Function that returns the bootstrapped samples of a set of parameters


def bootstrap_resample(X, samples, n=None):
    if n is None:
        n = np.size(X, axis=0)
    X_resample = np.empty((np.size(X, axis=0), 0), dtype=np.float16)
    for i in range(0, samples):
        resample_i = np.floor(
            np.random.rand(n) *
            np.size(
                X,
                axis=0)).astype(int)
        X_resample = np.append(
            X_resample, np.transpose([X[resample_i]]), axis=1)
    return X_resample

# Main function that writes LTS files from rain series file


def writeLTS(parameters, scriptFolder):
    # Open log file
    logFile = open(
        os.path.join(
            scriptFolder,
            'log.txt'),
        'w')
    global local_vars
    # Read old config file and write new config file
    if not __name__ == '__main__':
        configStr = ''
        parametersDict = {}
        configWrite = configparser.ConfigParser(allow_no_value=True)
        configWrite.add_section("ArcGIS input parameters")
        for i, par in enumerate(parameters):
            if par.value is None:
                par.value = 0
            parametersDict[par.name] = str(par.valueAsText)
            if i == 0:
                configWrite.set(
                    "ArcGIS input parameters",
                    "# " + par.displayName)
            else:
                configWrite.set(
                    "ArcGIS input parameters",
                    "\r\n# " + par.displayName)
            configWrite.set("ArcGIS input parameters", str(par.name), str(par.value))
            configStr += r"// " + par.displayName + \
                " = " + str(par.value) + "\n"
        with open(scriptFolder + r"\config.ini", "w") as config_file:
            configWrite.write(config_file)
    else:
        configStr = ""
        parametersDict = parameters

    rain_events = []
    class RainEvent:
        def __init__(self):
            start_i = None
            stop_i = None
            statistics = None
            accumulated_rain = None
            dts = None

        @property
        def reduce_timestep(self):


    if ".dfs0" in parametersDict["input_file"]:
        gaugetime, gaugeint = readDFS0(
            parametersDict["input_file"], scriptFolder)
    else:
        km2 = rainreader.KM2(parametersDict["input_file"])
        gaugetime, gaugeint = km2.gaugetime, km2.gaugeint

    if (float(parametersDict["time_series_duration"].replace(",",".")) == 0):
        # Calculate total data period for time series
        dataperiod = np.float((tminutes[-1] - tminutes[0]) / 60 / 24)
    else:
        dataperiod = float(parametersDict["time_series_duration"]) * 365

    time_aggregate_periods = list(map(int, parametersDict["time_aggregate_periods"].split(';'))) if parametersDict["time_aggregate_enable"] == "true" else []
    merge_period = max(
        time_aggregate_periods + [float(parametersDict["rain_event_merge_duration"])] + [5])
    rain_statistics, event_time = km2.rainStatistics(time_aggregate_periods, merge_period)

    for event_i in range(rain_statistics.shape[0]):
        rain_event = RainEvent()
        rain_event.start_i, rain_event.stop_i = event_time[event_i,:]
        rain_event.statistics = rain_statistics[event_i, :-1]
        rain_event.accumulated_rain = rain_statistics[event_i, -1]
        rain_event.time_aggregate_periods = time_aggregate_periods


    event_reduce_timestep = []
    while j < np.size(tminutes) - 1:
        # End of each event is when there's a dry period of xxx minutes
        jend = np.argmax(tdiff[j:] > mergePeriod) + j
        # Initialize time aggregate set for this event
        RDAgg = np.append(
            RDAgg, np.zeros(
                (1, len(dts) + 1), dtype=np.float), axis=0)
        # Start time of this event
        startj = np.append(startj, np.int(j))
        # Calculate total rain depth for event
        RDAgg[eventidx, -1] = np.sum(gaugeint[j:jend]) / 1000 * 60

        reduce_timestep = False
        if RDAgg[eventidx, -1] > 5:
            reduce_timestep = True
        else:
            for j in range(j, jend):
                for limit, dt in zip([4.5, 9.5, 11.5],[int(7), int(30), int(60)]):
                    mm = np.sum(
                        gaugeint[j:j + bisect.bisect_left((tminutes[j:j + dt] - tminutes[j]), dt)]) / 1000 * 60
                    if mm > limit:
                        reduce_timestep = True
                        break
                if reduce_timestep:
                    break

        event_reduce_timestep.append(reduce_timestep)

        # Loop over all intensities in event
        if strtobool(parametersDict["time_aggregate_enable"]):
            for j in range(j, jend):
                # Loop over all time aggregate periods
                for i, dt in enumerate(dts):
                    # Calculate total rain depth over aggregate period
                    dt = int(dt)
                    mm = np.sum(
                        gaugeint[j:j + bisect.bisect_left((tminutes[j:j + dt] - tminutes[j]), dt)]) / 1000 * 60
                    if (mm > RDAgg[eventidx, i]):
                        RDAgg[eventidx, i] = mm
        else:
            j = jend
        # End time of this event
        endj = np.append(endj, np.int(jend))
        j += 1
        eventidx += 1  # Change event index
    logFile.write(str(dtnow.now()) +
                  ": Initializing time series duration\n")

    if (float(parametersDict["time_series_duration"].replace(",",".")) == 0):
        # Calculate total data period for time series
        dataperiod = np.float((tminutes[-1] - tminutes[0]) / 60 / 24)
    else:
        dataperiod = float(parametersDict["time_series_duration"]) * 365

    # Sort the time aggregates according to rain depth
    sortidx = np.argsort(-RDAgg, axis=0)
    RDAggSort = np.flipud(np.sort(RDAgg, axis=0))

    logFile.write(str(dtnow.now(
    )) + ": Calculating total number of events for time aggregate and total rain depth\n")
    # Calculate number events of time aggregate from return period or just
    # number of events
    if strtobool(parametersDict["time_aggregate_enable"]) == False:
        eventsAgg = 0
        eventsAggByRP = False
    elif not (int(parametersDict["time_aggregate_number_events"]) == 0):
        eventsAgg = float(parametersDict["time_aggregate_number_events"])
        eventsAggByRP = False
    elif not (float(parametersDict["time_aggregate_return_period"].split("+")[0]) == 0):
        eventsAgg = np.int(np.floor(
            dataperiod / (float(parametersDict["time_aggregate_return_period"].split("+")[0]) * 365)))
        if "+" in parametersDict["time_aggregate_return_period"]:
            eventsAgg = eventsAgg + \
                float(parametersDict["time_aggregate_return_period"].split("+")[1])
        eventsAggByRP = True
    else:
        eventsAgg = 0
        eventsAggByRP = False

    # Calculate number events of total rain depth from return period or
    # total rain depth
    if not (float(parametersDict["include_events_total_rain_depth"]) == 0):
        eventsSum = sum(
            RDAgg[:, -1] > float(parametersDict["include_events_total_rain_depth"]))
    elif not (float(parametersDict["include_events_return_period"].split("+")[0]) == 0):
        eventsSum = np.int(np.floor(
            dataperiod / (float(parametersDict["include_events_return_period"].split("+")[0]) * 365)))
        if "+" in parametersDict["include_events_return_period"]:
            eventsSum = eventsSum + \
                float(parametersDict["include_events_return_period"].split("+")[1])
    else:
        eventsSum = 0
    eventsSum = int(eventsSum)

    # Parse date criteria as number

    if parametersDict["date_criteria"] and not parametersDict["date_criteria"] == "0":
        date_criteria = parametersDict["date_criteria"].split(' - ')
        date_criteriaNum = []
        for i, d in enumerate(date_criteria):
            date_criteriaNum.append(
                dates.date2num(
                    dateutil.parser.parse(
                        d,
                        dayfirst=True,
                        default=datetime.datetime(
                            datetime.datetime.now().year,
                            1,
                            1))))
    else:  # If it's empty, just set it to this
        date_criteriaNum = [0, 1e10]

    # Calculate 95% confidence interval for return period of rain event
    logFile.write(str(dtnow.now(
    )) + ": Bootstrapping and calculating 95% confidence interval for return period of rain event\n")

    # Initialize event information vectors
    eventstarttime = eventstoptime = eventstarttimeStr = []
    eventstoptimeStr = []
    durHour = []
    accrain = []
    durTotal = 0
    eventdts = []

    logFile.write(str(dtnow.now()) +
                  ": Assembling selected events to LTS data set by total rain depth\n")
    # Start selecting and including events by total rain depth
    for eventidx in range(0, np.int(eventsSum)):
        # Event start time and stop time for the included event
        eventstarttime = dates.num2date(
            gaugetime[startj[sortidx[eventidx, -1]]] - float(parametersDict["soft_start_time"]) / 60 / 24)
        eventstoptime = dates.num2date(
            gaugetime[endj[sortidx[eventidx, -1]]] + float(parametersDict["soft_stop_time"]) / 60 / 24)
        # If the event starts within the selected dates, include it.
        if (date_criteriaNum[0] < gaugetime[startj[sortidx[eventidx, -1]]]
                < date_criteriaNum[1]):
            # If this event is not already included in the set, include it
            if not (eventstarttime.strftime(
                    '%Y-%m-%d %H:%M:00') in eventstarttimeStr):
                # Write the times to strings
                eventstarttimeStr.append(
                    eventstarttime.strftime('%Y-%m-%d %H:%M:00'))
                eventstoptimeStr.append(
                    eventstoptime.strftime('%Y-%m-%d %H:%M:00'))
                eventdts.append(['total event'])
                # Calculate the duration of the event
                dur = (dates.date2num(eventstoptime) -
                       dates.date2num(eventstarttime)) * 24 * 3600
                durTotal += float(dur) / 3600
                durHour.append(float(dur) / 3600)

                # Calculate the rain depth of the event
                accrain.append("%1.1f" %
                               (RDAgg[sortidx[eventidx, -1], -1]))
            else:
                eventdts[eventstarttimeStr.index(eventstarttime.strftime(
                    '%Y-%m-%d %H:%M:00'))].append('total event')

    logFile.write(str(dtnow.now()) +
                  ": Assembling selected events to LTS data set by time aggregate\n")

    # Start selecting and including events by time aggregate
    eventidx = 0  # Index of event
    # If select by return period is enabled, continue until all events with return period above criteria have been tested
    # If select by number of events is enabled, continue until the number
    # of included events has reached this number
    # Amount of events that should be included for each aggregate period
    eventsAggByTA = np.zeros((len(dts))) + eventsAgg
    # Amount of events concluded for each aggregate period
    eventsAggByTAIncl = np.zeros((len(dts)))
    while ((not eventsAggByRP) * sum(eventsAggByTAIncl)
           ) < np.sum(eventsAggByTA) and eventsAggByRP * eventidx < eventsAgg:
        for i in range(0, np.size(dts)):
            # Event start time and stop time for the included event
            eventstarttime = dates.num2date(
                gaugetime[startj[sortidx[eventidx, i]]] - float(parametersDict["soft_start_time"]) / 60 / 24)
            eventstoptime = dates.num2date(
                gaugetime[endj[sortidx[eventidx, i]]] + float(parametersDict["soft_stop_time"]) / 60 / 24)
            # If the event starts within the selected dates, include it. If
            # there's already enough events for that aggregate period,
            # exclude it.
            if (date_criteriaNum[0] <= gaugetime[startj[sortidx[eventidx, i]]] <=
                    date_criteriaNum[1]) and eventsAggByTAIncl[i] < eventsAggByTA[i]:
                # If this event is not already included in the set, include
                # it
                if not (eventstarttime.strftime(
                        '%Y-%m-%d %H:%M:00') in eventstarttimeStr):
                    # Write the times to strings
                    eventstarttimeStr.append(
                        eventstarttime.strftime('%Y-%m-%d %H:%M:00'))
                    eventstoptimeStr.append(
                        eventstoptime.strftime('%Y-%m-%d %H:%M:00'))
                    eventdts.append(["%d min" % (dts[i])])
                    # Calculate the duration of the event
                    dur = (dates.date2num(eventstoptime) -
                           dates.date2num(eventstarttime)) * 24 * 3600
                    durTotal += float(dur) / 3600
                    durHour.append(float(dur) / 3600)

                    # Calculate the rain depth of the event
                    accrain.append("%1.1f" %
                                   (RDAgg[sortidx[eventidx, i], -1]))
                else:
                    eventdts[eventstarttimeStr.index(eventstarttime.strftime(
                        '%Y-%m-%d %H:%M:00'))].append("%d min" % (dts[i]))
                eventsAggByTAIncl[i] += 1
        eventidx += 1

    # Write strings of duration and data period for insertion into LTS.MJL
    # file
    durTotalDay, remainder = divmod(durTotal, 24)
    durTotalDay = int(durTotalDay)
    _, durTotalHour = divmod(remainder, 60)
    durTotalHour = int(durTotalHour)
    durTotalStr = "%d days, %dh" % (durTotalDay, durTotalHour)
    dataperiodYears, remainder = divmod(dataperiod, 365)
    dataperiodMonths = int(remainder / 365 * 12)
    dataperiodStr = "%d years, %d months" % (
        dataperiodYears, dataperiodMonths)

    dur_time_str = []
    for dur in durHour:
        dur_day, remainder = divmod(dur, 24)
        dur_day = int(dur_day)
        _, dur_hour = divmod(remainder, 60)
        dur_hour = int(dur_hour)
        if dur_day > 0:
            dur_time_str.append("%d days, %d hours" % (dur_day, dur_hour))
        else:
            dur_time_str.append("%d hours" % (dur_hour))

    # Sort rain events by start time
    zipped = sorted(zip(eventstarttimeStr, eventstoptimeStr,
                        durHour, accrain, eventdts, dur_time_str))
    eventstarttimeStr, eventstoptimeStr, durHour, accrain, eventdts, dur_time_str = zip(
        *zipped)

    #####
    # Write LTS file through jinga2 module
    # Jinga2 reads a text file and uses it as a template for creating LTS files
    #####
    logFile.write(str(dtnow.now()) + ": Writing LTS file using Jinga2\n")

    Environment(loader=FileSystemLoader(''))
    eventlist = range(0, np.size(eventstarttime))
    # Load LTS template file
    templateFile = open(
        os.path.join(
            scriptFolder,
            'LTSTemplate.MJL'),
        'r')
    templateFileStr = templateFile.read()

    hotstart_time = []
    getTimeRE = re.compile("(\d+:\d+:\d+)")
    for i in range(len(eventstarttimeStr)):
        hotstart_time.append(getTimeRE.findall(eventstarttimeStr[i])[0])

    # Write LTS file with event information vectors
    fout = open(parametersDict["output_mjl"], 'w+')
    fout.write(
        Environment().from_string(templateFileStr).render(
            inputfile=parametersDict["input_file"],
            eventlist=eventlist,
            simulation_start=eventstarttimeStr,
            simulation_stop=eventstoptimeStr,
            job_number=range(
                1,
                len(eventstarttimeStr) + 1),
            dur_time=dur_time_str,
            jobs=len(eventstarttimeStr),
            total_dur_time=durTotalStr,
            dataperiod=dataperiodStr,
            accumulated_rain=accrain,
            eventdts=eventdts,
            date_criteria=parametersDict["date_criteria"],
            configStr=configStr,
            hotstart_param = parametersDict["hotstart_param"],
            hotstart_name = parametersDict["hotstart_name"],
            hotstart_date = parametersDict["hotstart_date"],
            hotstart_time = hotstart_time,
            event_reduce_timestep = event_reduce_timestep,
            event_increase_timestep = [False] + event_reduce_timestep))
    fout.close()

    # Create csv-file with LTS events
    if not True:
        logFile.write(str(dtnow.now()) +
                      ": Creating csv-file with LTS-events\n")
        csvDict = OrderedDict()
        csvDict["Simulation start"] = eventstarttimeStr
        csvDict["Simulation stop"] = eventstoptimeStr
        csvDict["Accumulated rain [mm]"] = accrain
        with open('LTSJobList.csv', 'w') as csvfile:
            csvfile.write("Job")
            [csvfile.write(',{0}'.format(key))
             for key, _ in csvDict.items()]
            for i in range(0, len(csvDict["Simulation start"])):
                csvfile.write("\n{0}".format(i + 1))
                [csvfile.write(',{0}'.format(value[i]))
                 for _, value in csvDict.items()]

    logFile.write(str(dtnow.now()) + ": Succesful run\n")
    logFile.close()
    local_vars = inspect.currentframe().f_locals
    return

# Second function that combines LTS files to one total LTS file


def combineLTS(parameters, scriptFolder):
    lts_files = str(parameters[0].valueAsText).split(";")
    starttime = []
    stoptime = []
    starttimeRE = re.compile(r"Simulation_start {0,}= {0,}'([^']+)'")
    stoptimeRE = re.compile(r"Simulation_end {0,}= {0,}'([^']+)'")
    hotstartFileRE = re.compile(r"Hotstart_file {0,}= {0,}'([^']+)'")
    hotstartTimeRE = re.compile(r"Hotstart_time {0,}= {0,}'([^']+)'")
    hotstart = "false"
    hotstartFile = []
    hotstartTime = []
    input_files = []
    for lts_file in lts_files:
        with open(lts_file.replace('\'', ''), 'r') as lts_file_open:
            text = lts_file_open.read()
            starttime.extend(starttimeRE.findall(text))
            stoptime.extend(stoptimeRE.findall(text))
            input_files.extend([lts_file] * len(starttimeRE.findall(text)))
            if "Hotstart_file" in text:
                hotstart = "true"
                hotstartFile.extend(hotstartFileRE.findall(text))
                hotstartTime.extend(hotstartTimeRE.findall(text))

    starttime = np.array(starttime)
    stoptime = np.array(stoptime)

    for j in range(100):
        idx,starttime = [np.argsort(starttime), np.sort(starttime)]
        stoptime = stoptime[idx]
        deleteA = []
        deleteB = []
        for i in range(len(starttime)-1):
            if stoptime[i]>=stoptime[i+1]:
                print(stoptime[i],stoptime[i+1])
                deleteA.append(i)
                deleteB.append(i+1)
        starttime = np.delete(starttime,deleteA)
        stoptime = np.delete(stoptime,deleteB)

    idx,starttime = [np.argsort(starttime), np.sort(starttime)]
    stoptime = stoptime[idx]
    
    dur = []
    durTotal = 0
    for i in range(0, len(starttime)):
        dur.append(
            '%d' %
            ((dates.date2num(
                datetime.datetime.strptime(
                    stoptime[i],
                    "%Y-%m-%d %H:%M:00")) -
                dates.date2num(
                datetime.datetime.strptime(
                    starttime[i],
                    "%Y-%m-%d %H:%M:00"))) *
                24))
        durTotal += (
            dates.date2num(
                datetime.datetime.strptime(
                    stoptime[i],
                    "%Y-%m-%d %H:%M:00")) - dates.date2num(
                datetime.datetime.strptime(
                    starttime[i],
                    "%Y-%m-%d %H:%M:00"))) * 24

    durTotalDay, remainder = divmod(durTotal, 24)
    durTotalDay = int(durTotalDay)
    _, durTotalHour = divmod(remainder, 60)
    durTotalHour = int(durTotalHour)
    durTotalStr = "%d days, %dh" % (durTotalDay, durTotalHour)

    job_number = range(1, len(starttime) + 1)

    Environment(loader=FileSystemLoader(''))
    # Load LTS template file
    templateFile = open(
        os.path.join(
            scriptFolder.encode(
                'ascii',
                'ignore'),
            'LTSTemplateCombined.MJL'),
        'r')
    templateFileStr = templateFile.read()

    # Write LTS file with event information vectors
    fout = open(parameters[1].valueAsText, 'w+')
    fout.write(Environment().from_string(templateFileStr).render(inputfile=lts_files,
                                                                 simulation_start=starttime,
                                                                 simulation_stop=stoptime,
                                                                 job_number=job_number,
                                                                 dur_time=dur,
                                                                 total_dur_time=durTotalStr,
                                                                 jobs=job_number[-1],
                                                                 hotstart_param = hotstart,
                                                                 hotstart_name = hotstartFile,
                                                                 hotstart_time = hotstartTime))
    fout.close()
    return


# If this .py-file is run as a stand-alone, write an LTS-file with
# parameters from .ini-file (generated from previous run from ArcGIS)
if __name__ == '__main__':
    # Read config file
    config = configparser.ConfigParser()
    config.read(os.path.dirname(os.path.realpath(__file__)) + "\\config.ini")
    parametersDict = {}
    for option in config.options("ArcGIS input parameters"):
        parametersDict[option] = config.get("ArcGIS input parameters", option)
#
    scriptFolder = os.path.dirname(os.path.realpath(__file__))
#    scriptFolder = r"C:\Users\Eniel\Documents\enielsen93\LTS\Long Term Statistics\scripts\\"
#    # Write LTS using parameters from config file
    writeLTS(parametersDict, scriptFolder)
