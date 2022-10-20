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
            self.start_i = None
            self.stop_i = None
            self.statistics = None
            self.accumulated_rain = None
            self.dts = None
            self.include = False
            self.reduce_timestep = False
            self.selected_because_of = []

        @property
        def event_start_time(self):
            return dates.num2date(
                gaugetime[self.start_i] - float(parametersDict["soft_start_time"]) / 60 / 24)

        @property
        def event_stop_time(self):
            return dates.num2date(
                gaugetime[min(self.stop_i, len(gaugetime)-1)] + float(parametersDict["soft_stop_time"]) / 60 / 24)

        @property
        def duration(self):
            return gaugetime[min(self.stop_i, len(gaugetime)-1)]-gaugetime[self.start_i] + float(parametersDict["soft_stop_time"]) / 60 / 24 + float(parametersDict["soft_start_time"]) / 60 / 24


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

    time_aggregate_periods = list(map(int, parametersDict["time_aggregate_periods"].split(';'))) if parametersDict["time_aggregate_enable"].lower() == "true" else [7, 30, 60]
    merge_period = max(
        time_aggregate_periods + [float(parametersDict["rain_event_merge_duration"])] + [5])
    rain_statistics, event_time = km2.rainStatistics(time_aggregate_periods, merge_period)
    rain_statistics_sort = np.flipud(np.sort(rain_statistics, axis = 0))

    # Create rain events
    for event_i in range(rain_statistics.shape[0]):
        rain_event = RainEvent()
        rain_event.start_i, rain_event.stop_i = event_time[event_i,:]
        rain_event.statistics = rain_statistics[event_i, :-1]
        rain_event.accumulated_rain = rain_statistics[event_i, -1]
        rain_event.time_aggregate_periods = time_aggregate_periods
        rain_events.append(rain_event)

    time_aggregate_return_period = float(parametersDict["time_aggregate_return_period"].split("+")[0])
    time_aggregate_extra_events = parametersDict["time_aggregate_return_period"].split("+")[1] if "+" in parametersDict["time_aggregate_return_period"] else 0

    # Filter rain events based on input
    include_events_total_rain_depth = float(parametersDict["include_events_total_rain_depth"])
    for rain_event in rain_events:
        if include_events_total_rain_depth > 0 and rain_event.accumulated_rain > include_events_total_rain_depth:
            rain_event.include = True
            if rain_event.accumulated_rain > 15:
                rain_event.reduce_timestep = True

        if parametersDict["time_aggregate_enable"].lower() == "true":
            for period_i in range(len(time_aggregate_periods)-1):
                if rain_event.statistics[period_i] > rain_statistics_sort[int(dataperiod/365.0/time_aggregate_return_period)-time_aggregate_extra_events, period_i]:
                    rain_event.include = True
                    rain_event.selected_because_of.append(str(time_aggregate_periods[period_i]))

        if rain_event.include:
            for period_i in range(len(time_aggregate_periods)-1):
                if rain_event.statistics[period_i] > rain_statistics_sort[int(dataperiod/365.0/2), period_i]:
                    rain_event.reduce_timestep = True

    rain_events_included = [rain_event for rain_event in rain_events if rain_event.include]

    total_duration = np.sum([rain_event.duration for rain_event in rain_events if rain_event.include])


    def duration_to_string(duration):
        label = ""
        years, remainder = divmod(duration, 365)
        months, remainder = divmod(remainder, 365.0/12)
        days, remainder = divmod(remainder, 1)
        hours = remainder/ (1.0/24)

        label += "%d years " % years if years else ""
        label += "%d months " % months if months else ""
        label += "%d days " % days if days else ""
        label += "%1.1f hours " % hours if hours else ""

        return label

    #####
    # Write LTS file through jinga2 module
    # Jinga2 reads a text file and uses it as a template for creating LTS files
    #####
    logFile.write(str(dtnow.now()) + ": Writing LTS file using Jinga2\n")

    Environment(loader=FileSystemLoader(''))
    # Load LTS template file
    templateFile = open(
        os.path.join(
            scriptFolder,
            'LTSTemplate.MJL'),
        'r')
    templateFileStr = templateFile.read()

    getTimeRE = re.compile("(\d+:\d+:\d+)")

    # Write LTS file with event information vectors
    fout = open(parametersDict["output_mjl"], 'w+')
    fout.write(
        Environment().from_string(templateFileStr).render(
            inputfile=parametersDict["input_file"],
            simulation_start=[rain_event.event_start_time.strftime('%Y-%m-%d %H:%M:00') for rain_event in rain_events_included],
            simulation_stop=[rain_event.event_stop_time.strftime('%Y-%m-%d %H:%M:00') for rain_event in rain_events_included],
            job_number=range(
                1,
                len(rain_events_included) + 1),
            dur_time=[duration_to_string(rain_event.duration) for rain_event in rain_events_included],
            jobs=len(rain_events_included),
            total_dur_time=duration_to_string(total_duration),
            dataperiod=duration_to_string(dataperiod),
            accumulated_rain=[round(rain_event.accumulated_rain,2) for rain_event in rain_events_included],
            eventdts=[[", ".join(rain_event.selected_because_of)] for rain_event in rain_events_included],
            date_criteria=parametersDict["date_criteria"],
            configStr=configStr,
            hotstart_param = parametersDict["hotstart_param"],
            hotstart_name = parametersDict["hotstart_name"],
            hotstart_date = parametersDict["hotstart_date"],
            hotstart_time = [getTimeRE.findall(rain_event.event_start_time.strftime('%Y-%m-%d %H:%M:00'))[0] for rain_event in rain_events_included],
            event_reduce_timestep = [rain_event.reduce_timestep for rain_event in rain_events_included],
            event_increase_timestep = [False] + [rain_event.reduce_timestep for rain_event in rain_events_included]))
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
