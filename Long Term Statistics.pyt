# Tool for reading DFS0 or KM2 files and creating LTS files from it
# Created by Emil Nielsen
# Contact: 
# E-mail: enielsen93@hotmail.com

import arcpy
import os
import sys
import numpy as np
thisFolder = os.path.dirname(__file__)
scriptFolder = os.path.join(thisFolder, r"scripts")
sys.path.append(scriptFolder)
import generate_lts
from shutil import copyfile
import re
from datetime import datetime
import matplotlib.dates as dates
import datetime
import re
import bisect  # math bisection method
from time import sleep
import subprocess
import copy
import time

def run_mex(mouse_sim_launch, mex_file):
    run_cmd = r'"%s" "%s" "HD" "Run" "Close" "NoPrompt"' % (mouse_sim_launch, mex_file)
    subprocess.Popen(run_cmd)

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

class ToolParameters:
    def __init__(self, parameters):
        self.parameters = parameters
        self.values = {}

        for parameter in parameters:
            setattr(self, parameter.name, parameter)
            if " " in parameter.name:
                arcpy.AddError("Error: Space in parameter %s" % parameter.name)
            if type(parameter.value) is float or parameter.value is None or parameter.value is bool:
                self.values[parameter.name] = parameter.Value
            else:
                self.values[parameter.name] = parameter.ValueAsText

class Toolbox(object):
    def __init__(self):
        """Define the toolbox (the name of the toolbox is the name of the
        .pyt file)."""
        self.label = "Long Term Statistics"
        self.alias = ""

        # List of tool classes associated with this toolbox
        self.tools = [LTSGenerator,LTSCombiner, LTSSplitter, LTSExtractor, DFS0Reducer, LTSSplitterMex]


class LTSGenerator(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Create LTS-file"
        self.description = "Generate a Mike Urban LTS (Long Term Statistics) file from a DFS0 or KM2 file for a rain gauge. \n\nCreated by: Emil Nielsen \nContact: enielsen93@hotmail.com"
        self.canRunInBackground = True

    def getParameterInfo(self):
        #Define parameter definitions

        # First parameter
        input_file = arcpy.Parameter(
            displayName="Input DFS0 or KM2 file (Rain series)",
            name="input_file",
            datatype="File",
            parameterType="Required",
            direction="Input")
        input_file.filter.list=["dfs0","km2", "kmd", "txt"]
            
        output_mjl = arcpy.Parameter(
            displayName="Output MOUSE job list (MJL) file",
            name="output_mjl",
            datatype="File",
            parameterType="Required",
            direction="Output")
        output_mjl.filter.list=["MJL"]
        
        time_aggregate_enable = arcpy.Parameter(
            displayName="Use time aggregates",
            name="time_aggregate_enable",
            datatype="Boolean",
            parameterType="optional",
            direction="Output")
        time_aggregate_enable.value = False
        time_aggregate_enable.category = "Time Aggregate"
            
        time_aggregate_periods = arcpy.Parameter(
            displayName="Time aggregates: [min]",
            name="time_aggregate_periods",
            datatype="Long",
            parameterType="Optional",
            direction="Input",
            multiValue=True)
        time_aggregate_periods.value = [10,30,60,180,360]
        time_aggregate_periods.enabled = False
        time_aggregate_periods.filter.type = "Range"
        time_aggregate_periods.filter.list = [1, 100000]
        time_aggregate_periods.category = "Time Aggregate"
        
        include_events_total_rain_depth = arcpy.Parameter(
            displayName="Include all rain events with total rain depth above: [mm]",
            name="include_events_total_rain_depth",
            datatype="double",
            parameterType="Optional",
            direction="Input")
        
        soft_start_time = arcpy.Parameter(
            displayName="Soft start time [min]",
            name="soft_start_time",
            datatype="double",
            parameterType="Optional",
            direction="Input")
        soft_start_time.value = 0
        
        soft_stop_time = arcpy.Parameter(
            displayName="Soft stop time [min]",
            name="soft_stop_time",
            datatype="double",
            parameterType="Optional",
            direction="Input")
        soft_stop_time.value = 0
        
        include_events_return_period = arcpy.Parameter(
            displayName="Include all rain events with return period above: [year]",
            name="include_events_return_period",
            datatype="String",
            parameterType="Optional",
            direction="Input")
        include_events_return_period.enabled = False
            
        time_aggregate_return_period = arcpy.Parameter(
            displayName="Return period above: [year]",
            name="time_aggregate_return_period",
            datatype="String",
            parameterType="Optional",
            direction="Input")
        time_aggregate_return_period.category = "Time Aggregate"
        
        time_aggregate_number_events = arcpy.Parameter(
            displayName="Number of included events for each aggregate:",
            name="time_aggregate_number_events",
            datatype="long",
            parameterType="Optional",
            direction="Input")
        time_aggregate_number_events.category = "Time Aggregate"
        time_aggregate_number_events.enabled = False
        
        time_series_duration = arcpy.Parameter(
            displayName="Duration of time series: [years]",
            name="time_series_duration",
            datatype="double",
            parameterType="Optional",
            direction="Input")
            
        date_criteria = arcpy.Parameter(
            displayName="Select only rain events between these two dates:",
            name="date_criteria",
            datatype="String",
            parameterType="Optional",
            direction="Input")
        date_criteria.enabled = True
        
        dfs0_output_enable = arcpy.Parameter(
            displayName="Save DFS0 file",
            name="dfs0_output_enable",
            datatype="Boolean",
            parameterType="Optional",
            direction="Input")
        dfs0_output_enable.category = "Save DFS0 file"

        dfs_output = arcpy.Parameter(
            displayName="Output DFS0 file",
            name="dfs_output",
            datatype="File",
            parameterType="Optional",
            direction="Output")
        # dfs_output.filter = ["dfs0"]
        dfs_output.category = "Save DFS0 file"
        dfs_output.enabled = False

        input_file = arcpy.Parameter(
            displayName="Input DFS0 or KM2 file (Rain series)",
            name="input_file",
            datatype="File",
            parameterType="Required",
            direction="Input")
        input_file.filter.list = ["dfs0", "km2", "kmd", "txt"]
        
        rain_event_merge = arcpy.Parameter(
            displayName="Merge rain events over dry periods",
            name="rain_event_merge",
            datatype="Boolean",
            parameterType="Optional",
            direction="input")
        rain_event_merge.category = "Merge rain events"
        
        rain_event_merge_duration = arcpy.Parameter(
            displayName="Merge rain events with dry periods shorter than [min]",
            name="rain_event_merge_duration",
            datatype="double",
            parameterType="Optional",
            direction="input")
        rain_event_merge_duration.category = "Merge rain events"
        rain_event_merge_duration.enabled = False
        rain_event_merge_duration.value = 60
        
        hotstart_param = arcpy.Parameter(
            displayName="Include hotstart",
            name="hotstart_param",
            datatype="Boolean",
            parameterType="Optional",
            direction="Input")    
        hotstart_param.category = "Include hotstart"
        # hotstart_param = False
        
        hotstart_name = arcpy.Parameter(
            displayName="Name of hotstart-file",
            name="hotstart_name",
            datatype="String",
            parameterType="Optional",
            direction="Input")    
        hotstart_name.category = "Include hotstart"
        # hotstart_name.filter.list = ["HotstartBase.PRF"]
        
        hotstart_date = arcpy.Parameter(
            displayName="Date of hotstart-file (format: YYYY-MM-DD)",
            name="hotstart_date",
            datatype="String",
            parameterType="Optional",
            direction="Input")
        hotstart_date.category = "Include hotstart"
        # hotstart_date.filter.list = ["2010-01-09"]

        enable_dtmin = arcpy.Parameter(
            displayName="Reduce timestep for extreme rain events",
            name="enable_dtmin",
            datatype="Boolean",
            parameterType="Optional",
            direction="input")
        enable_dtmin.category = "Additional Settings"
        
        
        
        # date_criteria = arcpy.Parameter(
            # displayName="Shorten DFS0 file",
            # name="dfs0_shorten_enable",
            # datatype="Boolean",
            # parameterType="Optional",
            # direction="Input")
    
        #    0    input_file: Input DFS0-file (Rain series) 
        #    1    output_mjl: Output MOUSE job list (MJL) file
        #    10    time_aggregate_enable: Use time aggregates
        #    9    time_aggregate_periods: Time aggregates: [min]
        #    3    include_events_total_rain_depth: Include all rain events with total rain depth above: [mm]
        #    5    soft_start_time: Soft start time [min]
        #    6    soft_stop_time: Soft stop time [min]
        #    4    include_events_return_period: Include all rain events with return period for total rain depth above: [year]
        #    7    time_aggregate_return_period: Return period above: [year]
        #    8    time_aggregate_number_events: Number of included events for each aggregate:
        #    2    time_series_duration: Duration of time series: [years]


        params = [input_file, output_mjl, time_series_duration, date_criteria, include_events_total_rain_depth, include_events_return_period, soft_start_time, soft_stop_time, time_aggregate_enable, time_aggregate_return_period, time_aggregate_number_events, time_aggregate_periods, rain_event_merge, rain_event_merge_duration, hotstart_param, hotstart_name, hotstart_date, dfs0_output_enable, dfs_output, enable_dtmin] #, param13, param14, date_criteria

        return params

    def isLicensed(self):
        return True

    def updateParameters(self, parameters):
        # Set the default distance threshold to 1/100 of the larger of the width
        #  or height of the extent of the input features.  Do not set if there is no 
        #  input dataset yet, or the user has set a specific distance (Altered is true).
        #
        tool_parameters = ToolParameters(parameters)

        if tool_parameters.input_file.altered:
            if not tool_parameters.output_mjl.valueAsText:
                filename, _ = os.path.splitext(tool_parameters.input_file.valueAsText)
                tool_parameters.output_mjl.value = filename +".MJL"
            if not parameters[2].valueAsText or not parameters[3].valueAsText:
                dataperiod,daterange = generate_lts.getDataPeriod(parameters,scriptFolder)
                if not parameters[2].valueAsText:
                    parameters[2].value = round(float(dataperiod),1)
                if not parameters[3].valueAsText:
                    parameters[3].value = daterange
            if not tool_parameters.dfs_output.valueAsText:
                if ".km2" in str(tool_parameters.input_file.valueAsText) and tool_parameters.dfs0_output_enable.value == True:
                    filename, _ = os.path.splitext(tool_parameters.input_file.valueAsText)
                    tool_parameters.dfs_output.value = filename + "reduced.dfs0"
        if parameters[8].value == True:
            parameters[9].enabled = True
            parameters[11].enabled = True
            parameters[10].enabled = True
        else:
            parameters[9].enabled = False
            parameters[11].enabled = False
            # parameters[10].enabled = False
        if tool_parameters.dfs0_output_enable.value == True:
            tool_parameters.dfs_output.enabled = True
        else:
            parameters[13].enabled = False
        if parameters[12].value == True:
            parameters[13].enabled = True
        else:
            parameters[13].enabled = False
        if parameters[14].value == True:
            parameters[15].enabled = True
            parameters[16].enabled = True
        else:
            parameters[15].enabled = False
            parameters[16].enabled = False

        return

    def updateMessages(self, parameters):
        if True:
            if parameters[4].value and parameters[4].value > 0:
                if parameters[5].value and parameters[5].value > 0:
                    parameters[5].setErrorMessage("Can't include by both total rain depth and return period.")
                    parameters[4].setErrorMessage("Can't include by both total rain depth and return period.")
            if parameters[8].value and parameters[8].value == True:
                if parameters[10].value and parameters[10].value > 0:
                    if parameters[9].value > 0:
                        parameters[9].setErrorMessage("Can't include by both return period and number of events.")
                        parameters[10].setErrorMessage("Can't include by both return period and number of events.")
            if parameters[3].value and parameters[3].value:
                if not " - " in parameters[3].value:
                    parameters[3].setErrorMessage("Date range must contain \" - \" without the quotation marks.")
                    return
            if parameters[12].value == True:
                if not parameters[13].value > 0:
                    parameters[13].setErrorMessage("Must select duration of dry period")            
    def execute(self, parameters, messages):
        # tool_parameters = ToolParameters(parameters)
        generate_lts.writeLTS(parameters, scriptFolder)
        return
        
class LTSCombiner(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Combine LTS-files"
        self.description = "Combine two LTS-files. \n\nCreated by: Emil Nielsen \nContact: enielsen93@hotmail.com"
        self.canRunInBackground = False

    def getParameterInfo(self):
        #Define parameter definitions

        # First parameter
        input_file = arcpy.Parameter(
            displayName="Input LTS files",
            name="LTS_input_files",
            datatype="File",
            parameterType="Required",
            direction="Input",
            multiValue=True)
        input_file.filter.list=["mjl"]
        
        param1 = arcpy.Parameter(
            displayName="Output combined LTS file",
            name="LTS_output_file",
            datatype="File",
            parameterType="Required",
            direction="Output")
        param1.filter.list=["mjl"]
        
        params = [input_file, param1]

        return params

    def isLicensed(self):
        return True

    def updateParameters(self, parameters):        
        # Set the default distance threshold to 1/100 of the larger of the width
        #  or height of the extent of the input features.  Do not set if there is no 
        #  input dataset yet, or the user has set a specific distance (Altered is true).
        #
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        generate_lts.combineLTS(parameters,scriptFolder)
        return       
        
class LTSSplitter(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Split LTS file"
        self.description = "Split LTS-file. \n\nCreated by: Emil Nielsen \nContact: enielsen93@hotmail.com"
        self.canRunInBackground = True

    def getParameterInfo(self):
        #Define parameter definitions

        # First parameter
        LTSFile = arcpy.Parameter(
            displayName="Input LTS file",
            name="LTSFile",
            datatype="File",
            parameterType="Optional",
            direction="Input",)
        LTSFile.filter.list=["mjl"]
        
        LTSCount = arcpy.Parameter(
            displayName="Split into how many files?",
            name="LTSCount",
            datatype="Long",
            parameterType="Required",
            direction="Input")
            
        RunoffFile = arcpy.Parameter(
            displayName="Runoff file (for making dupilcates)",
            name="RunoffFile",
            datatype="File",
            parameterType="Optional",
            direction="Input",)
        RunoffFile.filter.list=["crf"]
        
        MU_database = arcpy.Parameter(
            displayName="Mike Urban database (for duplicating network simulations)",
            name="database",
            datatype="DEWorkspace",
            parameterType="Optional",
            direction="Input")
            
        network_simulation = arcpy.Parameter(
            displayName= "Network Simulation (for duplicating network simulations)",  
            name="network_simulation",  
            datatype="GPString",  
            parameterType="Optional",  
            direction="Input") 
        
        params = [LTSFile, LTSCount, RunoffFile, MU_database, network_simulation]

        return params

    def isLicensed(self):
        return True

    def updateParameters(self, parameters):  
        MU_Database = parameters[3].ValueAsText
        if parameters[3].altered:
                parameters[4].filter.list = [row[0] for row in arcpy.da.SearchCursor(MU_Database + "\msm_Project","MUID")]
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        LTSFile = parameters[0].ValueAsText
        splitCount = int(parameters[1].ValueAsText)
        RunoffFile = parameters[2].ValueAsText
        MU_Database = parameters[3].ValueAsText
        networkSimulation = parameters[4].ValueAsText
        
        if LTSFile:
            with open(LTSFile,'r') as f:
                LTSTxt = f.read().split("\n")
            LTSTxt = np.array(LTSTxt)

            LTSStartLine = [i for i,a in enumerate(LTSTxt) if r"[SIMULATION_EVENT]" in a]
            LTSEndLine = [i for i,a in enumerate(LTSTxt) if r"EndSect  // SIMULATION_EVENT" in a]
            
            #LTSNewTxt = LTSTxt[0:LTSStartLine[0]]
            jobsPerLTS = int(np.ceil(float(len(LTSStartLine))/splitCount))
            LTSNewTxt = [[] for i in range(splitCount)]
            job = 0
            for i in np.arange(0,len(LTSStartLine),jobsPerLTS):
                LTSNewTxt[job] = LTSTxt[0:LTSStartLine[0]]
                for j in np.arange(i,min(i+jobsPerLTS,len(LTSStartLine))):
                    LTSNewTxt[job] = np.concatenate((LTSNewTxt[job],LTSTxt[LTSStartLine[j]:LTSEndLine[j]+1]))
                LTSNewTxt[job] = np.concatenate((LTSNewTxt[job],LTSTxt[LTSEndLine[-1]+1:len(LTSTxt)]))
                
                with open(LTSFile[0:-4] + "_Split%d.mjl" % (job+1),'w') as f:
                    f.write("\n".join(LTSNewTxt[job]))
                job += 1
        
        if MU_Database:
            networkSims = [row[0] for row in arcpy.da.SearchCursor(MU_Database + "\msm_Project","MUID")]
            fields = [a.name for a in arcpy.ListFields(MU_Database + r"\msm_Project")][1:]
            fieldValues = []

            with arcpy.da.SearchCursor(MU_Database + "\msm_Project", fields, where_clause = "MUID = '%s'" % networkSimulation) as cursor:
                for row in cursor:
                    fieldValues = []
                    for f in row:
                        fieldValues.append(f)
            MUIDs = [row[0] for row in arcpy.da.SearchCursor(MU_Database + "\msm_Project","MUID")]
            
            for job in range(splitCount):
                splitName = networkSimulation + "_Split%d" % (job+1)
                if not splitName in MUIDs:
                    with arcpy.da.InsertCursor(MU_Database + "\msm_Project",fields) as cursor:
                        fieldValues[0] = splitName
                        fieldValues[[i for i,a in enumerate(fields) if a=="MJLFileName"][0]] = LTSFile[0:-4] + "_Split%d.mjl" % (job+1)
                        fieldValues[[i for i,a in enumerate(fields) if a=="CRFFileName"][0]] = RunoffFile[0:-4] + "_Split%d.CRF" % (job+1)
                        cursor.insertRow(fieldValues)
                else:
                    arcpy.AddWarning("Job %s already exists in msm_Project and will not be modified. Delete it from the Mike Urban model to update it." % (splitName))
        
        if RunoffFile:
            for job in range(splitCount):
                arcpy.AddMessage("Copying %s" % RunoffFile)
                copyfile(RunoffFile,RunoffFile[0:-4] + "_Split%d.CRF" % (job+1))
                
        return

class LTSSplitterMex(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Split LTS simulations and run simultaneously"
        self.description = "Split LTS-file. \n\nCreated by: Emil Nielsen \nContact: enielsen93@hotmail.com"
        self.canRunInBackground = True

    def getParameterInfo(self):
        #Define parameter definitions

        # First parameter
        LTSFile = arcpy.Parameter(
            displayName="Input LTS file",
            name="LTSFile",
            datatype="File",
            parameterType="Optional",
            direction="Input",)
        LTSFile.filter.list=["mjl"]
        
        LTSCount = arcpy.Parameter(
            displayName="Split into how many files?",
            name="LTSCount",
            datatype="Long",
            parameterType="Required",
            direction="Input")
            
        RunoffFile = arcpy.Parameter(
            displayName="Runoff file (for making dupilcates)",
            name="RunoffFile",
            datatype="File",
            parameterType="Optional",
            direction="Input",)
        RunoffFile.filter.list=["crf"]
        
        mex_file = arcpy.Parameter(
            displayName="Mex Network File",
            name="mex_file",
            datatype="File",
            parameterType="Required",
            direction="Input")
        mex_file.filter.list=["mex"]
        
        mouse_sim_launch = arcpy.Parameter(
            displayName="MOUSE Sim Launch Executable path",
            name="mouse_sim_launch",
            datatype="File",
            parameterType="Optional",
            direction="Input")
        mouse_sim_launch.filter.list=["exe"]
        
        mouse_sim_launch_paths = [r"C:\Program Files (x86)\DHI\2020\bin\x64\MOUSESimLaunch.exe"]
        for path in mouse_sim_launch_paths:
            for year in reversed(range(2010,2030)):
                if os.path.exists(path.replace("2020",str(year))):
                    mouse_sim_launch.value = path.replace("2020",str(year))
                    break
        
        
        params = [LTSFile, LTSCount, RunoffFile, mex_file, mouse_sim_launch]

        return params

    def isLicensed(self):
        return True

    def updateParameters(self, parameters):  
                
        
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        LTSFile = parameters[0].ValueAsText
        splitCount = int(parameters[1].ValueAsText)
        RunoffFile = parameters[2].ValueAsText
        mex_file = parameters[3].ValueAsText
        mouse_sim_launch = parameters[4].ValueAsText
        
        with open(LTSFile,'r') as f:
            LTSTxt = f.read().split("\n")
        LTSTxt = np.array(LTSTxt)

        LTSStartLine = [i for i,a in enumerate(LTSTxt) if r"[SIMULATION_EVENT]" in a]
        LTSEndLine = [i for i,a in enumerate(LTSTxt) if r"EndSect  // SIMULATION_EVENT" in a]
        
        #LTSNewTxt = LTSTxt[0:LTSStartLine[0]]
        jobsPerLTS = int(np.ceil(float(len(LTSStartLine))/splitCount))
        LTSNewTxt = [[] for i in range(splitCount)]
        
        LTS_files = []
        if LTSFile:
            job = 0
            for i in np.arange(0,len(LTSStartLine),jobsPerLTS):
                LTSNewTxt[job] = LTSTxt[0:LTSStartLine[0]]
                for j in np.arange(i,min(i+jobsPerLTS,len(LTSStartLine))):
                    LTSNewTxt[job] = np.concatenate((LTSNewTxt[job],LTSTxt[LTSStartLine[j]:LTSEndLine[j]+1]))
                LTSNewTxt[job] = np.concatenate((LTSNewTxt[job],LTSTxt[LTSEndLine[-1]+1:len(LTSTxt)]))
                
                LTS_files.append(LTSFile[0:-4] + "_Split%d.mjl" % (job+1))
                with open(LTS_files[-1],'w') as f:
                    f.write("\n".join(LTSNewTxt[job]))
                job += 1
        
        with open(mex_file,"r") as f:
            mex_file_text = f.readlines()
        mex_file_text_mjl_lineno = [lineno for lineno,line in enumerate(mex_file_text) if "MJL_file" in line][0]
        mex_file_text_CRF_lineno = [lineno for lineno,line in enumerate(mex_file_text) if "CRF_file" in line][0]
        
        for job in range(splitCount):
            
            if RunoffFile:
                arcpy.AddMessage("Copying %s" % RunoffFile)
                RunoffFileNew = RunoffFile[0:-4] + "_Split%d.CRF" % (job+1)
                if (not os.path.isfile(RunoffFileNew) or 
                    not os.path.getsize(RunoffFile)==os.path.getsize(RunoffFileNew) or 
                    os.path.getmtime(RunoffFile)>os.path.getmtime(RunoffFileNew)):
                    copyfile(RunoffFile, RunoffFileNew)
                else:
                    time.sleep(5)
    
            mex_file_new = copy.deepcopy(mex_file.replace(".mex","_Split%d.mex" % (job+1)))
            
            mex_file_new_text = copy.deepcopy(mex_file_text)
            if LTSFile:
                mex_file_new_text[mex_file_text_mjl_lineno] = re.sub("'[^']*'", "'%s'" % LTS_files[job], mex_file_new_text[mex_file_text_mjl_lineno])
            else:
                mex_file_new_text[mex_file_text_mjl_lineno] = mex_file_new_text[mex_file_text_mjl_lineno].replace(".MJL", "_Split%d.MJL" % (job+1))
            
            if RunoffFile:
                mex_file_new_text[mex_file_text_CRF_lineno] = re.sub("'[^']*'", "'%s'" % RunoffFileNew, mex_file_new_text[mex_file_text_CRF_lineno])
            else:
                mex_file_new_text[mex_file_text_CRF_lineno] = mex_file_new_text[mex_file_text_CRF_lineno].replace(".CRF", "_Split%d.CRF" % (job+1))
            
            # copyfile(mex_file, mex_file_new)
            with open(mex_file_new,'w+') as f:
                f.writelines(mex_file_new_text)
            if mouse_sim_launch:
                run_mex(mouse_sim_launch, mex_file_new)                
        return
        
class LTSExtractor(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Extract jobs from LTS file"
        self.description = "Extract jobs from LTS-file. \n\nCreated by: Emil Nielsen \nContact: enielsen93@hotmail.com"
        self.canRunInBackground = False

    def getParameterInfo(self):
        #Define parameter definitions

        # First parameter
        LTSFile = arcpy.Parameter(
            displayName="Input LTS file",
            name="LTSFile",
            datatype="File",
            parameterType="Required",
            direction="Input",)
        LTSFile.filter.list=["mjl"]
            
        dates = arcpy.Parameter(
            displayName="Date and time of jobs to extract from (same format as in MJL-file - split each date by a comma) Example: '2016-08-28 19:10:13',  '2014-09-06 18:17:13',  '2017-07-23 07:12:13'",
            name="dates",
            datatype="GPString",
            parameterType="Required",
            direction="Input",)
            
        LTSNewFile = arcpy.Parameter(
            displayName="Output LTS file",
            name="LTSNewFile",
            datatype="File",
            parameterType="Required",
            direction="Output",)
        LTSNewFile.filter.list=["mjl"]
        
        params = [LTSFile, dates, LTSNewFile]

        return params

    def isLicensed(self):
        return True

    def updateParameters(self, parameters):  
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        eventsTxt = parameters[1].ValueAsText
        LTSFile = parameters[0].ValueAsText
        LTSNewFile = parameters[2].ValueAsText
        with open(LTSFile,'r') as f:
            LTSTxt = f.read()
        simStarts = [datetime.datetime.strptime(a,"%Y-%m-%d %H:%M:%S") for a in re.findall(r"Simulation_start ?= ?'([^']+)'",LTSTxt)]
        simEnds = [datetime.datetime.strptime(a,"%Y-%m-%d %H:%M:%S") for a in re.findall(r"Simulation_end ?= ?'([^']+)'",LTSTxt)]
        
        events = re.findall(r"\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}",eventsTxt)
        matches = []
        for eventi in range(len(events)):
            above = [i for i,a in enumerate(simEnds) if a>=datetime.datetime.strptime(events[eventi],"%Y-%m-%d %H:%M:%S")]
            below = [i for i,a in enumerate(simStarts) if a<=datetime.datetime.strptime(events[eventi],"%Y-%m-%d %H:%M:%S")]
            match = [value for value in above if value in below]
            if len(match)>1:
                arcpy.AddWarning("Event %s is in multiple jobs" % events[eventi])
            elif len(match)==0:
                arcpy.AddWarning("Event %s is not included in the job list?" % events[eventi])
            for m in match:
                matches.append(m)
        matches = np.sort(matches)

        with open(LTSFile,'r') as f:
            LTSTxt = f.read().split("\n")
        LTSTxt = np.array(LTSTxt)

        LTSStartLine = [i for i,a in enumerate(LTSTxt) if r"[SIMULATION_EVENT]" in a]
        LTSEndLine = [i for i,a in enumerate(LTSTxt) if r"EndSect  // SIMULATION_EVENT" in a]
        LTSNewTxt = LTSTxt[0:LTSStartLine[0]]
        for m in matches:
            LTSNewTxt = np.concatenate((LTSNewTxt,LTSTxt[LTSStartLine[m]:LTSEndLine[m]+1]))
            
        LTSNewTxt = np.concatenate((LTSNewTxt,LTSTxt[LTSEndLine[-1]+1:len(LTSTxt)]))
        with open(LTSNewFile,'w') as f:
            f.write("\n".join(LTSNewTxt))
        arcpy.AddMessage("Extracted %d jobs from LTS file" % len(matches))
        return
        
class DFS0Reducer(object):
    def __init__(self):
        self.label       = "Reduce DFS0 to jobs in LTS"
        self.description = "Reduce DFS0 to jobs in LTS"
        self.canRunInBackground = True

    def getParameterInfo(self):
        #Define parameter definitions

        # Input Features parameter        
        dfs0File = arcpy.Parameter(
            displayName="Input DFS0 File",
            name="dfs0File",
            datatype="DEFile",
            parameterType="Required",
            direction="Input")
        dfs0File.filter.list=["dfs0"]
        
        ltsFile = arcpy.Parameter(
            displayName="LTS file",
            name="ltsFile",
            datatype="DEFile",
            parameterType="Required",
            direction="Input")
        ltsFile.filter.list=["mjl"]
        
        dfs0FileNew = arcpy.Parameter(
            displayName="Output DFS0 File",
            name="dfs0FileNew",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        dfs0FileNew.filter.list=["dfs0"]
            
        parameters = [dfs0File, ltsFile, dfs0FileNew]
        
        return parameters

    def isLicensed(self): #optional
        return True

    def updateParameters(self, parameters): #optional
        if parameters[1].altered:
            if not parameters[2].valueAsText:
                filename, _ = os.path.splitext(parameters[1].valueAsText)
                parameters[2].value = filename +"Reduced.dfs0"
        return

    def updateMessages(self, parameters): #optional
       
        return

    def execute(self, parameters, messages):
        dfs0File = parameters[0].ValueAsText
        ltsFile = parameters[1].ValueAsText
        dfs0FileNew = parameters[2].ValueAsText
        
        scriptFolder = os.path.dirname(__file__)

        [gaugetime,gaugeint] = readDFS0(dfs0File)

        with open(ltsFile,'r') as f:
            ltsFileTxt = f.read()

        getSimStart = re.compile(r"Simulation_start ?= ?'([^']+)'")
        getSimEnd = re.compile(r"Simulation_end ?= ?'([^']+)'")

        simStart = getSimStart.findall(ltsFileTxt)
        simEnd = getSimEnd.findall(ltsFileTxt)

        idx = np.array([],dtype=int)
        for i in range(len(simStart)):
            idx = np.concatenate((idx,np.arange(
                bisect.bisect_left(gaugetime,dates.date2num(
                    datetime.datetime.strptime(
                            simStart[i],
                        "%Y-%m-%d %H:%M:%S"))),
                bisect.bisect_right(gaugetime,dates.date2num(
                    datetime.datetime.strptime(
                        simEnd[i],
                        "%Y-%m-%d %H:%M:%S")))).tolist()))
        gaugetimeFilt = np.array([gaugetime[i] for i in idx])
        gaugeintFilt = np.array([gaugeint[i] for i in idx])
        
        gaugetimeFilt_gap_filled = list(gaugetimeFilt)
        gaugeintFilt_gap_filled = list(gaugeintFilt)
        for i in range(1,len(gaugetimeFilt)):
            if (gaugetimeFilt[i]-gaugetimeFilt[i-1])*24*60>30:
                gaugetimeFilt_gap_filled.extend([gaugetimeFilt[i-1] + 1.0/60/24, gaugetimeFilt[i] - 1.0/60/24])
        gaugetimeFilt_gap_filled.extend([gaugetimeFilt[-1] + 1.0/60/24])
        
        gaugeintFilt_gap_filled.extend(np.zeros(len(gaugetimeFilt_gap_filled)-len(gaugetimeFilt)))
        gaugetimeFilt_gap_filled = np.array(gaugetimeFilt_gap_filled)
        gaugeintFilt_gap_filled = np.array(gaugeintFilt_gap_filled)
        sort_idx = np.argsort(gaugetimeFilt_gap_filled)
        gaugetimeFilt_gap_filled = gaugetimeFilt_gap_filled[sort_idx]
        gaugeintFilt_gap_filled = gaugeintFilt_gap_filled[sort_idx]

        gaugetimeFilt_gap_filled = np.concatenate((gaugetime[0:3],gaugetimeFilt_gap_filled,gaugetime[-3:]))
        gaugeintFilt_gap_filled = np.concatenate(([0,0,0],gaugeintFilt_gap_filled,[0,0,0]))

        writeDFS0(gaugetimeFilt_gap_filled, gaugeintFilt_gap_filled, dfs0FileNew)
        return