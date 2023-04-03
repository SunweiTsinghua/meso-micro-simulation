#######################################################################################################################
#   the import of various system modules and packages
#######################################################################################################################

import json
import os
import shutil
import threading

from datetime import datetime, timedelta

from runningFoam.simulatingFoam import SimulatingFoam
from runningWRF.downloading import DownloadingMetData
from runningWRF.simulatingWPS import SimulatingWPS
from runningWRF.simulatingWRF import SimulatingWRF

#######################################################################################################################
#   run WPS
#######################################################################################################################
def execute_WPS(main_control_json):
    '''
    the function prepares and runs the WPS, which is a typical preprocessor to run WRF
    the running of WPS is controlled by the dict of wps_Dict with the following keys:
    1. start_date_time : the start date time (simulation time) to run the coupled simulation
    2. end_date_time: the end date time (simulation time) to terminate the coupled simulation
    3. domain_lonlat: the longitude/latitude pairs defining the corners of the coupling interface, which is always a box
    4. domain_coordinates: the coordinates (x, y coordinates in meters) pairs defining the corners of the coupling interface
    5. wps_directory: the directory containing the executables of the WPS, including geogrid.exe, ungrib.exe and metgrid.exe
    6. run_directory: the directory containing the running files and required files to initiate the WPS run
    7. template_directory: the directory containing the templates to run the WPS
    8. geo_directory: the directory containing the static geo files required to run WPS
    9. domain_dimension: the length and width of the box defning the coupled interface

    Input :
    main_control_json : a python dict to contain all the information required to run the coupled simulation

    Output :
    None
    '''

    start_date_time = main_control_json["start_date_time"]
    end_date_time = main_control_json["end_date_time"]

    domain_lonlat = main_control_json["coreLonsLats"]
    domain_coordinates = main_control_json["coreCoords"]

    wps_Dict = {
                    "wps_directory" : main_control_json["wps_directory"],
                    "run_directory" : os.path.join(main_control_json["case_directory"], "WPS_run"),
                    "template_directory" : os.path.join(main_control_json["template_directory"], "WRF"),
                    "met_data_directory" : os.path.join(main_control_json["case_directory"], "met_data"),
                    "geo_directory" : main_control_json["geo_directory"],

                    "start_date_time" : start_date_time,
                    "end_date_time" : end_date_time,
                    "center_latitude" : (domain_lonlat[0][0] + domain_lonlat[1][0]) / 2.0,
                    "center_longitude" : (domain_lonlat[0][1] + domain_lonlat[1][1]) / 2.0,
                    "domain_dimension" : [domain_coordinates[1][0]-domain_coordinates[0][0], domain_coordinates[1][1]-domain_coordinates[0][1]]
                }

    run_WPS = SimulatingWPS(wps_Dict)        
    run_WPS.run_geogrid()
    run_WPS.run_ungrib()
    run_WPS.run_metgrid()

#######################################################################################################################
#   run WRF
#######################################################################################################################
def execute_WRF(main_control_json):
    '''
    the function prepares and runs the WRF, which is the simulation tool at the meso-scale 
    the running of WRF is controlled by the dict of wrf_Dict with the following keys:
    1. start_date_time : the start date time (simulation time) to run the coupled simulation
    2. end_date_time: the end date time (simulation time) to terminate the coupled simulation
    3. wps_directory: the directory containing the running files of the WPS exectuable (run_directory in wps_Dict)
    4. run_directory: the directory containing the running files of the WRF exectuable
    5. wrf_directory: the directory containing the exectuable of WRF (wrf.exe and real.exe)
    6. template_directory: the directory containing the templates to run the WRF simulation
    7. coarse_grids: the number of the grids to run the first-round WRF simulation, the rest is for the fine-resolution simulation
    8. precice_xml_name: the name of the XML file controlling the precice coupling interface
    9. precice_participant: the name of the WRF part used in the precice coupling
    10. precice_mesh_name: the name of the mesh from the WRF part used in precice coupling
    11. precice_lats: the two latitudes gives the north and south boundaries for the precice interface (always a box)
    12. precice_lons: the two longitudes gives the east and west boundaries for the precice interface (always a box)
    13. precice_exchange_data: the name of the data used in the exchange information between WRF and OpenFOAM (U, k, epsilon etc)
    14. precice_cornerx/precice_cornery: the x and y coordinates defining the precice coupling interface
    15. precice_height: the height of the precice coupling interface
    16. precice_start_min: the minutes after which the WRF simualtion starts to port the simulation results to the OpenFOAM run
    17. precice_write_trans: the boolean to control wheter the WRF simulation results are written to the trans_dictory from which the OpenFOAM can read to initialize its run
    18. precice_trans_directory: the directory to store the WRF simulation results to initialize the OpenFOAM run
    19. precice_lid: wheter the lid of box defining the precice interface is included

    Input :
    main_control_json : a python dict to contain all the information required to run the coupled simulation

    Output :
    None
    '''

    start_date_time = main_control_json["start_date_time"]
    end_date_time = main_control_json["end_date_time"]

    precice_lonlat = main_control_json["coreLonsLats"]
    precice_coords = main_control_json["coreCoords"]
    precice_data = []

    for ii in range(len(main_control_json["precice_exchange_data"])) :
        precice_data.append(main_control_json["precice_exchange_data"][ii]["name"])

    wrf_Dict = {    
                    "wps_directory" : os.path.join(main_control_json["case_directory"], "WPS_run"),
                    "run_directory" : os.path.join(main_control_json["case_directory"], "WRF_run"), 
                    "wrf_directory" : main_control_json["wrf_directory"],
                    "template_directory" : os.path.join(main_control_json["template_directory"], "WRF"),
                    
                    "start_date_time" : start_date_time,
                    "end_date_time" : end_date_time,

                    "coarse_grids" : main_control_json["coarse_grids"],

                    "precice_xml_name": os.path.join(main_control_json["case_directory"], "precice-config.xml"),
                    "precice_participant_name": "WRF",
                    "precice_mesh_name" : "WRF-Mesh",
                    "precice_lats" : [precice_lonlat[0][0], precice_lonlat[1][0]],
                    "precice_lons" : [precice_lonlat[0][1], precice_lonlat[1][1]],
                    "precice_exchange_data" : precice_data,
                    "precice_Cmu" : 0.09,
                    "precice_cornerx": [precice_coords[0][0], precice_coords[1][0]],
                    "precice_cornery": [precice_coords[0][1], precice_coords[1][1]],
                    "precice_height" : main_control_json["precice_height"],
                    "precice_start_mins" : main_control_json["precice_wrf_start_mins"],
                    "precice_write_trans": bool(main_control_json["precice_write_init"]),
                    "precice_trans_directory" : os.path.join(main_control_json["case_directory"], "trans_data"),
                    "precice_lid" : False
               }

    run_WRF = SimulatingWRF(wrf_Dict)
    run_WRF.initialize()
    run_WRF.prepare_precice()
    run_WRF.run_real(kind="coarse")
    run_WRF.run_wrf(kind="coarse")
    
#######################################################################################################################
#   run OpenFOAM
#######################################################################################################################
def execute_foam(main_control_json):
    '''
    the function prepares and runs the openFOAM, which is the tool to present the micro-scale simulation
    the running of OpenFOAM is controlled by the dict of foam_Dict with the following keys:
    1. case_directory: the directory containing the running files of OpenFOAM
    2. mesh_directory: the directory containing the mesh for runing the OpenFOAM
    3. case_template: the directory containing the templates for running the OpenFOAM
    4. foam_mesh: the boolean to control whether the mesh is already in the OpenFOAM format
    5. min_mesh_size: the minimal size of the cell for meshing up the domain, not used in the case where foam_mesh=1
    6. case_time: the duration for runing the WRF-OpenFOAM coupled simulation 
    7. time_step: the time step for running the OpenFOAM simulation
    8. output_interval: the interval (a few seconds) for openFOAM to store the simulation
    9. number_process: the number of the preocess to run the parallel version of the OpenFOAM
    10. set_patches: the name of patches of the OpenFOAM domain to take the initial values from the WRF simulation
    11. precice_interface: the name of patches to exchange to precice and WRF
    12. trans_directory: the directory containing the file output from the WRF simualtion to initialize the openFOAM simulation

    Input :
    main_control_json : a python dict to contain all the information required to run the coupled simulation

    Output :
    None
    '''

    precice_data = []

    for ii in range(len(main_control_json["precice_exchange_data"])) :
        precice_data.append(main_control_json["precice_exchange_data"][ii]["name"])

    foam_Dict = {
                    "case_directory": os.path.join(main_control_json["case_directory"],"foam_run"),
                    "mesh_directory": main_control_json["mesh_directory"],
                    "case_template": os.path.join(main_control_json["template_directory"], "foam"),
                    "foam_mesh" : bool(main_control_json["foam_mesh"]),
                    "min_mesh_size" : main_control_json["min_mesh_size"],
                    "case_time" : main_control_json["exchange_maxTime"],
                    "time_step" : main_control_json["foam_time_step"],
                    "output_interval" : main_control_json["foam_output_interval"],

                    "number_process" : main_control_json["number_process"],

                    "set_patches" : main_control_json["foam_set_patches"],
                    "precice_interface" : main_control_json["foam_interface"],
                    "precice_data" : precice_data,
                    "trans_directory" : os.path.join(main_control_json["case_directory"], "trans_data")
                }

    run_Foam = SimulatingFoam(foam_Dict)
    run_Foam.prepare_case_directory()
    run_Foam.transMesh()
    run_Foam.setInitials()
    run_Foam.runPimpleFoam()

#######################################################################################################################
#   set up the precice xml control file
#######################################################################################################################
def setup_preciceXML(main_control_json) :
    '''
    the function write up a XML file to control the precice coupled run
    XML file is only file needed for the precice to run the WRF-OpenFOAM coupled run, the reference of the XML file can be found
    https://precice.org/configuration-xml-reference.html

    Input :
    main_control_json : a python dict to contain all the information required to run the coupled simulation

    Output :
    None
    '''

    import xml.etree.ElementTree as et
    import xml.dom.minidom as dm

    exchange_directory = main_control_json["case_directory"]

    if os.path.exists(os.path.join(exchange_directory, "precice-run")) :
        shutil.rmtree(os.path.join(exchange_directory, "precice-run"), ignore_errors=True)
        
    target_xml_name = os.path.join(main_control_json["case_directory"], "precice-config.xml")
    orig_xml_name = os.path.join(main_control_json["template_directory"], "precice.xml")

    xml_content = et.parse(orig_xml_name)
    xml_root = xml_content.getroot()

    solver_interface = xml_root[1]

    exchange_data = main_control_json["precice_exchange_data"]

    for ii in range(len(exchange_data)) :
        current_data = et.SubElement(solver_interface, "data:{:s}".format(exchange_data[ii]["type"]))
        current_data.set("name", exchange_data[ii]["name"])

    mesh_WRF = et.SubElement(solver_interface, "mesh")
    mesh_Foam = et.SubElement(solver_interface, "mesh")
    mesh_WRF.set("name", "WRF-Mesh")
    mesh_Foam.set("name", "OpenFOAM-Mesh")

    for ii in range(len(exchange_data)) :
        data_WRF = et.SubElement(mesh_WRF, "use-data")
        data_WRF.set("name", exchange_data[ii]["name"])

    for ii in range(len(exchange_data)) :
        data_Foam = et.SubElement(mesh_Foam, "use-data")
        data_Foam.set("name", exchange_data[ii]["name"])

    participant_WRF = et.SubElement(solver_interface, "participant")
    participant_WRF.set("name", "WRF")

    participant_WRF_mesh = et.SubElement(participant_WRF, "use-mesh")
    participant_WRF_mesh.set("name", "WRF-Mesh")
    participant_WRF_mesh.set("provide", "yes")   

    for ii in range(len(exchange_data)) :
        participant_WRF_writeData = et.SubElement(participant_WRF, "write-data")
        participant_WRF_writeData.set("name", exchange_data[ii]["name"])
        participant_WRF_writeData.set("mesh", "WRF-Mesh")

    participant_WRF_export = et.SubElement(participant_WRF, "export:vtu")
    participant_WRF_export.set("directory", "precice-export")

    participant_Foam = et.SubElement(solver_interface, "participant")
    participant_Foam.set("name", "OpenFOAM")

    participant_Foam_meshes = [et.SubElement(participant_Foam, "use-mesh"), et.SubElement(participant_Foam, "use-mesh")]
    participant_Foam_meshes[0].set("name", "OpenFOAM-Mesh")
    participant_Foam_meshes[0].set("provide", "yes")
    participant_Foam_meshes[1].set("name", "WRF-Mesh")
    participant_Foam_meshes[1].set("from", "WRF")

    for ii in range(len(exchange_data)) :
        participant_Foam_readData = et.SubElement(participant_Foam, "read-data")
        participant_Foam_readData.set("name", exchange_data[ii]["name"])
        participant_Foam_readData.set("mesh", "OpenFOAM-Mesh")

    participant_Foam_mapping = et.SubElement(participant_Foam, "mapping:nearest-projection")
    participant_Foam_mapping.set("direction", "read")
    participant_Foam_mapping.set("from", "WRF-Mesh")
    participant_Foam_mapping.set("to", "OpenFOAM-Mesh")
    participant_Foam_mapping.set("constraint", "consistent")

    socket = et.SubElement(solver_interface, "m2n:sockets")
    socket.set("from", "WRF")
    socket.set("to", "OpenFOAM")
    socket.set("exchange-directory", exchange_directory)

    if main_control_json["exchange_timeStep"] == 0.0 :
        coupling_scheme = et.SubElement(solver_interface, "coupling-scheme:serial-explicit")

    else :
        coupling_scheme = et.SubElement(solver_interface, "coupling-scheme:parallel-explicit")

    coupling_participants = et.SubElement(coupling_scheme, "participants")
    coupling_participants.set("first", "WRF")
    coupling_participants.set("second", "OpenFOAM")

    coupling_max_time = et.SubElement(coupling_scheme, "max-time")
    coupling_max_time.set("value", str(main_control_json["exchange_maxTime"]))

    coupling_time_step = et.SubElement(coupling_scheme, "time-window-size")
    
    if main_control_json["exchange_timeStep"] == 0.0 :
        coupling_time_step.set("value", "-1")
        coupling_time_step.set("method", "first-participant")

    else :
        coupling_time_step.set("value", str(main_control_json["exchange_timeStep"]))

    for ii in range(len(exchange_data)) :
        coupling_exchange_data = et.SubElement(coupling_scheme, "exchange")
        coupling_exchange_data.set("data", exchange_data[ii]["name"])
        coupling_exchange_data.set("mesh", "WRF-Mesh")
        coupling_exchange_data.set("from", "WRF")
        coupling_exchange_data.set("to", "OpenFOAM") 

    if os.path.isfile(target_xml_name):
        os.remove(target_xml_name)

    xml_content.write(target_xml_name)

#######################################################################################################################
#   execute the mesoscale simulation
#######################################################################################################################

def execute_MesoSimulation(main_control_json, meso_log_name):
    '''
    this function execute the meso-scale simulation, which includes downloading the ERA5 dataset from the website, running the WPS 
    to prepare the main run of the WRF and running WRF
    Input :
    main_control_json : a python dict to contain all the information required to run the coupled simulation
    meso_log_name: the name of the log file containing the top-level output messages from the meso-scale simualtion, these messages
    are mainly output from this python script
    '''
    from copy import deepcopy

    target_date_time = main_control_json["start_date_time"]
    target_place = main_control_json["target_place"]

    with open(os.path.join(main_control_json["template_directory"], target_place+".json"), "r") as json_file : 
        place_dict = json.load(json_file)

    meso_dict = deepcopy(main_control_json)
    meso_dict.update(place_dict)

    with open(meso_log_name, "a", buffering=1) as log_file:

        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S.%f")        
        log_file.write("{0:s}: Start the meso simulation for {1:s}-{2:s}\n".format(current_time, target_place, target_date_time))
        print("{0:s}: Start the meso simulation for {1:s}-{2:s}".format(current_time, target_place, target_date_time))

        # preparing the directory for downloading ERA5, running WPS and WRF
        prepare_WRF_directories(meso_dict)
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S.%f")       
        print("{0:s}: Prepared the WRF directories for {1:s}-{2:s}".format(current_time, target_place, target_date_time))
        log_file.write("{0:s}: Prepared the WRF directories for {1:s}-{2:s}\n".format(current_time, target_place, target_date_time))

        # download the ERA5 data if not already downloaded, whether the ERA5 files are already downloaded is indicated by "downloaded_era5" key in the main control dict
        if not main_control_json["downloaded_era5"] :
            downloading_era5(meso_dict)
            current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S.%f")  
            print("{0:s}: Downloaded the met data for {1:s}-{2:s}".format(current_time, target_place, target_date_time))
            log_file.write("{0:s}: Downloaded the met data for : {1:s}-{2:s}\n".format(current_time, target_place, target_date_time))

        # running WPS, which prepare necessary files for running WRF
        execute_WPS(meso_dict)
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S.%f")  
        print("{0:s}: Executed the WPS for {1:s}-{2:s}".format(current_time, target_place, target_date_time))
        log_file.write("{0:s}: Executed the WPS for {1:s}-{2:s}\n".format(current_time, target_place, target_date_time))

        # running WRF itself, which certainly links to the PreCICE to present the coupled run
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S.%f")  
        print("{0:s}: Executing the WRF for {1:s}-{2:s}\n".format(current_time, target_place, target_date_time))
        log_file.write("{0:s}: Executing the WRF for {1:s}-{2:s}\n".format(current_time, target_place, target_date_time))
        execute_WRF(meso_dict)

#######################################################################################################################
#   execute the microscale simulation
#######################################################################################################################

def execute_MicroSimulation(main_control_json, micro_log_name):
    '''
    this function execute the micro-scale simulation, which mainly just the OpenFOAM, it is noted that the running of OpenFOAM includes
    setting up the initial fields from the results of the WRF, and the simulation coupled WRF and OpenFOAM
    Input :
    main_control_json : a python dict to contain all the information required to run the coupled simulation
    micro_log_name: the name of the log file containing the top-level output messages from the micro-scale simualtion, these messages
    are mainly output from this python script
    '''

    from copy import deepcopy

    target_date_time = main_control_json["start_date_time"]
    target_place = main_control_json["target_place"]

    with open(os.path.join(main_control_json["template_directory"], target_place+".json"), "r") as json_file :
        place_dict = json.load(json_file)

    micro_dict = deepcopy(main_control_json)
    micro_dict.update(place_dict)

    with open(micro_log_name, "a", buffering=1) as log_file:

        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S.%f")  
        print("{0:s}: Start the micro simulation for {1:s}-{2:s}".format(current_time, target_place, target_date_time))        
        log_file.write("{0:s}: Start the micro simulation for {1:s}-{2:s}\n".format(current_time, target_place, target_date_time))

        # preparing the directory for uning the OpenFOAM simulation
        prepare_Foam_directories(micro_dict)
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S.%f")  
        print("{0:s}: Prepared the FOAM directories for {1:s}-{2:s}".format(current_time, target_place, target_date_time))
        log_file.write("{0:s}: Prepared the FOAM directories for {1:s}-{2:s}\n".format(current_time, target_place, target_date_time))
        
        # running OpenFOAM simulation, including take the WRF simulation results to set up the initial fields and to run the coupled simulation
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S.%f")  
        print("{0:s}: Executing the OpenFOAM for {1:s}-{2:s}".format(current_time, target_place, target_date_time))
        log_file.write("{0:s}: Executing the OpenFOAM for {1:s}-{2:s}\n".format(current_time, target_place, target_date_time))     
        execute_foam(micro_dict)   

#######################################################################################################################
#   the main script to do the precice coupled simulation (Meso-to-Micro scales)
#######################################################################################################################

tut_json_name = "tut.json"

with open(tut_json_name, "r") as tut_json_file :
    tut_control = json.load(tut_json_file)

setup_preciceXML(tut_control)

micro_log_name = os.path.join(tut_control["case_directory"], "micro.log")
meso_log_name = os.path.join(tut_control["case_directory"], "meso.log")

WRF_thread = threading.Thread(target=execute_MesoSimulation, args=(tut_control, meso_log_name))
Foam_thread = threading.Thread(target=execute_MicroSimulation, args=(tut_control, micro_log_name))

WRF_thread.start()
Foam_thread.start()

WRF_thread.join()
Foam_thread.join()
