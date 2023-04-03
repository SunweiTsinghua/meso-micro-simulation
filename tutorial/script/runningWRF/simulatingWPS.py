import os
import shutil

import f90nml
import subprocess

from netCDF4 import Dataset
from datetime import datetime

# adapted from the WRFpy suite to run WRF simulation controlled via the python interface

#################################################################################
# the class to use the WRFpy functionality to run the WPS simulation cases. The
# configurations are adjusted according to the json file generated from the web-frontend
# this class works with the simulateWRF class to run the whole set of the WRF simulation,
# this class is responsible for preprocessing
#################################################################################

class SimulatingWPS:

    def __init__(self, json_data):

        self.wps_directory = json_data["wps_directory"]
        self.run_directory = json_data["run_directory"]
        self.template_directory = json_data["template_directory"]

        self.met_data_directory = json_data["met_data_directory"]
        self.geo_directory = json_data["geo_directory"]

        self.start_dateTime = datetime.strptime(json_data["start_date_time"], "%Y-%m-%d_%H:%M:%S")
        self.end_dateTime = datetime.strptime(json_data["end_date_time"], "%Y-%m-%d_%H:%M:%S")

        self.center_longitude = json_data["center_longitude"]
        self.center_latitude = json_data["center_latitude"]
        self.domain_dimension = json_data["domain_dimension"]

        self._initialize()

    def _initialize(self):
        '''
        this function prepare and init the run of the WPS exectuable, which includes 
        1) removes the previous run files for WPS
        2) prepare the namelist for the WPS
        3) make soft links connecting the downloaded data to the boundary files known by WPS
        4) make soft links for the variable table
        5) make soft links for other tables
        '''
        self.clean_wps()
        self.prepare_namelist()
        self.link_boundary_files()
        self.link_vtable()
        self.link_tbl_files()

    def clean_wps(self):
        '''
        this function cleans up the directory containing files from the previous WPS run
        '''
        run_directory_list = os.listdir(self.run_directory)

        for current_item in run_directory_list:

            current_file = os.path.join(self.run_directory, current_item)

            if os.path.isfile(current_file) and current_file.find("log"):
                os.remove(current_file)

            if os.path.isfile(current_file) and current_file.find("err"):
                os.remove(current_file)

            if os.path.islink(current_file) and current_file.find("GRIBFILE"):
                os.remove(current_file)

            if os.path.isfile(current_file) and current_file.find("FILE:"):
                os.remove(current_file)

            if os.path.isfile(current_file) and current_file.find("PFILE:"):
                os.remove(current_file)

            if os.path.isfile(current_file) and current_file.find("PRES:"):
                os.remove(current_file)

            if os.path.isfile(current_file) and current_file.find("geo_em"):
                os.remove(current_file)

    def prepare_namelist(self):
        '''
        this function prepare the namelist for running the WPS exectuable, the namelist the main
        control file for running WPS and WRF. Using the f90nml libnrary of python, the namelist is
        organized as the dict, the meaning of the setting-ups of the namelist can be found in 
        https://blog.csdn.net/qq_40505953/article/details/124768253
        '''
        from datetime import datetime

        # read WPS namelist in WPS work_dir
        wps_name_list = f90nml.read(os.path.join(self.template_directory, "namelist.wps"))
        number_domains = wps_name_list["share"]["max_dom"]

        # share section
        starts = datetime.strftime(self.start_dateTime, "%Y-%m-%d_%H:%M:%S")
        ends = datetime.strftime(self.end_dateTime, "%Y-%m-%d_%H:%M:%S")

        wps_name_list["share"]["start_date"] = [starts for ii in range(number_domains)]
        wps_name_list["share"]["end_date"] = [ends for ii in range(number_domains)]
        wps_name_list["share"]["interval_seconds"] = 3600

        # domain settings
        wps_name_list["geogrid"]["parent_id"] = [ii for ii in range(number_domains)]
        wps_name_list["geogrid"]["parent_grid_ratio"] = [1] + [3 for ii in range(number_domains-1)]
        wps_name_list["geogrid"]["geog_data_res"] = ["default" for ii in range(number_domains)]

        # domain locations
        wps_name_list["geogrid"]["ref_lat"] = self.center_latitude
        wps_name_list["geogrid"]["ref_lon"] = self.center_longitude
        wps_name_list["geogrid"]["truelat1"] = self.center_latitude - 4.0
        wps_name_list["geogrid"]["truelat2"] = self.center_latitude + 4.0
        wps_name_list["geogrid"]["stand_lon"] = self.center_longitude
        wps_name_list["geogrid"]["geog_data_path"] = self.geo_directory

        if os.path.isfile(os.path.join(self.run_directory, "namelist.wps")):
            os.remove(os.path.join(self.run_directory, "namelist.wps"))

        wps_name_list.write(os.path.join(self.run_directory, "namelist.wps"))

    def link_boundary_files(self):
        '''
        this function make soft links of the boundary files downloaded from the website to provide
        the boundary conditions to run the WPS
        '''

        temporary_list = os.listdir(self.met_data_directory)
        boundary_file_list = []

        for current_item in temporary_list:
            current_file = os.path.join(self.met_data_directory, current_item)

            if os.path.isfile(current_file) :
                boundary_file_list.append(current_item)

        link_extensions = self._get_ext_list(len(boundary_file_list))

        for ii in range(len(boundary_file_list)):
            current_file = os.path.join(self.met_data_directory, boundary_file_list[ii])
            os.symlink(current_file, os.path.join(self.run_directory, "GRIBFILE." + link_extensions[ii]))

    def _get_ext_list(self, num):
        '''
        auxiliary function to generate the required the file name extensions, intented to be used only internally
        in this module
        '''

        from string import ascii_uppercase
        # create list of uppercase letters used linkname extension
        ext = [ascii_uppercase[idx] for idx in range(0, len(ascii_uppercase))]
        i1, i2, i3 = 0, 0, 0
        list_ext = []

        for filenum in range(num):  # loop over number of files
            # append extension to list (or create list for first iteration)

            try:
                list_ext.append(ext[i3] + ext[i2] + ext[i1])
            except NameError:
                list_ext = [ext[i3] + ext[i2] + ext[i1]]

            i1 += 1  # increment i1

            if i1 >= len(ascii_uppercase):
                i1 = 0
                i2 += 1  # increment i2
                if i2 >= len(ascii_uppercase):
                    i2 = 0
                    i3 += 1  # increment i3
                    if i3 >= len(ascii_uppercase):
                        message = 'Too many files to link'
                        # logger.error(message)
                        raise IOError(message)
        return list_ext

    def link_vtable(self):
        '''
        make soft links to the variable table in ungribing the met data downloaded
        '''
        if os.path.islink(os.path.join(self.run_directory, "Vtables")):
            os.remove(os.path.join(self.run_directory, "Vtables"))

        os.symlink(
                    os.path.join(self.template_directory, "Variable_Tables", "Vtable.ECMWF"), 
                    os.path.join(self.run_directory, "Vtable")
                  )

    def link_tbl_files(self):
        '''
        make soft links to other tables for running exectubales of geogrid, metgrid
        '''
        # geogrid
        if os.path.exists(os.path.join(self.run_directory, "geogrid")):
            shutil.rmtree(os.path.join(self.run_directory, "geogrid"), ignore_errors=True)

        os.makedirs(os.path.join(self.run_directory, "geogrid"))

        if not os.path.isfile(os.path.join(self.template_directory, "GeoGrid", "geogrid.tbl")):
            geogrid_table = os.path.join(self.template_directory, "GeoGridTbls", "GEOGRID.TBL.ARW")
        else:
            geogrid_table = os.path.join(self.template_directory, "GeoGridTbls", "geogrid.tbl")

        os.symlink(geogrid_table, os.path.join(self.run_directory, "geogrid", "GEOGRID.TBL"))

        # metgrid
        if os.path.exists(os.path.join(self.run_directory, "metgrid")):
            shutil.rmtree(os.path.join(self.run_directory, "metgrid"), ignore_errors=True)

        os.makedirs(os.path.join(self.run_directory, "metgrid"))

        if not os.path.isfile(os.path.join(self.template_directory, "MetGrid", "METGRID.TBL")):
            metgrid_table = os.path.join(self.template_directory, "MetGridTbls", "METGRID.TBL.ARW")
        else:
            metgrid_table = os.path.join(self.template_directory, "MetGridTbls", "metgrid.tbl")

        os.symlink(metgrid_table, os.path.join(self.run_directory, "metgrid", "METGRID.TBL"))

    def run_geogrid(self):
        '''
        this function uses the subprocess library of python to run the geogrid.exe
        '''
        geogrid_command = os.path.join(self.wps_directory, "geogrid", "geogrid.exe")

        try:
            geogrid_process = subprocess.run(
                                                geogrid_command,
                                                cwd=self.run_directory,
                                                check=True,
                                                stdout=subprocess.PIPE,
                                                stderr=subprocess.PIPE,
                                                universal_newlines=True
                                            )

            geogrid_out = geogrid_process.stdout
            geogrid_err = geogrid_process.stderr

            with open(os.path.join(self.run_directory, "wps.log"), "a") as log_file:
                log_file.write("--------------------------------------------------\n")
                log_file.write("Starts the output from the geogrid.exe\n")
                log_file.write("--------------------------------------------------\n")
                log_file.writelines(geogrid_out)
                log_file.write(" \n")

            with open(os.path.join(self.run_directory, "wps.err"), "a") as err_file:
                err_file.write("--------------------------------------------------\n")
                err_file.write("Starts the error logs from the geogrid.exe\n")
                err_file.write("--------------------------------------------------\n")
                err_file.writelines(geogrid_err)
                err_file.write(" \n")

        except subprocess.CalledProcessError as err:
            raise ValueError(
                "failed in running geogrid with the exit code of " + str(err.returncode))

    def run_ungrib(self):
        '''
        this function uses the subprocess library of python to run the ungrib.exe
        '''
        ungrib_command = os.path.join(self.wps_directory, "ungrib", "ungrib.exe")

        try:
            ungrib_process = subprocess.run(ungrib_command,
                                            cwd=self.run_directory,
                                            check=True,
                                            stdout=subprocess.PIPE,
                                            stderr=subprocess.PIPE,
                                            universal_newlines=True)

            ungrib_out = ungrib_process.stdout
            ungrib_err = ungrib_process.stderr

            with open(os.path.join(self.run_directory, "wps.log"), "a") as log_file:
                log_file.write("--------------------------------------------------\n")
                log_file.write("Starts the output from the ungrib.exe\n")
                log_file.write("--------------------------------------------------\n")
                log_file.writelines(ungrib_out)
                log_file.write(" \n")

            with open(os.path.join(self.run_directory, "wps.err"), "a") as err_file:
                err_file.write("--------------------------------------------------\n")
                err_file.write("Starts the error logs from the ungrib.exe\n")
                err_file.write("--------------------------------------------------\n")
                err_file.writelines(ungrib_err)
                err_file.write(" \n")

        except subprocess.CalledProcessError as err:
            raise ValueError("failed in running ungrib with the exit code of " + str(err.returncode))

    def run_metgrid(self):
        '''
        this function uses the subprocess library of python to run the metgrid.exe
        '''
        metgrid_command = os.path.join(self.wps_directory, "metgrid", "metgrid.exe")

        try:
            metgrid_process = subprocess.run(metgrid_command,
                                             cwd=self.run_directory,
                                             check=True,
                                             stdout=subprocess.PIPE,
                                             stderr=subprocess.PIPE,
                                             universal_newlines=True)

            metgrid_out = metgrid_process.stdout
            metgrid_err = metgrid_process.stderr

            with open(os.path.join(self.run_directory, "wps.log"), "a") as log_file:
                log_file.write("--------------------------------------------------\n")
                log_file.write("Starts the output from the metgrid.exe\n")
                log_file.write("--------------------------------------------------\n")
                log_file.writelines(metgrid_out)
                log_file.write(" \n")

            with open(os.path.join(self.run_directory, "wps.err"), "a") as err_file:
                err_file.write("--------------------------------------------------\n")
                err_file.write("Starts the error logs from the metgrid.exe\n")
                err_file.write("--------------------------------------------------\n")
                err_file.writelines(metgrid_err)
                err_file.write(" \n")

        except subprocess.CalledProcessError as err:
            raise ValueError("failed in running metgrid with exit code of " + str(err.returncode))


if __name__ == "__main__":

    # this part is for debugging purpose, not used in the actual run to prepare for the main WRF run

    test_dict = {
                    "wps_directory" : "/home/sunwei/Apps/WRF/WPS",
                    "run_directory" : "/home/sunwei/Work/weatherForecast/test/WPS-run",
                    "template_directory" : "/home/sunwei/Work/weatherForecast/etc/template/WRF",
                    "met_data_directory" : "/home/sunwei/Work/weatherForecast/test/met-data",
                    "geo_directory" : "/home/sunwei/Apps/WRF/Data/WPS_GEOG",

                    "start_date_time" : "2021-8-11_00:00:00",
                    "end_date_time" : "2021-08-11_01:00:00",
                    "center_longitude" : 88.875,
                    "center_latitude" : 42.591,
                    "domain_dimension" : [6000.0, 6000.0]
                }

    test_wps = SimulatingWPS(test_dict)
    test_wps.run_geogrid()
    test_wps.run_ungrib()
    test_wps.run_metgrid()
