import os
import f90nml
import subprocess

from datetime import datetime
# adapted from the WRFpy suite to run WRF simulation controlled via the python interface

#################################################################################
# the class to use the WRFpy functionality to run the WRF simulation cases. The
# configurations are adjusted according to the json file generated from the web-frontend
# the purpose is to run the WRF simulation focusing on the core area of the urban
# environment selected by the web user. The WRF source code is adjusted to output
# the simulation results at every time steps in a small area. The simulation of OpenFOAM
# is then use the output to specify the boundary and initial conditions
#################################################################################

class SimulatingWRF:

    def __init__(self, json_data):

        self.run_directory = json_data["run_directory"]
        self.wrf_directory = json_data["wrf_directory"]
        self.wps_directory = json_data["wps_directory"]
        self.template_directory = json_data["template_directory"]

        self.start_dateTime = datetime.strptime(json_data["start_date_time"], "%Y-%m-%d_%H:%M:%S")
        self.end_dateTime = datetime.strptime(json_data["end_date_time"], "%Y-%m-%d_%H:%M:%S")
        
        self.coarse_grids = json_data["coarse_grids"] 
        self.precice_xml_name = json_data["precice_xml_name"]
        self.participant_name = json_data["precice_participant_name"]
        self.mesh_name = json_data["precice_mesh_name"]

        self.exchange_data = json_data["precice_exchange_data"]
        self.precice_lats = json_data["precice_lats"]
        self.precice_lons = json_data["precice_lons"]
        self.precice_cornerx = json_data["precice_cornerx"]
        self.precice_cornery = json_data["precice_cornery"]
        self.precice_height = json_data["precice_height"]
        self.precice_start_mins = json_data["precice_start_mins"]
        self.precice_Cmu = json_data["precice_Cmu"]
        self.write_trans = json_data["precice_write_trans"]
        self.trans_directory = json_data["precice_trans_directory"]
        self.precice_lid = json_data["precice_lid"]

    def initialize(self):
        '''
        this function prepare and init the run of the WPS exectuable, which includes
        1) removes the files associated with a previous run
        2) prepare the namelists for the WRF run, which includes the namelist for a coarse grid run and fine grid run with LES
        3) prepare the run directory of the WRF simulation, which is the folder of runnning and result files
        '''
        self.cleanup_previous_runs()
        self.prepare_namelist()
        self.prepare_wrf_run_directory()

    def cleanup_previous_runs(self):
        '''
        this function is to remove the files in the run directory from previou srun
        '''
        run_directory_list = os.listdir(self.run_directory)

        for current_item in run_directory_list:
            current_file = os.path.join(self.run_directory, current_item)

            if os.path.islink(current_file) and "met_em" in current_file:
                os.remove(current_file)

            if os.path.isfile(current_file) and "wrfout" in current_file:
                os.remove(current_file)

            if os.path.isfile(current_file) and "wrfbdy" in current_file:
                os.remove(current_file)

            if os.path.isfile(current_file) and "wrfinput" in current_file:
                os.remove(current_file)

            if os.path.isfile(current_file) and "namelist" in current_file:
                os.remove(current_file)

            if os.path.isfile(current_file) and "precice" in current_file:
                os.remove(current_file)

            if os.path.isfile(current_file) and "rsl" in current_file:
                os.remove(current_file)

            if os.path.isfile(current_file) and "wrfndi" in current_file:
                os.remove(current_file)

    def namelist_time_settings(self, name_list, start_dateTime, end_dateTime, number_domains):
        '''
        this function set up the time controls in the namelist correspnding to both coarse gird
        and fine grid runs
        '''
        name_list["time_control"]["start_year"] = [start_dateTime.year for ii in range(number_domains)]
        name_list["time_control"]["start_month"] = [start_dateTime.month for ii in range(number_domains)]
        name_list["time_control"]["start_day"] = [start_dateTime.day for ii in range(number_domains)]
        name_list["time_control"]["start_hour"] = [start_dateTime.hour for ii in range(number_domains)]
        name_list["time_control"]["start_minute"] = [start_dateTime.minute for ii in range(number_domains)]
        name_list["time_control"]["start_second"] = [start_dateTime.second for ii in range(number_domains)]

        name_list["time_control"]["end_year"] = [end_dateTime.year for ii in range(number_domains)]
        name_list["time_control"]["end_month"] = [end_dateTime.month for ii in range(number_domains)]
        name_list["time_control"]["end_day"] = [end_dateTime.day for ii in range(number_domains)]
        name_list["time_control"]["end_hour"] = [end_dateTime.hour for ii in range(number_domains)]
        name_list["time_control"]["end_minute"] = [end_dateTime.minute for ii in range(number_domains)]
        name_list["time_control"]["end_second"] = [end_dateTime.second for ii in range(number_domains)]

        time_difference = end_dateTime - start_dateTime

        delta_days = time_difference.days
        delta_hours = int((time_difference.total_seconds() - delta_days * 24 * 3600) / 3600)
        delta_minutes = int((time_difference.total_seconds() - delta_days * 24 * 3600 - delta_hours * 3600) / 60)
        delta_seconds = int(time_difference.total_seconds() - delta_days * 24 * 3600 - delta_hours * 3600 - delta_minutes * 60)

        name_list["time_control"]["run_days"] = delta_days
        name_list["time_control"]["run_hours"] = delta_hours
        name_list["time_control"]["run_minutes"] = delta_minutes
        name_list["time_control"]["run_seconds"] = delta_seconds

    def prepare_namelist(self):
        '''
        this function prepare the namelist for the two sets of WRF runs, including a file for the coarse grid run with namelist.coarse
        and a file for the fine grid run with namelist.fine. in addition, the namelist governs the down scaling from the coarse grid
        to the fine grid is generated with namelist.ndown
        '''

        wrf_coarse_name_list = f90nml.read(os.path.join(self.template_directory, "namelist.coarse"))
        number_domains = wrf_coarse_name_list["domains"]["max_dom"]
        self.namelist_time_settings(wrf_coarse_name_list, self.start_dateTime, self.end_dateTime, number_domains)

        wrf_coarse_name_list["time_control"]["history_interval"][-1] = 1
        wrf_coarse_name_list["time_control"]["frames_per_outfile"][-1] = 6000

        wrf_coarse_name_list.write(os.path.join(self.run_directory, "namelist.coarse"), force=True)

        ndown_name_list = f90nml.read(os.path.join(self.template_directory, "namelist.ndown"))
        self.namelist_time_settings(ndown_name_list, self.start_dateTime, self.end_dateTime, 2)

        ndown_name_list["time_control"]["io_form_auxinput2"] = 2
        ndown_name_list["time_control"]["interval_seconds"] = 60
        ndown_name_list["domains"]["time_step"] = 5

        ndown_name_list.write(os.path.join(self.run_directory,"namelist.ndown"), force=True)

        wrf_fine_name_list = f90nml.read(os.path.join(self.template_directory, "namelist.fine"))
        number_domains = wrf_fine_name_list["domains"]["max_dom"]
        self.namelist_time_settings(wrf_fine_name_list, self.start_dateTime, self.end_dateTime, number_domains)

        wrf_fine_name_list["time_control"]["interval_seconds"] = 60
        wrf_fine_name_list.write(os.path.join(self.run_directory, "namelist.fine"), force=True)

    def prepare_precice(self):
        '''
        prepare the namelist for the precice coupling with the name of precice.nml. the setting of the precice
        namelist is organized as a dict
        '''

        if os.path.exists(os.path.join(self.run_directory, "precice.nml")):
            os.remove(os.path.join(self.run_directory, "precice.nml"))

        precice_ctrl = {
                            "precice_config" : 
                            {
                                "xml_name" : self.precice_xml_name,
                                "participant_name" : self.participant_name,
                                "mesh_name" : self.mesh_name,
                                "num_data" : len(self.exchange_data),
                                "write_init" : self.write_trans,
                                "init_directory" : self.trans_directory
                            }
                        }
        coupling_ctrl = {
                            "coupling_config":
                            {
                                "exchange_data" : self.exchange_data,
                                "lats" : self.precice_lats,
                                "lons" : self.precice_lons,
                                "cornerx" : self.precice_cornerx,
                                "cornery" : self.precice_cornery,
                                "height" : self.precice_height,
                                "start_mins" : self.precice_start_mins,
                                "Cmu" : self.precice_Cmu,
                                "has_lid" : self.precice_lid
                            }
                       }

        precice_namelist = f90nml.Namelist(precice_ctrl)
        coupling_namelist = f90nml.Namelist(coupling_ctrl)

        with open(os.path.join(self.run_directory, "precice.nml"), "w") as precice_file:
            precice_file.write("# --------- precice coupling control file --------------#\n")
            f90nml.write(precice_namelist, precice_file)
            f90nml.write(coupling_namelist, precice_file)

    def run_real(self, kind="coarse"):
        '''
        this function is used to run the real.exe to prepare the boundary and initial conditions for the run
        of WRF, the default is to run the real.exe for the coarse grid
        '''

        if os.path.islink(os.path.join(self.run_directory, "namelist.input")):
            os.remove(os.path.join(self.run_directory, "namelist.input"))

        target_name_list = os.path.join(self.run_directory, "namelist." + kind)
        os.symlink(target_name_list, os.path.join(self.run_directory, "namelist.input"))

        run_directory_list = os.listdir(self.run_directory)

        for current_item in run_directory_list:
            current_file = os.path.join(self.run_directory, current_item)

            if os.path.islink(current_file) and "met_em" in current_file:
                os.remove(current_file)

        wps_directory_list = os.listdir(self.wps_directory)

        for current_item in wps_directory_list:
            current_file = os.path.join(self.wps_directory, current_item)

            if os.path.isfile(current_file) and "met_em" in current_file:

                if kind == "coarse" and int(current_item[9]) <= self.coarse_grids:
                    os.symlink(current_file, os.path.join(self.run_directory, current_item))

                elif kind == "ndown" and int(current_item[9]) >= self.coarse_grids:
                    old_domain_id = current_item[7:10]
                    new_domain_id = "d{0:02d}".format(int(current_item[9]) - self.coarse_grids + 1)
                    symbol_item = current_item.replace(old_domain_id, new_domain_id)
                    os.symlink(current_file, os.path.join(self.run_directory, symbol_item))

                elif kind == "fine" and int(current_item[9]) > self.coarse_grids:
                    old_domain_id = current_item[7:10]
                    new_domain_id = "d{0:02d}".format(int(current_item[9]) - self.coarse_grids)
                    symbol_item = current_item.replace(old_domain_id, new_domain_id)
                    os.symlink(current_file, os.path.join(self.run_directory, symbol_item))

        real_command = os.path.join(self.wrf_directory, "main", "real.exe")

        try:
            real_process = subprocess.run(real_command,
                                          cwd=self.run_directory,
                                          check=True,
                                          stdout=subprocess.PIPE,
                                          stderr=subprocess.PIPE,
                                          universal_newlines=True)

        except subprocess.CalledProcessError as err:
            raise ValueError("failed in running real.exe with the exit code of " + str(err.returncode))

        os.rename(os.path.join(self.run_directory, "rsl.error.0000"), os.path.join(self.run_directory, "rsl.error.real." + kind))
        os.rename(os.path.join(self.run_directory, "rsl.out.0000"), os.path.join(self.run_directory, "rsl.out.real." + kind))

        if kind == "fine":
            os.remove(os.path.join(self.run_directory, "wrfinput_d01"))
            os.remove(os.path.join(self.run_directory, "wrfbdy_d01"))

            os.rename(os.path.join(self.run_directory, "wrfinput_d01.ndown"),
                      os.path.join(self.run_directory, "wrfinput_d01"))
            os.rename(os.path.join(self.run_directory, "wrfbdy_d01.ndown"),
                      os.path.join(self.run_directory, "wrfbdy_d01"))

    def link_run_file(self, current_file):
        '''
        make the soft links of various files in the run directory of the WRF installation to the current directory.
        this function is for internal use only
        '''

        if os.path.islink(os.path.join(self.run_directory, current_file)):
            os.remove(os.path.join(self.run_directory, current_file))

        os.symlink(os.path.join(self.template_directory, "run", current_file),
                   os.path.join(self.run_directory, current_file))

    def prepare_wrf_run_directory(self):
        '''
        make the soft links of various files in the run directory of the WRF installation to the current directory.
        this function is for internal use only
        '''
        self.link_run_file("LANDUSE.TBL")
        self.link_run_file("ozone.formatted")
        self.link_run_file("ozone_lat.formatted")
        self.link_run_file("ozone_plev.formatted")
        self.link_run_file("RRTMG_LW_DATA")
        self.link_run_file("RRTMG_LW_DATA_DBL")
        self.link_run_file("RRTMG_SW_DATA")
        self.link_run_file("RRTMG_SW_DATA_DBL")
        self.link_run_file("RRTM_DATA")
        self.link_run_file("RRTM_DATA_DBL")
        self.link_run_file("VEGPARM.TBL")
        self.link_run_file("SOILPARM.TBL")
        self.link_run_file("SOILPARM.TBL_Kishne_2017")
        self.link_run_file("GENPARM.TBL")
        self.link_run_file("URBPARM.TBL")
        self.link_run_file("CAMtr_volume_mixing_ratio")

    def run_wrf(self, kind="coarse"):
        '''
        this function is used to run the wrf.exe 
        the default is to run the real.exe for the coarse grid
        '''

        if os.path.islink(os.path.join(self.run_directory, "namelist.input")):
            os.remove(os.path.join(self.run_directory, "namelist.input"))

        os.symlink(os.path.join(self.run_directory, "namelist." + kind),
                   os.path.join(self.run_directory, "namelist.input"))

        wrf_command = os.path.join(self.wrf_directory, "main", "wrf.exe")

        try:
            wrf_process = subprocess.run(wrf_command,
                                         cwd=self.run_directory,
                                         check=True,
                                         stdout=subprocess.PIPE,
                                         stderr=subprocess.PIPE,
                                         shell=True,
                                         universal_newlines=True)

        except subprocess.CalledProcessError as err:
            raise ValueError("failed in running wrf.exe with the exit code of " + str(err.returncode))

        os.rename(os.path.join(self.run_directory, "rsl.error.0000"), os.path.join(self.run_directory, "rsl.error.wrf." + kind))
        os.rename(os.path.join(self.run_directory, "rsl.out.0000"), os.path.join(self.run_directory, "rsl.out.wrf." + kind))            

    def run_ndown(self):
        '''
        this function is to run the ndown.exe to port the results from a coarse grid of WRF run to the fine grid
        '''
        dn_out_name = "wrfout_d{0:02d}".format(self.coarse_grids)

        current_list = os.listdir(self.run_directory)

        for current_item in current_list:
            current_file = os.path.join(self.run_directory, current_item)

            if os.path.islink(current_file) and "met_em" in current_file:
                os.remove(current_file)

            if os.path.isfile(current_file) and "wrfout" in current_file:
                backup_file = current_file + ".coarse"
                os.rename(current_file, backup_file)

            if os.path.isfile(current_file) and "wrfinput" in current_file:
                backup_file = current_file + ".coarse"
                os.rename(current_file, backup_file)

            if os.path.isfile(current_file) and "wrfbdy" in current_file:
                backup_file = current_file + ".coarse"
                os.rename(current_file, backup_file)

        new_list = os.listdir(self.run_directory)

        for current_item in new_list:
            current_file = os.path.join(self.run_directory, current_item)

            if os.path.isfile(current_file) and dn_out_name in current_file:
                interim_file = current_file.replace("wrfout_d{0:02d}".format(self.coarse_grids), "wrfout_d01")
                replace_file = interim_file.replace(".coarse", "")
                os.rename(current_file, replace_file)

        self.run_real(kind="ndown")

        os.rename(os.path.join(self.run_directory, "wrfinput_d02"), os.path.join(self.run_directory, "wrfndi_d02"))

        if os.path.islink(os.path.join(self.run_directory, "namelist.input")):
            os.remove(os.path.join(self.run_directory, "namelist.input"))

        os.symlink(os.path.join(self.run_directory, "namelist.ndown.main"),
                   os.path.join(self.run_directory, "namelist.input"))

        ndown_command = os.path.join(self.wrf_directory, "main", "ndown.exe")

        try:
            ndown_process = subprocess.run(ndown_command,
                                           cwd=self.run_directory,
                                           check=True,
                                           stdout=subprocess.PIPE,
                                           stderr=subprocess.PIPE,
                                           universal_newlines=True)

            ndown_out = ndown_process.stdout
            ndown_err = ndown_process.stderr

            with open(self.run_directory + "/ndown.log", "w") as log_file:
                log_file.write("--------------------------------------------------\n")
                log_file.write("Starts the output from the ndown.exe\n")
                log_file.write("--------------------------------------------------\n")
                log_file.writelines(ndown_out)
                log_file.write(" \n")

            with open(self.run_directory + "/ndown.err", "w") as err_file:
                err_file.write("--------------------------------------------------\n")
                err_file.write("Starts the error logs from the ndown.exe\n")
                err_file.write("--------------------------------------------------\n")
                err_file.writelines(ndown_err)
                err_file.write(" \n")

        except subprocess.CalledProcessError as err:
            raise ValueError("failed in running ndown.exe with the exit code of " + str(err.returncode))

        os.rename(os.path.join(self.run_directory, "rsl.error.0000"), os.path.join(self.run_directory, "rsl.error.ndown"))
        os.rename(os.path.join(self.run_directory, "rsl.out.0000"), os.path.join(self.run_directory, "rsl.out.ndown"))            

        os.rename(os.path.join(self.run_directory, "wrfinput_d02"), os.path.join(self.run_directory, "wrfinput_d01.ndown"))
        os.rename(os.path.join(self.run_directory, "wrfbdy_d02"), os.path.join(self.run_directory, "wrfbdy_d01.ndown"))

if __name__ == "__main__":

    # this part is for debugging purpose, not used in the actual run of WRF
    
    test_dict = {
                    "run_directory" : "/home/sunwei/Work/weatherForecast/test/WRF-run",
                    "wrf_directory" : "/home/sunwei/Apps/WRF/WRF",
                    "wps_directory" : "/home/sunwei/Work/weatherForecast/test/WPS-run",
                    "template_directory" : "/home/sunwei/Work/weatherForecast/etc/template/WRF",

                    "coarse_grids" : 3,
                    
                    "start_date_time" : "2021-8-11_00:00:00",
                    "end_date_time" : "2021-08-11_01:00:00",

                    "precice_xml_name": "/home/sunwei/Work/weatherForecast/test/precice-config.xml",
                    "precice_participant_name": "WRF",
                    "precice_mesh_name" : "WRF-Mesh",
                    "precice_lats" : [42.560, 42.622],
                    "precice_lons" : [88.827, 88.923],
                    "precice_exchange_data" : ["pressure", "velocity", "TKE", "DissipationRate"],
                    "precice_cornerx": [-3000.0, 3000.0],
                    "precice_cornery": [-3000.0, 3000.0],
                    "precice_height" : 5000.0,
                    "precice_start_mins" : 10.0,
                    "precice_Cmu": 0.09,
                    "precice_write_trans": True,
                    "precice_trans_directory" : "/home/sunwei/Work/weatherForecast/test/foam-run/constant/domainData",
                    "precice_lid" : False
                }

    test_wrf = SimulatingWRF(test_dict)
    test_wrf.initialize()
    test_wrf.run_real(kind="coarse")
    test_wrf.run_wrf(kind="coarse")
    test_wrf.run_ndown()
    test_wrf.run_real(kind="fine")
    test_wrf.run_wrf(kind="fine")
