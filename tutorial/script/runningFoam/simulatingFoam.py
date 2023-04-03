import json
import os

# import from the PyFoam wrapper to run the OpenFOAM case

from PyFoam.Execution.ConvergenceRunner import ConvergenceRunner
from PyFoam.Execution.UtilityRunner import UtilityRunner
from PyFoam.Execution.ParallelExecution import LAMMachine
from PyFoam.LogAnalysis.BoundingLogAnalyzer import BoundingLogAnalyzer
from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile

#################################################################################
# The class use the pyFoam wrapper to run the CFD simulation with the mesh
# generated from preprocessed mesh utility. The simulation first copied the WRF
# output to set up the initial and boundary conditions, then run the simpleFoam/pisoFoam
# simulations to obtain the wind and temperature fields with finner meshes. Afterwards,
# the transportation of the scalars are simulated with the wind field produced by various
# openFOAM utilities
#################################################################################

class SimulatingFoam:

    def __init__(self, json_data):

        self.run_directory = json_data["case_directory"]
        self.template_directory = json_data["case_template"]
        self.mesh_directory = json_data["mesh_directory"]
        self.trans_directory = json_data["trans_directory"]

        self.set_patches = json_data["set_patches"]
        self.precice_interface = json_data["precice_interface"]
        self.number_process = json_data["number_process"]

        self.use_foam_mesh = json_data["foam_mesh"]
        self.min_mesh_size = json_data["min_mesh_size"]

        self.precice_data = json_data["precice_data"]
        self.simulation_time = json_data["case_time"]
        self.simulation_interval = json_data["time_step"]
        self.output_interval = json_data["output_interval"]

    def prepare_case_directory(self):
        '''
        this function is to prepare the run directory for the OpenFOAM simulation, which includes the sub-directories of
        system
        constant
        0
        the constant subdirctory contains a directory of domainData, which is used to store the results from WRF simulation
        to initialize the OpenFOAM run
        '''

        from shutil import rmtree, copy

        if os.path.exists(self.run_directory):
            rmtree(self.run_directory)

        os.makedirs(self.run_directory)
        os.makedirs(os.path.join(self.run_directory, "system"))
        os.makedirs(os.path.join(self.run_directory, "constant"))
        os.makedirs(os.path.join(self.run_directory, "constant", "domainData"))
        os.makedirs(os.path.join(self.run_directory, "0"))

        copy(os.path.join(self.template_directory, "controlDict.template"), os.path.join(self.run_directory, "system", "controlDict"))
        copy(os.path.join(self.template_directory, "fvSchemes.template"), os.path.join(self.run_directory, "system", "fvSchemes"))
        copy(os.path.join(self.template_directory, "fvSolution.template"), os.path.join(self.run_directory, "system", "fvSolution"))

        copy(os.path.join(self.template_directory, "transportProperties.template"), os.path.join(self.run_directory, "constant", "transportProperties"))
        copy(os.path.join(self.template_directory, "turbulenceProperties.template"), os.path.join(self.run_directory, "constant", "turbulenceProperties"))

        copy(os.path.join(self.template_directory, "U"), os.path.join(self.run_directory, "0"))
        copy(os.path.join(self.template_directory, "T"), os.path.join(self.run_directory, "0"))
        copy(os.path.join(self.template_directory, "p"), os.path.join(self.run_directory, "0"))
        copy(os.path.join(self.template_directory, "k"), os.path.join(self.run_directory, "0"))
        copy(os.path.join(self.template_directory, "epsilon"), os.path.join(self.run_directory, "0"))
        copy(os.path.join(self.template_directory, "nut"), os.path.join(self.run_directory, "0"))

    def setFvSchemes(self):
        '''
        this function is used to set up the fvScheme file in the system sub-direcotry from a template
        '''

        fvScheme_temp = os.path.join(self.template_directory, "fvSchemes.template")
        fvSchemes_dict = ParsedParameterFile(fvScheme_temp)

        ddtSchemes = fvSchemes_dict["ddtSchemes"]
        ddtSchemes["default"] = "CrankNicolson	0.8"

        divSchemes = fvSchemes_dict["divSchemes"]
        divSchemes["default"] = "Gauss upwind"
        divSchemes["div(phi,U)"] = "Gauss linearUpwind grad(U)"
        divSchemes["div(phi,T)"] = "Gauss linearUpwind grad(T)"
        divSchemes["div(phi,k)"] = "Gauss upwind"
        divSchemes["div(phi,epsilon)"] = "Gauss upwind"

        interpolateScheme = fvSchemes_dict["interpolationSchemes"]
        interpolateScheme["default"] = "linear"

        wallDist = fvSchemes_dict["wallDist"]
        wallDist["method"] = "none"

        fvSchemes_dict_name = os.path.join(self.run_directory, "system", "fvSchemes")
        fvSchemes_dict.writeFileAs(fvSchemes_dict_name)

    def setFvSolution(self):
       '''
        this function is used to set up the fvSolution file in the system sub-direcotry from a template
        '''
        fvSolution_temp = os.path.join(self.template_directory, "fvSolution.template")
        fvSolution_dict = ParsedParameterFile(fvSolution_temp)

        # ---------------- fvSolution solver settings -------------------

        solvers = fvSolution_dict["solvers"]

        p_solver = "GAMG"
        p_smoother = "GaussSeidel"

        U_etc_solver = "smoothSolver"
        U_etc_smoother = "symGaussSeidel"
        U_etc_preconditioner = "DILU"
        T_solver = "PBiCGStab"
        T_preconditioner = "DILU"
        T_tolerance = 1e-6
        T_relTol = 0
        U_etc_tolerance = 1e-5
        U_etc_relTol = 0
        p_tolerance = 1e-4
        p_relTol = 0.001

        solvers["p"] = {"solver": p_solver,
                        "smoother": p_smoother,
                        "tolerance": p_tolerance,
                        "relTol": p_relTol}

        solvers["pFinal"] = {"solver": p_solver,
                             "smoother": p_smoother,
                             "tolerance": p_tolerance*0.1,
                             "relTol": p_relTol*0.1}

        solvers["\"(U|k|omega|epsilon|R|nuTilda)\""] = {"solver": U_etc_solver,
                                                        "preconditioner": U_etc_preconditioner,
                                                        "smoother": U_etc_smoother,
                                                        "tolerance": U_etc_tolerance,
                                                        "relTol": U_etc_relTol}


        solvers["\"(U|k|omega|epsilon|R|nuTilda)Final\""] = {"solver": U_etc_solver,
                                                             "preconditioner": U_etc_preconditioner,
                                                             "smoother": U_etc_smoother,
                                                             "tolerance": U_etc_tolerance*0.1,
                                                             "relTol": 0.0}

        solvers["T"] = {"solver": T_solver,
                        "preconditioner": T_preconditioner,
                        "tolerance": T_tolerance,
                        "relTol": T_relTol}

    
        # ---------------- pimple algorithm settings  -------------------

        PIMPLE = fvSolution_dict["PIMPLE"]

        n_outter_corrects = 3
        n_corrects = 3
        n_nonorthogonal_corrects = 1

        PIMPLE["momentumPredictor"] = "yes"
        PIMPLE["correctPhi"] = "yes"
        PIMPLE["nCorrectors"] = n_corrects
        PIMPLE["nOuterCorrectors"] = n_outter_corrects
        PIMPLE["nNonOrthogonalCorrectors"] = n_nonorthogonal_corrects
        PIMPLE["consistent"] = "yes"

        # ---------------- relaxation settings  -------------------
        relaxation = fvSolution_dict["relaxationFactors"]
        p_relaxation = 0.5
        U_relaxation = 0.7
        turbulence_relaxation = 0.5

        relaxation["fields"] = {"p": p_relaxation}
        relaxation["equations"] = {"U": U_relaxation, "\"(k|omega|epsilon).*\"": turbulence_relaxation}

        fvSolution_dict_name = os.path.join(self.run_directory, "system", "fvSolution")
        fvSolution_dict.writeFileAs(fvSolution_dict_name)

    def setFvOptions(self, source_information):
       '''
        this function is used to set up the fvOption file in the system sub-direcotry from a template
        the fvOption is used to govern the scalarTransportFoam to run the pollutant disperion simulation
        '''
        fvOptions_dict_temp = os.path.join(self.case_dict_directory, "fvOptions.template")
        fvOptions_dict = ParsedParameterFile(fvOptions_dict_temp)

        pollution_source = fvOptions_dict["pollution_source"]

        pollution_source["type"] = "scalarFixedValueConstraint"
        pollution_source["active"] = True

        source_points = []
        source_intensities = []

        for ii in range(len(source_information)):
            ci = source_information[ii]
            source_points.append([ci[0], ci[1], ci[2]])
            source_intensities.append(("T", ci[3]))

        pollution_source["scalarFixedValueConstraintCoeffs"] = {
            "selectionMode": "points",
            "points": source_points,
            "mode": "uniform",
            "fieldValues": {"T": 1.0}
        }

        fvOptions_dict_name = os.path.join(self.case_directory, "constant", "fvOptions")
        fvOptions_dict.writeFileAs(fvOptions_dict_name)

    def setSamples(self, sample_locations):
        '''
        this function is used to set up the samples for measuring the simulated variables from a template
        '''

        sample_dict_temp = os.path.join(self.case_dict_directory, "sample.template")
        sample_dict = ParsedParameterFile(sample_dict_temp)

        sample_dict["type"] = "sets"
        sample_dict["libs"] = ["\"libsampling.so\""]
        sample_dict["interpolationScheme"] = "cellPoint"
        sample_dict["setFormat"] = "raw"
        sample_dict["fields"] = ["U", "T"]

        field_sample = {"type": "points",
                        "ordered": True,
                        "axis": "xyz",
                        "points": sample_locations}
        sample_dict["sets"] = [("sample_points", field_sample)]

        sample_dict_name = os.path.join(self.case_directory, "system", "sample")
        sample_dict.writeFileAs(sample_dict_name)

    def transMesh(self):
       '''
        this function is used to translate a mesh from other format (could be generated by Gmsh)
        '''

        if self.use_foam_mesh :
            from shutil import rmtree, copytree
            if os.path.exists(os.path.join(self.run_directory, "constant", "polyMesh")) :
                rmtree(os.path.join(self.run_directory, "constant", "polyMesh"))

            copytree(os.path.join(self.mesh_directory, "case.foam"), os.path.join(self.run_directory, "constant", "polyMesh"))

        else :

            trans_mesh_command = "gmshToFoam"
            mesh_name = os.path.join(self.mesh_directory, "case.msh")

            trans_mesh_logName = trans_mesh_command
            current_case = self.run_directory
            trans_mesh_run = UtilityRunner(argv=[trans_mesh_command, "-case", current_case, mesh_name], silent=True, logname=trans_mesh_logName)
            trans_mesh_run.start()

            transed_boundaries = ParsedParameterFile(
                                                        os.path.join(self.run_directory, "constant", "polyMesh", "boundary"), 
                                                        boundaryDict=True, 
                                                        treatBinaryAsASCII=True
                                                    )

            for ii in range(len(transed_boundaries.content)):
                current_boundary_name = transed_boundaries.content[ii]
                if current_boundary_name == "Bottom" or current_boundary_name == "terrains":
                    transed_boundaries[ii + 1]["type"] = "wall"

            transed_boundaries.writeFile()

    def setInitials(self):
        '''
        this function is used to set up the initial fields for the OpenFOAM run from the WRF simualted variables
        therefore, it constantly checks existances of the files generated from the WRF. Once the files containing 
        the WRF simulation results are presented in the target directory, this function use the data to set up the 
        initial condition of the OpenFOAM run
        '''
        from time import sleep
        from shlex import split
        from shutil import copy

        current_directory = self.run_directory
        trans_directory = self.trans_directory

        set_field_dict_temp = os.path.join(self.template_directory, "interpDict.template")
        set_field_dict = ParsedParameterFile(set_field_dict_temp)
        set_field_dict["nearestSwitch"] = False
        setFields = []

        if "Velocity" in  self.precice_data :
            setFields.append("U")

        if "TKE" in self.precice_data :
            setFields.append("k")

        if "DissipationRate" in self.precice_data :
            setFields.append("epsilon")

        if "Pressure" in self.precice_data :
            setFields.append("p")

        set_field_dict["setFields"] = setFields
        set_field_dict["setPatches"] = self.set_patches

        set_field_dict_name = os.path.join(current_directory, "system", "interpDict")
        set_field_dict.writeFileAs(set_field_dict_name)

        trans_point_name = os.path.join(trans_directory, "points")
        trans_index_name = os.path.join(trans_directory, "indices")
        trans_u_name = os.path.join(trans_directory, "inputU")
        trans_mx_name = os.path.join(trans_directory, "maxU")
        trans_p_name = os.path.join(trans_directory, "inputP")
        trans_k_name = os.path.join(trans_directory, "inputK")
        trans_epsilon_name = os.path.join(trans_directory, "inputEpsilon")

        input_point_name = os.path.join(current_directory, "constant", "domainData", "points")
        input_index_name = os.path.join(current_directory, "constant", "domainData", "indices")
        input_u_name = os.path.join(current_directory, "constant", "domainData", "inputU")
        input_p_name = os.path.join(current_directory, "constant", "domainData", "inputP")
        max_u_name = os.path.join(current_directory, "constant", "domainData", "maxU")
        input_k_name = os.path.join(current_directory, "constant", "domainData", "inputK")
        input_epsilon_name = os.path.join(current_directory, "constant", "domainData", "inputEpsilon")

        print("check for the existances of the initial field data in {0:s}".format(current_directory))

        while True :
            trans_exist = [os.path.exists(trans_point_name), os.path.exists(trans_index_name)]
    
            if "Velocity" in self.precice_data :
                trans_exist.append(os.path.exists(trans_u_name))
                trans_exist.append(os.path.expanduser(trans_mx_name))

            if "TKE" in self.precice_data :
                trans_exist.append(os.path.exists(trans_k_name))

            if "DissipationRate" in self.precice_data :
                trans_exist.append(os.path.exists(trans_epsilon_name))
            
            if "Pressure" in self.precice_data :
                trans_exist.append(os.path.exist(trans_p_name))

            if all(trans_exist) :
                
                copy(trans_point_name, input_point_name)
                copy(trans_index_name, input_index_name)
                
                if "Velocity" in self.precice_data :
                    copy(trans_u_name, input_u_name)
                    copy(trans_mx_name, max_u_name)

                if "TKE" in self.precice_data :
                    copy(trans_k_name, input_k_name)

                if "DissipationRate" in self.precice_data :
                    copy(trans_epsilon_name, input_epsilon_name)

                if "Pressure" in self.precice_data :
                    copy(trans_p_name, input_p_name)

                print("found the initial field data in {0:s}".format(current_directory))
                break

            else:
                       
                sleep(10)

        print("start to set initial fields")
        set_field_command = "trilinearInterpolateSetFields"
        set_field_logName = "setField"

        setField_run = UtilityRunner(argv=[set_field_command, "-case", current_directory], silent=True, logname=set_field_logName)
        setField_run.start()

        with open(max_u_name, "r") as max_u_file :
            u_data = split(max_u_file.readline(), ":")
            min_interval = 1.1 * (self.min_mesh_size / float(u_data[2]))
            
            if self.simulation_interval > min_interval :
                self.simulation_interval = min_interval

        return

    def setPrecice(self):
       '''
        this function set up the precice control file in the system sub-directory, the settings of the OpenFOAM precice
        can be found in https://precice.org/adapter-openfoam-config.html
        '''
        precice_dict_temp = os.path.join(self.template_directory, "preciceDict.template")
        precice_dict = ParsedParameterFile(precice_dict_temp)

        precice_dict["preciceConfig"] = os.path.join(self.run_directory, "..", "precice-config.xml")
        precice_dict["participant"] = "OpenFOAM"
        precice_dict["modules"] = ["FF"]

        coupling_interface = {
                                "WRF-OpenFOAM" : 
                                {
                                    "mesh" : "OpenFOAM-Mesh",
                                    "patches" : self.precice_interface,
                                    "locations" : "faceCenters",
                                    "readData" : self.precice_data,
                                    "writeData" : []
                                }
                             }
        precice_dict["interfaces"] = coupling_interface

        precice_dict_name = os.path.join(self.run_directory, "system", "preciceDict")
        precice_dict.writeFileAs(precice_dict_name)

    def runPimpleFoam(self):
        '''
        this function runs the pimpleFoam for the micro-scale simulation
        '''
        control_dict_temp = os.path.join(self.template_directory, "controlDict.template")
        control_dict = ParsedParameterFile(control_dict_temp)

        control_dict["application"] = "pimpleFoam"
        control_dict["startFrom"] = "startTime"
        control_dict["startTime"] = 0
        control_dict["stopAt"] = "endTime"
        control_dict["endTime"] = self.simulation_time
        control_dict["deltaT"] = self.simulation_interval
        control_dict["writeControl"] = "adjustableRunTime"
        control_dict["writeInterval"] = self.output_interval

        control_dict_name = os.path.join(self.run_directory, "system", "controlDict")
        control_dict.writeFileAs(control_dict_name)

        self.setFvSchemes()
        self.setFvSolution()
        self.setPrecice()

        pimpleFoam_command = "pimpleFoam"
        pimpleFoam_logName = "pimpleFoam"
        pimpleFoam_case = self.run_directory

        if self.number_process > 0:
            self.decomposeCase()
            pimpleFoam_lam = LAMMachine(nr=self.number_process)
            pimpleFoam_run = ConvergenceRunner(BoundingLogAnalyzer(),
                                               argv=[pimpleFoam_command, "-case", pimpleFoam_case], silent=True,
                                               logname=pimpleFoam_logName,
                                               lam=pimpleFoam_lam)
        else:
            pimpleFoam_run = ConvergenceRunner(BoundingLogAnalyzer(),
                                               argv=[pimpleFoam_command, "-case", pimpleFoam_case], silent=True,
                                               logname=pimpleFoam_logName)

        pimpleFoam_run.start()

    def decomposeCase(self):
        '''
        this function decompose the openFOAM case into parts for the parallel runs
        '''
        from math import sqrt
        decompose_dict_temp = os.path.join(self.template_directory, "decomposeParDict.template")
        decompose_dict = ParsedParameterFile(decompose_dict_temp)

        decompose_dict["numberOfSubdomains"] = self.number_process
        decompose_dict["method"] = "simple"
        decompose_dict["simpleCoeffs"] = {"n": [int(sqrt(self.number_process)), int(sqrt(self.number_process)), 1],
                                          "delta": 0.01}

        decompose_dict_name = os.path.join(self.run_directory, "system", "decomposeParDict")
        decompose_dict.writeFileAs(decompose_dict_name)

        decompose_command = "decomposePar"
        decompose_logName = "decompose"
        decompose_case = self.run_directory

        decompose_run = UtilityRunner(argv=[decompose_command, "-case", decompose_case], silent=True, logname=decompose_logName)
        decompose_run.start()

        return

    def reconstructCase(self):
        '''
        this function reconstruct the case from the decomposed case
        '''
        reconstruct_command = "reconstructPar"
        reconstruct_logName = "reconstruct"
        reconstruct_case = self.case_directory

        reconstruct_run = UtilityRunner(argv=[reconstruct_command, "-case", reconstruct_case], silent=True, logname=reconstruct_logName)
        reconstruct_run.start()

        return


if __name__ == "__main__":

    # this part is for debugging purpose, not used in the actual run of OpenFOAM

    test_dict = {
                    "case_directory": "/home/sunwei/Work/weatherForecast/test/foam_run",
                    "case_template": "/home/sunwei/Work/weatherForecast/etc/template/foam",
                    "mesh_directory": "/home/sunwei/Work/weatherForecast/data/geography/",
                    "set_patches": ["interface", "top"],
                    "precice_interface": ["interface"],
                    "precice_data" : ["pressure", "velocity"],
                    "time_step" : 0.001,
                    "case_time" : 30.0,
                    "number_process" : 4,
                    "output_interval" : 1.0
                }

    test_simulation = SimulatingFoam(test_dict)

    test_simulation.prepare_case_directory()
    test_simulation.transMesh(foamMesh=True)
    test_simulation.setPrecice()
    test_simulation.runPimpleFoam()
