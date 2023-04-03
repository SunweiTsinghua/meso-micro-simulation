import cdsapi
import os
from datetime import datetime
#################################################################################
# the class uses the python script suggested by the ERA5 website to download the
# meteorology data from the ERA5 database (hourly reanalysis data from 1979 to
# present) at the pressure levels below 100mb, the variables to download includes
# teh geopotential height, relative humidity, temperature and all three components
# of the wind velocity. The downloaded met data is primarily used to initialize
# a WRF run providing the boundary conditions for the CFD simulation
#################################################################################

class DownloadingMetData:

    def __init__(self, control_Dict):

        self.start_dateTime = datetime.strptime(control_Dict["start_date_time"], "%Y-%m-%d_%H:%M:%S")
        self.end_dateTime = datetime.strptime(control_Dict["end_date_time"], "%Y-%m-%d_%H:%M:%S")
        self.area = control_Dict["area"]
        self.output_directory = control_Dict["out_directory"]

    def retrieve_era5_data(self):
        '''
        the main function retrieve the ERA5 data from the website database, the code is rrevised from the snippet generated from
        the website, the data downloaded covers the duration with resultion of 1 hour and all available pressure levels, the data
        includes pressure levels (vertical profile) and single level (at the ground surface)

        pressure-level: geopotential(geopotential height), relative humidity, specific_humidity,
                        temperature, u_component_of_wind, v_component_of_wind, vertical_velocity(w_component_of_wind)

        single-level: 10m_u_component_of_wind, 10m_v_component_of_wind, 2m_dewpoint_temperature, 2m_temperature, land_sea_mask, 
                      mean_sea_level_pressure, sea_ice_cover, sea_surface_temperature, skin_temperature (gournd surface temperature),
                      snow_depth, surface_pressure, soil_temperature_level_1, soil_temperature_level_2, soil_temperature_level_3, 
                      soil_temperature_level_4 (these are soil temperatures at different depths),
                      volumetric_soil_water_layer_1, volumetric_soil_water_layer_2, volumetric_soil_water_layer_3, 
                      volumetric_soil_water_layer_4 (these are volumetric of soil waters at different depths)
        '''
        import calendar

        for current_year in range(self.start_dateTime.year, self.end_dateTime.year + 1):

            if current_year == self.start_dateTime.year:
                start_month = self.start_dateTime.month
            else:
                start_month = 1

            if current_year == self.end_dateTime.year:
                end_month = self.end_dateTime.month
            else:
                end_month = 12

            for current_month in range(start_month, end_month + 1):

                if current_month == self.start_dateTime.month and current_year == self.start_dateTime.year:
                    start_day = self.start_dateTime.day
                else:
                    start_day = 1

                if current_month == self.end_dateTime.month and current_year == self.end_dateTime.year:
                    end_day = self.end_dateTime.day
                else:
                    end_day = calendar.monthrange(current_year, current_month)[1]

                single_level_file = os.path.join(self.output_directory, "{0:04d}-{1:02d}-single_level.grib".format(current_year, current_month))
                pressure_level_file = os.path.join(self.output_directory, "{0:04d}-{1:02d}-pressure_levels.grib".format(current_year, current_month))

                era_client = cdsapi.Client()
                retrieve_data_dict = {
                                        'product_type': 'reanalysis', 'format': 'grib',
                                        'year': [str(current_year)], 'month': [str(current_month)],
                                        'day': [str(current_day) for current_day in range(start_day, end_day+1)],
                                        'time': ['00:00', '01:00', '02:00',
                                                 '03:00', '04:00', '05:00',
                                                 '06:00', '07:00', '08:00',
                                                 '09:00', '10:00', '11:00',
                                                 '12:00', '13:00', '14:00',
                                                 '15:00', '16:00', '17:00',
                                                 '18:00', '19:00', '20:00',
                                                 '21:00', '22:00', '23:00']
                                        }

                retrieve_data_dict['variable'] = ['10m_u_component_of_wind', '10m_v_component_of_wind',
                                                   '2m_dewpoint_temperature', '2m_temperature', 
                                                   'land_sea_mask', 'mean_sea_level_pressure',
                                                   'sea_ice_cover', 'sea_surface_temperature', 'skin_temperature',
                                                   'snow_depth',  'surface_pressure',
                                                   'soil_temperature_level_1', 'soil_temperature_level_2',
                                                   'soil_temperature_level_3', 'soil_temperature_level_4',
                                                   'volumetric_soil_water_layer_1', 'volumetric_soil_water_layer_2',
                                                   'volumetric_soil_water_layer_3', 'volumetric_soil_water_layer_4']

                retrieve_data_dict['area'] = self.area
                era_client.retrieve('reanalysis-era5-single-levels', retrieve_data_dict, single_level_file)

                era_client = cdsapi.Client()
                retrieve_data_dict = {
                                        'product_type': 'reanalysis', 'format': 'grib',
                                        'year': [str(current_year)], 'month': [str(current_month)],
                                        'day': [str(current_day) for current_day in range(start_day, end_day+1)],
                                        'time': ['00:00', '01:00', '02:00',
                                                 '03:00', '04:00', '05:00',
                                                 '06:00', '07:00', '08:00',
                                                 '09:00', '10:00', '11:00',
                                                 '12:00', '13:00', '14:00',
                                                 '15:00', '16:00', '17:00',
                                                 '18:00', '19:00', '20:00',
                                                 '21:00', '22:00', '23:00']
                                        }

                retrieve_data_dict['variable'] = ['geopotential', 'relative_humidity', 'specific_humidity',
                                                    'temperature', 'u_component_of_wind', 'v_component_of_wind', 'vertical_velocity']

                retrieve_data_dict['pressure_level'] = ['5', '7', '10', '20', '30', '50', '70', '100', '125',
                                                        '150', '175', '200', '225', '250', '300', '350', '400', '450',
                                                        '500', '550', '600', '650', '700', '750', '775', '800', '825',
                                                        '850', '875', '900', '925', '950', '975', '1000']

                retrieve_data_dict['area'] = self.area
                era_client.retrieve('reanalysis-era5-pressure-levels', retrieve_data_dict, pressure_level_file)


if __name__ == "__main__":

    # this part is for debugging purpose, not used in the actual run to download the data

    test_dict = {
                    "start_date_time" : "2021-8-11_00:00:00",
                    "end_date_time" : "2021-08-11_23:00:00",
                    "area" : [43.591,87.875,41.591,89.875],
                    "out_directory" : "/home/sunwei/Work/weatherForecast/test/met-data"
                }

    test_download = DownloadingMetData(test_dict)
    test_download.retrieve_era5_data()
