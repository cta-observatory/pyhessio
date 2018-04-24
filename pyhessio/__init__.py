# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""

"""

import numpy as np
import os
import ctypes
from contextlib import contextmanager

__all__ = ['HessioError', 'HessioChannelIndexError',
           'HessioTelescopeIndexError', 'HessioGeneralError',
           'HessioFile', 'open_hessio', 'close_file']

__version__ = '2.0.1'

TEL_INDEX_NOT_VALID = -2
PIXEL_INDEX_NOT_VALID = -3





class HessioError(Exception):
    pass


class HessioGeneralError(HessioError):
    pass


class HessioTelescopeIndexError(HessioError):
    pass


class HessioChannelIndexError(HessioError):
    pass


@contextmanager
def open_hessio(filename):
    """
    Context manager
    Parameters
    ----------
    filename: str
    Yields
    ------
        HessioFile instance with file opened
    Raises
    ------
    HessioError: when a file is already open. Use pyhessio.close_file()
    """
    hessfile = HessioFile(filename, enter_by_context_mng=True)
    try:
        yield hessfile
    finally:
        hessfile.close_file()

_path = os.path.dirname(__file__)
lib_close = np.ctypeslib.load_library('pyhessioc', _path)
lib_close.close_file.restype = None


def close_file():
    """
    Because C library could not open two file at the same time, we provide
    this method to force C library to properly close file a free memory
    """
    lib_close.close_file()


class HessioFile:
    """
    Represents a pyhessio file instance
    """
    def __init__(self, filename=None, enter_by_context_mng=False):
        if not enter_by_context_mng:
            raise HessioError('HessioFile can be only use with'
                              ' context manager thanks to pyhessio.open'
                              ' function: \n'
                              'with pyhessio.open(\'pyhessio-extra/datasets/'
                              'gamma_test.simtel.gz\')'
                              ' as f: \n f.fill_next_event()')

        self.__enter_by_context_mng = False  # private
        self.__opened_filename = None        # private
        self.lib = None
        self.init_lib()
        if filename:
            self.open_file(filename)

    def init_lib(self):
        lib_path = os.path.dirname(__file__)
        self.lib = np.ctypeslib.load_library('pyhessioc', lib_path)
        self.lib.close_file.restype = None
        self.lib.file_open.argtypes = [ctypes.c_char_p]
        self.lib.file_open.restype = ctypes.c_int
        self.lib.get_adc_sample.argtypes = [ctypes.c_int, ctypes.c_int,
                                       np.ctypeslib.ndpointer(ctypes.c_uint16,
                                                              flags="C_CONTIGUOUS")]
        self.lib.get_adc_sample.restype = ctypes.c_int
        self.lib.get_adc_sum.argtypes = [ctypes.c_int, ctypes.c_int,
                                    np.ctypeslib.ndpointer(ctypes.c_uint32,
                                                           flags="C_CONTIGUOUS")]
        self.lib.get_adc_sum.restype = ctypes.c_int
        self.lib.get_significant.argtypes = [ctypes.c_int,
                                        np.ctypeslib.ndpointer(ctypes.c_uint8,
                                                               flags="C_CONTIGUOUS")]
        self.lib.get_significant.restype = ctypes.c_int
        self.lib.get_calibration.argtypes = [ctypes.c_int,
                                        np.ctypeslib.ndpointer(ctypes.c_float,
                                                               flags="C_CONTIGUOUS")]
        self.lib.get_calibration.restype = ctypes.c_int
        self.lib.get_pedestal.argtypes = [ctypes.c_int,
                                     np.ctypeslib.ndpointer(ctypes.c_float,
                                                            flags="C_CONTIGUOUS")]
        self.lib.get_pedestal.restype = ctypes.c_int
        self.lib.get_global_event_count.restype = ctypes.c_int
        self.lib.get_mirror_area.argtypes = [ctypes.c_int,
                                        np.ctypeslib.ndpointer(ctypes.c_double,
                                                               flags="C_CONTIGUOUS")]
        self.lib.get_mirror_area.restype = ctypes.c_int
        self.lib.get_num_channel.argtypes = [ctypes.c_int]
        self.lib.get_num_channel.restype = ctypes.c_int
        self.lib.get_num_pixels.argtypes = [ctypes.c_int]
        self.lib.get_num_pixels.restype = ctypes.c_int
        self.lib.get_num_trig_pixels.argtypes = [ctypes.c_int]
        self.lib.get_num_trig_pixels.restype = ctypes.c_int
        self.lib.get_trig_pixels.argtypes = [ctypes.c_int,
                                                     np.ctypeslib.ndpointer(ctypes.c_int32,
                                                         flags="C_CONTIGUOUS")]
        self.lib.get_trig_pixels.restype = ctypes.c_int
        self.lib.get_event_num_samples.argtypes = [ctypes.c_int]
        self.lib.get_event_num_samples.restype = ctypes.c_int
        self.lib.get_zero_sup_mode.argtypes = [ctypes.c_int,
                                          np.ctypeslib.ndpointer(ctypes.c_int,
                                                                 flags="C_CONTIGUOUS")]
        self.lib.get_zero_sup_mode.restype = ctypes.c_int
        self.lib.get_data_red_mode.argtypes = [ctypes.c_int,
                                          np.ctypeslib.ndpointer(ctypes.c_int,
                                                                 flags="C_CONTIGUOUS")]
        self.lib.get_data_red_mode.restype = ctypes.c_int
        self.lib.get_num_teldata.restype = ctypes.c_int
        self.lib.get_num_telescope.restype = ctypes.c_int
        self.lib.get_num_tel_trig.restype = ctypes.c_int
        self.lib.get_pixel_timing_num_times_types.argtypes = [ctypes.c_int]
        self.lib.get_pixel_timing_num_times_types.restype = ctypes.c_int
        self.lib.get_pixel_position.argtypes = [ctypes.c_int, np.ctypeslib.ndpointer(
            ctypes.c_double, flags="C_CONTIGUOUS"),
                                           np.ctypeslib.ndpointer(
                                               ctypes.c_double,
                                               flags="C_CONTIGUOUS")]
        self.lib.get_pixel_position.restype = ctypes.c_int
        self.lib.get_pixel_timing_peak_global.argtypes = [ctypes.c_int,
                                                     np.ctypeslib.ndpointer(
                                                         ctypes.c_float,
                                                         flags="C_CONTIGUOUS")]
        self.lib.get_pixel_timing_peak_global.restype = ctypes.c_int
        self.lib.get_pixel_timing_threshold.argtypes = [ctypes.c_int,
                                                   np.ctypeslib.ndpointer(
                                                       ctypes.c_int,
                                                       flags="C_CONTIGUOUS")]
        self.lib.get_pixel_timing_threshold.restype = ctypes.c_int
        self.lib.get_pixel_timing_timval.argtypes = [ctypes.c_int,
                                                np.ctypeslib.ndpointer(
                                                    ctypes.c_float,
                                                    flags="C_CONTIGUOUS")]
        self.lib.get_pixel_timing_timval.restype = ctypes.c_int
        self.lib.get_pixel_shape.argtypes = [ctypes.c_int,
                                        np.ctypeslib.ndpointer(ctypes.c_double,
                                                               flags="C_CONTIGUOUS")]
        self.lib.get_pixel_shape.restype = ctypes.c_int
        self.lib.get_pixel_area.argtypes = [ctypes.c_int,
                                       np.ctypeslib.ndpointer(ctypes.c_double,
                                                              flags="C_CONTIGUOUS")]
        self.lib.get_pixel_area.restype = ctypes.c_int
        self.lib.get_run_number.restype = ctypes.c_int
        self.lib.get_telescope_with_data_list.argtypes = [
            np.ctypeslib.ndpointer(ctypes.c_int, flags="C_CONTIGUOUS")]
        self.lib.get_telescope_with_data_list.restype = ctypes.c_int
        self.lib.get_telescope_position.argtypes = [ctypes.c_int,
                                               np.ctypeslib.ndpointer(
                                                   ctypes.c_double,
                                                   flags="C_CONTIGUOUS")]
        self.lib.get_telescope_position.restype = ctypes.c_int
        self.lib.move_to_next_event.argtypes = [np.ctypeslib.ndpointer(ctypes.c_int)]
        self.lib.move_to_next_event.restype = ctypes.c_int
        self.lib.move_to_next_mc_event.argtypes = [
            np.ctypeslib.ndpointer(ctypes.c_int)]
        self.lib.move_to_next_mc_event.restype = ctypes.c_int
        self.lib.get_mc_event_xcore.restype = ctypes.c_double
        self.lib.get_mc_event_ycore.restype = ctypes.c_double
        self.lib.get_mc_run_array_direction.argtypes = [
            np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS")]
        self.lib.get_mc_run_array_direction.restype = ctypes.c_int
        self.lib.get_azimuth_raw.argtypes = [ctypes.c_int]
        self.lib.get_azimuth_raw.restype = ctypes.c_double
        self.lib.get_altitude_raw.argtypes = [ctypes.c_int]
        self.lib.get_altitude_raw.restype = ctypes.c_double
        self.lib.get_azimuth_cor.argtypes = [ctypes.c_int]
        self.lib.get_azimuth_cor.restype = ctypes.c_double
        self.lib.get_altitude_cor.argtypes = [ctypes.c_int]
        self.lib.get_altitude_cor.restype = ctypes.c_double
        self.lib.get_mc_event_offset_fov.argtypes = [
            np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS")]
        self.lib.get_mc_event_offset_fov.restype = ctypes.c_int
        self.lib.get_mc_number_photon_electron.argtypes = [ctypes.c_int,
                                                      np.ctypeslib.ndpointer(
                                                          ctypes.c_int,
                                                          flags="C_CONTIGUOUS")]
        self.lib.get_mc_number_photon_electron.restype = ctypes.c_int
        self.lib.get_mc_shower_energy.restype = ctypes.c_double
        self.lib.get_mc_shower_xmax.restype = ctypes.c_double
        self.lib.get_mc_shower_azimuth.restype = ctypes.c_double
        self.lib.get_mc_shower_altitude.restype = ctypes.c_double
        self.lib.get_mc_shower_primary_id.restype = ctypes.c_int
        self.lib.get_mc_shower_h_first_int.restype = ctypes.c_double
        self.lib.get_spectral_index.restype = ctypes.c_double
        self.lib.get_mc_obsheight.restype = ctypes.c_double
        self.lib.get_mc_num_showers.restype = ctypes.c_int
        self.lib.get_mc_num_use.restype = ctypes.c_int
        self.lib.get_mc_core_pos_mode.restype = ctypes.c_int
        self.lib.get_mc_core_range_X.restype = ctypes.c_double
        self.lib.get_mc_core_range_Y.restype = ctypes.c_double
        self.lib.get_mc_alt_range_Min.restype = ctypes.c_double
        self.lib.get_mc_alt_range_Max.restype = ctypes.c_double
        self.lib.get_mc_az_range_Min.restype = ctypes.c_double
        self.lib.get_mc_az_range_Max.restype = ctypes.c_double
        self.lib.get_mc_viewcone_Min.restype = ctypes.c_double
        self.lib.get_mc_viewcone_Max.restype = ctypes.c_double
        self.lib.get_mc_E_range_Min.restype = ctypes.c_double
        self.lib.get_mc_E_range_Max.restype = ctypes.c_double
        self.lib.get_mc_diffuse.restype = ctypes.c_int
        self.lib.get_mc_injection_height.restype = ctypes.c_double
        self.lib.get_B_total.restype = ctypes.c_double
        self.lib.get_B_inclination.restype = ctypes.c_double
        self.lib.get_B_declination.restype = ctypes.c_double
        self.lib.get_atmosphere.restype = ctypes.c_int
        self.lib.get_corsika_version.restype = ctypes.c_int
        self.lib.get_simtel_version.restype = ctypes.c_int
        self.lib.get_corsika_iact_options.restype = ctypes.c_int
        self.lib.get_corsika_low_E_model.restype = ctypes.c_int
        self.lib.get_corsika_high_E_model.restype = ctypes.c_int
        self.lib.get_corsika_bunchsize.restype = ctypes.c_double
        self.lib.get_corsika_wlen_min.restype = ctypes.c_double
        self.lib.get_corsika_wlen_max.restype = ctypes.c_double
        self.lib.get_corsika_low_E_detail.restype = ctypes.c_int
        self.lib.get_corsika_high_E_detail.restype = ctypes.c_int
        self.lib.get_adc_known.restype = ctypes.c_int
        self.lib.get_adc_known.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_int]
        self.lib.get_ref_shape.restype = ctypes.c_double
        self.lib.get_ref_shape.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_int]
        self.lib.get_time_slice.restype = ctypes.c_double
        self.lib.get_time_slice.argtypes = [ctypes.c_int]
        self.lib.get_ref_step.restype = ctypes.c_double
        self.lib.get_ref_step.argtypes = [ctypes.c_int]
        self.lib.get_ref_shapes.restypes = ctypes.c_int
        self.lib.get_ref_shapes.argtypes = [ctypes.c_int, ctypes.c_int,
                                       np.ctypeslib.ndpointer(ctypes.c_double,
                                                              flags="C_CONTIGUOUS")]
        self.lib.get_nrefshape.restypes = ctypes.c_int
        self.lib.get_nrefshape.argtypes = [ctypes.c_int]
        self.lib.get_lrefshape.restypes = ctypes.c_int
        self.lib.get_lrefshape.argtypes = [ctypes.c_int]
        self.lib.get_tel_event_gps_time.restype = ctypes.c_int
        self.lib.get_tel_event_gps_time.argtypes = [ctypes.c_int,
                                               np.ctypeslib.ndpointer(
                                                   ctypes.c_long,
                                                   flags="C_CONTIGUOUS"),
                                               np.ctypeslib.ndpointer(
                                                   ctypes.c_long,
                                                   flags="C_CONTIGUOUS")]
        self.lib.get_central_event_gps_time.restype = ctypes.c_int
        self.lib.get_central_event_gps_time.argtypes = [
            np.ctypeslib.ndpointer(ctypes.c_long, flags="C_CONTIGUOUS"),
            np.ctypeslib.ndpointer(ctypes.c_long, flags="C_CONTIGUOUS")]
        self.lib.get_central_event_teltrg_list.restype = ctypes.c_int
        self.lib.get_central_event_teltrg_list.argtypes = [
            np.ctypeslib.ndpointer(ctypes.c_int, flags="C_CONTIGUOUS")]
        self.lib.get_central_event_teltrg_time.restype = ctypes.c_int
        self.lib.get_central_event_teltrg_time.argtypes = [
            np.ctypeslib.ndpointer(ctypes.c_float, flags="C_CONTIGUOUS")]
        self.lib.get_camera_rotation_angle.argtypes = [ctypes.c_int]
        self.lib.get_camera_rotation_angle.restype = ctypes.c_double
        self.lib.get_mirror_number.restype = ctypes.c_int
        self.lib.get_mirror_number.argtypes = [ctypes.c_int]
        self.lib.get_optical_foclen.argtypes = [ctypes.c_int]
        self.lib.get_optical_foclen.restype = ctypes.c_double
        self.lib.get_telescope_ids.argtypes = [
            np.ctypeslib.ndpointer(ctypes.c_int, flags="C_CONTIGUOUS")]
        self.lib.get_telescope_ids.restype = ctypes.c_int

    def fill_next_event(self):
        """
        Fill corresponding container with the next event in file.
        Data can be then access with other available functions in this module
        Returns
        -------
        event_number
        Raises
        ------
        HessioError: when error occurs while reading next event
        """
        event_number = np.zeros(1, dtype=np.int32)
        run_id = self.lib.move_to_next_event(event_number)
        if run_id == -1 or event_number[0] == -1:
            raise HessioError("Error while reading next event")
        return event_number[0]

    def move_to_next_event(self, limit=0):
        """
        Read data from input file and fill corresponding container.
        Data can be then access with other available functions in
        this module.
        By default all events are computed


        Parameters
        ----------
        limit: int, optional
            limit the number of event generated
        Yields
        ------
        event id
        Raises
        ------
        HessioError: When input file is not open
        """
        if not self.__opened_filename:
            raise HessioError('input file is not open')
        result = np.zeros(1, dtype=np.int32)
        run_id = 0
        evt_num = 0
        while run_id >= 0 and (limit == 0 or evt_num < limit):
            run_id = self.lib.move_to_next_event(result)
            if run_id != -1:
                yield result[0]
                evt_num += 1

    def fill_next_mc_event(self):
        """
        Fill corresponding container with the next MC event in file.
        Data can be then access with other available functions in this module
        Returns
        -------
        event_number
        Raises
        ------
        HessioError: when error occurs while reading next event
        """
        event_number = np.zeros(1, dtype=np.int32)
        run_id = self.lib.move_to_next_mc_event(event_number)
        if run_id == -1 or event_number == -1:
            raise HessioError("Error while reading next event")
        return event_number[0]

    def move_to_next_mc_event(self, limit=0):
        """
        Read MC data form input file and fill corresponding container
        Data can be then access with other available functions in
        this module.
        This iterator scans all the simulated events, not only the triggered
        ones. By default all events are computed

        Parameters
        ----------
        limit: int, optional
            limit the number of event generated
        Yields
        ------
          event id
        Raises
        ------
        HessioError: When input file is not open
        """
        if not self.__opened_filename:
            raise HessioError('input file is not open')
        result = np.zeros(1, dtype=np.int32)
        run_id = 0
        sim_evt_num = 0
        while run_id >= 0 and (limit == 0 or sim_evt_num < limit):
            run_id = self.lib.move_to_next_mc_event(result)
            if run_id != -1:
                yield result[0]
                sim_evt_num += 1

    def open_file(self, filename):
        """
        Open input data file
        Parameters
        ----------
        filename: str
        :raises: HessioError: When another file is already open or if file can
        not be open
        """
        b_filename = filename.encode('utf-8')
        self.__opened_filename = filename
        res = self.lib.file_open(b_filename)
        if res == -1:
            raise HessioError('could not open file {}'.format(filename))
        elif res == -2:
            raise HessioError('a file is already open.'
                              ' Use pyhessio.close_file()')

    def close_file(self):
        """
        Close opened iobuf
        """
        self.lib.close_file()
        self.__opened_filename = None

    def get_event_id(self):
        """
        Returns
        -------
        int : current event id
        """
        return self._event_id

    def get_corsika_version(self):
        """
        Returns
        -------
        int : CORSIKA version (x1000)
        """
        return self.lib.get_corsika_version()

    def get_simtel_version(self):
        """
        Returns
        -------
        int : sim_telarray version (x1000)
        """
        return self.lib.get_simtel_version()

    def get_global_event_count(self):
        """
        Returns
        -------
        int : counter  for system trigger
        """
        return self.lib.get_global_event_count()

    def get_run_number(self):
        """
        Returns
        -------
        int : run number
        Raises
        ------
        HessioGeneralError: when hsdata->run_header.run is not available
        """
        run = self.lib.get_run_number()
        if run > 0:
            return run
        else:
            raise(HessioGeneralError("run number not available"))

    def get_num_telescope(self):
        """
        Returns
        -------
        int : the number of telescopes in current run.
        Raises
        ------
        HessioGeneralError: when hsdata->event.num_tel is not available
        """
        number = self.lib.get_num_telescope()
        if number > 0:
            return number
        else:
            raise(HessioGeneralError("number of telescopes in current run"
                                     " not available"))

    def get_num_tel_trig(self):
        """
        Returns
        -------
        int : How many telescopes triggered in Central Event
        Raises  HessioGeneralError: when hsdata is not available
        """
        number = self.lib.get_num_tel_trig()
        if number > 0:
            return number
        else:
            raise(HessioGeneralError("number of triggered telescopes in central"
                                     " event not available"))

    def get_mirror_area(self, telescope_id):
        """
        Parameters
        ----------
        telescope_id: int
            Telescope's id
        Returns
        -------
        int:
         total area of individual mirrors corrected for inclination [m^2].
        Raises
        ------
        HessioGeneralError: when hsdata->camera_set[itel].mirror_area not
         available
        HessioTelescopeIndexError: when no telescope exist with this id
        """

        data = np.zeros(1, dtype=np.double)
        result = self.lib.get_mirror_area(telescope_id,data)
        if result == 0:
            return data[0]
        elif result == TEL_INDEX_NOT_VALID:
            raise(HessioTelescopeIndexError("no telescope with id " +
                                            str(telescope_id)))
        raise(HessioGeneralError("hsdata->camera_set[itel].mirror_area not"
                                 " available"))

    def get_telescope_with_data_list(self):
        """
        Returns
        -------
        numpy.ndarray(num_teldata,dtype=np.int32)
        list of telescope with data for current event

        Raises
        ------
        HessioGeneralError: when information is not available
        """
        try:
            return self.get_teldata_list()
        except:
            raise(HessioGeneralError("hsdata->event.teldata_list is"
                                     " not available"))

    def get_teldata_list(self):
        """
        Returns
        -------
        numpy.ndarray(num_teldata,dtype=np.int32)
        list of IDs of telescopes with data for current event
        Raises
        ------
        HessioGeneralError: when information is not available
        """
        num_teldata= self.get_num_teldata()
        if num_teldata >= 0:
            array = np.zeros(num_teldata,dtype=np.int32)
            self.lib.get_telescope_with_data_list(array)
            return array
        else:
            raise(HessioGeneralError("hsdata->event.num_teldata is "
                                     "not available"))

    def get_telescope_position(self, telescope_id):
        """
        Parameters
        ----------
        int :
        telescope_id: The telescope id
        Returns
        -------
        numpy.ndarray(3,dtype=np.double)
        Telescope position for a telescope id.

        * x is counted from array reference position towards North
        * y towards West
        * z upwards
        Raises
        ------
        HessioGeneralError: when  telescope position not available for
         this telescope
        HessioTelescopeIndexError: when no telescope exist with this id
        """
        pos = np.zeros(3,dtype=np.double)

        result = self.lib.get_telescope_position(telescope_id,pos)
        if result == 0:
            return pos
        elif result == TEL_INDEX_NOT_VALID:
            raise(HessioTelescopeIndexError("no telescope with id " +
                                            str(telescope_id)))
        else:
            raise(HessioGeneralError("no telescope position for telescope " +
                                     str(telescope_id)))

    def get_num_teldata(self):
        """
        Returns
        -------
        int:
            number of telescopes for which we actually have data
        Raises
        ------
        HessioGeneralError: when hsdata->event.num_teldata is not available
        """
        number = self.lib.get_num_teldata()
        if number >= 0:
            return number
        else:
            raise(HessioGeneralError("hsdata->event.num_teldata"
                                     " is not available"))

    def get_num_channel(self, telescope_id):
        """
        Parameters
        ----------
        int : telescope_id: telescope's id

        Returns
        -------
        int : the number of different gains per pixel for a telscope id
        * HI_GAIN  0 Index to high-gain channels in adc_sum, adc_sample, ...
        * LO_GAIN  1 Index to low-gain channels in adc_sum, adc_sample, ...

        Raises
        ------
        HessioGeneralError: when hsdata->camera_org[itel].num_gains
        HessioTelescopeIndexError: when no telescope exist with this id
        """
        result = self.lib.get_num_channel(telescope_id)
        if result >= 0:
            return result
        elif result == TEL_INDEX_NOT_VALID:
            raise(HessioTelescopeIndexError("no telescope with id " +
                                            str(telescope_id)))
        else:
            raise(HessioGeneralError(" hsdata->event.teldata[itel].raw"
                                     " not available"))

    def get_num_pixels(self, telescope_id):
        """
        Parameters
        ----------
        int:
            telescope_id: telescope's id
        Returns
        -------
        int:
            the number of pixels in the camera (as in configuration)

        Raises
        ------
        HessioGeneralError: when hsdata->camera_set[itel].num_pixels
        HessioTelescopeIndexError: when no telescope exist with this id
        """
        result = self.lib.get_num_pixels(telescope_id)
        if result >= 0 :
            return result
        elif result == TEL_INDEX_NOT_VALID:
            raise(HessioTelescopeIndexError("no telescope with id " +
                                            str(telescope_id)))
        else:
            raise(HessioGeneralError("hsdata->camera_set[itel]."
                                     "num_pixels not available"))

    def get_num_trig_pixels(self, telescope_id):
        """
        Parameters
        ----------
        int:
            telescope_id: telescope's id
        Returns
        -------
        int:
            the number of trigger pixels in the camera

        Raises
        ------
        HessioGeneralError: when hsdata->camera_set[itel].num_pixels
        HessioTelescopeIndexError: when no telescope exist with this id
        """
        result = self.lib.get_num_trig_pixels(telescope_id)
        if result >= 0 :
            return result
        elif result == TEL_INDEX_NOT_VALID:
            raise(HessioTelescopeIndexError("no telescope with id " +
                                            str(telescope_id)))
        else:
            raise(HessioGeneralError("hsdata->camera_set[itel]."
                                     "num_pixels not available"))

    def get_trig_pixels(self, telescope_id):
        """
        Parameters
        ----------
        int:
            telescope_id: telescope's id
        Returns
        -------
        int:
            the list of trigger pixels in the camera

        Raises
        ------
        HessioGeneralError: when hsdata->camera_set[itel].num_pixels
        HessioTelescopeIndexError: when no telescope exist with this id
        """
        npix = self.get_num_trig_pixels(telescope_id)
        trig_pixels = np.zeros(npix, dtype=np.int32)
        result = self.lib.get_trig_pixels(telescope_id,trig_pixels)
        if result == 0:
            return trig_pixels
        elif result == TEL_INDEX_NOT_VALID:
            raise(HessioTelescopeIndexError("no telescope with id " +
                                            str(telescope_id)))
        else:
            raise(HessioGeneralError("no pixel trigger for telescope "
                                     + str(telescope_id)))

    def get_pixel_timing_threshold(self, telescope_id):
        """
        Parameters
        ----------
        int :
            telescope_id: telescope's id
        Returns
        -------
        numpy.int32:
            PixelTiming threshold: Minimum base-to-peak raw amplitude difference
            applied in pixel selection

        Raises
        ------
        HessioGeneralError: When hsdata->event.teldata[itel].pixtm
        HessioTelescopeIndexError: when  no telescope exist with this id
        """
        threshold = np.zeros(1,dtype=np.int32)
        result = self.lib.get_pixel_timing_threshold(telescope_id,threshold)
        if result == 0: return threshold[0]
        elif result == TEL_INDEX_NOT_VALID:
            raise(HessioTelescopeIndexError("no telescope with id " +
                                            str(telescope_id)))
        else:
            raise(HessioGeneralError("hsdata->event.teldata[itel].pixtm"
                                     " not available"))

    def get_pixel_timing_peak_global(self, telescope_id):

        """
        Parameters
        ----------
        int :
            telescope_id: telescope's id

        Returns
        -------
        numpy.float32:
            PixelTiming peak_global: Camera-wide (mean) peak position
             [time slices]

        Raises
        ------
        HessioGeneralError: when hsdata->event.teldata[itel].pixtm; not
         available
        HessioTelescopeIndexError when no telescope exist with this id
        """
        peak = np.zeros(1,dtype=np.float32)
        result = self.lib.get_pixel_timing_peak_global(telescope_id,peak)
        if result == 0: return peak[0]
        elif result == TEL_INDEX_NOT_VALID:
            raise(HessioTelescopeIndexError("no telescope with id " +
                                            str(telescope_id)))
        else:
            raise(HessioGeneralError("hsdata->event.teldata[itel].pixtm;"
                                     " not available"))

    def get_pixel_timing_num_times_types(self, telescope_id):
        """
        Parameters
        ----------
        telescope_id: int
            telescope's id

        Returns:
         int :
            how many different types of times can we store

        Raises
        ------
        HessioGeneralError: when hsdata->event.teldata[itel].pixtm->num_types
         not available
        HessioTelescopeIndexError: when no telescope exist with this id
        """
        result = self.lib.get_pixel_timing_num_times_types(telescope_id)
        if result >= 0:
            return result
        elif result == TEL_INDEX_NOT_VALID:
            raise(HessioTelescopeIndexError("no telescope with id " +
                                            str(telescope_id)))
        else:
            raise(HessioGeneralError("hsdata->event.teldata[itel]."
                                     "pixtm->num_types  not available"))

    def get_significant(self, telescope_id):
        """
        Returns
        -------
        Was amplitude large enough to record it? Bit 0: sum, 1: samples.
        Parameters
        ----------
        telescope_id: int
        Raises
        ------
        HessioTelescopeIndexError
        if no telescope exist with this id
        """
        try:
            npix = self.get_num_pixels(telescope_id)
            data = np.zeros(npix,dtype=np.uint8)
            result = self.lib.get_significant(telescope_id ,data)
            if result == 0:
                return data
            elif result == TEL_INDEX_NOT_VALID:
                raise(HessioTelescopeIndexError("no telescope with id " +
                                    str(telescope_id)))
            else:
                raise(HessioGeneralError("significant not available for"
                                         " telescope " + str(telescope_id)))

        except HessioTelescopeIndexError: raise(HessioTelescopeIndexError(
                                    "significant not available for telescope " +
                                    str(telescope_id)))

    def get_event_num_samples(self, telescope_id):
        """
        Returns the number of samples (time slices) recorded.
        Only contains the number of samples for the telescopes that have
        data in the current event.

        Parameters
        ----------
        telescope_id : int
         telescope's id
        Raises
        ------
        HessioGeneralError: when data->event.teldata[itel].raw->num->samples not
         available
        HessioTelescopeIndexError: when no telescope exist with this id
        """
        result = self.lib.get_event_num_samples(telescope_id)
        if result >= 0: return result
        elif result == TEL_INDEX_NOT_VALID:
            raise(HessioTelescopeIndexError("no telescope with id " +
                                            str(telescope_id)))
        else:
            raise(HessioGeneralError("ata->event.teldata[itel].raw->num->samples"
                                     " not available"))

    def get_zero_sup_mode(self, telescope_id):
        """
        Returns
        -------
        The desired or used zero suppression mode.

        Parameters
        ----------
        telescope_id: int

        Raises
        ------
        HessioGeneralError
        data->event.teldata[itel].raw->num->samples not available

        HessioTelescopeIndexError
        if no telescope exist with this id
        """
        mode  = np.zeros(1, dtype=np.int32)
        result = self.lib.get_zero_sup_mode(telescope_id,mode)
        if result == 0:
            return mode[0]
        elif result == TEL_INDEX_NOT_VALID:
            raise(HessioTelescopeIndexError("no telescope with id " +
                                            str(telescope_id)))
        else:
            raise(HessioGeneralError("get_zero_sup_mode not available"))

    def get_data_red_mode(self, telescope_id):
        """
        Returns
        -------
        The desired or used data reduction mode.

        Parameters
        ----------
        telescope_id: int

        Raises
        ------
        HessioGeneralError
        data->event.teldata[itel].raw->num->samples not available

        HessioTelescopeIndexError
        if no telescope exist with this id
        """
        mode = np.zeros(1, dtype=np.int32)
        result = self.lib.get_data_red_mode(telescope_id,mode)
        if result == 0:
            return mode[0]
        elif result == TEL_INDEX_NOT_VALID:
            raise(HessioTelescopeIndexError("no telescope with id " +
                                            str(telescope_id)))
        else:
            raise(HessioGeneralError("data_red_mode not available"))

    def get_adc_sample(self, telescope_id):
        """
        Parameters
        ----------
        telescope_id: int
            telescope's id
        Returns
        -------
        ndarray

        Raises
        ------
        HessioGeneralError: when information is not available
        HessioTelescopeIndexError when no telescope exist with this id
        """
        n_chan = self.get_num_channel(telescope_id)
        n_pix = self.get_num_pixels(telescope_id)
        n_samples = self.get_event_num_samples(telescope_id)

        if not n_samples > 0:
            return np.zeros(0)
        data = np.zeros((n_chan, n_pix, n_samples), dtype=np.uint16)
        try:
            for chan in range(n_chan):  # (0->HI_GAIN, 1->LOW_GAIN)
                result = self.lib.get_adc_sample(telescope_id, chan,
                                                 data[chan])
                if result == 0:
                    continue
                elif result == TEL_INDEX_NOT_VALID:
                    raise (HessioTelescopeIndexError("no telescope with id " +
                                                     str(telescope_id)))
                else:
                    raise (HessioGeneralError("adc sample not available"
                                              " for telescope " +
                                              str(telescope_id) +
                                              " and channel " + str(chan)))
        except HessioTelescopeIndexError:
            raise(HessioTelescopeIndexError("no telescope with id " +
                                            str(telescope_id)))
        except HessioGeneralError:
            raise (HessioGeneralError("adc sample not available for telescope "
                                      + str(telescope_id)))
        return data

    def get_adc_sum(self, telescope_id):
        """
        Parameters
        ----------
        telescope_id: int
            telescope's id
        Returns
        -------
        ndarray

        Raises
        ------
        HessioGeneralError: when information is not available
        HessioTelescopeIndexError when no telescope exist with this id
        """
        n_chan = self.get_num_channel(telescope_id)
        n_pix = self.get_num_pixels(telescope_id)

        data = np.zeros((n_chan, n_pix), dtype=np.uint32)
        try:
            for chan in range(n_chan):  # (0->HI_GAIN, 1->LOW_GAIN)
                result = self.lib.get_adc_sum(telescope_id, chan, data[chan])
                if result == 0:
                    continue
                elif result == TEL_INDEX_NOT_VALID:
                    raise (HessioTelescopeIndexError("no telescope with id " +
                                                     str(telescope_id)))
                else:
                    raise (HessioGeneralError("adc sample not available"
                                              " for telescope " +
                                              str(telescope_id) +
                                              " and channel " + str(chan)))
        except HessioTelescopeIndexError:
            raise(HessioTelescopeIndexError("no telescope with id " +
                                            str(telescope_id)))
        except HessioGeneralError:
            raise (HessioGeneralError("adc sample not available for telescope "
                                      + str(telescope_id)))
        return data

    def get_pixel_timing_timval(self, telescope_id):
        """
        Parameters
        ----------
        telescope_id: int
            telescope's id
        Returns
        -------
            PixelTiming.timval : numpy.array(npix,ntimes,dtype=np.float32)
        Raises
        ------
        HessioGeneralError: when hsdata->event.teldata[itel]->timval[ipix]
        [itimes] not available
        HessioTelescopeIndexError: when no telescope exist with this id
        """
        npix = self.get_num_pixels(telescope_id)
        ntimes = self.get_pixel_timing_num_times_types(telescope_id)
        data = np.zeros(npix*ntimes, dtype=np.float32)
        result = self.lib.get_pixel_timing_timval(telescope_id, data)
        if result == 0:
            d_data = data.reshape(npix, ntimes)
            return d_data
        elif result == TEL_INDEX_NOT_VALID:
            raise(HessioTelescopeIndexError("no telescope with id " +
                                            str(telescope_id)))
        else:
            raise(HessioGeneralError("no pixel timing timval for telescope "
                                     + str(telescope_id)))

    def get_calibration(self, telescope_id):
        """
        Parameters
        ----------
        telescope_id: int
            telescope's id
        Returns
        -------
        nupy.array(ngain,npix,dtype=np.double)
            calibration numpy array (ngain*npix dimention)
        Raises
        ------
        HessioGeneralError: when data not available for this telescope
        HessioTelescopeIndexError: when no telescope exist with this id
        """
        npix = self.get_num_pixels(telescope_id)
        ngain = self.get_num_channel(telescope_id)

        calibration = np.zeros((ngain, npix), dtype=np.float32)

        result = self.lib.get_calibration(telescope_id, calibration)
        if result == 0:
            return calibration
        elif result == TEL_INDEX_NOT_VALID:
            raise(HessioTelescopeIndexError("no telescope with id " +
                                            str(telescope_id)))
        else:
            raise(HessioGeneralError("no calibration data for telescope "
                                     + str(telescope_id)))

    def get_pedestal(self, telescope_id):
        """
        Parameters
        ----------
        telescope_id: int
            telescope's id
        Returns
        -------
        pedestal numpy.ndarray(ngain,npix,dtype=np.double)
            numpy array (ngain*npix dimension)
        Raises
        ------
        HessioGeneralError: when data not available for this telescope
        HessioTelescopeIndexError: when no telescope exist with this id
        """
        npix = self.get_num_pixels(telescope_id)
        ngain = self.get_num_channel(telescope_id)

        pedestal = np.zeros((ngain, npix), dtype=np.float32)

        result = self.lib.get_pedestal(telescope_id, pedestal)
        if result == 0:
            return pedestal
        elif result == TEL_INDEX_NOT_VALID:
            raise(HessioTelescopeIndexError("no telescope with id " +
                                            str(telescope_id)))
        else:
            raise(HessioGeneralError("no pedestal data for telescope "
                                     + str(telescope_id)))

    def get_pixel_position(self, telescope_id):
        """
        numpy.ndarray(ngain,npix,dtype=np.double)
        Parameters
        ----------
        telescope_id: int
            telescope's id
        Returns
        -------
        tuple(numpy.ndarray(npix,dtype=np.double),
        numpy.ndarray(npix,dtype=np.double))
        pixels position for a telescope id (pos_x,pos_y)

        Raises
        ------
        HessioGeneralError: when pixel position not available for this telescope
        HessioTelescopeIndexError: when no telescope exist with this id
        """
        npix = self.get_num_pixels(telescope_id)

        pos_x = np.zeros(npix,dtype=np.double)
        pos_y = np.zeros(npix,dtype=np.double)

        result = self.lib.get_pixel_position(telescope_id, pos_x, pos_y)
        if result == 0:
            return pos_x, pos_y
        elif result == TEL_INDEX_NOT_VALID:
            raise(HessioTelescopeIndexError("no telescope with id " +
                                            str(telescope_id)))
        else:
            raise(HessioGeneralError("no pixel position for telescope "
                                     + str(telescope_id)))

    def get_pixel_shape(self, telescope_id):
        """
        Parameters
        ----------
        telescope_id: int
            telescope's id
        Returns
        -------
        numpy.ndarray(npix,dtype=np.double)
            pixels shape for a telescope id
        Raises
        ------
        HessioGeneralError: when pixel shape not available for this telescope
        HessioTelescopeIndexError: when no telescope exist with this id
        """
        npix = self.get_num_pixels(telescope_id)
        shape = np.zeros(npix, dtype=np.double)
        result = self.lib.get_pixel_shape(telescope_id,shape)
        if result == 0:
            return shape
        elif result == TEL_INDEX_NOT_VALID:
            raise(HessioTelescopeIndexError("no telescope with id " +
                                            str(telescope_id)))
        else:
            raise(HessioGeneralError("no pixel position for telescope "
                                     + str(telescope_id)))

    def get_pixel_area(self, telescope_id):
        """
        Returns pixels area for a telescope id
        Parameters
        ----------
        telescope_id: int
            telescope's id

        Returns
        --------
        numpy.ndarray(npix,dtype=np.double)
        Raises
        ------
        HessioGeneralError: when pixel area not available for this telescope
        HessioTelescopeIndexError: when no telescope exist with this id
        """
        npix = self.get_num_pixels(telescope_id)

        area = np.zeros(npix,dtype=np.double)

        result = self.lib.get_pixel_area(telescope_id, area)
        if result == 0:
            return area
        elif result == TEL_INDEX_NOT_VALID:
            raise(HessioTelescopeIndexError("no telescope with id " +
                                            str(telescope_id)))
        else:
            raise(HessioGeneralError("no pixel position for telescope "
                                     + str(telescope_id)))


    def get_mc_event_xcore(self):
        """
        Returns
        -------
        float
            x core position w.r.t. array reference point [m], x -> N
        """
        return self.lib.get_mc_event_xcore()


    def get_mc_event_ycore(self):
        """
        Returns y core position w.r.t. array reference point [m],y -> W
        Returns
        -------
        float
        """
        return self.lib.get_mc_event_ycore()


    def get_mc_run_array_direction(self):
        """
        Returns the tracking/pointing direction in [radians]. Depending on
        "tracking_mode" this either contains:
        [0]=Azimuth, [1]=Altitude in mode 0,
            OR
        [0]=R.A., [1]=Declination in mode 1.

        Returns
        -------
        numpy.ndarray(2,dtype=np.double)
        Raises
        ------
        HessioGeneralError: when information is not available
        """
        direction = np.zeros(2,dtype=np.double)

        result = self.lib.get_mc_run_array_direction(direction)
        if result == 0:
            return direction
        else:
            raise(HessioGeneralError("hsdata is not available"))


    def get_azimuth_raw(self, telescope_id):
        """
        Returns the Raw azimuth angle [radians from N->E] for the telescope
         If telescope_id is not valid return 0.
        Parameters
        ----------
        telescope_id: int
            telescope's id
        Returns double
        """
        return self.lib.get_azimuth_raw(telescope_id)


    def get_altitude_raw(self, telescope_id):
        """
        Returns the Raw altitude angle [radians] for the telescope
         If telescope_id is not valid return 0.
        Parameters
        ----------
        telescope_id: int
            telescope's id
        Returns double
        """
        return self.lib.get_altitude_raw(telescope_id)


    def get_azimuth_cor(self, telescope_id):
        """
        Returns the tracking Azimuth corrected for pointing errors for the telescope
         If telescope_id is not valid return 0.
        Parameters
        ----------
        telescope_id: int
            telescope's id
        Returns double
        """
        return self.lib.get_azimuth_cor(telescope_id)


    def get_altitude_cor(self, telescope_id):
        """
        Returns the tracking Altitude corrected for pointing errors for the telescope
         If telescope_id is not valid return 0.
        Parameters
        ----------
        telescope_id: int
            telescope's id
        Returns double
        """
        return self.lib.get_altitude_cor(telescope_id)


    def get_mc_event_offset_fov(self):
        """
        Returns offset of pointing direction in camera f.o.v. divided by focal
         length, i.e. converted to radians: [0] = Camera x (downwards in normal
          pointing, i.e. increasing Alt) [1] = Camera y -> Az.
        Returns
        -------
        numpy.ndarray(2,dtype=np.double)
        Raises
        ------
        HessioGeneralError: when information is not available
        """
        offset = np.zeros(2,dtype=np.double)

        result = self.lib.get_mc_event_offset_fov(offset)
        if result == 0:
            return offset
        else:
            raise(HessioGeneralError("hsdata is not available"))

    def get_mc_number_photon_electron(self, telescope_id):
        """
        Returns numbers of photon electron
        Parameters
        ----------
         telescope_id: int
          telescope's id
        Raises
        ------
        HessioTelescopeIndexError when no telescope with this id exists
        HessioGeneralError if hsdata is not available
        """
        npix = self.get_num_pixels(telescope_id)
        pe = np.zeros(npix,dtype=np.int32)
        result = self.lib.get_mc_number_photon_electron(telescope_id, pe)
        if result == TEL_INDEX_NOT_VALID:
            raise(HessioTelescopeIndexError("no telescope with id " +
                                            str(telescope_id)))
        elif result < 0:
            raise(HessioGeneralError("numbers of photon electron not"
                                     " available"))
        return pe

    def get_mc_shower_energy(self):
        """
        Returns shower primary energy [TeV]
        Returns
        -------
        float
        """
        return self.lib.get_mc_shower_energy()

    def get_mc_shower_xmax(self):
        """
        Returns shower primary energy [TeV]
        Returns
        -------
        float
        """
        return self.lib.get_mc_shower_xmax()

    def get_mc_shower_azimuth(self):
        """
        Returns shower azimuth (N->E) [rad]
        Returns
        -------
        float
        """
        return self.lib.get_mc_shower_azimuth()

    def get_mc_shower_altitude(self):
        """
        Returns shower altitude [rad]
        Returns
        -------
        float
        """
        return self.lib.get_mc_shower_altitude()

    def get_mc_shower_primary_id(self):
        """
        Returns shower primary ID 0 (gamma), 1(e-), 2(mu-), 100*A+Z for
         nucleons and nuclei, negative for antimatter.
        Returns
        -------
        int
        """
        return self.lib.get_mc_shower_primary_id()

    def get_mc_shower_h_first_int(self):
        """
        Returns shower height of first interaction a.s.l. [m]
        Returns
        -------
        float
        """
        return self.lib.get_mc_shower_h_first_int()

    def get_spectral_index(self):
        """
        Returns spectral_index.
        Returns
        -------
        float
        """
        return self.lib.get_spectral_index()

    def get_mc_obsheight(self):
        """
        Returns mc_obsheight.
        Returns
        -------
        float
        """
        return self.lib.get_mc_obsheight()

    def get_mc_num_showers(self):
        """
        Returns mc_num_showers.
        Returns
        -------
        int
        """
        return self.lib.get_mc_num_showers()

    def get_mc_num_use(self):
        """
        Returns mc_num_use.
        Returns
        -------
        int
        """
        return self.lib.get_mc_num_use()

    def get_mc_core_pos_mode(self):
        """
        Returns mc_core_pos_mode.
        Returns
        -------
        int
        """
        return self.lib.get_mc_core_pos_mode()

    def get_mc_core_range_X(self):
        """
        Returns mc_core_range_X.
        Returns
        -------
        float
        """
        return self.lib.get_mc_core_range_X()

    def get_mc_core_range_Y(self):
        """
        Returns mc_core_range_Y.
        Returns
        -------
        float
        """
        return self.lib.get_mc_core_range_Y()

    def get_mc_alt_range_Min(self):
        """
        Returns mc_alt_range_Min.
        Returns
        -------
        float
        """
        return self.lib.get_mc_alt_range_Min()

    def get_mc_alt_range_Max(self):
        """
        Returns mc_alt_range_Max.
        Returns
        -------
        float
        """
        return self.lib.get_mc_alt_range_Max()

    def get_mc_az_range_Min(self):
        """
        Returns mc_az_range_Min.
        Returns
        -------
        float
        """
        return self.lib.get_mc_az_range_Min()

    def get_mc_az_range_Max(self):
        """
        Returns mc_az_range_Max.
        Returns
        -------
        float
        """
        return self.lib.get_mc_az_range_Max()

    def get_mc_viewcone_Min(self):
        """
        Returns mc_viewcone_Min.
        Returns
        -------
        float
        """
        return self.lib.get_mc_viewcone_Min()

    def get_mc_viewcone_Max(self):
        """
        Returns mc_viewcone_Max.
        Returns
        -------
        float
        """
        return self.lib.get_mc_viewcone_Max()

    def get_mc_E_range_Min(self):
        """
        Returns mc_E_range_Min.
        Returns
        -------
        float
        """
        return self.lib.get_mc_E_range_Min()

    def get_mc_E_range_Max(self):
        """
        Returns mc_E_range_Max.
        Returns
        -------
        float
        """
        return self.lib.get_mc_E_range_Max()

    def get_mc_diffuse(self):
        """
        Returns mc_diffuse.
        Returns
        -------
        int
        """
        return self.lib.get_mc_diffuse()

    def get_mc_injection_height(self):
        """
        Returns mc_injection_height.
        Returns
        -------
        float
        """
        return self.lib.get_mc_injection_height()

    def get_B_total(self):
        """
        Returns B_total.
        Returns
        -------
        float
        """
        return self.lib.get_B_total()

    def get_B_inclination(self):
        """
        Returns B_inclination.
        Returns
        -------
        float
        """
        return self.lib.get_B_inclination()

    def get_B_declination(self):
        """
        Returns B_declination.
        Returns
        -------
        float
        """
        return self.lib.get_B_declination()

    def get_atmosphere(self):
        """
        Returns atmosphere.
        Returns
        -------
        int
        """
        return self.lib.get_atmosphere()

    def get_corsika_iact_options(self):
        """
        Returns corsika_iact_options.
        Returns
        -------
        int
        """
        return self.lib.get_corsika_iact_options()

    def get_corsika_low_E_model(self):
        """
        Returns corsika_low_E_model.
        Returns
        -------
        int
        """
        return self.lib.get_corsika_low_E_model()

    def get_corsika_high_E_model(self):
        """
        Returns corsika_high_E_model.
        Returns
        -------
        int
        """
        return self.lib.get_corsika_high_E_model()

    def get_corsika_bunchsize(self):
        """
        Returns corsika_bunchsize.
        Returns
        -------
        float
        """
        return self.lib.get_corsika_bunchsize()

    def get_corsika_wlen_min(self):
        """
        Returns corsika_wlen_min.
        Returns
        -------
        float
        """
        return self.lib.get_corsika_wlen_min()

    def get_corsika_wlen_max(self):
        """
        Returns corsika_wlen_max.
        Returns
        -------
        float
        """
        return self.lib.get_corsika_wlen_max()

    def get_corsika_low_E_detail(self):
        """
        Returns corsika_low_E_detail.
        Returns
        -------
        int
        """
        return self.lib.get_corsika_low_E_detail()

    def get_corsika_high_E_detail(self):
        """
        Returns corsika_high_E_detail.
        Returns
        -------
        int
        """
        return self.lib.get_corsika_high_E_detail()

    def get_adc_known(self, telescope_id, channel, pixel_id):
        """
        Returns individual channel recorded information ? Bit 0: sum, 1:
         samples, 2: ADC was in saturation.
        Parameters
        ----------
        telescope_id: int
            telescope's id
        channel: int
            HI_GAIN, LOW_GAIN
        pixel_id: int
            pixel's id
        """
        return self.lib.get_adc_known(telescope_id, channel, pixel_id)

    def get_ref_shape(self, telescope_id, channel, fshape):
        """
        Returns  Reference pulse shape(s) If telescope_id, channel or fshape
         are not valid return 0.

        Parameters
        ----------
        telescope_id: int
            telescope's id
        Parameters channel: int
            HI_GAIN, LOW_GAIN
        fshape: int
        Returns
        -------
        float
        """
        return self.lib.get_ref_shape(telescope_id, channel, fshape)

    def get_ref_shapes(self, telescope_id):
        """
        Returns Array of Reference pulse shape(s).
        Parameters
        ----------
        telescope_id: int
            telescope's id
        Returns
        -------
        ndarray

        Raises
        ------
        HessioGeneralError: when information is not available
        """
        n_chan = self.get_nrefshape(telescope_id)
        n_samples = self.get_lrefshape(telescope_id)

        if not n_samples > 0:
            raise (HessioGeneralError("Pulse reference shape(s) not "
                                      "available"))
        data = np.zeros((n_chan, n_samples), dtype=np.double)
        for chan in range(n_chan):  # (0->HI_GAIN, 1->LOW_GAIN)
            self.lib.get_ref_shapes(telescope_id, chan, data[chan])
        return data

    def get_nrefshape(self, telescope_id):
        """
        Returns  Number of following reference pulse shapes (num_gains or 0)
        ----------
        telescope_id: int
          telescope's id
        Returns int
        """
        return self.lib.get_nrefshape(telescope_id)

    def get_lrefshape(self, telescope_id):
        """
        Returns Length of following reference pulse shape(s).
        Parameters
        ----------
        telescope_id: int
            telescope's id
        Returns
        -------
        lrefshape
        """
        return self.lib.get_lrefshape(telescope_id)


    def get_ref_step(self, telescope_id):
        """
        Returns  If telescope_id, channel or fshape are not valid return 0.
        Parameters
        ----------
        telescope_id: int
         telescope's id
        Returns
        -------
        int
        """
        return self.lib.get_ref_step(telescope_id)

    def get_time_slice(self, telescope_id):
        """
        Returns Width of readout time slice (i.e. one sample) [ns].
         If telescope_id is not valid return 0.
        Parameters
        ----------
        telescope_id: int
            telescope's id
        Returns float
        """
        return self.lib.get_time_slice(telescope_id)

    def get_tel_event_gps_time(self, telescope_id):
        """
        Returns telescope event gps tine in a 2D array: -seconds  -nonosecond
        Parameters
        ----------
        telescope_id: int
            telescope's id
        Returns
        -------
        numpy.ndarray(1,dtype=np.long)
        """
        seconds = np.zeros(1,dtype=np.long)
        nanoseconds = np.zeros(1,dtype=np.long)

        result = self.lib.get_tel_event_gps_time(telescope_id,
                                                 seconds,nanoseconds)
        if result == 0:
            return seconds[0], nanoseconds[0]
        else:
            raise(HessioGeneralError("no event gps time for telescope "))

    def get_central_event_gps_time(self):
        """
        Returns telescope central envent gps tine in a 2D array:
            -seconds  -nonosecond
        Returns
        -------
        numpy.ndarray(1, dtype=np.long)
        """
        seconds = np.zeros(1, dtype=np.long)
        nanoseconds = np.zeros(1, dtype=np.long)
        result = self.lib.get_central_event_gps_time(seconds, nanoseconds)
        if result == 0:
            return seconds[0], nanoseconds[0]
        else:
            raise(HessioGeneralError("no central event  gps time"))

    def get_central_event_teltrg_list(self):
        """
        Returns List of IDs of triggered telescopes
        Returns
        -------
        np.ndarray(num_teltrig,dtype=np.int32)
        Raises
        ------
        HessioGeneralError: when information is not available
        """
        num_teltrig= self.lib.get_num_tel_trig()
        if num_teltrig >= 0:
            array = np.zeros(num_teltrig,dtype=np.int32)
            self.lib.get_central_event_teltrg_list(array)
            return array
        else:
            raise(HessioGeneralError("hsdata is not available"))

    def get_central_event_teltrg_time(self):
        """
        Returns List of relative time of trigger signal after correction for
         nominal delay (in ns) for each triggered telescope
        Returns
        -------
        np.ndarray(num_teltrig,dtype=np.float32)
        Raises
        ------
        HessioGeneralError: when information is not available
        """
        num_teltrig= self.lib.get_num_tel_trig()
        if num_teltrig >= 0:
            array = np.zeros(num_teltrig,dtype=np.float32)
            self.lib.get_central_event_teltrg_time(array)
            return array
        else:
            raise(HessioGeneralError("hsdata is not available"))


    def get_camera_rotation_angle(self, telescope_id):
        """
        Returns rotation angle of the camera of a given telescope
         (counter-clock-wise from back side for prime focus camera)L
        Parameters
        ----------
        telescope_id: int
            telescope's id
        Returns
        -------
        float
        Raises
        ------
        HessioGeneralError: when hsdata->camera_set[itel].cam_rot not available
        HessioTelescopeIndexError: when no telescope exist with this id
        """
        result = self.lib.get_camera_rotation_angle(telescope_id)
        if result >= 0 :
            return result
        elif result == TEL_INDEX_NOT_VALID:
            raise(HessioTelescopeIndexError("no telescope with id " +
                                            str(telescope_id)))
        else:
            raise(HessioGeneralError("hsdata->camera_set[itel].cam_rot"
                                     " not available"))

    def get_mirror_number(self, telescope_id):
        """
        Returns total number of mirror tiles of a telescope
        Parameters
        ----------
        telescope_id: int
            telescope's id
        Returns
        -------
        int
        Raises
        ------
        HessioGeneralError: when hsdata->camera_set[itel].num_mirrors
         not available
        HessioTelescopeIndexError when no telescope exist with this id
        """
        result = self.lib.get_mirror_number(telescope_id)
        if result >= 0 : return result
        elif result == TEL_INDEX_NOT_VALID:
            raise(HessioTelescopeIndexError("no telescope with id " +
                                            str(telescope_id)))
        else:
            raise(HessioGeneralError("hsdata->camera_set[itel].num_mirrors"
                                     " not available"))

    def get_optical_foclen(self, telescope_id):
        """
        Returns focal length of optics of a telescope [m]
        Parameters:
          telescope_id: int
        telescope's id
        Returns
        -------
        float
        Raises
        ------
        HessioGeneralError: when hsdata->camera_set[itel].flen not available
        HessioTelescopeIndexError: when no telescope exist with this id
        """
        result = self.lib.get_optical_foclen(telescope_id)
        if result >= 0 : return result
        elif result == TEL_INDEX_NOT_VALID:
            raise(HessioTelescopeIndexError("no telescope with id " +
                                            str(telescope_id)))
        else:
            raise(HessioGeneralError("hsdata->camera_set[itel].flen not"
                                     " available"))
    def get_telescope_ids(self):
        """
        Returns list of IDs of telescopes used in the run
        Returns
        -------
        numpy.ndarray(num_tel,dtype=np.int32)
        Raises
        ------
        HessioGeneralError: when information is not available
        """
        num_tel = self.get_num_telescope()
        if num_tel >= 0:
            array = np.zeros(num_tel,dtype=np.int32)
            self.lib.get_telescope_ids(array)
            return array
        else:
            raise(HessioGeneralError("hsdata->run_header.tel_id is not"
                                     " available"))
