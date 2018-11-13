
import numpy as np

from pyhessio import *


def test_hessio():
    """
    v move_to_next_event(limit=0, event_type = EventType.CHERENKOV.value):
    v fill_next_event():
    v open_file(filename):
    v close_file():
    v get_global_event_count():
    v get_run_number():
    v get_num_telescope():
    v get_teldata_list():
    v get_num_teldata():
    v get_num_channel(telescope_id):
    v get_num_pixels(telescope_id):
    v get_event_num_samples(telescope_id):
    v get_zero_sup_mode(telescope_id):
    v get_data_red_mode(telescope_id):
    v get_adc_sample(telescope_id):
    v get_adc_sum(telescope_id):
    v get_significant(telescope_id):
    v get_data_for_calibration(telescope_id):
    v get_pixel_position(telescope_id):
    v get_pixel_shape(telescope_id)
    v get_pixel_area(telescope_id)
    v get_telescope_with_data_list()
    v get_pixel_timing_timval(telescope_id)
    v get_mirror_area(telescope_id)
    v get_pixel_timing_num_times_types(telescope_id)
    v get_pixel_timing_threshold(telescope_id)
    v get_pixel_timing_peak_global(telescope_id)
    v get_mirror_number(telescope_id)
    v get_optical_foclen(telescope_id)
    v get_mc_number_photon_electron(tel_id, pixel_id)
    v get_telescope_ids()
    v get_mc_shower_num()
    v get_mc_shower_primary_id()
    v get_mc_number_photon_eletron()
    v get_mc_shower_h_first_int()
    v get_mc_shower_xmax()
    v get_mc_shower_hmax()
    v get_telescope_position(telescope_id)
    v get_mc_event_offset_fov()
    """
    tel_id = 47


    with open_hessio('pyhessio-extra/datasets/gamma_test.simtel.gz') as hessio:
# test exception by using getter before read the first event
        try:
            hessio.get_num_pixels(1)
            raise
        except HessioGeneralError: pass

        try:
            hessio.open_file("pyhessio-extra/datasets/gamma_test.simtel.gz")
            raise HessioError('You cannot open 2 file simultaneity')
        except HessioError:
            pass
        event_id = next(hessio.move_to_next_event())

        assert hessio.get_run_number() == 31964
        assert hessio.get_global_event_count() == 408 == event_id
        assert hessio.get_num_telescope() == 126
        assert hessio.get_num_teldata() == 2
        assert len(hessio.get_telescope_ids()) == 126
        assert hessio.get_telescope_ids()[100] == 101

    with open_hessio('pyhessio-extra/datasets/gamma_test.simtel.gz') as hessio:
        event_id = hessio.fill_next_event()
        assert hessio.get_global_event_count() == 408 == event_id
        assert hessio.get_num_channel(tel_id) == 1
        try:
            hessio.get_num_channel(-1)
            raise
        except HessioTelescopeIndexError:
            pass

        assert set(hessio.get_teldata_list()) == set([38, 47])

        #get_num_pixels
        assert hessio.get_num_pixels(tel_id)== 2048
        try:
            hessio.get_num_pixels(4000)
            raise
        except HessioTelescopeIndexError:
            pass

        # test the trigger data
        assert (hessio.get_num_trig_pixels(tel_id) == 3)
        assert (np.array_equal(hessio.get_trig_pixels(tel_id), [68, 1242, 1338]) == True)

        #get_significant
        data_sig = hessio.get_significant(tel_id)
        model = np.ones(2048)
        assert np.array_equal(data_sig, model) is True

        #get_event_num_samples
        nb_sample = hessio.get_event_num_samples(tel_id)
        assert nb_sample == 25
        try:
            hessio.get_event_num_samples(70000)
            raise
        except HessioTelescopeIndexError:
            pass

        #get_zero_sup_mode
        mode = hessio.get_zero_sup_mode(tel_id)
        assert mode == 0
        try:
            hessio.get_zero_sup_mode(70000)
            raise
        except HessioTelescopeIndexError:pass

        #get_data_red_mode
        mode = hessio.get_data_red_mode(tel_id)
        assert mode == 0
        try:
            hessio.get_data_red_mode(70000)
            raise
        except HessioTelescopeIndexError:pass

        #get_calibration
        calibration = hessio.get_calibration(tel_id)
        assert np.isclose(calibration[0][2], 0.086124487221240997)
        try :
            hessio.get_calibration(0)
            raise
        except HessioTelescopeIndexError: pass

        #get_pedestal
        pedestal = hessio.get_pedestal(tel_id)
        assert np.isclose(pedestal[0][0], 457.36550903320312)
        try :
            hessio.get_pedestal(0)
            raise
        except HessioTelescopeIndexError: pass


        #get_pixel_position
        pos_x, pos_y = hessio.get_pixel_position(tel_id)
        assert np.isclose(pos_x[2], -0.085799999535083771)

        assert np.isclose(pos_y[2], -0.14880000054836273)
        try:
            hessio.get_pixel_position(0)
            raise
        except HessioTelescopeIndexError: pass

        assert(np.array_equal(hessio.get_telescope_with_data_list(), [38, 47]))

        #get_pixel_shape
        shape = hessio.get_pixel_shape(tel_id)
        assert shape[0] == -1.0
        assert len(shape) == 2048
        try:
            hessio.get_pixel_shape(0)
            raise
        except HessioTelescopeIndexError: pass

        #get_pixel_area
        p_area = hessio.get_pixel_area(tel_id)
        assert p_area[0] == 3.3640000765444711e-05
        try:
            hessio.get_pixel_area(0)
            raise
        except HessioTelescopeIndexError: pass

        #get_camera_rotation_angle
        assert(float(hessio.get_camera_rotation_angle(tel_id)) == 0.0)
        try:
            hessio.get_camera_rotation_angle(-1)
            raise
        except HessioTelescopeIndexError: pass

        #get_mirror_area
        assert(np.isclose(hessio.get_mirror_area(tel_id), 14.562566757202148))
        try:
            hessio.get_mirror_area(-1)
            raise
        except HessioTelescopeIndexError: pass

        #get_mirror_number
        assert(hessio.get_mirror_number(tel_id) == 2)
        try:
            hessio.get_mirror_number(-1)
            raise
        except HessioTelescopeIndexError: pass

        # get_optical_foclen
        assert(np.isclose((hessio.get_optical_foclen(tel_id)),2.1500000953674316))
        try:
            hessio.get_optical_foclen(-1)
            raise
        except HessioTelescopeIndexError: pass

        # get_pixel_timing_num_times_types
        assert(hessio.get_pixel_timing_num_times_types(tel_id) == 7)
        try:
            hessio.get_pixel_timing_num_times_types(4000)
            raise
        except HessioTelescopeIndexError: pass
        assert(hessio.get_pixel_timing_num_times_types(1) == 0)


        #get_pixel_threashold
        assert(hessio.get_pixel_timing_threshold(tel_id)== -6)
        try:
            hessio.get_pixel_timing_threshold(-1)
            raise
        except  HessioTelescopeIndexError: pass

        #get_pixel_timing_peak_global
        assert(float(hessio.get_pixel_timing_peak_global(tel_id)) == float(9.740449905395508))
        try:
            hessio.get_pixel_timing_peak_global(1000)
            raise
        except HessioTelescopeIndexError:
            pass

        assert(float(hessio.get_pixel_timing_timval(tel_id)[8][0]) == float(11.069999694824219) )
        try:
            hessio.get_pixel_timing_timval(-1)
            raise
        except HessioTelescopeIndexError:
            pass

        # Telescope position
        tel_x, tel_y, tel_z = hessio.get_telescope_position(tel_id)
        assert(float(tel_x) == float(1223.800048828125))
        assert(float(tel_y) == float(704.0999755859375))
        assert(float(tel_z) == float(5.0))
        try:
            hessio.get_telescope_position(-1)
            raise
        except HessioTelescopeIndexError:
            pass

        az, alt = hessio.get_mc_run_array_direction()
        assert(float(az) == float(0.0))
        assert(float(alt) == float(1.2217304706573486))
        assert(float(hessio.get_azimuth_raw(tel_id)) == float(0.0))
        assert(float(hessio.get_altitude_raw(tel_id)) == float(1.2217304706573486))
        assert(float(hessio.get_azimuth_cor(tel_id)) == float(0.0))
        assert(float(hessio.get_altitude_cor(tel_id)) == float(0.0))

        mc_offset_x, mc_offset_y = hessio.get_mc_event_offset_fov()
        assert(mc_offset_x == 0)
        assert(mc_offset_y == 0)

        assert hessio.get_mc_event_num() == 408
        assert(float(hessio.get_mc_event_xcore()) == float(1129.6055908203125))
        assert(float(hessio.get_mc_event_ycore()) == float(547.77001953125))
        assert(np.isclose(hessio.get_mc_event_shower_num(), 4))


        assert(float(hessio.get_mc_shower_energy()) == float(.3820943236351013))
        assert(hessio.get_mc_number_photon_electron(1)[0] == 0)
        assert(float(hessio.get_mc_shower_azimuth()) == float(6.283185005187988))
        assert(float(hessio.get_mc_shower_altitude()) == float(1.2217304706573486))
        assert np.isclose(hessio.get_mc_shower_xmax(), 339.1954)
        assert np.isclose(hessio.get_mc_shower_hmax(), 8872.70)
        assert np.isclose(hessio.get_mc_core_range_min(), 0.0)

        assert np.isclose(hessio.get_mc_core_range_max(), 2500.)
        assert np.isclose(hessio.get_mc_core_range_X(), 0.0)
        assert np.isclose(hessio.get_mc_core_range_Y(), 2500.)

        assert (hessio.get_mc_shower_num() == 4)
        assert(hessio.get_mc_shower_primary_id() == 0)
        assert(float(hessio.get_mc_shower_h_first_int()) == float(17846.654296875))

        assert(hessio.get_adc_known(38,0,1000)== 3 )

        assert(float(hessio.get_ref_shape(38,0,2)) == float(.0269622802734375))
        assert(float(hessio.get_ref_step(38)) == float(.3003003001213074))
        assert(float(hessio.get_time_slice(38))== float(3.0030031204223633))

        seconds, nanoseconds = hessio.get_tel_event_gps_time(38)
        assert(seconds == 0)
        assert(nanoseconds == 0)

        seconds, nanoseconds = hessio.get_central_event_gps_time()
        assert(seconds == 1408549473)
        assert(nanoseconds ==35597000)

        num_tel_trig = hessio.get_num_tel_trig()
        assert(num_tel_trig == 2 )

        assert(np.array_equal(hessio.get_central_event_teltrg_list(), [38, 47]) == True)
        assert(np.array_equal(hessio.get_central_event_teltrg_time(), np.array([3.75554633, 0.0], dtype=np.float32)) == True)

        """
        size 80 ref_shapes [  1.37252808e-02   1.89666748e-02   2.69622803e-02   3.62854004e-02
       4.95605469e-02   6.60400391e-02   8.57543945e-02   1.12304688e-01
       1.41845703e-01   1.79443359e-01   2.22412109e-01   2.70507812e-01
       3.27392578e-01   3.87207031e-01   4.54101562e-01   5.23437500e-01
       5.95703125e-01   6.67480469e-01   7.37792969e-01   8.03710938e-01
       8.62304688e-01   9.15039062e-01   9.52148438e-01   9.82421875e-01
       9.94140625e-01   9.93652344e-01   9.82910156e-01   9.55566406e-01
       9.22851562e-01   8.79882812e-01   8.32519531e-01   7.81738281e-01
       7.28027344e-01   6.73828125e-01   6.19140625e-01   5.66406250e-01
       5.14648438e-01   4.66552734e-01   4.19921875e-01   3.77197266e-01
       3.37402344e-01   3.00048828e-01   2.67089844e-01   2.35839844e-01
       2.08496094e-01   1.83593750e-01   1.60766602e-01   1.41357422e-01
       1.23107910e-01   1.07604980e-01   9.36889648e-02   8.11157227e-02
       7.06787109e-02   6.09741211e-02   5.28259277e-02   4.56237793e-02
       3.91540527e-02   3.38745117e-02   2.90222168e-02   2.49786377e-02
       2.14385986e-02   1.82800293e-02   1.57318115e-02   1.33972168e-02
       1.14746094e-02   9.80377197e-03   8.30841064e-03   7.12966919e-03
       6.03866577e-03   5.15365601e-03   4.38308716e-03   3.70216370e-03
       3.16429138e-03   2.67219543e-03   2.27165222e-03   1.92737579e-03
       1.62220001e-03   1.38282776e-03   1.16348267e-03   9.87052917e-04]
    nrefshapes 1
    lref_shape 80

        """



        hessio.close_file()
        hessio.open_file("pyhessio-extra/datasets/gamma_test.simtel.gz")
        event_id = hessio.fill_next_event( EventType.MC.value )
        run_number = hessio.get_run_number()

        assert run_number == 31964
        assert event_id == 100

        # Configuration is the same
        assert len(hessio.get_telescope_ids()) == 126
        assert hessio.get_telescope_ids()[100] == 101

        # now we don't have any telescope data
        assert(hessio.get_num_teldata() == 0)

        # but we have simulation data
        assert(hessio.get_mc_shower_primary_id() == 0)
        assert(float(hessio.get_mc_event_xcore()) == float(-591.63360595703125))
        assert(float(hessio.get_mc_event_ycore()) == float(1080.77392578125))

        close_file()

    # test calibration events
    with open_hessio('pyhessio-extra/datasets/calibevents_test.simtel.gz') as calib_hessio:
        try:
            calib_hessio.fill_next_event( EventType.PEDESTAL.value )
            calib_run_id = calib_hessio.get_run_number()
            calib_tel_with_data = len(calib_hessio.get_telescope_with_data_list())
            calib_numsamples = calib_hessio.get_event_num_samples(47)
            calib_sample = calib_hessio.get_adc_sample(47)
            calib_sample_v = calib_sample[0,455,20]

            assert calib_run_id == 22
            assert calib_tel_with_data == 99
            assert calib_numsamples == 50
            assert calib_sample_v == 39
        except HessioTelescopeIndexError:
            pass



# Test data on 2 channels event
def test_hessio_two_channels():
    tel_id = 2

    with open_hessio('pyhessio-extra/datasets/gamma_test_large.simtel.gz') as hessio:
        # read an event containing 2 channels telscopes

        for loop in range(4):
            next(hessio.move_to_next_event())

        # get_adc_sample
            data_ch = hessio.get_adc_sample(tel_id)
        #get_adc_sample channel 0
        channel = 0
        assert np.array_equal(data_ch[channel][10:11], [[292, 302, 295, 283, 256, 237, 249, 253, 250,
                                                291, 287, 307, 276, 274, 288, 351, 319, 300,
                                                294, 306, 365, 364, 357, 293, 267, 288, 294,
                                                345, 378, 365]]) == True
        #get_adc_sample channel 1
        channel = 1
        assert np.array_equal(data_ch[channel][10:11], [[300, 296, 298, 298, 296, 299, 298, 296, 299,
                                                298, 299, 302, 297, 296, 296, 303, 299, 299,
                                                299, 302, 301, 302, 304, 302, 300, 295, 300,
                                                300, 309, 301]]) == True

        #  get_adc_sample on wrong tel_id
        try:
            data_ch = hessio.get_adc_sample(-1)[channel]
            raise
        except HessioTelescopeIndexError:
            pass

        # get_adc_sum
        data_ch_sum = hessio.get_adc_sum(tel_id)

        #get_adc_sum channel 0
        channel = 0
        assert np.array_equal(data_ch_sum[channel][0:10], [8893, 8605, 8031, 9105, 9957,
                                                  8451, 8820, 9578, 8010, 9091]) is True
        #get_adc_sum channel 1
        channel = 1
        assert np.array_equal(data_ch_sum[channel][0:10], [8988, 8950, 8970, 8949, 9028,
                                                  8988, 8965, 8969, 8880, 8997]) is True

        try:
            data_ch_sum = hessio.get_adc_sum(-1)[channel]
            raise
        except HessioTelescopeIndexError:
            pass

        #get_ref_shapes channel 0
        channel = 0
        ref_shapes = hessio.get_ref_shapes(tel_id)[channel]
        assert (np.isclose(ref_shapes[0], 2.980232238769531e-07))
        assert (np.isclose(ref_shapes[79], 0.64013671875))

        # get_ref_shapes channel 1
        channel = 1
        ref_shapes = hessio.get_ref_shapes(tel_id)[channel]
        assert (np.isclose(ref_shapes[0], 2.980232238769531e-07))
        assert (np.isclose(ref_shapes[79], 0.64013671875))

        assert(hessio.get_nrefshape(tel_id) == 2)
        assert(hessio.get_lrefshape(tel_id) == 250)
