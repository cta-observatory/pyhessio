import pytest
import numpy as np
from pyhessio import HessioChannelIndexError

try:
    from pyhessio import *
except ImportError as err:
    print("the `hessio` python module is required to access MC data: {}"
                 .format(err))
    assert(err)

def test_hessio():
    """
    v move_to_next_event(limit=0):
    v move_to_next_mc_event(limit=0):
    v file_open(filename):
    close_file():
    v get_global_event_count():
    v get_run_number():
    v get_num_telescope():
    v get_teldata_list():
    v get_num_teldata():
    v get_num_channel(telescope_id):
    v get_num_pixels(telescope_id):
    v get_num_samples(telescope_id):
    v get_adc_sample(telescope_id,channel):
    v get_adc_sum(telescope_id,channel):
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
    v get_mc_shower_primary_id()
    v get_mc_number_photon_eletron()
    v get_mc_shower_h_first_int()
    v get_telescope_position(telescope_id)
    v get_mc_event_offset_fov()
    """
    tel_id = 47
    channel = 0

    # test exception by usging getter before read the first event
    try:
        print("DEBUG", get_num_pixels(1))
        assert()
    except HessioGeneralError: pass

    # test reading file
    assert file_open("pyhessio-extra/datasets/gamma_test.simtel.gz") == 0

    run_id, event_id = next(move_to_next_event())

    assert run_id == 31964
    assert event_id == 408

    assert get_run_number() == 31964
    assert get_global_event_count() == 408
    assert get_num_telescope() == 126
    assert get_num_teldata() == 2
    assert len(get_telescope_ids()) == 126
    assert get_telescope_ids()[100] == 101

    #get_num_channel
    assert get_num_channel(tel_id) == 1
    try:
        get_num_channel(-1)
        assert()
    except HessioTelescopeIndexError: pass
    try:
        get_num_channel(1)
        assert()
    except HessioGeneralError: pass

    assert set(get_teldata_list()) == set([38, 47])

    #get_num_pixels
    assert get_num_pixels(tel_id)== 2048
    try:
        get_num_pixels(4000)
        assert()
    except HessioTelescopeIndexError: pass

    #get_adc_sample
    data_ch = get_adc_sample(tel_id, channel)
    assert np.array_equal(data_ch[10:11],[[22,20,21,24,22,19,22,27,22,21,20,22,21,20,19,22,23,20,22,20,20,23,20,20,22]]) == True

    try:
        get_adc_sample(-1, 0)
        assert()
    except HessioTelescopeIndexError: pass

    try:
        data_ch = get_adc_sample(47, 5)
        assert()
    except HessioChannelIndexError: pass

    #get_adc_sum
    data_ch_sum = get_adc_sum(tel_id,channel)
    assert  np.array_equal(data_ch_sum[0:10], [451, 550,505,465,519,467,505,496,501,478]) == True

    try:
        data_ch_sum = get_adc_sum(-1,channel)
        assert()
    except HessioTelescopeIndexError: pass

    try:
        data_ch_sum = get_adc_sum(47,2)
        assert()
    except HessioChannelIndexError: pass



    #get_num_sample
    nb_sample = get_num_samples(tel_id)
    assert nb_sample == 25
    try:
        get_num_samples(70000)
        assert()
    except HessioTelescopeIndexError:pass


    #get_calibration
    calibration = get_calibration(tel_id)
    assert calibration[0][2] ==  0.086124487221240997
    try :
        get_calibration(0)
        assert()
    except HessioTelescopeIndexError: pass

    #get_pedestal
    pedestal = get_pedestal(tel_id)
    assert pedestal[0][0] == 457.36550903320312
    try :
        get_pedestal(0)
        assert()
    except HessioTelescopeIndexError: pass


    #get_pixel_position
    pos_x,pos_y = get_pixel_position(tel_id)
    assert pos_x[2] == -0.085799999535083771
    assert pos_y[2] == -0.14880000054836273
    try:
        get_pixel_position(0)
        assert()
    except HessioTelescopeIndexError: pass

    assert(np.array_equal(get_telescope_with_data_list() , [38, 47]) == True)

    #get_pixel_shape
    shape = get_pixel_shape(tel_id)
    assert shape[0] == -1.0
    assert len(shape) == 2048
    try:
        get_pixel_shape(0)
        assert()
    except HessioTelescopeIndexError: pass

    #get_pixel_area
    p_area = get_pixel_area(tel_id)
    assert p_area[0] == 3.3640000765444711e-05
    try:
        get_pixel_area(0)
        assert()
    except HessioTelescopeIndexError: pass

    #get_camera_rotation_angle
    assert(float(get_camera_rotation_angle(tel_id)) == 0.0)
    try:
        get_camera_rotation_angle(-1)
        assert()
    except HessioTelescopeIndexError: pass

    #get_mirror_area
    assert(get_mirror_area(tel_id) ==  14.562566757202148)
    try:
        get_mirror_area(-1)
        assert()
    except HessioTelescopeIndexError: pass

    #get_mirror_number
    assert(get_mirror_number(tel_id) == 2)
    try:
        get_mirror_number(-1)
        assert()
    except HessioTelescopeIndexError: pass

    # get_optical_foclen
    assert(float(get_optical_foclen(tel_id)) == float(2.1500000953674316))
    try:
        get_optical_foclen(-1)
        assert()
    except HessioTelescopeIndexError: pass

    # get_pixel_timing_num_times_types
    assert(get_pixel_timing_num_times_types(tel_id) == 7)
    try:
        get_pixel_timing_num_times_types(4000)
        assert()
    except HessioTelescopeIndexError: pass
    assert(get_pixel_timing_num_times_types(1) == 0)


    #get_pixel_threashold
    assert(get_pixel_timing_threshold(tel_id)== -6)
    try:
        get_pixel_timing_threshold(-1)
        assert()
    except  HessioTelescopeIndexError: pass

    #get_pixel_timing_peak_global
    assert(float(get_pixel_timing_peak_global(tel_id)) == float(9.740449905395508))
    try:
        get_pixel_timing_peak_global(1000)
        assert()
    except HessioTelescopeIndexError: pass

    assert(float(get_pixel_timing_timval(tel_id)[8][0]) ==  float(11.069999694824219) )
    try:
        get_pixel_timing_timval(-1)
        assert()
    except HessioTelescopeIndexError: pass

    # Telescope position
    tel_x, tel_y, tel_z = get_telescope_position(tel_id)
    assert(float( tel_x ) == float(1223.800048828125))
    assert(float( tel_y ) == float(704.0999755859375))
    assert(float( tel_z ) == float(5.0))
    try:
        get_telescope_position(-1)
        assert()
    except HessioTelescopeIndexError: pass

    mc_offset_x, mc_offset_y = get_mc_event_offset_fov()
    assert( mc_offset_x == 0)
    assert( mc_offset_y == 0)

    assert(get_mc_number_photon_electron(1,1) == 0)


    """
    xcode 1129.6055908203125
    ycode 547.77001953125
    shower energy 0.3820943236351013
    shower azimuth 6.283185005187988
    shower altitude 1.2217304706573486
    teltrig_list [38 47]
    get_adc_knows(38,0,1000) 1
    get_ref_shape(38,0,2) 0.0269622802734375
    get_ref_step(38) 0.3003003001213074
    get_time_slice(38) 3.0030031204223633
    tel_event gps seconds 0 nanoseconds 0
    central_event_gps_time 1319627141 nanoseconds 1275579058
    num tel trig 2
    teltrig_list [38 47]

    """


    assert(float(get_mc_event_xcore()) == float(1129.6055908203125))
    assert(float(get_mc_event_ycore()) == float(547.77001953125))
    assert(float(get_mc_shower_energy()) == float(.3820943236351013))
    assert(get_mc_number_photon_electron(1,1) == 0)
    assert(float(get_mc_shower_azimuth()) == float(6.283185005187988))
    assert(float(get_mc_shower_altitude()) == float( 1.2217304706573486))

    assert(get_mc_shower_primary_id() == 0)
    assert(float(get_mc_shower_h_first_int()) == float(17846.654296875))

    assert(get_adc_known(38,0,1000)== 1 )

    assert(float(get_ref_shape(38,0,2)) == float(.0269622802734375))
    assert(float(get_ref_step(38)) == float(.3003003001213074))
    assert(float(get_time_slice(38))== float(3.0030031204223633))

    seconds, nanoseconds = get_tel_event_gps_time(38)
    assert(seconds == 0)
    assert(nanoseconds == 0)

    seconds, nanoseconds = get_central_event_gps_time()
    assert(seconds == 1408549473)
    assert(nanoseconds ==35597000)

    num_tel_trig = get_num_tel_trig()
    assert(num_tel_trig == 2 )

    assert(np.array_equal(get_central_event_teltrg_list() ,[38, 47]) == True)

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

    ref_shapes = get_ref_shapes(38,0)
    assert(float(ref_shapes[0]) == float(0.01372528076171875))
    assert(float(ref_shapes[79]) == float(0.0009870529174804688))
    assert(get_nrefshape(38) == 1)
    assert(get_lrefshape(38) == 80)

    close_file()

    assert file_open("pyhessio-extra/datasets/gamma_test.simtel.gz") == 0

    # Testing move_to_next_mc_event iterator
    run_id, event_id = next(move_to_next_mc_event())

    assert run_id == 31964
    assert event_id == 100 # Different from before

    # Configuration is the same
    assert len(get_telescope_ids()) == 126
    assert get_telescope_ids()[100] == 101

    # now we don't have any telescope data
    assert(get_num_teldata() == 0)

    # but we have simulation data
    assert(get_mc_shower_primary_id() == 0)
    assert(float(get_mc_event_xcore()) == float(-591.63360595703125))
    assert(float(get_mc_event_ycore()) == float(1080.77392578125))

    close_file()

if __name__ == "__main__":
    test_hessio()
