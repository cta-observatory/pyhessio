/** @file pyhessio.c
 *  @short A wrapper program reading H.E.S.S. data with python
 *
 * Suggestions/complains to jacquem@lapp.in2p3.fr
 *
 */
#include "initial.h"
#include "io_basic.h"
#include "io_hess.h"
#include "io_history.h"
#include "io_histogram.h"
#include "fileopen.h"
#include "stdio.h"
void close_file (void);
int file_open (const char *filename);
void free_hsdata(void);
int fill_hsdata (int *event_id);
int get_adc_sample (int telescope_id, int channel, uint16_t * data);
int get_adc_sum (int telescope_id, int channel, uint32_t * data);
uint8_t get_significant (int telescope_id, uint8_t * data);
int get_pedestal (int telescope_id, float *pedestal);
int get_calibration (int telescope_id, float *calib);
int get_global_event_count (void);
int get_mirror_area (int telescope_id, double *mirror_area);
int get_num_channel (int telescope_id);
int get_num_pixels (int telescope_id);
int get_num_trig_pixels (int telescope_id);
int get_trig_pixels (int telescope_id, int *trigpix);
int get_event_num_samples (int telescope_id);
int get_zero_sup_mode(int telescope_id,int* result);
int get_data_red_mode(int telescope_id,int* result);
int get_num_teldata (void);
int get_num_telescope (void);
int get_pixel_timing_num_times_types (int telescope_id);
int get_pixel_position (int telescope_id, double *xpos, double *ypos);
int get_pixel_timing_threshold (int telescope_id, int *result);
int get_pixel_timing_timval (int telescope_id, float *data);
int get_pixel_timine_peak_global (int telescope_id, float *peak);
int get_pixel_shape(int telescope_id, double *pixel_shape);
int get_pixel_area(int telescope_id, double *pixel_area);
int get_corsika_version (void);
int get_simtel_version (void);
int get_run_number (void);
int get_telescope_with_data_list (int *list);
int get_telescope_position (int telescope_id, double *pos);
int get_telescope_index (int telescope_id);
int move_to_next();
int move_to_next_event (int *event_id, int event_type );
int move_to_next_mc_event (int *event_id);
int move_to_next_calib_event (int *event_id);
double get_mc_event_xcore (void);
double get_mc_event_ycore (void);
long get_mc_num_generated_events(const char *filename);
int get_mc_run_array_direction (double *dir);
double get_azimuth_raw (int telescope_id);
double get_altitude_raw (int telescope_id);
double get_azimuth_cor (int telescope_id);
double get_altitude_cor (int telescope_id);
int get_mc_event_offset_fov (double *off);
int get_mc_shower_num (void);
double get_mc_shower_energy (void);
double get_mc_shower_xmax (void);
double get_mc_shower_hmax (void);
double get_mc_shower_azimuth (void);
double get_mc_shower_altitude (void);
int get_mc_shower_primary_id(void);
double get_mc_shower_h_first_int(void);
int get_mc_number_photon_electron(int telescope_id, int* pe);
double get_spectral_index(void);
double get_mc_obsheight(void);
int get_mc_num_showers(void);
int get_mc_num_use(void);
int get_mc_core_pos_mode(void);
double get_mc_core_range_X(void);
double get_mc_core_range_Y(void);
double get_mc_alt_range_Min(void);
double get_mc_alt_range_Max(void);
double get_mc_az_range_Min(void);
double get_mc_az_range_Max(void);
double get_mc_viewcone_Min(void);
double get_mc_viewcone_Max(void);
double get_mc_E_range_Min(void);
double get_mc_E_range_Max(void);
double get_B_total(void);
double get_B_inclination(void);
double get_B_declination(void);
double get_atmosphere(void);
double get_corsika_iact_options(void);
double get_corsika_low_E_model(void);
double get_corsika_high_E_model(void);
double get_corsika_bunchsize(void);
double get_corsika_wlen_min(void);
double get_corsika_wlen_max(void);
double get_corsika_low_E_detail(void);
double get_corsika_high_E_detail(void);
int get_mc_diffuse(void);
double get_mc_injection_height(void);
uint8_t get_adc_known (int telescope_id, int channel, int pixel_id);
double get_ref_shape (int telescope_id, int channel, int fshape);
double get_ref_step (int telescope_id);
double get_time_slice (int telescope_id);
int get_tel_event_gps_time (int telescope_id, long *seconds,
			    long *nanoseconds);
int get_central_event_gps_time (long *seconds, long *nanoseconds);
int get_central_event_teltrg_list (int *tel_list);
int get_central_event_teltrg_time (float *teltrg_time);
int get_num_tel_trig (void);
int get_ref_shapes (int telescope_id, int channel, double *ref_shapes);
int get_nrefshape (int telescope_id);
int get_lrefshape (int telescope_id);
int get_mirror_number(int telescope_id);
double get_camera_rotation_angle(int telescope_id);
double get_optical_foclen(int telescope_id);
int get_telescope_ids(int* list);
int show_history(void);

static AllHessData *hsdata = NULL;
static IO_ITEM_HEADER item_header;
static IO_BUFFER *iobuf = NULL;
static int file_is_opened = 0;
static int showhistory = 0;
#define TEL_INDEX_NOT_VALID -2
#define PIXEL_INDEX_NOT_VALID -3
//-----------------------------------
// Returns array index for specific id
// Returns TEL_INDEX_NOT_VALID if telescope index is not valid
//-----------------------------------
int get_telescope_index (int telescope_id){
	int itel = 0;
	for (itel = 0; itel < hsdata->run_header.ntel; itel++){
		if ((hsdata)->run_header.tel_id[itel] == telescope_id)
		return itel;
	}
	return TEL_INDEX_NOT_VALID;
}

//----------------------------------
//Read input file and fill hsdata
// and item_header global var
// Return :
//    -1 if filename is does not exit
//    -2 if a file is already opened
//----------------------------------
int file_open (const char *filename){
	if (filename)
	{
		if (file_is_opened){
			return -2;
		}
		/* Check assumed limits with the ones compiled into the library. */
		H_CHECK_MAX();

		if ((iobuf = allocate_io_buffer (1000000L)) == NULL){
			Error ("Cannot allocate I/O buffer");
			exit (1);
		}
		iobuf->max_length = 100000000L;
		if ((iobuf->input_file = fileopen (filename, READ_BINARY)) == NULL){
			perror (filename);
			Error ("Cannot open input file.");
			return -1;
		}
		file_is_opened = 1;
	}
	return 0;
}

//----------------------------------
//Read input file and fill hsdata
// and item_header global var
// return item type
//----------------------------------
int move_to_next (){
	if (!file_is_opened)
		return -1;
	int foo = 0;
	return  fill_hsdata (&foo);
}


//----------------------------------
//Read input file and fill hsdata and and new event is found
// and item_header global var
//----------------------------------
int move_to_next_event (int *event_id, int event_type ){
	if (!file_is_opened)
		return -1;
	int rc = 0;
	while (rc != event_type){
		rc = fill_hsdata (event_id);
		if (rc < 0){
			return -1;
		}
	} 
	return get_run_number ();
}
//----------------------------------
//Read input file and fill hsdata
// and item_header global var
//Scan all simulated events
//----------------------------------
int move_to_next_mc_event (int *event_id){
	if (!file_is_opened) return -1;
	int rc = 0;
	while (rc != IO_TYPE_HESS_MC_EVENT)	{
		rc = fill_hsdata (event_id);
		if (rc < 0){
			return -1;
		}
	}
	return get_run_number ();
}
//----------------------------------
// Scan all calibration events
//----------------------------------
int move_to_next_calib_event (int *event_id){
	if (!file_is_opened) return -1;
	int rc = 0;
	while (rc != IO_TYPE_HESS_CALIBEVENT)	{
		rc = fill_hsdata (event_id);
		if (rc < 0){
			return -1;
		}
	}
	return get_run_number ();
}
/* 
 * show how sim_telarray was run and configured
*/
int show_history(){
    if (!file_is_opened) return -1;
    int rc = 0;
    int *event_id = 0;
    showhistory = 1;
    while (rc >=0 ) {
        rc = fill_hsdata (event_id);
    }
    showhistory = 0;
    return 1;
}

/*--------------------------------*/
//  Cleanly close iobuf
//----------------------------------
void close_file(){
    if (iobuf == NULL || file_is_opened == 0) return;

	if (iobuf->input_file != NULL && iobuf->input_file != stdin && file_is_opened)
	{
		fileclose (iobuf->input_file);
		if (hsdata != NULL) free_hsdata();
		iobuf->input_file = NULL;
		reset_io_block (iobuf);
		free_io_buffer(iobuf);
	}
	if (iobuf->output_file != NULL) fileclose(iobuf->output_file);
	iobuf = NULL;
	file_is_opened = 0;

}
//------------------------------------------
//  return run number from last readed event
//------------------------------------------
int get_run_number (void){
	if (hsdata != NULL){
		return hsdata->run_header.run;
	}
	return -1;
}
//------------------------------------
// Returns number of telescopes in run.
//------------------------------------
int get_num_telescope (void){
	if (hsdata != NULL){
		return hsdata->event.num_tel;
	}
return -1;
}
//------------------------------------------------------------
// Returns number of telescopes for which we actually have data
//------------------------------------------------------------
int get_num_teldata (void)
{
	if (hsdata != NULL){
		return hsdata->event.num_teldata;
		}
	return -1;
}
//-------------------------------------------
// Get list of IDs of telescopes with data
//-------------------------------------------
int get_telescope_with_data_list (int *list){
	if (hsdata != NULL){
		int num_teldata = get_num_teldata ();
		int loop = 0;
		for (loop = 0; loop < num_teldata; loop++){
			*list++ = hsdata->event.teldata_list[loop];
		}
		return 0;
	}
	return -1;
}
//----------------------------------------------------------------
// Returns x,y,z positions of the telescopes [m].
//   x is counted from array reference position towards North,
//   y towards West,
//   z upwards.
// Returns TEL_INDEX_NOT_VALID if telescope index is not valid
// -1 if hsdata == NULL
//----------------------------------------------------------------
int get_telescope_position (int telescope_id, double *pos){
	if (hsdata != NULL){
		int itel = get_telescope_index (telescope_id);
		if (itel == TEL_INDEX_NOT_VALID) return TEL_INDEX_NOT_VALID;
		int loop = 0;
		for (loop = 0; loop < 3; ++loop)  // loop over coordinates
			*pos++ = hsdata->run_header.tel_pos[itel][loop];
		return 0;
	}
	return -1;
}
//-------------------------------------------
// Get number of triggered telescope.
//-------------------------------------------
int get_num_tel_trig (){
	if (hsdata != NULL){
		return hsdata->event.central.num_teltrg;
	}
	else
		return -1;
}
//-------------------------------------------
// Get List of IDs of triggered telescopes.
//-------------------------------------------
int get_central_event_teltrg_list (int *tel_list){
	if (hsdata != NULL)
	{
		int num_teltrig = get_num_tel_trig ();
		int loop = 0;
		for (loop = 0; loop < num_teltrig; loop++){
			*tel_list++ = hsdata->event.central.teltrg_list[loop];
		}
		return 0;
	}
	return -1;
}
//-------------------------------------------
// Get List of relative time of trigger signal
// after correction for nominal delay (in ns) for each
// triggered telescope
//-------------------------------------------
int get_central_event_teltrg_time (float *teltrg_time){
	if (hsdata != NULL)
	{
		int num_teltrig = get_num_tel_trig ();
		int loop = 0;
		for (loop = 0; loop < num_teltrig; loop++){
			*teltrg_time++ = hsdata->event.central.teltrg_time[loop];
		}
		return 0;
	}
	return -1;
}
//-------------------------------------------
// Returns  Global event count
//-------------------------------------------
int get_global_event_count (void){
	if (hsdata != NULL)
		{
		return hsdata->event.central.glob_count;
		}
	return -1;
}
//-------------------------------------------
// Set seconds and nanosecond parameter with
// the central trigger time
//-------------------------------------------
int get_central_event_gps_time (long *seconds, long *nanoseconds){
	if (hsdata != NULL)
		{
		if (seconds != NULL)
			*seconds = hsdata->event.central.gps_time.seconds;
		if (nanoseconds != NULL)
			*nanoseconds = hsdata->event.central.gps_time.nanoseconds;
		return 0;
	}
	return -1;
}
//----------------------------------------------------------------
// Returns the number of different gains per pixel for a telscope id
// Returns TEL_INDEX_NOT_VALID if telescope index is not valid
//----------------------------------------------------------------
int get_num_channel (int telescope_id)
{
	if (hsdata != NULL)
		{
		int itel = get_telescope_index (telescope_id);
		if (itel == TEL_INDEX_NOT_VALID)
		return TEL_INDEX_NOT_VALID;
		return hsdata->camera_org[itel].num_gains;
		}
	return -1;
}
//----------------------------------------------------------------
// Returns Width of readout time slice (i.e. one sample) [ns].
// Returns TEL_INDEX_NOT_VALID if telescope index is not valid
//----------------------------------------------------------------
double get_time_slice (int telescope_id) {
	if (hsdata != NULL)
	{
		int itel = get_telescope_index (telescope_id);
		if (itel == TEL_INDEX_NOT_VALID)
			return TEL_INDEX_NOT_VALID;
		PixelSetting setting = hsdata->pixel_set[itel];
		return setting.time_slice;
	}
	return 0.;
}
//----------------------------------------------------------------
// Fill ref_shapes array with lrefshape values from
// PixelSetting.refshape[channel]
// Returns -1 if channel or fshape are not valid
// Returns TEL_INDEX_NOT_VALID if telescope index is not valid
//----------------------------------------------------------------
int get_ref_shapes (int telescope_id, int channel, double *ref_shapes){
	if (hsdata != NULL && channel < H_MAX_GAINS && ref_shapes != NULL)
	{
		int itel = get_telescope_index (telescope_id);
		if (itel == TEL_INDEX_NOT_VALID)
			return TEL_INDEX_NOT_VALID;
		PixelSetting setting = hsdata->pixel_set[itel];
		int i=0;
		for (i=0; i <= setting.lrefshape; ++i){
			ref_shapes[i] = setting.refshape[channel][i];
		}
		return 0;
	}
	return -1;
}
//----------------------------------------------------------------
// Returns  Reference pulse shape(s)
// refshape[H_MAX_GAINS][H_MAX_FSHAPE];
// If   channel or fshape are not valid return 0.
// Returns TEL_INDEX_NOT_VALID if telescope index is not valid
//----------------------------------------------------------------
double get_ref_shape (int telescope_id, int channel, int fshape){
	if (hsdata != NULL && channel < H_MAX_GAINS && fshape < H_MAX_FSHAPE){
		int itel = get_telescope_index (telescope_id);
		if (itel == TEL_INDEX_NOT_VALID)
			return TEL_INDEX_NOT_VALID;
		PixelSetting setting = hsdata->pixel_set[itel];
		return setting.refshape[channel][fshape];
	}
	return 0.;
}
//----------------------------------------------------------------
// Returns Number of following reference pulse shapes (num_gains or 0)
// Returns TEL_INDEX_NOT_VALID if telescope index is not valid
// Returns -1 if data is not accessible
//----------------------------------------------------------------
int get_nrefshape (int telescope_id){
	if (hsdata != NULL){
		int itel = get_telescope_index (telescope_id);
		if (itel == TEL_INDEX_NOT_VALID)
			return TEL_INDEX_NOT_VALID;
		PixelSetting setting = hsdata->pixel_set[itel];
		return setting.nrefshape;
		}
	return -1;
}
//----------------------------------------------------------------
// Returns Length of following reference pulse shape(s).
// Returns TEL_INDEX_NOT_VALID if telescope index is not valid
// Returns -1 if data is not accessible
//----------------------------------------------------------------
//----------------------------------------------------------------
int get_lrefshape (int telescope_id){
	if (hsdata != NULL)	{
		int itel = get_telescope_index (telescope_id);
		if (itel == TEL_INDEX_NOT_VALID)
			return TEL_INDEX_NOT_VALID;
		PixelSetting setting = hsdata->pixel_set[itel];
		return setting.lrefshape;
	}
	return -1;
}
//----------------------------------------------------------------
// Returns  Time step between refshape entries [ns]
// Returns TEL_INDEX_NOT_VALID if telescope index is not valid
//----------------------------------------------------------------
double get_ref_step (int telescope_id){
	if (hsdata != NULL){
		int itel = get_telescope_index (telescope_id);
		if (itel == TEL_INDEX_NOT_VALID)
			return TEL_INDEX_NOT_VALID;
		PixelSetting setting = hsdata->pixel_set[itel];
		return setting.ref_step;
		}
	return -0.;
}
//----------------------------------------------------------------
// Returns individual channel recorded information ?
// Bit 0: sum, 1: samples, 2: ADC was in saturation.
// Returns TEL_INDEX_NOT_VALID if telescope index is not valid
//----------------------------------------------------------------
uint8_t get_adc_known (int telescope_id, int channel, int pixel_id){
	if (hsdata != NULL && channel < H_MAX_GAINS && pixel_id < H_MAX_PIX){
		int itel = get_telescope_index (telescope_id);
		if (itel == TEL_INDEX_NOT_VALID)
		return TEL_INDEX_NOT_VALID;
			AdcData *raw = hsdata->event.teldata[itel].raw;
		if (raw != NULL && raw->known){
			return raw->adc_known[channel][pixel_id];
		}
	}
	return -1;
}
//----------------------------------------------------------------
// Returns shower primary ID
// 0 (gamma), 1(e-), 2(mu-), 100*A+Z for nucleons and nuclei,
// negative for antimatter.
// Returns -1 if data is not accessible
//----------------------------------------------------------------
int get_mc_shower_primary_id(){
	if ( hsdata != NULL){
		return hsdata->mc_shower.primary_id;
		}
	return -1;
}
//----------------------------------------------------------------
// pe is a output parameter, fill with numbers of photon electron
// Returns  0 on success otherwise -1
//----------------------------------------------------------------
int get_mc_number_photon_electron(int telescope_id, int* pe){


    if (hsdata != NULL){
        int itel = get_telescope_index (telescope_id);
		if (itel == TEL_INDEX_NOT_VALID)
			return TEL_INDEX_NOT_VALID;
		AdcData *raw = hsdata->event.teldata[itel].raw;
		if (raw != NULL && raw->known){	// If triggered telescopes
			int ipix = 0.;
			for (ipix = 0.; ipix < raw->num_pixels; ipix++){ 	//  loop over pixels
				if (raw->significant[ipix]){
					 *pe++ = hsdata->mc_event.mc_pe_list[itel].pe_count[ipix];
				}		// end if raw->significant[ipix]
			}			// end of   loop over pixels
		}			// end if triggered telescopes
		return 0;
	}
	return -1;
}

//----------------------------------------------------------------
// Returns shower number
//----------------------------------------------------------------
int get_mc_shower_num (void){
if ( hsdata != NULL)
	{
		return hsdata->mc_shower.shower_num;
		}
	return -0;
}

//----------------------------------------------------------------
// Returns shower height of first interaction a.s.l. [m]
//----------------------------------------------------------------
double get_mc_shower_h_first_int(){
	if ( hsdata != NULL)
		{
		return hsdata->mc_shower.h_first_int;
		}
	return -0.;
}
//----------------------------------------------------------------
// Returns CORSIKA version  *1000
// Returns -1 if data is not accessible
//----------------------------------------------------------------
int get_corsika_version() {
  if (hsdata != NULL) {
    MCRunHeader *mcrh = &(hsdata)->mc_run_header;
    return mcrh->shower_prog_vers;
  }
  return -1.0;
}
//----------------------------------------------------------------
// Returns sim_telarray version  *1000
// Returns -1 if data is not accessible
//----------------------------------------------------------------
int get_simtel_version() {
  if (hsdata != NULL) {
    MCRunHeader *mcrh = &(hsdata)->mc_run_header;
    return mcrh->detector_prog_vers;
  }
  return -1.0;
}
//----------------------------------------------------------------
// Returns shower altitude [rad]
//----------------------------------------------------------------
double get_mc_shower_altitude (){
	if (hsdata != NULL)
		{
		return hsdata->mc_shower.altitude;
		}
	return -0.;
}
//----------------------------------------------------------------
// Returns shower azimuth (N->E) [rad]
//----------------------------------------------------------------
double get_mc_shower_azimuth (){
	if (hsdata != NULL){
		return hsdata->mc_shower.azimuth;
	}
	return -0.;
}
//----------------------------------------------------------------
// Returns shower primary energy [TeV]
//----------------------------------------------------------------
double get_mc_shower_energy () {
	if (hsdata != NULL){
		return hsdata->mc_shower.energy;
	}
	return -0.;
}
//----------------------------------------------------------------
// Returns shower Xmax
//----------------------------------------------------------------
double get_mc_shower_xmax () {
	if (hsdata != NULL){
		return hsdata->mc_shower.xmax;
	}
	return -0.;
}

//----------------------------------------------------------------
// Returns shower Height of shower maximum [m] in xmax.
//----------------------------------------------------------------
double get_mc_shower_hmax (void){
	if (hsdata != NULL){
		return hsdata->mc_shower.hmax;
	}
	return -0.;
}

//----------------------------------------------------------------
// Returns  Shower number as in shower structure.
//
//----------------------------------------------------------------
int get_mc_event_shower_num(void){
	if (hsdata != NULL){
		return hsdata->mc_event.shower_num;
	}
	return -0.;
}


//----------------------------------------------------------------
// Returns  mc event number -> global counter
//----------------------------------------------------------------
int get_mc_event_num (){
	if (hsdata != NULL){
		return hsdata->mc_event.event;
	}
	return -0.;
}
//----------------------------------------------------------------
// Returns  x core position w.r.t. array reference point [m],
//  x -> N
//----------------------------------------------------------------
double get_mc_event_xcore (){
	if (hsdata != NULL){
		return hsdata->mc_event.xcore;
	}
	return -0.;
}
//----------------------------------------------------------------
// Returns  y core position w.r.t. array reference point [m],
//  y -> W
//----------------------------------------------------------------
double get_mc_event_ycore ()
{
	if (hsdata != NULL){
		return hsdata->mc_event.ycore;
	}
	return -0.;
}

//----------------------------------------------------------------
// Returns number of entries of histogram #6, corresponding to the
// total number of simulated events, including re-uses of showers
// If histogram is not found, it returns -1
// If more than one histogram block is found, the entries of the histogram
// in the different blocks are added. This is the case e.g. for files
// resulting from the merging of several simtelarray outputs
//
long get_mc_num_generated_events(const char *filename)
{
    const long excluded[9] = {1, 2, 3, 4, 7, 11, 12, 21, 22};

    HISTOGRAM* hist = get_histogram_by_ident(6);
    if (hist)
        free_histogram(hist);  // In case it exists from a previous call!

    // read histograms #6 from file, if more than one, add contents:
    read_histogram_file_x(filename, 1, excluded, 9);

    hist = get_histogram_by_ident(6);

    if (!hist){
        Error("Could not find histogram #6 to get number of generated MC events!");
        return -1;
    }
    struct Histogram_Extension *he = hist->extension;

    // Get number of entries of the histogram "Events, without weights (Ra3d, log10(E))"
    // This is the real number of Corsika events processed bysimtelarray (_including_ re-uses of each shower!)
    long numsimushowers = 0;
    int ibin;
    for (ibin = 0; ibin < hist->nbins * hist->nbins_2d; ibin++)
         numsimushowers += he->fdata[ibin];

    return numsimushowers;
}

//----------------------------------------------------------------
// Returns B_total
// Returns -1 if data is not accessible
//----------------------------------------------------------------
double get_B_total() {
  if (hsdata != NULL) {
    MCRunHeader *mcrh = &(hsdata)->mc_run_header;
    return mcrh->B_total;
  }
  return -1.0;
}
//----------------------------------------------------------------
// Returns B_inclination
// Returns -1 if data is not accessible
//----------------------------------------------------------------
double get_B_inclination() {
  if (hsdata != NULL) {
    MCRunHeader *mcrh = &(hsdata)->mc_run_header;
    return mcrh->B_inclination;
  }
  return -1.0;
}
//----------------------------------------------------------------
// Returns B_declination
// Returns -1 if data is not accessible
//----------------------------------------------------------------
double get_B_declination() {
  if (hsdata != NULL) {
    MCRunHeader *mcrh = &(hsdata)->mc_run_header;
    return mcrh->B_declination;
  }
  return -1.0;
}
//----------------------------------------------------------------
// Returns atmosphere
// Returns -1 if data is not accessible
//----------------------------------------------------------------
double get_atmosphere() {
  if (hsdata != NULL) {
    MCRunHeader *mcrh = &(hsdata)->mc_run_header;
    return mcrh->atmosphere;
  }
  return -1;
}
//----------------------------------------------------------------
// Returns corsika_iact_options
// Returns -1 if data is not accessible
//----------------------------------------------------------------
double get_corsika_iact_options() {
  if (hsdata != NULL) {
    MCRunHeader *mcrh = &(hsdata)->mc_run_header;
    return mcrh->corsika_iact_options;
  }
  return -1;
}
//----------------------------------------------------------------
// Returns corsika_low_E_model
// Returns -1 if data is not accessible
//----------------------------------------------------------------
double get_corsika_low_E_model() {
  if (hsdata != NULL) {
    MCRunHeader *mcrh = &(hsdata)->mc_run_header;
    return mcrh->corsika_low_E_model;
  }
  return -1;
}
//----------------------------------------------------------------
// Returns corsika_high_E_model
// Returns -1 if data is not accessible
//----------------------------------------------------------------
double get_corsika_high_E_model() {
  if (hsdata != NULL) {
    MCRunHeader *mcrh = &(hsdata)->mc_run_header;
    return mcrh->corsika_high_E_model;
  }
  return -1;
}
//----------------------------------------------------------------
// Returns corsika_bunchsize
// Returns -1 if data is not accessible
//----------------------------------------------------------------
double get_corsika_bunchsize() {
  if (hsdata != NULL) {
    MCRunHeader *mcrh = &(hsdata)->mc_run_header;
    return mcrh->corsika_bunchsize;
  }
  return -1.0;
}
//----------------------------------------------------------------
// Returns corsika_wlen_min
// Returns -1 if data is not accessible
//----------------------------------------------------------------
double get_corsika_wlen_min() {
  if (hsdata != NULL) {
    MCRunHeader *mcrh = &(hsdata)->mc_run_header;
    return mcrh->corsika_wlen_min;
  }
  return -1.0;
}
//----------------------------------------------------------------
// Returns corsika_wlen_max
// Returns -1 if data is not accessible
//----------------------------------------------------------------
double get_corsika_wlen_max() {
  if (hsdata != NULL) {
    MCRunHeader *mcrh = &(hsdata)->mc_run_header;
    return mcrh->corsika_wlen_max;
  }
  return -1.0;
}
//----------------------------------------------------------------
// Returns corsika_low_E_detail
// Returns -1 if data is not accessible
//----------------------------------------------------------------
double get_corsika_low_E_detail() {
  if (hsdata != NULL) {
    MCRunHeader *mcrh = &(hsdata)->mc_run_header;
    return mcrh->corsika_low_E_detail;
  }
  return -1;
}
//----------------------------------------------------------------
// Returns corsika_high_E_detail
// Returns -1 if data is not accessible
//----------------------------------------------------------------
double get_corsika_high_E_detail() {
  if (hsdata != NULL) {
    MCRunHeader *mcrh = &(hsdata)->mc_run_header;
    return mcrh->corsika_high_E_detail;
  }
  return -1;
}

//----------------------------------------------------------------
// Returns spectral index
// normaly -2
//----------------------------------------------------------------
double get_spectral_index() {
  if (hsdata != NULL) {
    MCRunHeader *mcrh = &(hsdata)->mc_run_header;
    return mcrh->spectral_index;
  }
  return 0.0;
}

//----------------------------------------------------------------
// Returns height
// Returns -1 if data is not accessible
//----------------------------------------------------------------
double get_mc_obsheight() {
  if (hsdata != NULL) {
    MCRunHeader *mcrh = &(hsdata)->mc_run_header;
    return mcrh->obsheight;
  }
  return -1.0;
}

//----------------------------------------------------------------
// Returns mc_num_showers
// Returns -1 if data is not accessible
//----------------------------------------------------------------
int get_mc_num_showers() {
  if (hsdata != NULL) {
    MCRunHeader *mcrh = &(hsdata)->mc_run_header;
    return mcrh->num_showers;
  }
  return -1;
}

//----------------------------------------------------------------
// Returns mc_num_use
// Returns -1 if data is not accessible
//----------------------------------------------------------------
int get_mc_num_use() {
  if (hsdata != NULL) {
    MCRunHeader *mcrh = &(hsdata)->mc_run_header;
    return mcrh->num_use;
  }
  return -1;
}

//----------------------------------------------------------------
// Returns mc_core_pos_mode
// Returns -1 if data is not accessible
//----------------------------------------------------------------
int get_mc_core_pos_mode() {
  if (hsdata != NULL) {
    MCRunHeader *mcrh = &(hsdata)->mc_run_header;
    return mcrh->core_pos_mode;
  }
  return -1;
}

//----------------------------------------------------------------
// Returns mc_core_range_X bzw mc_core_range_min
// Returns -1 if data is not accessible
//----------------------------------------------------------------
double get_mc_core_range_X() {
  if (hsdata != NULL) {
    MCRunHeader *mcrh = &(hsdata)->mc_run_header;
    return mcrh->core_range[0];
  }
  return -1.0;
}
//----------------------------------------------------------------
// Returns mc_core_range_Y bzw. mc_core_range_max
// Returns -1 if data is not accessible
//----------------------------------------------------------------
double get_mc_core_range_Y() {
  if (hsdata != NULL) {
    MCRunHeader *mcrh = &(hsdata)->mc_run_header;
    return mcrh->core_range[1];
  }
  return -1.0;
}

//----------------------------------------------------------------
// Returns mc_alt_range_Min [TeV]
// Returns -1 if data is not accessible
//----------------------------------------------------------------
double get_mc_alt_range_Min() {
  if (hsdata != NULL) {
    MCRunHeader *mcrh = &(hsdata)->mc_run_header;
    return mcrh->alt_range[0];
  }
  return -1.0;
}
//----------------------------------------------------------------
// Returns mc_alt_range_Max [TeV]
// Returns -1 if data is not accessible
//----------------------------------------------------------------
double get_mc_alt_range_Max() {
  if (hsdata != NULL) {
    MCRunHeader *mcrh = &(hsdata)->mc_run_header;
    return mcrh->alt_range[1];
  }
  return -1.0;
}

//----------------------------------------------------------------
// Returns mc_az_range_Min
// Returns -1 if data is not accessible
//----------------------------------------------------------------
double get_mc_az_range_Min() {
  if (hsdata != NULL) {
    MCRunHeader *mcrh = &(hsdata)->mc_run_header;
    return mcrh->az_range[0];
  }
  return -1.0;
}
//----------------------------------------------------------------
// Returns mc_az_range_Max
// Returns -1 if data is not accessible
//----------------------------------------------------------------
double get_mc_az_range_Max() {
  if (hsdata != NULL) {
    MCRunHeader *mcrh = &(hsdata)->mc_run_header;
    return mcrh->az_range[1];
  }
  return -1.0;
}

//----------------------------------------------------------------
// Returns mc_viewcone_Min
// Returns -1 if data is not accessible
//----------------------------------------------------------------
double get_mc_viewcone_Min() {
  if (hsdata != NULL) {
    MCRunHeader *mcrh = &(hsdata)->mc_run_header;
    return mcrh->viewcone[0];
  }
  return -1.0;
}
//----------------------------------------------------------------
// Returns mc_viewcone_Max
// Returns -1 if data is not accessible
//----------------------------------------------------------------
double get_mc_viewcone_Max() {
  if (hsdata != NULL) {
    MCRunHeader *mcrh = &(hsdata)->mc_run_header;
    return mcrh->viewcone[1];
  }
  return -1.0;
}

//----------------------------------------------------------------
// Returns mc_E_range_Min
// Returns -1 if data is not accessible
//----------------------------------------------------------------
double get_mc_E_range_Min() {
  if (hsdata != NULL) {
    MCRunHeader *mcrh = &(hsdata)->mc_run_header;
    return mcrh->E_range[0];
  }
  return -1.0;
}
//----------------------------------------------------------------
// Returns mc_E_range_Max
// Returns -1 if data is not accessible
//----------------------------------------------------------------
double get_mc_E_range_Max() {
  if (hsdata != NULL) {
    MCRunHeader *mcrh = &(hsdata)->mc_run_header;
    return mcrh->E_range[1];
  }
  return -1.0;
}

//----------------------------------------------------------------
// Returns mc_diffuse
// Returns -1 if data is not accessible
//----------------------------------------------------------------
int get_mc_diffuse() {
  if (hsdata != NULL) {
    MCRunHeader *mcrh = &(hsdata)->mc_run_header;
    return mcrh->diffuse;
  }
  return -1;
}

//----------------------------------------------------------------
// Returns mc_injection_height
// Returns -1 if data is not accessible
//----------------------------------------------------------------
double get_mc_injection_height() {
  if (hsdata != NULL) {
    MCRunHeader *mcrh = &(hsdata)->mc_run_header;
    return mcrh->injection_height;
  }
  return -1.0;
}
//----------------------------------------------------------------
// Returns the tracking/pointing direction in [radians]. Depending on
// "tracking_mode" this either contains:
//   [0]=Azimuth, [1]=Altitude in mode 0,
//   OR
//   [0]=R.A., [1]=Declination in mode 1.
// -1 if hsdata == NULL
//----------------------------------------------------------------
int get_mc_run_array_direction (double *dir){
	if (hsdata != NULL){
		int loop = 0;
		for (loop = 0; loop < 2; ++loop)  // loop over coordinates
			*dir++ = hsdata->run_header.direction[loop];
		return 0;
	}
	return -1;
}
//----------------------------------------------------------------
// Returns the Raw azimuth angle [radians from N->E] for the telescope
// Returns TEL_INDEX_NOT_VALID if telescope index is not valid
//----------------------------------------------------------------
double get_azimuth_raw (int telescope_id)
{
	if (hsdata != NULL)
		{
		int itel = get_telescope_index (telescope_id);
        if (itel == TEL_INDEX_NOT_VALID)
            return TEL_INDEX_NOT_VALID;
        return hsdata->event.trackdata[itel].azimuth_raw;
		}
	return -1;
}
//----------------------------------------------------------------
// Returns the Raw altitude angle [radians] for the telescope
// Returns TEL_INDEX_NOT_VALID if telescope index is not valid
//----------------------------------------------------------------
double get_altitude_raw (int telescope_id)
{
	if (hsdata != NULL)
		{
		int itel = get_telescope_index (telescope_id);
        if (itel == TEL_INDEX_NOT_VALID)
            return TEL_INDEX_NOT_VALID;
        return hsdata->event.trackdata[itel].altitude_raw;
		}
	return -1;
}
//----------------------------------------------------------------
// Returns the tracking Azimuth corrected for pointing errors for the telescope
// Returns TEL_INDEX_NOT_VALID if telescope index is not valid
//----------------------------------------------------------------
double get_azimuth_cor (int telescope_id)
{
	if (hsdata != NULL)
		{
		int itel = get_telescope_index (telescope_id);
        if (itel == TEL_INDEX_NOT_VALID)
            return TEL_INDEX_NOT_VALID;
        return hsdata->event.trackdata[itel].azimuth_cor;
		}
	return -1;
}
//----------------------------------------------------------------
// Returns the tracking Altitude corrected for pointing errors for the telescope
// Returns TEL_INDEX_NOT_VALID if telescope index is not valid
//----------------------------------------------------------------
double get_altitude_cor (int telescope_id)
{
	if (hsdata != NULL)
		{
		int itel = get_telescope_index (telescope_id);
        if (itel == TEL_INDEX_NOT_VALID)
            return TEL_INDEX_NOT_VALID;
        return hsdata->event.trackdata[itel].altitude_cor;
		}
	return -1;
}
//----------------------------------------------------------------
// Returns the offset of pointing direction in camera f.o.v.
// divided by focal length, i.e. converted to radians:
//   [0] = Camera x (downwards in normal pointing, i.e. increasing Alt)
//   [1] = Camera y -> Az.
// -1 if hsdata == NULL
//----------------------------------------------------------------
int get_mc_event_offset_fov (double *off){
	if (hsdata != NULL){
		int loop = 0;
		for (loop = 0; loop < 2; ++loop)  // loop over coordinates
			*off++ = hsdata->run_header.offset_fov[loop];
		return 0;
	}
	return -1;
}
//-------------------------------------------
// Returns  PixelTiming.timval[H_MAX_PIX][H_MAX_PIX_TIMES]
// Returns TEL_INDEX_NOT_VALID if telescope index is not valid
//-------------------------------------------
int get_pixel_timing_timval (int telescope_id, float *data){
	if (hsdata != NULL){
		int itel = get_telescope_index (telescope_id);
		if (itel == TEL_INDEX_NOT_VALID)
			return TEL_INDEX_NOT_VALID;
		PixelTiming *pt = hsdata->event.teldata[itel].pixtm;
		if (pt != NULL){
			int ipix = 0;
			for (ipix = 0; ipix < pt->num_pixels; ipix++){
				int itimes = 0;
				for (itimes = 0;
				itimes < pt->num_types && itimes < H_MAX_PIX_TIMES;
				itimes++){
					*data++ = pt->timval[ipix][itimes];
				}
			}			// end for ipix
		}			// end if pt != NULL
		return 0;
	}				// end if hsdata
	return -1;
}
//----------------------------------------------------------------
// Returns Pulses sampled
// Returns TEL_INDEX_NOT_VALID if telescope index is not valid
//----------------------------------------------------------------
int get_adc_sample (int telescope_id, int channel, uint16_t * data){
	if (hsdata != NULL){
		int itel = get_telescope_index (telescope_id);
		if (itel == TEL_INDEX_NOT_VALID)
			return TEL_INDEX_NOT_VALID;
		AdcData *raw = hsdata->event.teldata[itel].raw;
		if (raw != NULL && raw->known){	// If triggered telescopes
			int ipix = 0.;
			for (ipix = 0.; ipix < raw->num_pixels; ipix++){ 	//  loop over pixels
				if (raw->significant[ipix]){
					int isamp = 0.;
					for (isamp = 0.; isamp < raw->num_samples; isamp++){
						*data++ = raw->adc_sample[channel][ipix][isamp];
					}
				}		// end if raw->significant[ipix]
			}			// end of   loop over pixels
		}			// end if triggered telescopes
		return 0;
	}
	return -1;
}

//----------------------------------------------------------------
// Returns Was amplitude large enough to record it? Bit 0: sum, 1: samples.
// Returns TEL_INDEX_NOT_VALID if telescope index is not valid
//----------------------------------------------------------------
uint8_t get_significant (int telescope_id, uint8_t * data){
	if (hsdata != NULL){
		int itel = get_telescope_index (telescope_id);
		if (itel == TEL_INDEX_NOT_VALID)
			return TEL_INDEX_NOT_VALID;
		AdcData *raw = hsdata->event.teldata[itel].raw;
		if (raw != NULL ){
			int ipix = 0;
			for (ipix = 0; ipix < raw->num_pixels; ipix++){ 	//  loop over pixels
				*data++ = raw->significant[ipix];
			}			// end of   loop over pixels
			return 0;
		}

	}
	return -1;
}


//----------------------------------------------------------------
// Return adc sum for corresponding telescope and channel (HI_GAIN/LOW_GAIN)
// Returns TEL_INDEX_NOT_VALID if telescope index is not valid
//----------------------------------------------------------------
int get_adc_sum (int telescope_id, int channel, uint32_t * data){
	if (hsdata != NULL){
		int itel = get_telescope_index (telescope_id);
		if (itel == TEL_INDEX_NOT_VALID)
			return TEL_INDEX_NOT_VALID;
		AdcData *raw = hsdata->event.teldata[itel].raw;
		if (raw != NULL && raw->known){	// If triggered telescopes
			int ipix = 0.;
			for (ipix = 0.; ipix < raw->num_pixels; ipix++){	//  loop over pixels
				*data++ = raw->adc_sum[channel][ipix];
			}			// end of   loop over pixels
			return 0;
		}			// end if triggered telescopes
	}
	return -1;
}
//----------------------------------------------------------------
// Fill calibration data
//  double calib[H_MAX_GAINS][H_MAX_PIX]; /**< ADC to laser/LED p.e. conversion,
// Returns  0 for success,  TEL_INDEX_NOT_VALID if telescope index is not valid
//
int get_calibration (int telescope_id, float *calib)
//----------------------------------------------------------------
{
	if (hsdata != NULL)
	{
		int itel = get_telescope_index (telescope_id);
		if (itel == TEL_INDEX_NOT_VALID)
		return TEL_INDEX_NOT_VALID;
		LasCalData calibration = hsdata->tel_lascal[itel];
		unsigned int num_gain = calibration.num_gains;
		unsigned int num_pixels = hsdata->camera_set[itel].num_pixels;
		unsigned int igain = 0;
		for (igain = 0; igain < num_gain; igain++) {
			unsigned int ipix = 0.;
			for (ipix = 0.; ipix < num_pixels; ipix++) {	// loop over pixels
				*calib++ = (float)calibration.calib[igain][ipix];
			}// end loop pixels
		}// end of   loop over gains
		return 0;
	}
	return -1;
}
//----------------------------------------------------------------
// Fill pedestal data
//  double pedestal[H_MAX_GAINS][H_MAX_PIX];  ///< Average pedestal on ADC sums
// Returns 0 for success TEL_INDEX_NOT_VALID if telescope index is not valid
//
int get_pedestal (int telescope_id, float *pedestal)
//----------------------------------------------------------------
{
	if (hsdata != NULL)
		{
		int itel = get_telescope_index (telescope_id);
		if (itel == TEL_INDEX_NOT_VALID)
		return TEL_INDEX_NOT_VALID;
		TelMoniData monitor = hsdata->tel_moni[itel];
		unsigned int num_gain = monitor.num_gains;
		unsigned int num_pixels = hsdata->camera_set[itel].num_pixels;
		unsigned int igain = 0.;
        for (igain = 0; igain < num_gain; igain++) // loop over channel
		{
    		unsigned int ipix = 0;
    		for (ipix = 0.; ipix < num_pixels; ipix++)
    			{
    			*pedestal++ = (float)monitor.pedestal[igain][ipix];
    			}
		}
		return 0;
		}
	return -1;
}
//----------------------------------------------------------------
// Returns pixel position information
// Returns TEL_INDEX_NOT_VALID if telescope index is not valid
// -1 if hsdata == NULL
//----------------------------------------------------------------
int get_pixel_position (int telescope_id, double *xpos, double *ypos)
{
	if (hsdata != NULL)
		{
		int itel = get_telescope_index (telescope_id);
		if (itel == TEL_INDEX_NOT_VALID)
		return TEL_INDEX_NOT_VALID;
		int ipix = 0.;
		int num_pixels = hsdata->camera_set[itel].num_pixels;
		for (ipix = 0.; ipix < num_pixels; ipix++)	// loop over pixels
		{
		*xpos++ = hsdata->camera_set[itel].xpix[ipix];
		*ypos++ = hsdata->camera_set[itel].ypix[ipix];
		}
		return 0;
		}
	return -1;
}
//----------------------------------------------------------------
// Returns the shape of the pixels (0: circ., 1,3: hex, 2: square, -1: unknown)
// Returns TEL_INDEX_NOT_VALID if telescope index is not valid
// -1 if hsdata == NULL
//---------------------------------------------------------------
int get_pixel_shape (int telescope_id, double *pixel_shape)
{
        if (hsdata != NULL)
		{
		int itel = get_telescope_index (telescope_id);
		if (itel == TEL_INDEX_NOT_VALID)
		return TEL_INDEX_NOT_VALID;
		int ipix = 0.;
		int num_pixels = hsdata->camera_set[itel].num_pixels;
		for (ipix = 0.; ipix < num_pixels; ipix++)	// loop over pixels
		{
		*pixel_shape++ = hsdata->camera_set[itel].pixel_shape[ipix];
		}
		return 0;
		}
	return -1;
}
//----------------------------------------------------------------
// Returns the area of the pixels ([m^2])
// Returns TEL_INDEX_NOT_VALID if telescope index is not valid
// -1 if hsdata == NULL
//---------------------------------------------------------------
int get_pixel_area (int telescope_id, double *pixel_area)
{
        if (hsdata != NULL)
		{
		int itel = get_telescope_index (telescope_id);
		if (itel == TEL_INDEX_NOT_VALID)
		return TEL_INDEX_NOT_VALID;
		int ipix = 0.;
		int num_pixels = hsdata->camera_set[itel].num_pixels;
		for (ipix = 0.; ipix < num_pixels; ipix++)	// loop over pixels
		{
		*pixel_area++ = hsdata->camera_set[itel].area[ipix];
		}
		return 0;
		}
	return -1;
}
//----------------------------------------------------------------
// Returns the number of pixels in the camera (as in configuration)
// Returns TEL_INDEX_NOT_VALID if telescope index is not valid
//----------------------------------------------------------------
int get_num_pixels (int telescope_id)
{
	if (hsdata != NULL)
		{
		int itel = get_telescope_index (telescope_id);
		if (itel == TEL_INDEX_NOT_VALID)
		return TEL_INDEX_NOT_VALID;
		return hsdata->camera_set[itel].num_pixels;
		}
	return -1;
}

//----------------------------------------------------------------
// Returns the number of pixels used in the camera trigger
// Returns TEL_INDEX_NOT_VALID if telescope index is not valid
//----------------------------------------------------------------
int get_num_trig_pixels (int telescope_id)
{
	if (hsdata != NULL)
		{
		    int itel = get_telescope_index (telescope_id);
		    if (itel == TEL_INDEX_NOT_VALID)
                return TEL_INDEX_NOT_VALID;
            return hsdata->event.teldata[itel].trigger_pixels.pixels;
		}
    return -1;
}

//----------------------------------------------------------------
// Returns a list of pixels used in the camera trigger
// Returns TEL_INDEX_NOT_VALID if telescope index is not valid
//----------------------------------------------------------------
int get_trig_pixels (int telescope_id, int *trigpix)
{
	if (hsdata != NULL)
		{
		    int itel = get_telescope_index (telescope_id);
		    if (itel == TEL_INDEX_NOT_VALID)
                return TEL_INDEX_NOT_VALID;
            int npix = hsdata->event.teldata[itel].trigger_pixels.pixels;
            int ipix = 0;
            for (ipix=0.; ipix < npix; ipix++)
            {
                *trigpix++ = hsdata->event.teldata[itel].trigger_pixels.pixel_list[ipix];
            }
            return 0;
        }
    return -1;
}

//----------------------------------------------------------------
// Returns total area of individual mirrors corrected   for inclination [m^2].
// Returns TEL_INDEX_NOT_VALID if telescope index is not valid
//----------------------------------------------------------------
int  get_mirror_area (int telescope_id, double *result)
{
	if (hsdata != NULL && result != NULL)
		{
		int itel = get_telescope_index (telescope_id);
		if (itel == TEL_INDEX_NOT_VALID)
		return TEL_INDEX_NOT_VALID;
		*result = hsdata->camera_set[itel].mirror_area;
		return 0;
		}
	return -1.;
}
//-----------------------------------------------------
// Returns the number of samples (time slices) recorded.
// Only contains the number of samples for the telescopes that have
// data in the current event.
// Returns TEL_INDEX_NOT_VALID if telescope index is not valid
// Otherwise 0
//-----------------------------------------------------
int get_event_num_samples (int telescope_id)
{
	if (hsdata != NULL)
		{
		int itel = get_telescope_index (telescope_id);
		if (itel == TEL_INDEX_NOT_VALID)
		return TEL_INDEX_NOT_VALID;
		AdcData *raw = hsdata->event.teldata[itel].raw;
		if (raw != NULL)		//&& raw->known   )
		{
		return raw->num_samples;
		}
		}
	return -1;
}

//-----------------------------------------------------
// Returns the desired or used zero suppression mode.
// Returns TEL_INDEX_NOT_VALID if telescope index is not valid
//-----------------------------------------------------
int get_zero_sup_mode(int telescope_id,int* mode){
	if (hsdata != NULL)
	{
		int itel = get_telescope_index (telescope_id);
		if (itel == TEL_INDEX_NOT_VALID)
			return TEL_INDEX_NOT_VALID;
		AdcData *raw = hsdata->event.teldata[itel].raw;
		if (raw != NULL)
		{
			*mode = raw->zero_sup_mode;
			return 0;
		}
	}
	return -1;
}

//-----------------------------------------------------
// Returns the desired or used zero suppression mode.
// Returns TEL_INDEX_NOT_VALID if telescope index is not valid
// Otherwise 0
//-----------------------------------------------------
int get_data_red_mode(int telescope_id,int* mode){
	if (hsdata != NULL)
	{
		int itel = get_telescope_index (telescope_id);
		if (itel == TEL_INDEX_NOT_VALID)
			return TEL_INDEX_NOT_VALID;
		AdcData *raw = hsdata->event.teldata[itel].raw;
		if (raw != NULL)
		{
			*mode = raw->data_red_mode;
			return 0;
		}
	}
	return -1;
}


//-----------------------------------------------------
// Returns the number of different types of times can we store
// Returns TEL_INDEX_NOT_VALID if telescope index is not valid
//-----------------------------------------------------
int get_pixel_timing_num_times_types (int telescope_id)
{
	if (hsdata != NULL)
	{
		int itel = get_telescope_index (telescope_id);
		if (itel == TEL_INDEX_NOT_VALID)
			return TEL_INDEX_NOT_VALID;
		PixelTiming *pt = hsdata->event.teldata[itel].pixtm;
		if (pt != NULL)
		{
			return pt->num_types;
		}
	}
	return -1;
}
//-----------------------------------------------------
// Set the local telescope trigger time second and nanoseconds
// to function paramameters
// returns 0 if set, otherwise returns -1
// Returns TEL_INDEX_NOT_VALID if telescope index is not valid
//-----------------------------------------------------
int get_tel_event_gps_time (int telescope_id, long *seconds, long *nanoseconds)
{
	if (hsdata != NULL && seconds != NULL && nanoseconds != NULL)
		{
		int itel = get_telescope_index (telescope_id);
		if (itel == TEL_INDEX_NOT_VALID)
		return TEL_INDEX_NOT_VALID;
		*seconds = hsdata->event.teldata[itel].gps_time.seconds;
		*nanoseconds = hsdata->event.teldata[itel].gps_time.nanoseconds;
		return 0;
		}
	return -1;
}
//---------------------------------------------
// Returns PixelTiming threshold:
//  - Minimum base-to-peak raw amplitude difference applied in pixel selection
// Returns TEL_INDEX_NOT_VALID if telescope index is not valid
//-----------------------------------------------------
int get_pixel_timing_threshold (int telescope_id, int *result)
{
	if (hsdata != NULL && result != NULL)
		{
		int itel = get_telescope_index (telescope_id);
		if (itel == TEL_INDEX_NOT_VALID)
		return TEL_INDEX_NOT_VALID;
		PixelTiming *pt = hsdata->event.teldata[itel].pixtm;
		if (pt != NULL)
		*result = pt->threshold;
		return 0;
		}
	return -1;
}
//---------------------------------------------
// Returns PixelTiming peak_global:
//  Camera-wide (mean) peak position [time slices]
// Returns TEL_INDEX_NOT_VALID if telescope index is not valid
//-----------------------------------------------------
int get_pixel_timing_peak_global (int telescope_id, float *result)
{
	if (hsdata != NULL && result != NULL)
		{
		int itel = get_telescope_index (telescope_id);
		if (itel == TEL_INDEX_NOT_VALID)
		return TEL_INDEX_NOT_VALID;
		PixelTiming *pt = hsdata->event.teldata[itel].pixtm;
		if (pt != NULL)
		{
		*result = pt->peak_global;
		}
		return 0;
		}
	return -1;
}
//--------------------------------------------------
// fill hsdata global variable by decoding data file
//--------------------------------------------------
int fill_hsdata (int *event_id)	//,int *header_readed)
{
	int itel;
	int rc = 0;
	int ignore = 0;
	int tel_id;
	/* Find and read the next block of data. */
	/* In case of problems with the data, just give up. */
	if (find_io_block (iobuf, &item_header) != 0){
		return -1;
	}
	if (read_io_block (iobuf, &item_header) != 0){
		return -1;
	}

	// if ( ( !header_readed) &&
	if (hsdata == NULL &&
		item_header.type > IO_TYPE_HESS_RUNHEADER &&
		item_header.type < IO_TYPE_HESS_RUNHEADER + 200){
			fprintf (stderr, "Trying to read event data before run header.\n");
			fprintf (stderr, "Skipping this data block.\n");
			return (int) item_header.type;
	}
	switch ((int) item_header.type){
		/* =================================================== */
		case IO_TYPE_HESS_RUNHEADER:
			/* Structures might be allocated from previous run */
			if (hsdata != NULL){
				/* Free memory allocated inside ... */
				free_hsdata();
			}
			hsdata = (AllHessData *) calloc (1, sizeof (AllHessData));
			if (hsdata == NULL) {
			        Warning ("Memory Allocation Failed. Please free up RAM");
				exit (1);
			}
			if ((rc = read_hess_runheader (iobuf, &(hsdata)->run_header)) < 0){
				Warning ("Reading run header failed.");
				exit (1);
			}
			for (itel = 0; itel < (hsdata)->run_header.ntel; itel++){
				tel_id = (hsdata)->run_header.tel_id[itel];
				(hsdata)->camera_set[itel].tel_id = tel_id;
				(hsdata)->camera_org[itel].tel_id = tel_id;
				(hsdata)->pixel_set[itel].tel_id = tel_id;
				(hsdata)->pixel_disabled[itel].tel_id = tel_id;
				(hsdata)->cam_soft_set[itel].tel_id = tel_id;
				(hsdata)->tracking_set[itel].tel_id = tel_id;
				(hsdata)->point_cor[itel].tel_id = tel_id;
				(hsdata)->event.num_tel = (hsdata)->run_header.ntel;
				(hsdata)->event.teldata[itel].tel_id = tel_id;
				(hsdata)->event.trackdata[itel].tel_id = tel_id;
				if (((hsdata)->event.teldata[itel].raw =
					(AdcData *) calloc (1, sizeof (AdcData))) == NULL){
					Warning ("Not enough memory");
					exit (1);
					}
				(hsdata)->event.teldata[itel].raw->tel_id = tel_id;
				if (((hsdata)->event.teldata[itel].pixtm =
					(PixelTiming *) calloc (1, sizeof (PixelTiming))) == NULL){
					Warning ("Not enough memory");
					exit (1);
				}
				(hsdata)->event.teldata[itel].pixtm->tel_id = tel_id;
				if (((hsdata)->event.teldata[itel].img =
					(ImgData *) calloc (2, sizeof (ImgData))) == NULL){
					Warning ("Not enough memory");
					exit (1);
				}
				(hsdata)->event.teldata[itel].max_image_sets = 2;
				(hsdata)->event.teldata[itel].img[0].tel_id = tel_id;
				(hsdata)->event.teldata[itel].img[1].tel_id = tel_id;
				(hsdata)->tel_moni[itel].tel_id = tel_id;
				(hsdata)->tel_lascal[itel].tel_id = tel_id;
			}
			break;
		// end case IO_TYPE_HESS_RUNHEADER:
		/* =================================================== */
		case IO_TYPE_HESS_MCRUNHEADER:
			rc = read_hess_mcrunheader (iobuf, &(hsdata)->mc_run_header);
			break;
		/* =================================================== */
		case IO_TYPE_MC_INPUTCFG:
			break;
		/* =================================================== */
		case 70:			/* How sim_hessarray was run and how it was configured. */
                        if ( showhistory ) 
                            list_history(iobuf,NULL);
                            rc = 70;
			break;
		/* =================================================== */
		case IO_TYPE_HESS_CAMSETTINGS:
			tel_id = item_header.ident;	// Telescope ID is in the header
			if ((itel = find_tel_idx (tel_id)) < 0){
			char msg[256];
			snprintf (msg, sizeof (msg) - 1,
					"Camera settings for unknown telescope %d.", tel_id);
			Warning (msg);
			exit (1);
			}
			rc = read_hess_camsettings (iobuf, &(hsdata)->camera_set[itel]);
			break;
		/* =================================================== */
		case IO_TYPE_HESS_CAMORGAN:
			tel_id = item_header.ident;	// Telescope ID is in the header
			if ((itel = find_tel_idx (tel_id)) < 0){
				char msg[256];
				snprintf (msg, sizeof (msg) - 1,
						"Camera organisation for unknown telescope %d.", tel_id);
				Warning (msg);
				exit (1);
			}
			rc = read_hess_camorgan (iobuf, &(hsdata)->camera_org[itel]);
			break;
		/* =================================================== */
		case IO_TYPE_HESS_PIXELSET:
			tel_id = item_header.ident;	// Telescope ID is in the header
			if ((itel = find_tel_idx (tel_id)) < 0){
				char msg[256];
				snprintf (msg, sizeof (msg) - 1,
						"Pixel settings for unknown telescope %d.", tel_id);
				Warning (msg);
				exit (1);
			}
			rc = read_hess_pixelset (iobuf, &(hsdata)->pixel_set[itel]);
			break;
		/* =================================================== */
		case IO_TYPE_HESS_PIXELDISABLE:
			tel_id = item_header.ident;	// Telescope ID is in the header
			if ((itel = find_tel_idx (tel_id)) < 0){
				char msg[256];
				snprintf (msg, sizeof (msg) - 1,
						"Pixel disable block for unknown telescope %d.", tel_id);
				Warning (msg);
				exit (1);
			}
			rc = read_hess_pixeldis (iobuf, &(hsdata)->pixel_disabled[itel]);
			break;
		/* =================================================== */
		case IO_TYPE_HESS_CAMSOFTSET:
			tel_id = item_header.ident;	// Telescope ID is in the header
			if ((itel = find_tel_idx (tel_id)) < 0){
				char msg[256];
				snprintf (msg, sizeof (msg) - 1,
						"Camera software settings for unknown telescope %d.",
						tel_id);
				Warning (msg);
				exit (1);
			}
			rc = read_hess_camsoftset (iobuf, &(hsdata)->cam_soft_set[itel]);
			break;
		/* =================================================== */
		case IO_TYPE_HESS_POINTINGCOR:
			tel_id = item_header.ident;	// Telescope ID is in the header
			if ((itel = find_tel_idx (tel_id)) < 0){
				char msg[256];
				snprintf (msg, sizeof (msg) - 1,
						"Pointing correction for unknown telescope %d.", tel_id);
				Warning (msg);
				exit (1);
			}
			rc = read_hess_pointingcor (iobuf, &(hsdata)->point_cor[itel]);
			break;
		/* =================================================== */
		case IO_TYPE_HESS_TRACKSET:
			tel_id = item_header.ident;	// Telescope ID is in the header
			if ((itel = find_tel_idx (tel_id)) < 0){
				char msg[256];
				snprintf (msg, sizeof (msg) - 1,
						"Tracking settings for unknown telescope %d.", tel_id);
				Warning (msg);
				exit (1);
			}
			rc = read_hess_trackset (iobuf, &(hsdata)->tracking_set[itel]);
			break;
		/* =================================================== */
		/* =============   IO_TYPE_HESS_EVENT  =============== */
		/* =================================================== */
		case IO_TYPE_HESS_EVENT:
			rc = read_hess_event (iobuf, &(hsdata)->event, -1);
			*event_id = item_header.ident;
			break;
		/* =================================================== */
		case IO_TYPE_HESS_CALIBEVENT:
		{
                        int type = -1;
			rc = read_hess_calib_event(iobuf, &(hsdata)->event, -1, &type);
			*event_id = item_header.ident;
		}
		break;
		/* =================================================== */
		case IO_TYPE_HESS_MC_SHOWER:
			rc = read_hess_mc_shower (iobuf, &(hsdata)->mc_shower);
			break;
		/* =================================================== */
		case IO_TYPE_HESS_MC_EVENT:
			rc = read_hess_mc_event (iobuf, &(hsdata)->mc_event);
			*event_id = item_header.ident;
			break;
		/* =================================================== */
		case IO_TYPE_MC_TELARRAY:
			if (hsdata && (hsdata)->run_header.ntel > 0){
				rc = read_hess_mc_phot (iobuf, &(hsdata)->mc_event);
			}
			break;
		/* =================================================== */
		/* With extended output option activated, the particles
			arriving at ground level would be stored as seemingly
			stray photon bunch block. */
		case IO_TYPE_MC_PHOTONS:
			break;
		/* =================================================== */
		case IO_TYPE_MC_RUNH:
		case IO_TYPE_MC_EVTH:
		case IO_TYPE_MC_EVTE:
		case IO_TYPE_MC_RUNE:
			break;
		/* =================================================== */
		case IO_TYPE_HESS_MC_PE_SUM:
			rc = read_hess_mc_pe_sum (iobuf, &(hsdata)->mc_event.mc_pesum);
			break;
		/* =================================================== */
		case IO_TYPE_HESS_TEL_MONI:
			// Telescope ID among others in the header
			tel_id = (item_header.ident & 0xff) |
			((item_header.ident & 0x3f000000) >> 16);
			if ((itel = find_tel_idx (tel_id)) < 0){
				char msg[256];
				snprintf (msg, sizeof (msg) - 1,
						"Telescope monitor block for unknown telescope %d.",
						tel_id);
				Warning (msg);
				exit (1);
			}
			rc = read_hess_tel_monitor (iobuf, &(hsdata)->tel_moni[itel]);
			break;
		/* =================================================== */
		case IO_TYPE_HESS_LASCAL:
			tel_id = item_header.ident;	// Telescope ID is in the header
			if ((itel = find_tel_idx (tel_id)) < 0){
				char msg[256];
				snprintf (msg, sizeof (msg) - 1,
						"Laser/LED calibration for unknown telescope %d.",
						tel_id);
				Warning (msg);
				exit (1);
			}
			rc = read_hess_laser_calib (iobuf, &hsdata->tel_lascal[itel]);
			break;
		/* =================================================== */
		case IO_TYPE_HESS_RUNSTAT:
		rc = read_hess_run_stat (iobuf, &hsdata->run_stat);
		break;
		/* =================================================== */
		case IO_TYPE_HESS_MC_RUNSTAT:
			rc = read_hess_mc_run_stat (iobuf, &hsdata->mc_run_stat);
			break;
		/* (End-of-job or DST) histograms */
		case 100:
			break;
		default:
		if (!ignore) fprintf (stderr, "WARNING: Ignoring unknown data block type %ld\n",
			item_header.type);
		}				// end switch item_header.type
	/* What did we actually get? */
	return (int) item_header.type;
}
//-----------------------------------
// Free hsdata structure
//-----------------------------------
void free_hsdata(void)
{
	int itel=0;
	/* Structures might be allocated from previous run */
	if ( hsdata != NULL )
	{
		/* Free memory allocated inside ... */
		for (itel=0; itel<hsdata->run_header.ntel; itel++){
			if ( hsdata->event.teldata[itel].raw != NULL ){
				free(hsdata->event.teldata[itel].raw);
				hsdata->event.teldata[itel].raw = NULL;
			}
			if ( hsdata->event.teldata[itel].pixtm != NULL ){
				free(hsdata->event.teldata[itel].pixtm);
				hsdata->event.teldata[itel].pixtm = NULL;
			}
			if ( hsdata->event.teldata[itel].img != NULL ){
				free(hsdata->event.teldata[itel].img);
				hsdata->event.teldata[itel].img = NULL;
			}
			if ( hsdata->event.teldata[itel].pixcal != NULL ){
				free(hsdata->event.teldata[itel].pixcal);
				hsdata->event.teldata[itel].pixcal = NULL;
			}
		}
		if ( hsdata->run_header.target != NULL )
		{
			free(hsdata->run_header.target);
			hsdata->run_header.target = NULL;
		}
		if ( hsdata->run_header.observer != NULL )
		{
			free(hsdata->run_header.observer);
			hsdata->run_header.observer = NULL;
		}
		int j;
		for ( j=0; j<H_MAX_PROFILE; j++)
		{
			if ( hsdata->mc_shower.profile[j].content != NULL )
			{
				free(hsdata->mc_shower.profile[j].content);
				hsdata->mc_shower.profile[j].content = NULL;
			}
		}

		/* Free main structure */
		free(hsdata);
		hsdata = NULL;
	}
}

//-----------------------------------
// Returns Camera rotation angle (counter-clock-wise from back side for prime focus camera)
//-----------------------------------
double get_camera_rotation_angle(int telescope_id)
{
	if ( hsdata != NULL ){
		int itel = get_telescope_index(telescope_id);
		if (itel == TEL_INDEX_NOT_VALID) return TEL_INDEX_NOT_VALID;
		return hsdata->camera_set[itel].cam_rot;
		}
	return -1.;
}

//----------------------------------------------------------------
// Returns total number of mirror tiles.
// Returns TEL_INDEX_NOT_VALID if telescope index is not valid
//----------------------------------------------------------------
int get_mirror_number(int telescope_id)
{
	if ( hsdata != NULL ){
		int itel = get_telescope_index(telescope_id);
		if (itel == TEL_INDEX_NOT_VALID) return TEL_INDEX_NOT_VALID;
		return hsdata->camera_set[itel].num_mirrors;
	}
	return -1.;
}

//----------------------------------------------------------------
// Returns focal length of optics [m].
// Returns TEL_INDEX_NOT_VALID if telescope index is not valid
//----------------------------------------------------------------

double get_optical_foclen(int telescope_id)
{
	if ( hsdata != NULL ){
		int itel = get_telescope_index(telescope_id);
		if (itel == TEL_INDEX_NOT_VALID) return TEL_INDEX_NOT_VALID;
		return hsdata->camera_set[itel].flen;
		}
	return -1.;
}

//-----------------------------------
// Returns IDs of used telescope in the run
//-----------------------------------

int get_telescope_ids(int* list)
{
	if ( hsdata != NULL){
		int num_tel = get_num_telescope();
		int loop=0;
		for (loop=0; loop < num_tel; loop++)
		{
			*list++ =hsdata->run_header.tel_id[loop];
		}
		return 0;
		}
	return -1;
}
