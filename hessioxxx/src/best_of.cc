/// @file best_of.cc
/// @short Tool for extracting best values from listings of 'rh3' sensitivity
/// evaluations. Three versions of the 'rh3' output format are supported.
/// All of the input (from standard input) should be in the same format type.

/// One type is before the addition of 68% and 80% angular resolution values.
/// Another one is after addition of angular resolution but before addition
/// of the energy resolution, and the third one is after the energy
/// resolution got added to the output.

/// The different formats are recognized by the presence and position of
/// the histogram number (12056 to 12064 normally) on which the sensitivity
/// evaluation is mainly based.

/// @defgroup best_of_cc The best_of program
/// @{

#include "initial.h"
#include "straux.h"
#include "fileopen.h"
#include <vector>
#include <map>
#include <iostream>
#include <cstdio>
#include <cstring>

// using namespace stdtools;
using std::vector;
using std::map;
using std::string;
using std::cout;
using std::cerr;

enum SpecType { SPEC_NONE = -1, 
                SPEC_GAMMA = 0, SPEC_ELECTRON = 1,
                SPEC_PROTON = 101, SPEC_HE = 402,
                SPEC_CNO = 1407, SPEC_SI = 2814,
                SPEC_IRON = 5626 };

string particle_type(SpecType sp)
{
   switch ( sp )
   {
      case SPEC_GAMMA:
         return "Gamma";
         break;
      case SPEC_ELECTRON:
         return "Electron";
         break;
      case SPEC_PROTON:
         return "Proton";
         break;
      case SPEC_HE:
         return "Helium";
         break;
      case SPEC_CNO:
         return "Nitrogen (CNO)";
         break;
      case SPEC_SI:
         return "Silicon (heavy)";
         break;
      case SPEC_IRON:
         return "Iron (very heavy)";
         break;
      case SPEC_NONE:
         return "Invalid";
         break;
   }

   return "Unknown";
}

// TeV to erg conversion factor:
static double sce = 1.6022;
// m^-2 to cm^-2 conversion factor:
static double sca = 1e-4;
// Total conversion factor
static double sc = sce*sca;

double Crab_Unit(double E)
{
   return 2.79e-7*pow(E,-2.57);  // [1/(m^2 s TeV)]
}
static double cu(double x) {return Crab_Unit(x); }
double Crab_Unit_int(double E)
{
   return 2.79e-7*pow(E,-1.57)/1.57;  // [1/(m^2 s TeV)]
}

double ergs(double E)
{
   return sce*E;
}

// South, 50 h:
static double f50(double x)
{ 
   return 10.00*pow(x/1600.,1.50) +
          0.18*pow(0.01/x,1.3) +
          0.0028*pow(x,0.30) +
          0.08*pow(0.021/x,6.) +
          0.04*pow(0.016/x,14.) +
          0.02*pow(0.016/x,18.);
}
static double fsp50(double x) /* Note: required only from 20 GeV to 300 TeV */
{
   return 1.25*f50(x*0.94);
}
double Flux_req50_south(double E /* TeV */) 
{
   return fsp50(E)*cu(E); // [1/(m^2 s TeV)]
}
double Flux_req50_E2erg_south(double E /* TeV */) 
{
   return fsp50(E)*cu(E)*E*E*sc; // [erg/(cm^2 s)]
}
double Flux_req50_CU_south(double E /* TeV */) 
{
   return fsp50(E); // [Crab Units]
}

// North, 50 h:
static double fn50(double x)
{
   return 150.00*pow(x/1600.,1.74) + 
          0.18*pow(0.01/x,1.3) +
          0.0028*pow(x,0.30) +
          0.08*pow(0.021/x,6.) +
          0.04*pow(0.016/x,14.) + 
          0.02*pow(0.016/x,18.);
}
static double fnsp50(double x) /* Note: required only from 20 GeV to 20 TeV */
{
   return 1.25*fn50(x*0.94);
}
double Flux_req50_north(double E /* TeV */) 
{
   return fnsp50(E)*cu(E); // [1/(m^2 s TeV)]
}
double Flux_req50_E2erg_north(double E /* TeV */) 
{
   return fnsp50(E)*cu(E)*E*E*sc; // [erg/(cm^2 s)]
}
double Flux_req50_CU_north(double E /* TeV */) 
{
   return fnsp50(E); // [Crab Units]
}

// South, 5h:
static double f5(double x)
{
   return 100.0*pow(x/1500.,1.54) +
          0.50*pow(0.01/x,1.3) + 
          0.0120*pow(x,0.25) +
          0.25*pow(0.02/x,6.) +
          0.06*pow(0.016/x,14.) +
          0.03*pow(0.016/x,18.);
}   
static double fsp5(double x)
{
   return 1.25*f5(x*0.94);
}
double Flux_req5_south(double E /* TeV */) 
{
   return fsp5(E)*cu(E); // [1/(m^2 s TeV)]
}
double Flux_req5_E2erg_south(double E /* TeV */) 
{
   return fsp5(E)*cu(E)*E*E*sc; // [erg/(cm^2 s)]
}
double Flux_req5_CU_south(double E /* TeV */) 
{
   return fsp5(E); // [Crab Units]
}

// North, 5h:
static double fn5(double x)
{
   return 1300.0*pow(x/1500.,1.71) +
          0.50*pow(0.01/x,1.3) +
          0.0120*pow(x,0.25) +
          0.25*pow(0.02/x,6) +
          0.06*pow(0.016/x,14.) +
          0.03*pow(0.016/x,18.);
}
static double fnsp5(double x)
{
   return 1.25*fn5(x*0.94);
}
double Flux_req5_north(double E /* TeV */) 
{
   return fnsp5(E)*cu(E); // [1/(m^2 s TeV)]
}
double Flux_req5_E2erg_north(double E /* TeV */) 
{
   return fnsp5(E)*cu(E)*E*E*sc; // [erg/(cm^2 s)]
}
double Flux_req5_CU_north(double E /* TeV */) 
{
   return fnsp5(E); // [Crab Units]
}

// South, 0.5 h:
static double f05(double x)
{
   return 1000.*pow(x/1400.,1.58) +
          1.30*pow(0.01/x,1.3) +
          0.0600*pow(x,0.20) +
          0.50*pow(0.02/x,6.) +
          0.07*pow(0.016/x,14.) +
          0.035*pow(0.016/x,18.);
}
static double fsp05(double x)
{
   return 1.25*f05(x*0.94);
}
double Flux_req05_south(double E /* TeV */) 
{
   return fsp05(E)*cu(E); // [1/(m^2 s TeV)]
}
double Flux_req05_E2erg_south(double E /* TeV */) 
{
   return fsp05(E)*cu(E)*E*E*sc; // [erg/(cm^2 s)]
}
double Flux_req05_CU_south(double E /* TeV */) 
{
   return fsp05(E); // [Crab Units]
}

// North, 0.5 h:
static double fn05(double x)
{
   return 11000.*pow(x/1400.,1.68) +
          1.30*pow(0.01/x,1.3) +
          0.0600*pow(x,0.20) +
          0.50*pow(0.02/x,6.) +
          0.07*pow(0.016/x,14.) +
          0.035*pow(0.016/x,18.);
}
static double fnsp05(double x)
{
   return 1.25*fn05(x*0.94);
}
double Flux_req05_north(double E /* TeV */) 
{
   return fnsp05(E)*cu(E); // [1/(m^2 s TeV)]
}
double Flux_req05_E2erg_north(double E /* TeV */) 
{
   return fnsp05(E)*cu(E)*E*E*sc; // [erg/(cm^2 s)]
}
double Flux_req05_CU_north(double E /* TeV */) 
{
   return fnsp05(E); // [Crab Units]
}

// Goal sensitivities, only derived for 50 hours observation time:

// Goal South, 50 h:
static double fd50(double x)
{
   return 6.50*pow(x/1600.,1.50) +
          0.06*pow(0.01/x,1.3) +
          0.0013*pow(x,0.30) +
          0.05*pow(0.021/x,6.) +
          0.022*pow(0.016/x,14.) +
          0.02*pow(0.016/x,18.);
}
static double fdes50(double x)
{
   return 1.25*fd50(x*0.94);
}
double Flux_goal50_south(double E /* TeV */) 
{
   return fdes50(E)*cu(E); // [1/(m^2 s TeV)]
}
double Flux_goal50_E2erg_south(double E /* TeV */) 
{
   return fdes50(E)*cu(E)*E*E*sc; // [erg/(cm^2 s)]
}
double Flux_goal50_CU_south(double E /* TeV */) 
{
   return fdes50(E); // [Crab Units]
}

// Goal North, 50 h:
static double fnd50(double x)
{
   return 100.0*pow(x/1600.,1.74) +
          0.06*pow(0.01/x,1.3) +
          0.0016*pow(x,0.30) +
          0.05*pow(0.021/x,6.) +
          0.022*pow(0.016/x,14.) +
          0.02*pow(0.016/x,18.);
}
static double fndes50(double x)
{
   return 1.25*fnd50(x*0.94);
}
double Flux_goal50_north(double E /* TeV */) 
{
   return fndes50(E)*cu(E); // [1/(m^2 s TeV)]
}
double Flux_goal50_E2erg_north(double E /* TeV */) 
{
   return fndes50(E)*cu(E)*E*E*sc; // [erg/(cm^2 s)]
}
double Flux_goal50_CU_north(double E /* TeV */) 
{
   return fndes50(E); // [Crab Units]
}

double Angular_resolution_req(double E /* TeV */)
{
   return 1.1*sqrt(
      0.31*0.31 * pow(0.02/E,4.) +
      0.185*0.185 * pow(0.03/E,2.) +
      0.04*0.04 * (0.1/E) +
      0.035*0.035 * sqrt(5./E) +
      0.016*0.016 ); // [degrees]
}
double Angular_resolution_goal(double E /* TeV */)
{
   return sqrt(
      0.25*0.25 * pow(0.02/E,4.) +
      0.15*0.15 * pow(0.03/E,2.) +
      0.030*0.030 * (0.1/E) +
      0.02*0.02 * sqrt(5./E) +
      0.008*0.008 ); // [degrees]
}

static double eresb(double E /* TeV */)
{
   return sqrt(
      0.30*0.30 * pow(0.03/E,2.) +
      0.20*0.20 * (0.1/E) +
      0.05*0.05 * sqrt(5./E) ); // Fraction of E, r.m.s.,
}
double Energy_resolution_req(double E /* TeV */)
{
   double de = eresb(E);
   if ( de < 0.1 )
      return 0.1;
   else
      return de;
}

static double eresdb(double E)
{
   return 0.764*sqrt(
      0.12*0.12 * pow(0.03/E,2.) +
      0.16*0.16 * (0.1/E) +
      0.035*0.035 * sqrt(5./E) ); // Fraction of E, r.m.s.,
}
double Energy_resolution_goal(double E /* TeV */)
{
   double de = eresdb(E);
   if ( de < 0.05 )
      return 0.05;
   else
      return de;
}

enum espec_t { OLD_E_POWERLAW = 1, NEW_E_POWERLAW = 2, 
      NEW_E_PL_LGN1 = 3, NEW_E_PL_LGN2 = 4 };
espec_t espec_type = OLD_E_POWERLAW;

double flux_int (SpecType sp, double E1, double E2)
{   
   double sp_idx = -2.5, sp_norm = 1.0;
   double e_peak = 0., amp_peak = 0., width_peak=0.;
   double fi;
   switch ( sp )
   {
      case SPEC_GAMMA: // HEGRA(1) Crab power-law
         sp_idx = -2.57;
         sp_norm = 2.79e-7; // [TeV^-1 s^-1 m^-2]
         break;
      case SPEC_ELECTRON:
         if ( espec_type == OLD_E_POWERLAW )
         {
         sp_idx = -3.30;
         sp_norm = 4.75e-5;// [TeV^-1 s^-1 m^-2 sr^-1] Multiply with solid angle 
         }
         else if ( espec_type == NEW_E_POWERLAW )
         {
         // Not including the Fermi shoulder in the spectrum.
         // (Not to speak of the claimed ATIC peak).
         sp_idx = -3.21;
         sp_norm = 6.85e-5;// [TeV^-1 s^-1 m^-2 sr^-1] Multiply with solid angle 
         }
         else if ( espec_type == NEW_E_PL_LGN1 )
         {
         // Including the Fermi shoulder in the spectrum (power law + log-normal)
         sp_idx = -3.21;
         sp_norm = 6.85e-5; // [TeV^-1 s^-1 m^-2 sr^-1] Multiply with solid angle
         // F_peak(E) = Amp/(E*width*sqrt(2pi))*exp(-0.5*((log(E)-log(Epeak))/width)**2)
         e_peak = 0.107;     // TeV, peak position in log-normal distribution
         amp_peak = 3.186e-3;// [TeV^-1 s^-1 m^-2 sr^-1]
         width_peak = 0.776; // width in log E
         }
         else
         {
         // Including the Fermi shoulder in the spectrum (power law + log-normal)
         sp_idx = -3.201;
         sp_norm = 6.98e-5; // [TeV^-1 s^-1 m^-2 sr^-1] Multiply with solid angle
         // F_peak(E) = Amp/(E*width*sqrt(2pi))*exp(-0.5*((log(E)-log(Epeak))/width)**2)
         e_peak = 0.0180;     // TeV, peak position in log-normal distribution
         amp_peak = 2.927e-2;// [TeV^-1 s^-1 m^-2 sr^-1]
         width_peak = 1.171; // width in log E
         }
         break;
      case SPEC_PROTON: // BESS power-law without cut-off
         sp_idx = -2.70;
         sp_norm = 0.096; // [TeV^-1 s^-1 m^-2 sr^-1] Multiply with solid angle later.
         break;
      case SPEC_HE:
         sp_idx = -2.64;
         sp_norm = 0.0719; // [TeV^-1 s^-1 m^-2 sr^-1] Multiply with solid angle 
         break;
      case SPEC_CNO:
         sp_idx = -2.67;
         sp_norm = 0.0321; // [TeV^-1 s^-1 m^-2 sr^-1] Multiply with solid angle 
         break;
      case SPEC_SI: // Not just silicon but heavy element group around Si
         sp_idx = -2.66;
         sp_norm = 0.0284;
         break;
      case SPEC_IRON: // Group of elements around iron
         sp_idx = -2.63;
         sp_norm = 0.0134;
         break;
      default:
         sp_idx = -2.;
         sp_norm = 0.;
         break;
   }

   // dF/dE = N*E^(-x) -> dF/dlog10E = ln(10.)*N*E^(-x+1)
   // Int_E1^E2(dF/dE dE) = -N/(-x+1)*(E1^(-x+1)-E2^(-x+1))
   
   fi = -sp_norm / (sp_idx+1) * (pow(E1,sp_idx+1) - pow(E2,sp_idx+1));
   if ( e_peak > 0. && width_peak > 0. )
   {
      double int_lgn = 0.5 * (
          erf((log(E2)-log(e_peak))/(width_peak*sqrt(2.))) -
          erf((log(E1)-log(e_peak))/(width_peak*sqrt(2.))) );
      fi += amp_peak * int_lgn;
   }
   
   return fi;
}

double lima17(double on, double off, double alpha)
{
   double q = 
        on  * log( (1.+alpha)/alpha * (on /(on+off)) ) +
        off * log( (1.+alpha)       * (off/(on+off)) );
   if ( q <= 0. )
      q = 0.;
   double s_17 = sqrt(2.0) * sqrt(q);
   return s_17;
}

bool matching_required_diffsens(int calc_pput, bool with_flux, double E, double diff_sens)
{
   double diff_sens_req = 0.;
   switch ( calc_pput )
   {
      case -1:
      case 1:
         diff_sens_req =  (with_flux ? Flux_req50_south(E) : Flux_req50_CU_south(E));
         break;
      case -2:
      case 2:
         diff_sens_req =  (with_flux ? Flux_req5_south(E) : Flux_req5_CU_south(E));
         break;
      case -3:
      case 3:
         diff_sens_req =  (with_flux ? Flux_req05_south(E) : Flux_req05_CU_south(E));
         break;
      case -11:
      case 11:
         diff_sens_req =  (with_flux ? Flux_req50_north(E) : Flux_req50_CU_north(E));
         break;
      case -12:
      case 12:
         diff_sens_req =  (with_flux ? Flux_req5_north(E) : Flux_req5_CU_north(E));
         break;
      case -13:
      case 13:
         diff_sens_req =  (with_flux ? Flux_req05_north(E) : Flux_req05_CU_north(E));
         break;
   }
   if ( diff_sens > diff_sens_req )
      return false;
   return true;
}

bool matching_required_performance(int calc_pput, bool with_flux, double E, double diff_sens, double angres, double eres)
{
   if ( ! matching_required_diffsens(calc_pput,with_flux,E,diff_sens) )
      return false;
   if ( angres > Angular_resolution_req(E) )
      return false;
   if ( eres > Energy_resolution_req(E) )
      return false;
   return true;
}

bool matching_required_angres(double E, double angres)
{
   if ( angres > Angular_resolution_req(E) )
      return false;
   return true;
}

bool matching_required_eres(double E, double eres)
{
   if ( eres > Energy_resolution_req(E) )
      return false;
   return true;
}

struct best_value
{
   int kbin;
   double best;
   int q;
   string text;
   double A; ///< effective area (for gammas)
   double lgE, lgE1, lgE2;
   double diff_sens;
   double bg_rate;
   double gamma_rate;
   double angres;
   double eres;
   double ebias;
   double n_gamma_cu, nint_gamma_cu;
   double n_bg, nint_bg; // sum of electrons, protons, helium, CNO, heavy, iron nuclei, in source region
   best_value() : kbin(0), best(9999.), text(""), A(0.), 
      lgE(-12.), lgE1(-12.), lgE2(12.), diff_sens(0.), bg_rate(0.), gamma_rate(0.), angres(0.), eres(0.), ebias(0.),
      n_gamma_cu(0.), nint_gamma_cu(0.), n_bg(0.), nint_bg(0.) {}
   best_value(int k, double v, int qtr, const string& t, double aeff, 
      double vlgE, double vlgE1, double vlgE2, double vds, double vbr=0., double vgr=0., double var=0., double ver=0., double veb=0.,
      double ng=0., double nb=0.) : /* 16 parameters, 9 required */
         kbin(k), best(v), q(qtr), text(t), A(aeff), 
         lgE(vlgE), lgE1(vlgE1), lgE2(vlgE2), diff_sens(vds), bg_rate(vbr), gamma_rate(vgr), angres(var), eres(ver), ebias(veb),
         n_gamma_cu(ng), nint_gamma_cu(0.), n_bg(nb), nint_bg(0.) {}
   ~best_value() {}
};

enum BestChoice { BestDiff=1, BestIntegral=2, BestAngle=3, BestEres=4, BestRate=5, BestCombined=6, BestAll=7 };

int main(int argc, char **argv)
{
   map<int,best_value> all_best;
   // bool integral = false, best_angle = false, best_eres = false, best_rate = false;
   int max_qtr = 9, max_qtr2 = 8;
   bool showed_header = false, best_angle_ok = false;
   vector<string> fnames;
   vector<double> dv(22,0.);
   double min_lg_e = -3., max_lg_e = 5.;
   double req_diff_sens = 0.;
   int h_low = 12056, h_high = 12064;
   int t_low = 0, t_high = 0;
   BestChoice select = BestDiff;
   double hours = 0.;
   bool show_eff_area = 0;
   int calc_pput = 0;
   double s_pput = 0.;
   int n_pput = 0;
   bool with_flux = false;
   double comb_wt = 0.8; /* 80% by sensitivity, 20% by rate */
   bool with_zero_protons = true;
   bool with_extrapolation = true;
   bool match_required = false;
   bool match_required_angres = false;
   bool match_required_eres = false;
   bool match_required_diffsens = false;
   double req_upper_lgE = 2.0, req_lower_lgE=-1.5;
   bool no_fail = false, no_penalty = false;
   bool only_even = false, simple_overlap = false;
   double max_ebias = 1000.;

   std::cout << "# Cmdline: ";
   for ( int iarg=0; iarg<argc; ++iarg )
      std::cout << " " << argv[iarg];
   std::cout << "\n";

   int iarg = 1;
   for ( iarg=1; iarg<argc; ++iarg )
   {
      if ( strcmp(argv[iarg],"-i") == 0 )
      {
         // integral = true;
         select = BestIntegral;
      }
      else if ( strcmp(argv[iarg],"-a") == 0 )
      {
         // best_angle = true;
         select = BestAngle;
      }
      else if ( strcmp(argv[iarg],"-e") == 0 )
      {
         // best_eres = true;
         select = BestEres;
      }
      else if ( strcmp(argv[iarg],"-f") == 0 )
      {
         with_flux = true;
      }
      else if ( strcmp(argv[iarg],"-r") == 0 )
      {
         // best_rate = true;
         select = BestRate;
      }
      else if ( strcmp(argv[iarg],"-d") == 0 )
      {
         select = BestDiff;
      }
      else if ( strcmp(argv[iarg],"-c") == 0 )
      {
         select = BestCombined;
      }
      else if ( strcmp(argv[iarg],"-C") == 0 )
      {
         select = BestAll;
      }
      else if ( strcmp(argv[iarg],"--combine") == 0 && iarg+1<argc )
      {
         select = BestCombined;
         comb_wt = atof(argv[++iarg]);
         if ( comb_wt < 0 )
         {
            select = BestAll;
            comb_wt *= -1.;
         }
      }
      else if ( strcmp(argv[iarg],"--even") == 0 )
      {
         only_even = true;
      }
      else if ( strcmp(argv[iarg],"--simple-overlap") == 0 )
      {
         simple_overlap = true;
      }
      else if ( strcmp(argv[iarg],"--no-fail") == 0 )
      {
         no_fail = true;
      }
      else if ( strcmp(argv[iarg],"--no-penalty") == 0 )
      {
         no_penalty = true;
      }
      else if ( strcmp(argv[iarg],"--no-zero-protons") == 0 )
         with_zero_protons = false;
      else if ( strcmp(argv[iarg],"--no-extrapolation") == 0 )
         with_extrapolation = false;
      else if ( strcmp(argv[iarg],"-A") == 0 )
         show_eff_area = true;
      else if ( strcmp(argv[iarg],"-h") == 0 && iarg+1<argc )
         hours = atof(argv[++iarg]);
      else if ( strcmp(argv[iarg],"-t") == 0 && iarg+1<argc )
         max_qtr = atoi(argv[++iarg]);
      else if ( strcmp(argv[iarg],"-T") == 0 && iarg+1<argc )
         max_qtr2 = atoi(argv[++iarg]);
      else if ( strcmp(argv[iarg],"--min-lgE") == 0 && iarg+1<argc )
         min_lg_e = atof(argv[++iarg]);
      else if ( strcmp(argv[iarg],"--max-lgE") == 0 && iarg+1<argc )
         max_lg_e = atof(argv[++iarg]);
      else if ( strcmp(argv[iarg],"--max-Ebias") == 0 && iarg+1<argc )
         max_ebias = atof(argv[++iarg]);
      else if ( strcmp(argv[iarg],"--min-tel") == 0 && iarg+1<argc )
      {
         iarg++;
         int i1=t_low, i2=t_high;
         int rc = sscanf(argv[iarg],"%d,%d",&i1,&i2);
         if ( rc >= 1 && i1 > 0 )
            t_low = t_high = i1;
         if ( rc >= 2 && i2 > 0 )
            t_high = i2;
      }
      else if ( strcmp(argv[iarg],"--max-sens") == 0 && iarg+1<argc )
         req_diff_sens = atof(argv[++iarg]);
      else if ( strcmp(argv[iarg],"--match-required") == 0 )
         match_required = true;
      else if ( strcmp(argv[iarg],"--match-required-angres") == 0 )
         match_required_angres = true;
      else if ( strcmp(argv[iarg],"--match-required-eres") == 0 )
         match_required_eres = true;
      else if ( strcmp(argv[iarg],"--match-required-diffsens") == 0 )
         match_required_diffsens = true;
      else if ( strcmp(argv[iarg],"--pput") == 0 && iarg+1<argc )
      {
         if ( strcasecmp(argv[iarg+1],"S50") == 0 )
         {   calc_pput = 1; req_upper_lgE = 2.0; }
         else if ( strcasecmp(argv[iarg+1],"S5") == 0 )
         {   calc_pput = 2; req_upper_lgE = 2.0; }
         else if ( strcasecmp(argv[iarg+1],"S0.5") == 0 )
         {   calc_pput = 3; req_upper_lgE = 2.0; }
         else if ( strcasecmp(argv[iarg+1],"N50") == 0 )
         {   calc_pput = 11; req_upper_lgE = 1.3; }
         else if ( strcasecmp(argv[iarg+1],"N5") == 0 )
         {   calc_pput = 12; req_upper_lgE = 1.3; }
         else if ( strcasecmp(argv[iarg+1],"N0.5") == 0 )
         {   calc_pput = 13; req_upper_lgE = 1.3; }
         else
         {
            fprintf(stderr,"Invalid PPUT mode. Use one of: S50 S5 S0.5 N50 N5 N0.5\n");
            exit(1);
         }
         iarg++;
      }
      else if ( strcmp(argv[iarg],"--hist-range") == 0 && iarg+1<argc )
      {
         iarg++;
         int i1=h_low, i2=h_high;
         int rc = sscanf(argv[iarg],"%d,%d",&i1,&i2);
         if ( rc >= 1 && i1 > 0 )
            h_low = i1;
         if ( rc >= 2 && i2 > 0 )
            h_high = i2;
      }
      else if ( argv[iarg][0] == '-' )
      {
         std::cerr << "Best_of: post-optimization of read_simtel (read_cta) analysis results\n"
                   << "as derived with the rh3 utility.\n\n";
         std::cerr << "Syntax: best_of [ -i | -a | -e ] [ -t n ] ... [ files ]\n";
         std::cerr << "    -a : optimise by angular resolution.\n";
         std::cerr << "    -e : optimise by energy resolution.\n";
         std::cerr << "    -i : optimise by integral sensitivity.\n";
         std::cerr << "    -r : optimise by gamma-ray rate.\n";
         std::cerr << "    -d : optimise by diff. sensitivity (default).\n";
         std::cerr << "    -c : combine by diff. sens. with gamma-ray rate (4:1).\n";
         std::cerr << "    -C : comb. by diff. sens./rate/Eres/angres (4:1:2:2).\n";
         std::cerr << "    --combine <wt_sens> : similar but user-defined weighting\n";
         std::cerr << "    Default: optimisation by diff. sensitivity.\n";
         std::cerr << "    -A : add effective area for gamma rays to output.\n";
         std::cerr << "    -h <t> : hours observation time (default: 0=automatic).\n";
         std::cerr << "    -t <n> : Maximum threshold index value accepted, e.g.:\n";
         std::cerr << "      n=0 (well above peak in dN/d logE)\n";
         std::cerr << "      n=2 (just above peak in dN/d logE)\n";
         std::cerr << "      n=3 (just at peak in dN/d logE)\n";
         std::cerr << "      n=4 (just below peak in dN/d logE)\n";
         std::cerr << "      n=6 (well below peak in dN/d logE but not too steep)\n";
         std::cerr << "      n=9 (well below peak in dN/d logE and very steep rise)\n";
         std::cerr << "    -T <n> : Similar to '-t' but less strict.\n";
         std::cerr << "    -f : sensitivity in input is in flux units.\n";
         std::cerr << "    --min-lgE ... (minimum lg(E) to be included)\n";
         std::cerr << "    --max-lgE ... (maximum lg(E) to be included)\n";
         std::cerr << "    --min-tel ... (only entries with matching min. tel.)\n";
         std::cerr << "    --max-sens ... (maximum sensitivity of accepted rows)\n";
         std::cerr << "    --max-Ebias ... (maximum energy bias fraction)\n";
         std::cerr << "    --match-required : (full set of requirements is applied,\n";
         std::cerr << "                   with requirement set specified by --pput option)\n";
         std::cerr << "    --match-required-diffsens : (diff.sensitivity req. satisfied)\n";
         std::cerr << "    --match-required-angres : (angular resolution req. satisfied)\n";
         std::cerr << "    --match-required-eres : (energy resolution req. satisfied)\n";
         std::cerr << "    --hist-range h1,h2 : (use only histogram IDs in limited range)\n";
         std::cerr << "    --even : (only use intervals starting at even fraction of decade)\n";
         std::cerr << "    --simple-overlap : (integral rates with overlapping bins use even share).\n";
         std::cerr << "    --no-zero-protons ... (ignore cut combinations with zero proton bg.)\n";
         std::cerr << "    --no-extrapolation ... (ignore lines with extrapolated cut efficiencies.)\n";
         std::cerr << "    --no-fail : Do not mark cases which failed to cover the PPUT energy range.)\n";
         std::cerr << "    --no-penalty : No penalty for energy bins not covered.)\n";
         std::cerr << "    --pput { S50 | S5 | S0.5 | N50 | N5 | N0.5 } ...\n"
                   << "      (divide sensitivity by required one, and then the inverse geometric mean)\n";
         exit(1);
      }
      else
         fnames.push_back(argv[iarg]);
   }
   
   if ( match_required && calc_pput == 0 )
   {
      std::cerr << "Assuming requirements for CTA South, 50 h. Use '--pput' option to change.\n";
      calc_pput = -1;
   }

   // if ( integral+best_angle+best_eres+best_rate > 1 )
   // {
   //    std::cerr << "Please use only one of '-i', '-e', and '-a'.\n";
   //    exit(1);
   // }
   if ( fnames.size() == 0 )
      fnames.push_back("-");

   for ( size_t kn=0; kn<fnames.size(); ++kn )
   {
      FILE *input = 0;
      if ( fnames[kn] == "-" )
         input = stdin;
      else
         input = fileopen(fnames[kn].c_str(),"r");
      if ( input == 0 )
      {
         perror(fnames[kn].c_str());
         exit(1);
      }
      if ( strstr(fnames[kn].c_str(),"flux") != NULL )
         with_flux = true;

      char line[1024], line2[1024];
      size_t nline = 0;
      while ( fgets(line,sizeof(line)-1,input) != 0 )
      {
         nline++;
         if ( line[0] == '\n' || line[0] == '\0' )
            continue;
         if ( line[0] == '#' )
         {
            if ( strncmp(line,"# Cmdline:",10) == 0 )
               std::cout << line;
            if ( strncmp(line,
                 "##+ Differential and integral limits in units of 1/(m^2 s [TeV])",
                 64) == 0 )
               with_flux = 1;
            if ( !showed_header && 
               (strncmp(line,"##+lgE1",7) == 0 ||
                strncmp(line,"#lgE1 ",6) == 0) )
            {
               std::cout << line;
               showed_header = true;
            }
            if ( hours == 0. )
            {
               char *c = strstr(line," after ");
               if ( c != 0 )
                  if ( strstr(c," hours ") != 0 )
                     hours = atof(c+7);
            }
      	    continue;
         }

if ( line[0] < 32 || line[0] >= 127 )
{
   std::cerr << "Invalid character '" << line[0] 
      << "' (" << size_t((unsigned char)(line[0])) << 
      ") at the beginning of line " << nline << 
      " of file " << fnames[kn].c_str() << "\n";
   exit(1);
}

         if ( ! with_extrapolation )
         {
            if ( strstr(line,"\te\n") != NULL )
               continue;
         }
         
         strcpy(line2,line);
         char *trm;
         if ( (trm = strstr(line2," e")) != NULL )
            *trm = '\0';
         if ( (trm = strstr(line2,"\te")) != NULL )
            *trm = '\0';
         if ( (trm = strchr(line2,'(')) != NULL )
            *trm = '\0';
         if ( (trm = strchr(line2,'[')) != NULL )
            *trm = '\0';
         if ( (trm = strchr(line2,'{')) != NULL )
            *trm = '\0';

         bool bad = false;
         char word[512];
         int ipos=0;
         size_t nv=13, kv=0;
         int this_hist = 0;
         int this_min_tel = 0;
         for ( size_t i=0; i<dv.size(); i++ )
            dv[i] = 0.;
         for ( size_t i=0; i<nv+3 && i<dv.size(); i++ )
         {
            dv[i] = 0.;
      	    if ( getword(line2,&ipos,word,sizeof(word)-1,'\t','\n') <= 0 )
	    {
               if ( i<nv )
	          bad = true;
	       break;
	    }
	    dv[i] = atof(word);
            if ( i==10 && dv[i] > 12000. && dv[i] < 12100. ) // New format with 6 more values
            {
               nv = 19;
               kv = 6;
               best_angle_ok = true;
               this_hist = dv[i];
               this_min_tel = int(dv[3]+0.1);
            }
            else if ( i==9 && dv[i] > 12000. && dv[i] < 12100. ) // New format with 5 more values
            {
               nv = 18;
               kv = 5;
               best_angle_ok = true;
               this_hist = dv[i];
               this_min_tel = int(dv[3]+0.1);
            }
            else if ( i==8 && dv[i] > 12000. && dv[i] < 12100. ) // New format with 4 more values
            {
               nv = 17;
               kv = 4;
               best_angle_ok = true;
               this_hist = dv[i];
               this_min_tel = int(dv[2]+0.1);
            }
            else if ( i==7 && dv[i] > 12000. && dv[i] < 12100. ) // New format with 3 more values
            {
               nv = 16;
               kv = 3;
               if ( select == BestEres )
                  select = BestDiff; // Energy resolution not available, cannot select on it.
               // best_eres = false; // Not available.
               best_angle_ok = true;
               this_hist = dv[i];
               this_min_tel = int(dv[2]+0.1);
            }
         }
         if ( bad )
      	    continue;
         if ( this_hist < h_low || this_hist > h_high )
            continue;
         if ( t_low > 0 && this_min_tel < t_low )
            continue;
         if ( t_high > 0 && this_min_tel > t_high )
            continue;
         if ( (dv[0] == 0. && dv[1] == 0.) || 
              dv[0] < -3. || dv[0] > 5. ||
              dv[1]<=dv[0] || dv[1]-dv[0] > 0.5 )
         {
            std::cerr << "Bad contents in line " << nline << " of file "
                << fnames[kn] << "\n";
            exit(1);
         }
         double gamma_events = dv[nv-1];
         double proton_events = dv[nv];
         double electron_events = dv[nv+1];
         double nuclei_events = dv[nv+2];
         double bg_events = proton_events+electron_events+nuclei_events;
         double gamma_rate = 0., gamma_eff_area = 0.;
         double bg_rate = 0.;
         double angres = ( (kv>=3) ? dv[5+kv] : 0. );
         double eres = ( (kv>=6) ? dv[2+kv] : ( (kv>=4) ? dv[3+kv] : 0. ) );
         double ebias = ( (kv>=6) ? dv[3+kv] : 0. );
         
         if ( fabs(ebias) > max_ebias )
            continue;

         if ( hours > 0. )
         {
            double elow = pow(10.,dv[0]);
            double ehigh = pow(10.,dv[1]);
            gamma_rate = gamma_events / hours / 3600.;
            gamma_eff_area = gamma_rate / flux_int(SPEC_GAMMA, elow, ehigh);
            bg_rate = bg_events / hours / 3600.;
// std::cout << "E = " << elow << " to " << ehigh << " TeV: " << gamma_events << " gammas in " 
//   << hours << " hours -> rate of " << gamma_rate << "/s, eff.area=" << gamma_eff_area <<" (nv=" << nv << ")\n";
         }
         if ( !best_angle_ok )
         {
            // best_angle = best_eres = false; // Not available.
            if ( select == BestEres || select == BestAngle )
               select = BestDiff;
         }
         if ( dv[0] < min_lg_e || dv[1] > max_lg_e )
            continue;
         int k = int(0.5*(dv[0]+dv[1])*20+60+0.1);
         int qtr = int(dv[7+kv]+0.5);
         if ( qtr > max_qtr )
            continue;
         // int best_col = ( best_rate ) ? 10+kv :
         //                ( best_eres&&kv>=4 ) ? 3+kv :
         //                ( ( best_angle&&kv>=3 ) ? 5+kv : 
         //                  ( integral ? 11+kv : 9+kv ) );
         //    continue;
         int best_col = ( select == BestRate ) ? 12+kv :
                        ( (select==BestEres)&&(kv>=4) ) ? (kv>=6 ? 2+kv : 3+kv) :
                        ( ( (select==BestAngle)&&(kv>=3) ) ? 5+kv : 
                          ( (select==BestIntegral) ? 11+kv : 9+kv ) );
         int best_col2 = ( select == BestCombined || select == BestAll ) ? 12+kv : 0;
         int diff_sens_col = 9 + kv;
         double lgE = 0.5*(dv[0]+dv[1]);
         if ( kv >= 4 && dv[2] >= dv[0] && dv[2] <= dv[1] )
            lgE = dv[2];
         double diff_sens = dv[diff_sens_col];
         if ( match_required )
            if ( !matching_required_performance(calc_pput,with_flux,pow(10.,lgE),diff_sens,angres,eres) )
               continue;
         if ( match_required_diffsens )
            if ( !matching_required_diffsens(calc_pput,with_flux,pow(10.,lgE),diff_sens) )
               continue;
         if ( match_required_angres )
            if ( !matching_required_angres(pow(10.,lgE),angres) )
               continue;
         if ( match_required_eres )
            if ( !matching_required_eres(pow(10.,lgE),eres) )
               continue;
         if ( req_diff_sens != 0. && diff_sens > req_diff_sens )
            continue;
         if ( kv>=4 && (kv>=6 ? dv[2+kv] : dv[3+kv]) == 0. ) // Don't trust lines with zero energy resolution
            continue;
         if ( kv >= 3 && dv[5+kv] == 0. )  // Don't trust lines with zero angular resolution either
            continue;
         if ( ! with_zero_protons )
         {
            if ( 13+kv < dv.size() )
               if ( dv[13+kv] <= 0. )
                  continue;
         }
         if ( size_t(best_col) >= dv.size() )
         {
            std::cerr << "Bad best column in line " << nline << ".\n";
            exit(1);
         }
         double current = dv[best_col];
         if ( current == 0. )
            continue;
         // Note optimizing by best rate (and sensitivity combined with rate)
         // is generally not a good idea for arrays of multiple telescope types
         // because combining the good quality data from the larger type with
         // the poor quality data of the smaller type does not work out well.
         if ( select == BestRate )
         {
            if ( current > 0. )
               current = 1./current; // Rate wanted as large as possible.
         }
         if ( best_col2 > 0 && size_t(best_col2) < dv.size() ) // BestCombined
         {
            if ( dv[best_col] > 0. && dv[best_col2] > 0. )
               current = comb_wt * log(dv[best_col]) - (1.-comb_wt) * log(dv[best_col2]);
            else
               current = -1e99;
            if ( select == BestAll )
            {
               double w1 = comb_wt;
               double w2 = (1.-comb_wt);
               double w = w1 + 5.*w2;
               w1 /= w;
               w2 /= w;
               current = w1 * log(dv[best_col]) - w2 * log(dv[best_col2]) + w2 * log(angres*angres) + w2 * log(eres*eres);
            }
         }
         map<int,best_value>::iterator p = all_best.find(k);
         char *c = strchr(line,'\n');
         if ( c != 0 )
            *c = '\0';
         char *s = NULL;
         if ( show_eff_area && hours > 0. && (s=strstr(line,"(A=")) != NULL )
            *s = '\0';
         double lgE1 = dv[0];
         double lgE2 = dv[1];
         if ( only_even )
         {
            if ( int(10.*(lgE1+6.)+0.001)%2 != 0 )
               continue;
         }
         if ( p == all_best.end() )
      	    all_best[k] = best_value(k,current,qtr,line,gamma_eff_area,
                  lgE,lgE1,lgE2,diff_sens,bg_rate,gamma_rate,angres,eres,ebias,
                  gamma_events, bg_events);
         else if ( (*p).second.best > current &&
              (qtr <= max_qtr2 || qtr <= (*p).second.q) )
         {
      	    (*p).second.best = current;
	    (*p).second.text = line;
            (*p).second.q    = qtr;
            (*p).second.A    = gamma_eff_area;
            (*p).second.lgE  = lgE;
            (*p).second.lgE1  = lgE1;
            (*p).second.lgE2  = lgE2;
            (*p).second.diff_sens = diff_sens;
            (*p).second.bg_rate = bg_rate;
            (*p).second.gamma_rate = gamma_rate;
            (*p).second.angres = angres;
            (*p).second.eres = eres;
            (*p).second.ebias = ebias;
            (*p).second.n_gamma_cu = gamma_events;
            (*p).second.n_bg = bg_events;
         }
         else if ( (*p).second.q > max_qtr2 && qtr <= max_qtr2 )
         {
      	    (*p).second.best = current;
	    (*p).second.text = line;
            (*p).second.q    = qtr;
            (*p).second.A    = gamma_eff_area;
            (*p).second.lgE  = lgE;
            (*p).second.lgE1  = lgE1;
            (*p).second.lgE2  = lgE2;
            (*p).second.diff_sens = diff_sens;
            (*p).second.bg_rate = bg_rate;
            (*p).second.gamma_rate = gamma_rate;
            (*p).second.angres = angres;
            (*p).second.eres = eres;
            (*p).second.ebias = ebias;
            (*p).second.n_gamma_cu = gamma_events;
            (*p).second.n_bg = bg_events;
         }
      }

      if ( input != stdin )
         fileclose(input);
      // input = 0;
   }
   
   double lowest_lgE = 3., highest_lgE = -3.;
   size_t n_lgE = 0;
   double overlapping_correction_g = 1.;
   double overlapping_correction_p = 1.;
   map<int,best_value>::iterator p1 = all_best.begin();
   map<int,best_value>::reverse_iterator p2 = all_best.rbegin();
   size_t nbest = all_best.size();
   if ( p1 != all_best.end() && p2 != all_best.rend() && nbest > 1 )
   {
      // If two bins overlap, the correction factor for integral event counts
      // is not simply 0.5 because due to the spectrum more than half of the
      // events are in the lower half, for constant effective area.
      // Due to changes in cuts and effective area from bin to bin, it is not
      // actually that simple. Generally the assumption of constant effective
      // area gives a too large correction factor at low energies.
      double stepping = ((*p2).second.lgE2-(*p1).second.lgE1) / (nbest*((*p1).second.lgE2-(*p1).second.lgE1));
      
      if ( simple_overlap )
      {
         // An even share of all overlapping bins is assumed. This may still be too large
         // below threshold but is too small over most of the energy range.
         overlapping_correction_g = overlapping_correction_p = stepping;
      }
      else
      {
         // Constant effective area expectations:
         overlapping_correction_g = (1.-pow(pow(10.,stepping*0.2),-1.57))/(1.-pow(pow(10.,0.2),-1.57));
         overlapping_correction_p = (1.-pow(pow(10.,stepping*0.2),-1.7))/(1.-pow(pow(10.,0.2),-1.7));
         // Closer to even share between bins is a better match over the whole energy range:
         overlapping_correction_g = 0.5*(stepping + overlapping_correction_g);
         overlapping_correction_p = 0.5*(stepping + overlapping_correction_p);
      }
      if ( int(overlapping_correction_g*1000) != 1000 )
      {
         std::cout << "# Overlapping bins: correction factors for integral rate: " 
                   << overlapping_correction_g << ", " << overlapping_correction_p << "\n";
         std::cout << "# for " << nbest << " intervals of " << (*p1).second.lgE2-(*p1).second.lgE1
                   << " decades each covering a range of " << (*p2).second.lgE2-(*p1).second.lgE1 << " decades\n";
      }
   }

   // Iterate from highest energies to lowest energies, integrating event counts
   double nint_gamma_cu = 0., nint_bg = 0.;
   for ( map<int,best_value>::reverse_iterator p = all_best.rbegin(); p != all_best.rend(); ++p )
   {
      bool first = (p == all_best.rbegin());
      double E1 = pow(10.,(*p).second.lgE1);
      double E2 = pow(10.,(*p).second.lgE2);
      // For the highest energy bin we extrapolate to higher energies
      nint_gamma_cu += (*p).second.n_gamma_cu * (first ? 
         flux_int(SPEC_GAMMA,E1,10.*E1)/flux_int(SPEC_GAMMA,E1,E2) : overlapping_correction_g);
      nint_bg += (*p).second.n_bg * (first ? 
         flux_int(SPEC_PROTON,E1,10.*E1)/flux_int(SPEC_PROTON,E1,E2) : overlapping_correction_p);
      (*p).second.nint_gamma_cu = nint_gamma_cu;
      (*p).second.nint_bg = nint_bg;
   }
   
   double plgE2 = -999.;

   // Now iterate from lowest energies upward
   for ( map<int,best_value>::iterator p = all_best.begin(); p != all_best.end(); ++p )
   {
      double ngamma = (*p).second.nint_gamma_cu;
      double nbg = (*p).second.nint_bg;
      double alpha = 0.2;
      double f_5s = 999.;
      if ( ngamma > 0. && nbg > 0. )
      {
         // First guess: sigma = sqrt(ngamma/nbg) = 5
         f_5s = 25. * nbg / (ngamma+0.05*nbg) * sqrt(1.2); // Fraction of CU resulting in 5 sigma
         for (size_t k=0; k<10; k++ )
         {
            double on = f_5s*ngamma + nbg;
            double off = nbg/alpha;
            double sig = lima17(on,off,alpha);
            if ( sig <= 0. )
               f_5s *= 2.;
            else
               f_5s *= pow(5./sig,1.12);
         }
         if ( f_5s*ngamma < 10. )
            f_5s = 10./ngamma;
         if ( f_5s*ngamma < 0.05*nbg )
            f_5s = 0.05*nbg/ngamma;
      }

      // Since we re-calculate the integral sensitivity anyway, drop any such part
      string::size_type idx = (*p).second.text.find("\t[IS: ");
      if ( idx != string::npos )
         (*p).second.text.resize(idx);

      // Write an empty line if there are gaps in the energy coverage (useful with gnuplot etc.)
      if ( plgE2 > -10. && (*p).second.lgE1 > plgE2+0.05 )
         std::cout << "\n";
      plgE2 = (*p).second.lgE2;

      std::cout << (*p).second.text;

      if ( show_eff_area && hours > 0. )
         std::cout << "\t(A=" << (*p).second.A 
                   << " m^2, Bg=" << (*p).second.bg_rate*3600.
                   << "/h, Gr=" << (*p).second.gamma_rate*3600.
         << "/h/CU)";
      // std::cout << "\t{lgE=" << (*p).second.lgE << ", dS=" <<(*p).second.diff_sens << "}";


      if ( (*p).second.lgE1 < lowest_lgE )
         lowest_lgE = (*p).second.lgE1;
      if ( (*p).second.lgE2 > highest_lgE )
         highest_lgE = (*p).second.lgE2;
      n_lgE++;

      if ( calc_pput > 0 && (*p).second.diff_sens > 0. )
      {
         double E = pow(10.,(*p).second.lgE);
         double diff_sens = (*p).second.diff_sens;
         double diff_sens_req = 0.;
         switch ( calc_pput )
         {
            case 1:
               diff_sens_req =  (with_flux ? Flux_req50_south(E) : Flux_req50_CU_south(E));
               break;
            case 2:
               diff_sens_req =  (with_flux ? Flux_req5_south(E) : Flux_req5_CU_south(E));
               break;
            case 3:
               diff_sens_req =  (with_flux ? Flux_req05_south(E) : Flux_req05_CU_south(E));
               break;
            case 11:
               diff_sens_req =  (with_flux ? Flux_req50_north(E) : Flux_req50_CU_north(E));
               break;
            case 12:
               diff_sens_req =  (with_flux ? Flux_req5_north(E) : Flux_req5_CU_north(E));
               break;
            case 13:
               diff_sens_req =  (with_flux ? Flux_req05_north(E) : Flux_req05_CU_north(E));
               break;
         }
         if ( diff_sens > 0. && diff_sens_req > 0. )
         {
            std::cout << "\t[P: " << diff_sens/diff_sens_req << " ]";
            s_pput += log(diff_sens/diff_sens_req);
            n_pput++;
         }
      }
      // To skip over irrelevant parts in output use: sed 's/\s*\(\S*\s*\S*\s*\S*\).*IS:/\1     /'
      std::cout << "\t[IS: " << f_5s << " " 
         << f_5s*Crab_Unit_int(pow(10.,(*p).second.lgE1)) 
         << " " << ngamma << " " << nbg << " ]";
      std::cout << "\n";
   }

   if ( n_pput > 0 )
   {
      // Penalize and mark if the desired energy range is not fully covered:
      bool incomplete = false;
      double stepsize = (n_lgE > 1) ? (highest_lgE - lowest_lgE) / (n_lgE-1) : 5.0;
      if ( ! no_penalty )
      {
         if ( lowest_lgE-0.05 > req_lower_lgE ) // Each missing point assumed two orders of magnitude below requirement
         {
            s_pput += 2.*log(10.)*(lowest_lgE - req_lower_lgE) / stepsize;
            n_pput += (lowest_lgE - req_lower_lgE) / stepsize;
//            std::cerr << "\nlowest_lgE = " << lowest_lgE << ", req_lower_lgE = " << req_lower_lgE << "\n";
            incomplete = true;
         }
         if ( highest_lgE+0.05 < req_upper_lgE ) // Each missing point assumed one order of magnitude below requirement
         {
            s_pput += 1.*log(10.)*(req_upper_lgE-highest_lgE) / stepsize;
            n_pput += (req_upper_lgE-highest_lgE) / stepsize;
//            std::cerr << "\nhighest_lgE = " << highest_lgE << ", req_upper_lgE = " << req_upper_lgE << "\n";
            incomplete = true;
         }
      }

      if ( incomplete && !no_fail )
         std::cout << "\nPPUT = " << exp(-1.*s_pput/n_pput) << (incomplete ? "x" : "") 
         << " [" << lowest_lgE << ":" << highest_lgE << "]"
         << "\n";
      else
         std::cout << "\nPPUT = " << exp(-1.*s_pput/n_pput) << "\n";
   }

   return 0;
}

/// @}
