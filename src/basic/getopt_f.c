/*
	Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch

	This program is free software; you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation; either version 2, or (at your option)
	any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with this program; if not, write to the Free Software
	Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
	02111-1307, USA.

	$Id: getopt_f.c 2516 2006-10-24 21:31:59Z acastro $
*/

#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#if defined(HAVE_GETOPT_LONG)
#include <getopt.h>
#else
#include <unistd.h>
#endif
#include <string.h>
#include "string_f.h"

/* GENERAL FUNCTIONS AND VARIABLES */

char **argv;
int argc;

void FC_FUNC_(set_number_clarg, SET_NUMBER_CLP)(int *nargc)
{
  argc = *nargc+1;
  argv = (char **)malloc(argc*sizeof(char *));
}

void FC_FUNC_(set_clarg, SET_CLARG)(int *i, STR_F_TYPE arg STR_ARG1)
{
  char *c;
  TO_C_STR1(arg, c)
  argv[*i] = c;
}

void FC_FUNC_(clean_clarg, CLEAR_CLARG)()
{
  int i;
  for(i=0; i<argc; i++)
    free(argv[i]);
  free(argv);
}

void print_config(){
#ifdef HAVE_OPENMP
  printf("openmp");
#endif
#ifdef HAVE_MPI
  printf("mpi ");
#endif
#ifdef HAVE_OPENCL
  printf("opencl ");
#endif
#ifdef HAVE_M128D
  printf("sse2 ");
#endif
#ifdef HAVE_M256D
  printf("avx ");
#endif
#ifdef HAVE_BLUE_GENE
  printf("bluegene ");
#endif
#ifdef HAVE_MPI2
  printf("mpi2 ");
#endif
#ifdef HAVE_SCALAPACK
  printf("scalapack ");
#endif
#ifdef HAVE_NETCDF
  printf("netcdf ");
#endif
#ifdef HAVE_METIS
  printf("metis ");
#endif
#ifdef HAVE_GDLIB
  printf("gdlib ");
#endif
#ifdef HAVE_PAPI
  printf("papi ");
#endif
#ifdef HAVE_SPARSKIT
  printf("sparskit ");
#endif
#ifdef HAVE_ETSF_IO
  printf("etsf_io ");
#endif
#ifdef HAVE_PFFT
  printf("pfft ");
#endif
  printf("\n");
}


/* FUNCTIONS TO BE USED BY THE PROGRAM oct-oscillator-strength */
void oscillator_strength_help(){
  printf("Usage: oct-oscillator_strength [OPTIONS] [w]\n");
  printf("\n");
  printf("Options:\n");
  printf("  -h              Print this help and exits.\n");
  printf("  -m <mode>       Select the run mode:\n");
  printf("                    1 (default) analyzes the signal present in an 'ot' file.\n");
  printf("                      This should have been generated by this same utility\n");
  printf("                      (run mode 2).\n");
  printf("                    2 Reads a number of 'multipoles' files, which should be\n");
  printf("                      present in the working directory, and be called\n");
  printf("                      'multipoles.1', 'multipoles.2', ..., and generate an 'ot'\n");
  printf("                      file with the k-th order response of a given operator O.\n");
  printf("                      The order k is decided by the '-k' option. The operator\n");
  printf("                      is decided by the '-O' option.\n");
  printf("                    3 Peforms an analysis of the second-order response of an\n");
  printf("                      operator O, present in the working directory, and\n");
  printf("                      previously generated with run mode 2. It also reads a\n");
  printf("                      file with a list of frequecies around which the search\n");
  printf("                      for resonances is performed.\n");
  printf("                    4 Reads an 'ot' file, and generates an 'omega' file with\n");
  printf("                      either the sine or cosine Fourier transform of the\n");
  printf("                      signal present in 'ot'.\n");
  printf("  -O <operator>   Selects the operator to be analyzed:\n");
  printf("                    o If <operatot> is a pair of integers in the form '(l,m)'\n");
  printf("                      then the operator will be the (l,m) multipole.\n");
  printf("                    o If <operatot> is x, y, or z, then the response operator\n");
  printf("                      to be analyzed will be the dipole in the given direction.\n");
  printf("                    o If the -O option is not given in the command line, then\n");
  printf("                      the observation operator O will be the same as the\n");
  printf("                      perturbation operator that defines the initial kick.\n");
  printf("  -f <file>       This is the file where the frequencies needed in run mode\n");
  printf("                  3 are stored.\n");
  printf("  -d <gamma>      gamma is the damping factor used in the SOS formulae that");
  printf("                  produce (hyper)-polarizabilities.");
  printf("  -s <dw>         Limits of the search interval: [w-dw,w+dw]\n");
  printf("  -r <r>          Number of resonances to search for.\n");
  printf("  -n <N>          Number of frequencies in which the search interval\n");
  printf("                    is discretized (default 1000)\n");
  printf("  -k <k>          Process, or generate, the k-th order response.\n");
  printf("  -t <time>       The signal analysis will be done by integrating in the \n");
  printf("                    time interval [0, <time>]. If this argument is absent,\n");
  printf("                    it makes use of all the time-window present in the\n");
  printf("                    multipoles files.\n");
  exit(-1);
}

void FC_FUNC_(getopt_oscillator_strength, GETOPT_OSCILLATOR_STRENGTH)
  (int *mode, double *omega, double *searchinterval, int *order, 
   int *nresonances, int *nfrequencies, double *time, int *l, int *m, double *damping,
   STR_F_TYPE ffile STR_ARG1)
{
  int c;

  /* This line would be present if we wanted to make the omega a 
     mandatory argument. But for the moment I think it should not be mandatory.
     if(argc==1) oscillator_strength_help(); */

  while (1) {
    c = getopt(argc, argv, "hm:s:k:O:r:n:t:d:f:");
    if (c == -1) break;
    switch (c) {

    case 'h':
      oscillator_strength_help();
      break;

    case 'm':
      *mode = (int)atoi(optarg);
      break;

    case 's':
      *searchinterval = (double)atof(optarg);
      break;

    case 'O':
      c = sscanf(optarg, "(%d,%d)", l, m);
      if(c != 2){
        switch (optarg[0]){
        case 'x':
          *l = 0;
          *m = 1;
          break;
        case 'y':
          *l = 0;
          *m = 2;
          break;
        case 'z':
          *l = 0;
          *m = 3;
          break;
        default:
          printf("Problem reading the -O option value.\n\n");
          oscillator_strength_help();
        }
      }
      break;

    case 'k':
      *order = (int)atoi(optarg);
      break;

    case 'r':
      *nresonances = (int)atoi(optarg);
      break;

    case 'n':
      *nfrequencies = (int)atoi(optarg);
      break;

    case 't':
      *time = (double)atof(optarg);
      break;

    case 'f':
      TO_F_STR1(optarg, ffile);
      break;

    case 'd':
      *damping = (double)atof(optarg);
      break;

    case '?':
      oscillator_strength_help();
      break;

    }
  }
  if (optind < argc) {
    while (optind < argc) *omega = (double)atof(argv[optind++]);
  }

}
/***************************************************************/


/* FUNCTIONS TO BE USED BY THE PROGRAM oct-harmonic-spectrum */
void harmonic_spectrum_help(){
  printf("Usage: oct-harmonic-spectrum [OPTIONS] \n");
  printf("\n");
  printf("Options:\n");
  printf("  -h, --help      Prints this help and exits.\n");
  printf("  -v, --version   Prints octopus version.\n");
  printf("  -w, --freq=freq Specifies the fundamental frequency.\n");
  printf("  -p, --pol=pol   Specifies the direction of the light polarization.\n");
  printf("                  The oct-harmonic-spectrum utility program needs to know\n");
  printf("                  the direction along which the emission radiation is\n");
  printf("                  considered to be polarized. It may be linearly polarized\n");
  printf("                  or circularly polarized. The valid options are:\n");
  printf("                     'x' : Linearly polarized field in the x direction.\n");
  printf("                     'y' : Linearly polarized field in the x direction.\n");
  printf("                     'z' : Linearly polarized field in the x direction.\n");
  printf("                     '+' : Circularly polarized field, counterclockwise.\n");
  printf("                     '-' : Circularly polarized field, clockwise.\n");
  printf("                     'v' : Along a direction specified by -x X -y Y -z Z.\n");
  printf("                  The default is 'x'\n");
  printf("  -a, --ar        Calculates the angle-resolved harmonic-spectrum along a\n");
  printf("                  direction (X,Y,Z) specified by  by -x X -y Y -z Z.\n");
  printf("  -m, --mode=mode Whether the harmonic spectrum is computed by taking the\n");
  printf("                  second derivative of the dipole moment numerically, or by\n");
  printf("                  making use of the acceleration operator, stored in the\n:");
  printf("                  'acceleration' file. The options are:\n");
  printf("                     '1' : use the dipole, take second derivative numerically.\n");
  printf("                     '2' : use the acceleration file.\n");
  printf("                  The default is '1'\n");
  exit(-1);
}

void FC_FUNC_(getopt_harmonic_spectrum, GETOPT_HARMONIC_SPECTRUM)
     (double *w0, int *m, int *ar,  double *x, double *y, double *z, STR_F_TYPE pol STR_ARG1)
{
  int c;

#if defined(HAVE_GETOPT_LONG)
  static struct option long_options[] =
    {
      {"help", no_argument, 0, 'h'},
      {"version", no_argument, 0, 'v'},
      {"freq", required_argument, 0, 'w'},
      {"pol", required_argument, 0, 'p'},
      {"mode", required_argument, 0, 'm'},
      {"ar", required_argument, 0, 'a'},
      {"x", required_argument, 0, 'x'},
      {"y", required_argument, 0, 'y'},
      {"z", required_argument, 0, 'z'},
      {0, 0, 0, 0}
    };
#endif

  while (1) {
    int option_index = 0;
#if defined(HAVE_GETOPT_LONG)
    c = getopt_long(argc, argv, "hvw:p:m:x:y:z:a", long_options, &option_index);
#else
    c = getopt(argc, argv, "hvw:p:m:x:y:z:a");
#endif
    if (c == -1) break;
    switch (c) {

    case 'h':
      harmonic_spectrum_help();
      break;

    case 'v':
      printf("octopus %s (svn version %s)\n", PACKAGE_VERSION, LATEST_SVN);
      exit(0);

    case 'w':
      *w0 = (double)atof(optarg);
      break;

    case 'p':
      TO_F_STR1(optarg, pol);
      break;

    case 'm':
      *m = (int)atoi(optarg);
      break;

    case 'a':
      *ar = 1;
      break;

    case 'x':
      *x = (double)atof(optarg);
      break;

    case 'y':
      *y = (double)atof(optarg);
      break;

    case 'z':
      *z = (double)atof(optarg);
      break;

    }
  }

}
/***************************************************************/


/* FUNCTIONS TO BE USED BY THE PROGRAM oct-help */
void help_help(){
  printf("Usage: oct-help [options] \n");
  printf("\n");
  printf("Options:\n");
  printf("  -h, --help            Prints this help and exits.\n");
  printf("  -v, --version         Prints octopus version.\n");
  printf("  -s, --search=STRING   Search variables whose names contain string 'STRING'.\n");
  printf("  -p, --print=VARNAME   Prints description of variable 'VARNAME'.\n");
  exit(-1);
}


void FC_FUNC_(getopt_help, GETOPT_HELP)
     (STR_F_TYPE mode, STR_F_TYPE name STR_ARG2)
{
  int c;

#if defined(HAVE_GETOPT_LONG)
  static struct option long_options[] =
    {
      {"help", no_argument, 0, 'h'},
      {"version", no_argument, 0, 'v'},
      {"list", no_argument, 0, 'l'},
      {"search", required_argument, 0, 's'},
      {"print", required_argument, 0, 'p'},
      {0, 0, 0, 0}
    };
#endif


  while (1) {
    int option_index = 0;
#if defined(HAVE_GETOPT_LONG)
    c = getopt_long(argc, argv, "hvls:p:", long_options, &option_index);
#else
    c = getopt(argc, argv, "hvls:p:");
#endif
    if(argc==1) help_help();
    if (c == -1) break;
    switch (c) {

    case 'h':
      help_help();
      break;

    case 'v':
      printf("octopus %s (svn version %s)\n", PACKAGE_VERSION, LATEST_SVN);
      exit(0);

    case 'l':
      TO_F_STR1("list", mode);
      return;

    case 's':
      TO_F_STR1("search", mode);
      TO_F_STR2(optarg, name);
      return;

    case 'p':
      TO_F_STR1("print", mode);
      TO_F_STR2(optarg, name);
      return;

    }
  }
  if (optind < argc) help_help();

}


/***************************************************************/


/* FUNCTIONS TO BE USED BY THE PROGRAM octopus */
void octopus_help(){
  printf("Usage: octopus [options] \n");
  printf("\n");
  printf("Options:\n");
  printf("  -h, --help            Prints this help and exits.\n");
  printf("  -v, --version         Prints octopus version.\n");
  printf("  -c, --config          Prints compilation configuration options.\n");
  exit(-1);
}

void FC_FUNC_(getopt_octopus, GETOPT_OCTOPUS)
     ()
{
  int c;
#if defined(HAVE_GETOPT_LONG)
  static struct option long_options[] =
    {
      {"help", no_argument, 0, 'h'},
      {"version", no_argument, 0, 'v'},
      {"config", no_argument, 0, 'c'},
      {0, 0, 0, 0}
    };
#endif

  while (1) {
    int option_index = 0;
#if defined(HAVE_GETOPT_LONG)
    c = getopt_long(argc, argv, "hvc", long_options, &option_index);
#else
    c = getopt(argc, argv, "hvc");
#endif
    if (c == -1) break;
    switch (c) {
    case 'h':
      octopus_help();
      break;
    case 'v':
      printf("octopus %s (svn version %s)\n", PACKAGE_VERSION, LATEST_SVN);
      exit(0);
      break;
    case 'c':
      print_config();
      exit(0);
      break;
    }
  }
  if (optind < argc) octopus_help();

}

/***************************************************************/


/* FUNCTIONS TO BE USED BY THE PROGRAM oct-casida_spectrum */
void casida_spectrum_help(){
  printf("Usage: oct-casida_spectrum [options] \n");
  printf("\n");
  printf("Options:\n");
  printf("  -h, --help            Prints this help and exits.\n");
  printf("  -v, --version         Prints octopus version.\n");
  exit(-1);
}

void FC_FUNC_(getopt_casida_spectrum, GETOPT_CASIDA_SPECTRUM)
     ()
{
  int c;
#if defined(HAVE_GETOPT_LONG)
  static struct option long_options[] =
    {
      {"help", no_argument, 0, 'h'},
      {"version", no_argument, 0, 'v'},
      {0, 0, 0, 0}
    };
#endif

  while (1) {
    int option_index = 0;
#if defined(HAVE_GETOPT_LONG)
    c = getopt_long(argc, argv, "hv", long_options, &option_index);
#else
    c = getopt(argc, argv, "hv");
#endif
    if (c == -1) break;
    switch (c) {
    case 'h':
      casida_spectrum_help();
      break;
    case 'v':
      printf("octopus %s (svn version %s)\n", PACKAGE_VERSION, LATEST_SVN);
      exit(0);
    }
  }
  if (optind < argc) casida_spectrum_help();

}

/***************************************************************/


/* FUNCTIONS TO BE USED BY THE PROGRAM oct-center-geom */
void center_geom_help(){
  printf("Usage: oct-center-geom [options] \n");
  printf("\n");
  printf("Options:\n");
  printf("  -h, --help            Prints this help and exits.\n");
  printf("  -v, --version         Prints octopus version.\n");
  exit(-1);
}

void FC_FUNC_(getopt_center_geom, GETOPT_CENTER_GEOM)
     ()
{
  int c;
#if defined(HAVE_GETOPT_LONG)
  static struct option long_options[] =
    {
      {"help", no_argument, 0, 'h'},
      {"version", no_argument, 0, 'v'},
      {0, 0, 0, 0}
    };
#endif

  while (1) {
    int option_index = 0;
#if defined(HAVE_GETOPT_LONG)
    c = getopt_long(argc, argv, "hv", long_options, &option_index);
#else
    c = getopt(argc, argv, "hv");
#endif
    if (c == -1) break;
    switch (c) {
    case 'h':
      center_geom_help();
      break;
    case 'v':
      printf("octopus %s (svn version %s)\n", PACKAGE_VERSION, LATEST_SVN);
      exit(0);
    }
  }
  if (optind < argc) center_geom_help();

}

/***************************************************************/


/* FUNCTIONS TO BE USED BY THE PROGRAM oct-center-geom */
void dielectric_function_help(){
  printf("Usage: oct-dielectric-function [options] \n");
  printf("\n");
  printf("Options:\n");
  printf("  -h, --help            Prints this help and exits.\n");
  printf("  -v, --version         Prints octopus version.\n");
  exit(-1);
}

void FC_FUNC_(getopt_dielectric_function, GETOPT_DIELECTRIC_FUNCTION)
     ()
{
  int c;
#if defined(HAVE_GETOPT_LONG)
  static struct option long_options[] =
    {
      {"help", no_argument, 0, 'h'},
      {"version", no_argument, 0, 'v'},
      {0, 0, 0, 0}
    };
#endif

  while (1) {
    int option_index = 0;
#if defined(HAVE_GETOPT_LONG)
    c = getopt_long(argc, argv, "hv", long_options, &option_index);
#else
    c = getopt(argc, argv, "hv");
#endif
    if (c == -1) break;
    switch (c) {
    case 'h':
      dielectric_function_help();
      break;
    case 'v':
      printf("octopus %s (svn version %s)\n", PACKAGE_VERSION, LATEST_SVN);
      exit(0);
    }
  }
  if (optind < argc) dielectric_function_help();

}

/***************************************************************/


/* FUNCTIONS TO BE USED BY THE PROGRAM oct-propagation_spectrum */
void propagation_spectrum_help(){
  printf("Usage: oct-propagation_spectrum [options] \n");
  printf("\n");
  printf("Options:\n");
  printf("  -h, --help            Prints this help and exits.\n");
  printf("  -v, --version         Prints octopus version.\n");
  printf("  -r, --reference REF   REF should be the name of the 'reference' multipoles file,.\n");
  printf("                        whenever you want to do a time-dependent response function.\n");
  printf("                        calculation.\n");
  exit(-1);
}

void FC_FUNC_(getopt_propagation_spectrum, GETOPT_PROPAGATION_SPECTRUM)
     (STR_F_TYPE fname STR_ARG1)
{
  int c;
#if defined(HAVE_GETOPT_LONG)
  static struct option long_options[] =
    {
      {"help", no_argument, 0, 'h'},
      {"version", no_argument, 0, 'v'},
      {"reference", required_argument, 0, 'r'},
      {0, 0, 0, 0}
    };
#endif

  while (1) {
    int option_index = 0;
#if defined(HAVE_GETOPT_LONG)
    c = getopt_long(argc, argv, "hvr:", long_options, &option_index);
#else
    c = getopt(argc, argv, "hvr:");
#endif
    if (c == -1) break;
    switch (c) {
    case 'h':
      propagation_spectrum_help();
      break;
    case 'v':
      printf("octopus %s (svn version %s)\n", PACKAGE_VERSION, LATEST_SVN);
      exit(0);
    case 'r':
      TO_F_STR1(optarg, fname);
      break;
    }
  }
  if (optind < argc) propagation_spectrum_help();

}

/***************************************************************/


/* FUNCTIONS TO BE USED BY THE PROGRAM oct-rotatory_strength */
void rotatory_strength_help(){
  printf("Usage: oct-rotatory_strength [options] \n");
  printf("\n");
  printf("Options:\n");
  printf("  -h, --help            Prints this help and exits.\n");
  printf("  -v, --version         Prints octopus version.\n");
  exit(-1);
}

void FC_FUNC_(getopt_rotatory_strength, GETOPT_ROTATORY_STRENGTH)
     ()
{
  int c;
#if defined(HAVE_GETOPT_LONG)
  static struct option long_options[] =
    {
      {"help", no_argument, 0, 'h'},
      {"version", no_argument, 0, 'v'},
      {0, 0, 0, 0}
    };
#endif

  while (1) {
    int option_index = 0;
#if defined(HAVE_GETOPT_LONG)
    c = getopt_long(argc, argv, "hv", long_options, &option_index);
#else
    c = getopt(argc, argv, "hv");
#endif
    if (c == -1) break;
    switch (c) {
    case 'h':
      rotatory_strength_help();
      break;
    case 'v':
      printf("octopus %s (svn version %s)\n", PACKAGE_VERSION, LATEST_SVN);
      exit(0);
    }
  }
  if (optind < argc) rotatory_strength_help();

}

/***************************************************************/


/* FUNCTIONS TO BE USED BY THE PROGRAM oct-vibrational */
void vibrational_help(){
  printf("Usage: oct-vibrational [options] \n");
  printf("\n");
  printf("Options:\n");
  printf("  -h, --help            Prints this help and exits.\n");
  printf("  -v, --version         Prints octopus version.\n");
  printf("  -m, --mode            Run mode: 1 to obtain the vibrational spectrum, .\n");
  printf("                        through the velocity autocorrelation function, 2 to\n");
  printf("                        obtain the infrared spectrum through the dipole.\n");
  exit(-1);
}

void FC_FUNC_(getopt_vibrational, GETOPT_VIBRATIONAL)
     (int *mode)
{
  int c;
#if defined(HAVE_GETOPT_LONG)
  static struct option long_options[] =
    {
      {"help", no_argument, 0, 'h'},
      {"version", no_argument, 0, 'v'},
      {"mode", required_argument, 0, 'm'},
      {0, 0, 0, 0}
    };
#endif

  while (1) {
    int option_index = 0;
#if defined(HAVE_GETOPT_LONG)
    c = getopt_long(argc, argv, "hvm:", long_options, &option_index);
#else
    c = getopt(argc, argv, "hvm:");
#endif
    if (c == -1) break;
    switch (c) {
    case 'h':
      vibrational_help();
      break;
    case 'v':
      printf("octopus %s (svn version %s)\n", PACKAGE_VERSION, LATEST_SVN);
      exit(0);
    case 'm':
      *mode = (int)atoi(optarg);
      break;
    }
  }
  if (optind < argc) vibrational_help();

}

/***************************************************************/


/* FUNCTIONS TO BE USED BY THE PROGRAM oct-xyz-anim */
void xyz_anim_help(){
  printf("Usage: oct-xyz-anim [options] \n");
  printf("\n");
  printf("Options:\n");
  printf("  -h, --help            Prints this help and exits.\n");
  printf("  -v, --version         Prints octopus version.\n");
  exit(-1);
}

void FC_FUNC_(getopt_xyz_anim, GETOPT_XYZ_ANIM)
     ()
{
  int c;
#if defined(HAVE_GETOPT_LONG)
  static struct option long_options[] =
    {
      {"help", no_argument, 0, 'h'},
      {"version", no_argument, 0, 'v'},
      {0, 0, 0, 0}
    };
#endif

  while (1) {
    int option_index = 0;
#if defined(HAVE_GETOPT_LONG)
    c = getopt_long(argc, argv, "hv", long_options, &option_index);
#else
    c = getopt(argc, argv, "hv");
#endif
    if (c == -1) break;
    switch (c) {
    case 'h':
      xyz_anim_help();
      break;
    case 'v':
      printf("octopus %s (svn version %s)\n", PACKAGE_VERSION, LATEST_SVN);
      exit(0);
    }
  }
  if (optind < argc) xyz_anim_help();

}

/***************************************************************/


/* FUNCTIONS TO BE USED BY THE PROGRAM oct-photoelectron-spectrum */
void photoelectron_spectrum_help(){
  printf("Usage: oct-photoelectron-spectrum [OPTIONS] \n");
  printf("\n");
  printf("Options:\n");
  printf("  -h, --help          Prints this help and exits.\n");
  printf("  -v, --version       Prints octopus version.\n");
  printf("  -m, --mode=mode     Whether we want the angle- or energy-resolved photoelectron\n");
  printf("                      spectrum. The options are:\n");
  printf("                         '1' : energy-resolved.\n");
  printf("                         '2' : angle and energy resolved.\n");
  printf("                         '3' : velocity map on a plane.\n");
  printf("                         '4' : angle and energy resolved on the azimuthal plane.\n");
  printf("                      The default is '1'\n");
  printf("  -i, --int=Y/N       Interpolate the output. Default is Yes.\n");
  printf("  -p, --pol=x:y:z     The polarization axis direction. Default: 1:0:0. \n");
  printf("  -e, --de            The resolution in energy.\n");
  printf("  -E,                 Maximum and mium energy colon separated values. \n");
  printf("   --espan=Emin:Emax                                                  \n");
  exit(-1);
}

void FC_FUNC_(getopt_photoelectron_spectrum, GETOPT_PHOTOELECTRON_SPECTRUM)
     ( int *m, int *interp, double *estep, double *espan, double *pol)
{
  int c;
  char delims[] = ":";
  char *tok = NULL;


  *interp = 1;
  *estep  = -1.;
  espan[0]  = -1.;
  espan[1]  = -1.;
  pol[0] = 1;
  pol[1] = 0;
  pol[2] = 0;

#if defined(HAVE_GETOPT_LONG)
  static struct option long_options[] =
    {
      {"help", no_argument, 0, 'h'},
      {"version", no_argument, 0, 'v'},
      {"mode", required_argument, 0, 'm'},
      {"int", required_argument, 0, 'i'},
      {"de", required_argument, 0, 'e'},
      {"espan", required_argument, 0, 'E'},
      {"pol", required_argument, 0, 'p'},
      {0, 0, 0, 0}
    };
#endif

  while (1) {
    int option_index = 0;
#if defined(HAVE_GETOPT_LONG)
    c = getopt_long(argc, argv, "hvm:i:e:E:p:", long_options, &option_index);
#else
    c = getopt(argc, argv, "hvm:i:e:E:p:");
#endif
    if (c == -1) break;
    switch (c) {

    case 'h':
      photoelectron_spectrum_help();
      break;

    case 'v':
      printf("octopus %s (svn version %s)\n", PACKAGE_VERSION, LATEST_SVN);
      exit(0);
    break;


    case 'm':
      *m = (int)atoi(optarg);
      break;

    case 'i':
      if ((strcmp(optarg,"y") == 0)|| (strcmp(optarg,"yes")== 0))  *interp = 1;
      if ((strcmp(optarg,"n") == 0)|| (strcmp(optarg,"no") == 0))  *interp = 0;
    break;
   
    case 'e':
      *estep = atof(optarg);
    break;

    case 'E':
      tok = strtok(optarg,delims);
      espan[0] = atof(tok);
      tok = strtok(NULL,delims);
      espan[1] = atof(tok);     

    break;

    case 'p':
      tok = strtok(optarg,delims);
      pol[0] = atof(tok);
      tok = strtok(NULL,delims);
      pol[1] = atof(tok);     
      tok = strtok(NULL,delims);
      pol[2] = atof(tok);     

    break;

     
    }
  }

}


