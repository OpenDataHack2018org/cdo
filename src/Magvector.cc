/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2018 Uwe Schulzweida, <uwe.schulzweida AT mpimet.mpg.de>
  See COPYING file for copying and redistribution conditions.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/
#ifdef HAVE_CONFIG_H
#include "config.h" /* HAVE_LIBMAGICS */
#endif

#include <cdi.h>

#include "cdo_int.h"
#include "grid.h"
#include "pstream_int.h"

#ifdef HAVE_LIBMAGICS

#include "magics_api.h"

#include "magics_template_parser.h"
#include "results_template_parser.h"
#include "StringUtilities.h"

#define DBG 0

int VECTOR, STREAM;
const char *vector_params[] = { "thin_fac", "unit_vec", "device", "step_freq" };
int vector_param_count = sizeof(vector_params) / sizeof(char *);

/* Default Magics Values */
double THIN_FAC = 2.0, UNIT_VEC = 25.0;
extern int ANIM_FLAG, STEP_FREQ;

extern int checkdevice(char *device_in);

extern char *DEVICE;
extern char *DEVICE_TABLE;
extern int DEVICE_COUNT;

static void
magvector(const char *plotfile, int operatorID, const char *varname, long nlon, long nlat, double *grid_center_lon,
          double *grid_center_lat, double *uarray, double *varray, int nparam, char **params, char *datetime)

{
  long i;
  double dlon = 0, dlat = 0;
  int split_str_count;
  char plotfilename[4096];
  const char *sep_char = "=";
  char **split_str = NULL;
  char *temp_str = NULL;
  char *titlename;

  (void) varname;

  if (uarray == NULL && varray == NULL)
    {
      fprintf(stderr, " No Velocity Components in input file, cannot creaate "
                      "Vector PLOT!\n");
      return;
    }

  if (uarray == NULL || varray == NULL)
    {
      fprintf(stderr, " Found only one Velocity Component in input file, "
                      "cannot create Vector PLOT!\n");
      return;
    }

  if (DBG)
    {
      fprintf(stderr, "Num params %d\n", nparam);

      for (i = 0; i < nparam; i++) fprintf(stderr, "Param %s\n", params[i]);
      fflush(stderr);
    }

  for (i = 0; i < nparam; ++i)
    {
      split_str_count = 0;
      sep_char = "=";
      split_str_count = StringSplitWithSeperator(params[i], sep_char, &split_str);
      (void) split_str_count;

      if (!strcmp(split_str[0], "thin_fac"))
        {
          THIN_FAC = atof(split_str[1]);
          if (DBG) fprintf(stderr, "THIN FACTOR %g\n", THIN_FAC);
        }

      if (!strcmp(split_str[0], "unit_vec"))
        {
          UNIT_VEC = atof(split_str[1]);
          if (DBG) fprintf(stderr, "UNIT VECTOR %g\n", UNIT_VEC);
        }

      if (!strcmp(split_str[0], "device"))
        {
          temp_str = strdup(split_str[1]);
          StrToUpperCase(temp_str);
          DEVICE = temp_str;
          if (DBG) fprintf(stderr, "DEVICE %s\n", DEVICE);

          mag_setc("output_format", DEVICE);
        }

      if (!strcmp(split_str[0], "step_freq"))
        {
          STEP_FREQ = atoi(split_str[1]);
          if (DBG) fprintf(stderr, "STEP FREQ %d\n", STEP_FREQ);
        }

      Free(split_str);
    }

  if (nlon > 1)
    {
      for (i = 1; i < nlon; ++i) dlon += (grid_center_lon[i] - grid_center_lon[i - 1]);
      dlon /= (nlon - 1);
    }

  if (nlat > 1)
    {
      for (i = 1; i < nlat; ++i) dlat += (grid_center_lat[nlon * i] - grid_center_lat[nlon * (i - 1)]);
      dlat /= (nlat - 1);
    }

  /* magics_template_parser( magics_node ); */

  /* results_template_parser(results_node, varname ); */

  sprintf(plotfilename, "Velocity Vectors %s", datetime);
  titlename = strdup(plotfilename);
  sprintf(plotfilename, "%s", plotfile);

  mag_setc("output_name", plotfilename);
  mag_new("page");

  /* Set the input data */
  mag_setr("input_field_initial_latitude", grid_center_lat[0]);
  mag_setr("input_field_latitude_step", dlat);

  mag_setr("input_field_initial_longitude", grid_center_lon[0]);
  mag_setr("input_field_longitude_step", dlon);

  mag_set2r("input_wind_u_component", uarray, nlon, nlat);
  mag_set2r("input_wind_v_component", varray, nlon, nlat);

  mag_seti("map_label_latitude_frequency", 2);
  mag_seti("map_label_longitude_frequency", 2);
  /*mag_setr ("map_label_height",0.5);*/
  mag_setr("map_label_height", 0.4);

  if (operatorID == VECTOR)
    {
      /* Magics functions for performing vector operation */
      /*
        mag_setc("wind_legend_only", "on" );
        mag_setc("wind_legend_text", "on" );
      */

      mag_setc("legend", "on");
      mag_setc("wind_flag_cross_boundary", "on");
      mag_seti("wind_arrow_thickness", 1);
      mag_coast();

      if (IS_NOT_EQUAL(THIN_FAC, 2.0f)) mag_setr("wind_thinning_factor", THIN_FAC);

      /*wind_arrow_unit_velocity */
      if (IS_NOT_EQUAL(UNIT_VEC, 25.0f)) mag_setr("wind_arrow_unit_velocity", UNIT_VEC);

      mag_wind();

      mag_set1c("text_lines", (const char **) &titlename, 1);
      mag_setc("text_colour", "black");
      mag_setc("text_justification", "centre");
      mag_text();
    }
}

static void
init_MAGICS()

{

  setenv("MAGPLUS_QUIET", "1", 1); /* To suppress magics messages */
  mag_open();

  /* Some standard parameters affectng the magics environment, moved from the
   * xml file  ** begin ** */
  mag_setc("page_id_line", "off");
}

static void
quit_MAGICS()

{

  mag_close();
  if (DBG) fprintf(stdout, "Exiting From MAGICS\n");
}

static void
VerifyVectorParameters(int num_param, char **param_names, int opID)

{

  int i, j;
  int found = FALSE, syntax = TRUE, halt_flag = FALSE, split_str_count;
  int param_count = 0;
  const char **params = NULL;
  char **split_str = NULL;
  const char *sep_char = "=";

  /* char  *vector_params[] =
   * {"min","max","count","interval","list","colour","thickness","style","RGB"};
   */

  for (i = 0; i < num_param; ++i)
    {
      split_str_count = 0;
      found = FALSE;
      syntax = TRUE;
      split_str_count = StringSplitWithSeperator(param_names[i], sep_char, &split_str);

      if (DBG) fprintf(stderr, "Verifying params!\n");

      if (split_str_count > 1)
        {

          if (opID == VECTOR)
            {
              param_count = vector_param_count;
              params = vector_params;
            }

          for (j = 0; j < param_count; ++j)
            {
              if (!strcmp(split_str[0], params[j]))
                {
                  found = TRUE;

                  if (!strcmp(split_str[0], "thin_fac") || !strcmp(split_str[0], "unit_vec") || !strcmp(split_str[0], "step_freq"))
                    {
                      if (!IsNumeric(split_str[1])) syntax = FALSE;
                    }

                  if (!strcmp(split_str[0], "device"))
                    {
                      if (IsNumeric(split_str[1]))
                        syntax = FALSE;
                      else
                        {
                          if (!strcmp(split_str[0], "device"))
                            {
                              if (DBG) fprintf(stderr, "Parameter value '%s'\n", split_str[1]);
                              if (checkdevice(split_str[1])) syntax = FALSE;

                              /* Vector not supported in google earth format */
                              if (!strcmp(split_str[1], "KML") || !strcmp(split_str[1], "kml"))
                                {
                                  syntax = FALSE;
                                  if (DBG) fprintf(stderr, "Parameter value '%s'\n", split_str[1]);
                                }
                            }
                        }
                    }
                }
            }
        }
      else
        {
          syntax = FALSE;
        }

      if (found == FALSE)
        {
          halt_flag = TRUE;
          fprintf(stderr, "Invalid parameter  '%s'\n", param_names[i]);
        }
      if (found == TRUE && syntax == FALSE)
        {
          halt_flag = TRUE;
          fprintf(stderr, "Invalid parameter specification  '%s'\n", param_names[i]);
        }

      if (split_str) Free(split_str);
    }

  if (halt_flag == TRUE)
    {
      exit(0);
    }
}
#endif

void *
Magvector(void *process)
{
  cdoInitialize(process);

#ifdef HAVE_LIBMAGICS
  int nrecs;
  int levelID;
  size_t nmiss;
  char varname[CDI_MAX_NAME];
  char units[CDI_MAX_NAME];
  char vdatestr[32], vtimestr[32], datetimestr[64];

  int nparam = operatorArgc();
  char **pnames = operatorArgv();

  int VECTOR = cdoOperatorAdd("vector", 0, 0, NULL);
  int STREAM = cdoOperatorAdd("stream", 0, 0, NULL);

  int operatorID = cdoOperatorID();

  if (nparam)
    {
      if (DBG)
        {
          for (int i = 0; i < nparam; i++) fprintf(stderr, "Param %d is %s!\n", i + 1, pnames[i]);
        }

      VerifyVectorParameters(nparam, pnames, operatorID);
    }

  int streamID = cdoStreamOpenRead(cdoStreamName(0));

  int vlistID = cdoStreamInqVlist(streamID);
  int taxisID = vlistInqTaxis(vlistID);

  int found = 0;
  int varID = 0;
  int gridID = vlistInqVarGrid(vlistID, varID);
  // int zaxisID = vlistInqVarZaxis(vlistID, varID);

  if (gridInqType(gridID) == GRID_GME) cdoAbort("GME grid unspported!");
  if (gridInqType(gridID) == GRID_UNSTRUCTURED) cdoAbort("Unstructured grid unspported!");

  if (gridInqType(gridID) != GRID_CURVILINEAR) gridID = gridToCurvilinear(gridID, 1);

  size_t gridsize = gridInqSize(gridID);
  int nlon = gridInqXsize(gridID);
  int nlat = gridInqYsize(gridID);
  // int nlev     = zaxisInqSize(zaxisID);

  double *uarray = (double *) Malloc(gridsize * sizeof(double));
  double *varray = (double *) Malloc(gridsize * sizeof(double));
  double *grid_center_lat = (double *) Malloc(gridsize * sizeof(double));
  double *grid_center_lon = (double *) Malloc(gridsize * sizeof(double));

  gridInqYvals(gridID, grid_center_lat);
  gridInqXvals(gridID, grid_center_lon);

  /* Convert lat/lon units if required */
  gridInqXunits(gridID, units);
  grid_to_degree(units, gridsize, grid_center_lon, "grid center lon");
  gridInqYunits(gridID, units);
  grid_to_degree(units, gridsize, grid_center_lat, "grid center lat");

  int tsID = 0;

  /* HARDCODED THE FILE NAME .. TO BE SENT AS COMMAND LINE ARGUMENT FOR THE
   * MAGICS OPERATOR */
  /*
  init_XMLtemplate_parser( Filename );
  updatemagics_and_results_nodes( );
  */

  init_MAGICS();

  while ((nrecs = cdoStreamInqTimestep(streamID, tsID)))
    {
      if (ANIM_FLAG)
        {
          if (tsID % STEP_FREQ)
            {
              tsID++;
              continue;
            }
        }
      else
        {
          if (!STEP_FREQ && tsID)
            {
              cdoWarning("File has values at more than one time step! Image "
                         "created for first time step!!!");
              break;
            }
        }

      int64_t vdate = taxisInqVdate(taxisID);
      int vtime = taxisInqVtime(taxisID);

      date2str(vdate, vdatestr, sizeof(vdatestr));
      time2str(vtime, vtimestr, sizeof(vtimestr));
      sprintf(datetimestr, "%s %s", vdatestr, vtimestr);

      for (int recID = 0; recID < nrecs; recID++)
        {
          pstreamInqRecord(streamID, &varID, &levelID);

          vlistInqVarName(vlistID, varID, varname);

          if (operatorID == VECTOR)
            {
              if (!strcmp(varname, "var131") || !strcmp(varname, "u")) /* U Velocity as per GRIB is
                                                                          'var131, as per NC 'u' */
                {
                  if (DBG) fprintf(stderr, "Found U VEL in Varname %s\n", varname);
                  pstreamReadRecord(streamID, uarray, &nmiss);
                  if (nmiss) cdoSetNAN(vlistInqVarMissval(vlistID, varID), gridsize, uarray);
                  found++;
                }
              if (!strcmp(varname, "var132") || !strcmp(varname, "v")) /* V Velocity as per GRIB  is
                                                                          'var132, as per NC 'v'*/
                {
                  if (DBG) fprintf(stderr, "Found V VEL in Varname %s\n", varname);
                  pstreamReadRecord(streamID, varray, &nmiss);
                  if (nmiss) cdoSetNAN(vlistInqVarMissval(vlistID, varID), gridsize, varray);
                  found++;
                }
              if (found == 2) break;
            }
          else if (operatorID == STREAM)
            fprintf(stderr, " Stream Operator Un-Supported!\n");
          else
            fprintf(stderr, " Operator Un-Supported!\n");
        }

      if (operatorID == VECTOR)
        {
          if (found == 2)
            {
              if (DBG) fprintf(stderr, "Found Both U & V VEL, Creating vector fields! \n");
              magvector(cdoGetStreamName(1).c_str(), operatorID, varname, nlon, nlat, grid_center_lon, grid_center_lat, uarray,
                        varray, nparam, pnames, datetimestr);
            }
          else if (found == 1)
            {
              fprintf(stderr, "Found only one Velocity Component in input "
                              "file, cannot create Vector PLOT!\n");
              break;
            }
          else if (found == 0)
            {
              fprintf(stderr, "No Velocity Components in input file, cannot "
                              "create Vector PLOT!\n");
              break;
            }
        }

      tsID++;

      /*
      if( ANIM_FLAG )
        tsID++;
      else
        {
           cdoWarning("File has values at more than one time step! Image created
      for first time step!!!"); if( STEP_FREQ > 1 ) cdoWarning("Step frequency
      parameter ignored!!!"); break;
        }
      */
    }

  pstreamClose(streamID);

  if (uarray) Free(uarray);
  if (varray) Free(varray);
  if (grid_center_lon) Free(grid_center_lon);
  if (grid_center_lat) Free(grid_center_lat);

  /*   quit_XMLtemplate_parser( ); */

  quit_MAGICS();

#else

  cdoAbort("MAGICS support not compiled in!");

#endif

  cdoFinish();

  return 0;
}
