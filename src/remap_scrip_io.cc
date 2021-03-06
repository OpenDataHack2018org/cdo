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
#include "config.h"
#endif

#ifdef HAVE_LIBNETCDF
#include "netcdf.h"
#endif

#include <time.h>

#include <cdi.h>

#include "cdo_int.h"
#include "grid.h"
#include "remap.h"
#include "commandline.h"

#ifdef HAVE_LIBNETCDF
static void
nce(int istat)
{
  // This routine provides a simple interface to NetCDF error message routine.
  if (istat != NC_NOERR) cdoAbort(nc_strerror(istat));
}

static void
writeLinks(int nc_file_id, int nc_add_id, nc_type sizetype, size_t num_links, size_t *cell_add)
{
  if (num_links == 0) return;

  if (sizetype == NC_INT)
    {
      std::vector<int> intadd(num_links);
      for (size_t i = 0; i < num_links; ++i) intadd[i] = (int) cell_add[i];
      nce(nc_put_var_int(nc_file_id, nc_add_id, &intadd[0]));
    }
#ifdef HAVE_NETCDF4
  else
    {
      nce(nc_put_var_ulonglong(nc_file_id, nc_add_id, (unsigned long long *) cell_add));
    }
#endif
}

static void
readLinks(int nc_file_id, int nc_add_id, size_t num_links, size_t *cell_add)
{
  if (num_links < 0x7FFFFC00)  // 2GB
    {
      std::vector<int> intadd(num_links);
      nce(nc_get_var_int(nc_file_id, nc_add_id, &intadd[0]));
      for (size_t i = 0; i < num_links; ++i) cell_add[i] = (size_t) intadd[i];
    }
#ifdef HAVE_NETCDF4
  else
    {
      nce(nc_get_var_ulonglong(nc_file_id, nc_add_id, (unsigned long long *) cell_add));
    }
#endif
}
#endif

void
remapWriteDataScrip(const char *interp_file, RemapMethod mapType, SubmapType submapType, int numNeighbors, int remapOrder,
                    RemapGrid &src_grid, RemapGrid &tgt_grid, RemapVars &rv)
{
// Writes remap data to a NetCDF file using SCRIP conventions

#ifdef HAVE_LIBNETCDF

  // Local variables

  int nc_file_id;           /* id for NetCDF file                       */
  int nc_srcgrdsize_id;     /* id for source grid size                  */
  int nc_dstgrdsize_id;     /* id for destination grid size             */
  int nc_srcgrdcorn_id = 0; /* id for number of source grid corners     */
  int nc_dstgrdcorn_id = 0; /* id for number of dest grid corners       */
  int nc_srcgrdrank_id;     /* id for source grid rank                  */
  int nc_dstgrdrank_id;     /* id for dest grid rank                    */
  int nc_numlinks_id;       /* id for number of links in mapping        */
  int nc_numwgts_id;        /* id for number of weights for mapping     */
  int nc_srcgrddims_id;     /* id for source grid dimensions            */
  int nc_dstgrddims_id;     /* id for dest grid dimensions              */
  int nc_srcgrdcntrlat_id;  /* id for source grid center latitude       */
  int nc_dstgrdcntrlat_id;  /* id for dest grid center latitude         */
  int nc_srcgrdcntrlon_id;  /* id for source grid center longitude      */
  int nc_dstgrdcntrlon_id;  /* id for dest grid center longitude        */
  int nc_srcgrdimask_id;    /* id for source grid mask                  */
  int nc_dstgrdimask_id;    /* id for dest grid mask                    */
  int nc_srcgrdcrnrlat_id;  /* id for latitude of source grid corners   */
  int nc_srcgrdcrnrlon_id;  /* id for longitude of source grid corners  */
  int nc_dstgrdcrnrlat_id;  /* id for latitude of dest grid corners     */
  int nc_dstgrdcrnrlon_id;  /* id for longitude of dest grid corners    */
  int nc_srcgrdarea_id;     /* id for area of source grid cells         */
  int nc_dstgrdarea_id;     /* id for area of dest grid cells           */
  int nc_srcgrdfrac_id;     /* id for area fraction on source grid      */
  int nc_dstgrdfrac_id;     /* id for area fraction on dest grid        */
  int nc_srcadd_id;         /* id for map source address                */
  int nc_dstadd_id;         /* id for map destination address           */
  int nc_rmpmatrix_id;      /* id for remapping matrix                  */

  int nc_dims2_id[2]; /* NetCDF ids for 2d array dims             */

  const char *map_name = "SCRIP remapping with CDO";
  char normalize_opt[64] = "unknown";
  char map_method[64] = "unknown";
  char tmp_string[64] = "unknown";
  char src_grid_name[64] = "source grid";
  char tgt_grid_name[64] = "dest grid";
  const char *src_grid_units = "radians";
  const char *tgt_grid_units = "radians";
  bool lgridarea = false;
  int writemode = NC_CLOBBER;
  nc_type sizetype = NC_INT;

  switch (rv.normOpt)
    {
    case NormOpt::NONE: strcpy(normalize_opt, "none"); break;
    case NormOpt::FRACAREA: strcpy(normalize_opt, "fracarea"); break;
    case NormOpt::DESTAREA: strcpy(normalize_opt, "destarea"); break;
    }

  switch (mapType)
    {
    case RemapMethod::CONSERV:
      lgridarea = true;
      if (submapType == SubmapType::LAF)
        {
          strcpy(map_method, "Largest area fraction");
          break;
        }
      else
        {
          strcpy(map_method, "Conservative remapping");
          break;
        }
    case RemapMethod::CONSERV_YAC:
      lgridarea = true;
      /*
      if ( submapType == SubmapType::LAF )
        {
          strcpy(map_method, "Largest area fraction");
          break;
        }
      else
      */
      {
        strcpy(map_method, "Conservative remapping using clipping on sphere");
        break;
      }
    case RemapMethod::BILINEAR: strcpy(map_method, "Bilinear remapping"); break;
    case RemapMethod::BICUBIC: strcpy(map_method, "Bicubic remapping"); break;
    case RemapMethod::DISTWGT:
      if (numNeighbors == 1)
        strcpy(map_method, "Nearest neighbor");
      else
        strcpy(map_method, "Distance weighted avg of nearest neighbors");
      break;
    case RemapMethod::UNDEF: break;
    }

  /*
  if ( rv.num_links == 0 )
    cdoAbort("Number of remap links is 0, no remap weights found!");
  */
  {
    size_t nlinks = rv.num_links;
    size_t nele1 = 4 * 8 + 4;
    size_t nele2 = 4 * 8 + 4;
    if (src_grid.lneed_cell_corners) nele1 += src_grid.num_cell_corners * 2 * 8;
    if (tgt_grid.lneed_cell_corners) nele2 += tgt_grid.num_cell_corners * 2 * 8;
    size_t filesize = src_grid.size * nele1 + tgt_grid.size * nele2 + nlinks * (4 + 4 + rv.num_wts * 8);

    if (cdoVerbose)
      {
        cdoPrint("Number of remap links:       %zu", nlinks);
        cdoPrint("Filesize for remap weights: ~%zu", filesize);
      }

    if (filesize > 0x7FFFFC00)  // 2**31 - 1024 (<2GB)
      {
        size_t maxlinks = 0x3FFFFFFF;  // 1GB
        size_t gridsize_max = (src_grid.size > tgt_grid.size) ? src_grid.size : tgt_grid.size;
        if (nlinks > maxlinks || filesize > 8 * maxlinks || gridsize_max > 0x7FFFFC00)
          {
#ifdef HAVE_NETCDF4
            if (cdoVerbose) cdoPrint("Store weights and links to NetCDF4!");
            writemode |= NC_NETCDF4;
            if (gridsize_max > 0x7FFFFC00)
              sizetype = NC_UINT64;
            else
              writemode |= NC_CLASSIC_MODEL;
#else
            cdoPrint("Number of remap links %zu exceeds maximum of %zu and NetCDF 4 is not available!", nlinks, maxlinks);
#endif
          }
        else
          {
#if defined(NC_64BIT_OFFSET)
            writemode |= NC_64BIT_OFFSET;
            if (cdoVerbose) cdoPrint("Store weights and links to NetCDF2!");
#else
            cdoPrint("Filesize for remap weights maybe too large!");
#endif
          }
      }
  }

  // Create NetCDF file for mapping and define some global attributes
  nce(nc_create(interp_file, writemode, &nc_file_id));

  // Map name
  nce(nc_put_att_text(nc_file_id, NC_GLOBAL, "title", strlen(map_name), map_name));

  // Normalization option
  nce(nc_put_att_text(nc_file_id, NC_GLOBAL, "normalization", strlen(normalize_opt), normalize_opt));

  // Map method
  nce(nc_put_att_text(nc_file_id, NC_GLOBAL, "map_method", strlen(map_method), map_method));

  // Remap order
  if (mapType == RemapMethod::CONSERV && submapType == SubmapType::NONE)
    nce(nc_put_att_int(nc_file_id, NC_GLOBAL, "remap_order", NC_INT, 1L, &remapOrder));

  // File convention
  strcpy(tmp_string, "SCRIP");
  nce(nc_put_att_text(nc_file_id, NC_GLOBAL, "conventions", strlen(tmp_string), tmp_string));

  // Source and destination grid names
  gridName(gridInqType(src_grid.gridID), src_grid_name);
  nce(nc_put_att_text(nc_file_id, NC_GLOBAL, "source_grid", strlen(src_grid_name), src_grid_name));

  gridName(gridInqType(tgt_grid.gridID), tgt_grid_name);
  nce(nc_put_att_text(nc_file_id, NC_GLOBAL, "dest_grid", strlen(tgt_grid_name), tgt_grid_name));

  // History
  time_t date_and_time_in_sec = time(NULL);
  if (date_and_time_in_sec != -1)
    {
      char history[1024] = "date and time";
      struct tm *date_and_time = localtime(&date_and_time_in_sec);
      (void) strftime(history, 1024, "%d %b %Y : ", date_and_time);
      strcat(history, commandLine());
      nce(nc_put_att_text(nc_file_id, NC_GLOBAL, "history", strlen(history), history));
    }

  if (CDO_Version_Info) nce(nc_put_att_text(nc_file_id, NC_GLOBAL, "CDO", (int) strlen(cdoComment()) + 1, cdoComment()));

  // Prepare NetCDF dimension info

  // Define grid size dimensions
  nce(nc_def_dim(nc_file_id, "src_grid_size", src_grid.size, &nc_srcgrdsize_id));
  nce(nc_def_dim(nc_file_id, "dst_grid_size", tgt_grid.size, &nc_dstgrdsize_id));

  // Define grid corner dimension
  if (src_grid.lneed_cell_corners) nce(nc_def_dim(nc_file_id, "src_grid_corners", src_grid.num_cell_corners, &nc_srcgrdcorn_id));
  if (tgt_grid.lneed_cell_corners) nce(nc_def_dim(nc_file_id, "dst_grid_corners", tgt_grid.num_cell_corners, &nc_dstgrdcorn_id));

  // Define grid rank dimension
  nce(nc_def_dim(nc_file_id, "src_grid_rank", src_grid.rank, &nc_srcgrdrank_id));
  nce(nc_def_dim(nc_file_id, "dst_grid_rank", tgt_grid.rank, &nc_dstgrdrank_id));

  // Define map size dimensions
  nce(nc_def_dim(nc_file_id, "num_links", rv.num_links, &nc_numlinks_id));
  nce(nc_def_dim(nc_file_id, "num_wgts", rv.num_wts, &nc_numwgts_id));

  // Define grid dimensions
  nce(nc_def_var(nc_file_id, "src_grid_dims", sizetype, 1, &nc_srcgrdrank_id, &nc_srcgrddims_id));
  nce(nc_def_var(nc_file_id, "dst_grid_dims", sizetype, 1, &nc_dstgrdrank_id, &nc_dstgrddims_id));

  // Define all arrays for NetCDF descriptors

  // Define grid center latitude array
  nce(nc_def_var(nc_file_id, "src_grid_center_lat", NC_DOUBLE, 1, &nc_srcgrdsize_id, &nc_srcgrdcntrlat_id));
  nce(nc_def_var(nc_file_id, "dst_grid_center_lat", NC_DOUBLE, 1, &nc_dstgrdsize_id, &nc_dstgrdcntrlat_id));

  // Define grid center longitude array
  nce(nc_def_var(nc_file_id, "src_grid_center_lon", NC_DOUBLE, 1, &nc_srcgrdsize_id, &nc_srcgrdcntrlon_id));
  nce(nc_def_var(nc_file_id, "dst_grid_center_lon", NC_DOUBLE, 1, &nc_dstgrdsize_id, &nc_dstgrdcntrlon_id));

  // Define grid corner lat/lon arrays

  nc_dims2_id[0] = nc_srcgrdsize_id;
  nc_dims2_id[1] = nc_srcgrdcorn_id;

  if (src_grid.lneed_cell_corners)
    {
      nce(nc_def_var(nc_file_id, "src_grid_corner_lat", NC_DOUBLE, 2, nc_dims2_id, &nc_srcgrdcrnrlat_id));
      nce(nc_def_var(nc_file_id, "src_grid_corner_lon", NC_DOUBLE, 2, nc_dims2_id, &nc_srcgrdcrnrlon_id));
    }

  nc_dims2_id[0] = nc_dstgrdsize_id;
  nc_dims2_id[1] = nc_dstgrdcorn_id;

  if (tgt_grid.lneed_cell_corners)
    {
      nce(nc_def_var(nc_file_id, "dst_grid_corner_lat", NC_DOUBLE, 2, nc_dims2_id, &nc_dstgrdcrnrlat_id));
      nce(nc_def_var(nc_file_id, "dst_grid_corner_lon", NC_DOUBLE, 2, nc_dims2_id, &nc_dstgrdcrnrlon_id));
    }

  // Define units for all coordinate arrays

  nce(nc_put_att_text(nc_file_id, nc_srcgrdcntrlat_id, "units", strlen(src_grid_units), src_grid_units));
  nce(nc_put_att_text(nc_file_id, nc_dstgrdcntrlat_id, "units", strlen(tgt_grid_units), tgt_grid_units));
  nce(nc_put_att_text(nc_file_id, nc_srcgrdcntrlon_id, "units", strlen(src_grid_units), src_grid_units));
  nce(nc_put_att_text(nc_file_id, nc_dstgrdcntrlon_id, "units", strlen(tgt_grid_units), tgt_grid_units));
  if (src_grid.lneed_cell_corners)
    {
      nce(nc_put_att_text(nc_file_id, nc_srcgrdcrnrlat_id, "units", strlen(src_grid_units), src_grid_units));
      nce(nc_put_att_text(nc_file_id, nc_srcgrdcrnrlon_id, "units", strlen(src_grid_units), src_grid_units));
    }
  if (tgt_grid.lneed_cell_corners)
    {
      nce(nc_put_att_text(nc_file_id, nc_dstgrdcrnrlat_id, "units", strlen(tgt_grid_units), tgt_grid_units));
      nce(nc_put_att_text(nc_file_id, nc_dstgrdcrnrlon_id, "units", strlen(tgt_grid_units), tgt_grid_units));
    }

  // Define grid mask

  nce(nc_def_var(nc_file_id, "src_grid_imask", NC_INT, 1, &nc_srcgrdsize_id, &nc_srcgrdimask_id));
  nce(nc_put_att_text(nc_file_id, nc_srcgrdimask_id, "units", 8, "unitless"));

  nce(nc_def_var(nc_file_id, "dst_grid_imask", NC_INT, 1, &nc_dstgrdsize_id, &nc_dstgrdimask_id));
  nce(nc_put_att_text(nc_file_id, nc_dstgrdimask_id, "units", 8, "unitless"));

  // Define grid area arrays

  if (lgridarea)
    {
      nce(nc_def_var(nc_file_id, "src_grid_area", NC_DOUBLE, 1, &nc_srcgrdsize_id, &nc_srcgrdarea_id));
      nce(nc_put_att_text(nc_file_id, nc_srcgrdarea_id, "units", 14, "square radians"));

      nce(nc_def_var(nc_file_id, "dst_grid_area", NC_DOUBLE, 1, &nc_dstgrdsize_id, &nc_dstgrdarea_id));
      nce(nc_put_att_text(nc_file_id, nc_dstgrdarea_id, "units", 14, "square radians"));
    }

  // Define grid fraction arrays

  nce(nc_def_var(nc_file_id, "src_grid_frac", NC_DOUBLE, 1, &nc_srcgrdsize_id, &nc_srcgrdfrac_id));
  nce(nc_put_att_text(nc_file_id, nc_srcgrdfrac_id, "units", 8, "unitless"));

  nce(nc_def_var(nc_file_id, "dst_grid_frac", NC_DOUBLE, 1, &nc_dstgrdsize_id, &nc_dstgrdfrac_id));
  nce(nc_put_att_text(nc_file_id, nc_dstgrdfrac_id, "units", 8, "unitless"));

  // Define mapping arrays

  nce(nc_def_var(nc_file_id, "src_address", sizetype, 1, &nc_numlinks_id, &nc_srcadd_id));
  nce(nc_def_var(nc_file_id, "dst_address", sizetype, 1, &nc_numlinks_id, &nc_dstadd_id));

  nc_dims2_id[0] = nc_numlinks_id;
  nc_dims2_id[1] = nc_numwgts_id;

  nce(nc_def_var(nc_file_id, "remap_matrix", NC_DOUBLE, 2, nc_dims2_id, &nc_rmpmatrix_id));

  // End definition stage

  nce(nc_enddef(nc_file_id));

  // Write mapping data

  int dims[2];
  dims[0] = (int) src_grid.dims[0];
  dims[1] = (int) src_grid.dims[1];
  nce(nc_put_var_int(nc_file_id, nc_srcgrddims_id, dims));
  dims[0] = (int) tgt_grid.dims[0];
  dims[1] = (int) tgt_grid.dims[1];
  nce(nc_put_var_int(nc_file_id, nc_dstgrddims_id, dims));

  nce(nc_put_var_int(nc_file_id, nc_srcgrdimask_id, &src_grid.mask[0]));
  nce(nc_put_var_int(nc_file_id, nc_dstgrdimask_id, &tgt_grid.mask[0]));

  if (src_grid.cell_center_lat) nce(nc_put_var_double(nc_file_id, nc_srcgrdcntrlat_id, src_grid.cell_center_lat));
  if (src_grid.cell_center_lon) nce(nc_put_var_double(nc_file_id, nc_srcgrdcntrlon_id, src_grid.cell_center_lon));

  if (src_grid.lneed_cell_corners)
    {
      nce(nc_put_var_double(nc_file_id, nc_srcgrdcrnrlat_id, src_grid.cell_corner_lat));
      nce(nc_put_var_double(nc_file_id, nc_srcgrdcrnrlon_id, src_grid.cell_corner_lon));
    }

  if (tgt_grid.cell_center_lat) nce(nc_put_var_double(nc_file_id, nc_dstgrdcntrlat_id, tgt_grid.cell_center_lat));
  if (tgt_grid.cell_center_lon) nce(nc_put_var_double(nc_file_id, nc_dstgrdcntrlon_id, tgt_grid.cell_center_lon));

  if (tgt_grid.lneed_cell_corners)
    {
      nce(nc_put_var_double(nc_file_id, nc_dstgrdcrnrlat_id, tgt_grid.cell_corner_lat));
      nce(nc_put_var_double(nc_file_id, nc_dstgrdcrnrlon_id, tgt_grid.cell_corner_lon));
    }

  if (lgridarea) nce(nc_put_var_double(nc_file_id, nc_srcgrdarea_id, &src_grid.cell_area[0]));

  nce(nc_put_var_double(nc_file_id, nc_srcgrdfrac_id, &src_grid.cell_frac[0]));

  /*
  if ( luse_cell_area )
    nce(nc_put_var_double(nc_file_id, nc_dstgrdarea_id, tgt_grid.cell_area_in));
  else
  */
  if (lgridarea) nce(nc_put_var_double(nc_file_id, nc_dstgrdarea_id, &tgt_grid.cell_area[0]));

  nce(nc_put_var_double(nc_file_id, nc_dstgrdfrac_id, &tgt_grid.cell_frac[0]));

  for (size_t i = 0; i < rv.num_links; i++)
    {
      rv.src_cell_add[i]++;
      rv.tgt_cell_add[i]++;
    }

  writeLinks(nc_file_id, nc_srcadd_id, sizetype, rv.num_links, &rv.src_cell_add[0]);
  writeLinks(nc_file_id, nc_dstadd_id, sizetype, rv.num_links, &rv.tgt_cell_add[0]);
  nce(nc_put_var_double(nc_file_id, nc_rmpmatrix_id, &rv.wts[0]));

  nce(nc_close(nc_file_id));

#else
  cdoAbort("NetCDF support not compiled in!");
#endif

}  // remapWriteDataScrip

/*****************************************************************************/

#ifdef HAVE_LIBNETCDF
static RemapMethod
getMapType(int nc_file_id, SubmapType *submapType, int *numNeighbors, int *remapOrder)
{
  // Map method
  size_t attlen;
  char map_method[64];
  nce(nc_get_att_text(nc_file_id, NC_GLOBAL, "map_method", map_method));
  nce(nc_inq_attlen(nc_file_id, NC_GLOBAL, "map_method", &attlen));
  map_method[attlen] = 0;

  *submapType = SubmapType::NONE;
  *remapOrder = 1;

  RemapMethod mapType = RemapMethod::UNDEF;
  if (cmpstr(map_method, "Conservative") == 0)
    {
      if (cmpstr(map_method, "Conservative remapping using clipping on sphere") == 0)
        mapType = RemapMethod::CONSERV_YAC;
      else
        mapType = RemapMethod::CONSERV;

      int iatt;
      int status = nc_get_att_int(nc_file_id, NC_GLOBAL, "remap_order", &iatt);
      if (status == NC_NOERR) *remapOrder = iatt;
    }
  else if (cmpstr(map_method, "Bilinear") == 0)
    mapType = RemapMethod::BILINEAR;
  else if (cmpstr(map_method, "Bicubic") == 0)
    mapType = RemapMethod::BICUBIC;
  else if (cmpstr(map_method, "Distance") == 0)
    {
      mapType = RemapMethod::DISTWGT;
      *numNeighbors = 4;
    }
  else if (cmpstr(map_method, "Nearest") == 0)
    {
      mapType = RemapMethod::DISTWGT;
      *numNeighbors = 1;
    }
  else if (cmpstr(map_method, "Largest") == 0)
    {
      mapType = RemapMethod::CONSERV;
      *submapType = SubmapType::LAF;
    }
  else
    {
      cdoPrint("mapType = %s", map_method);
      cdoAbort("Invalid Map Type");
    }

  if (cdoVerbose) cdoPrint("mapType = %s", map_method);

  return mapType;
}
#endif

void
remapReadDataScrip(const char *interp_file, int gridID1, int gridID2, RemapMethod *mapType, SubmapType *submapType,
                   int *numNeighbors, int *remapOrder, RemapGrid &src_grid, RemapGrid &tgt_grid, RemapVars &rv)
{
// The routine reads a NetCDF file to extract remapping info in SCRIP format

#ifdef HAVE_LIBNETCDF

  // Local variables

  bool lgridarea = false;
  int status;
  int nc_srcgrdsize_id;    /* id for source grid size                  */
  int nc_dstgrdsize_id;    /* id for destination grid size             */
  int nc_srcgrdcorn_id;    /* id for number of source grid corners     */
  int nc_dstgrdcorn_id;    /* id for number of dest grid corners       */
  int nc_srcgrdrank_id;    /* id for source grid rank                  */
  int nc_dstgrdrank_id;    /* id for dest grid rank                    */
  int nc_numlinks_id;      /* id for number of links in mapping        */
  int nc_numwgts_id;       /* id for number of weights for mapping     */
  int nc_srcgrddims_id;    /* id for source grid dimensions            */
  int nc_dstgrddims_id;    /* id for dest grid dimensions              */
  int nc_srcgrdcntrlat_id; /* id for source grid center latitude       */
  int nc_dstgrdcntrlat_id; /* id for dest grid center latitude         */
  int nc_srcgrdcntrlon_id; /* id for source grid center longitude      */
  int nc_dstgrdcntrlon_id; /* id for dest grid center longitude        */
  int nc_srcgrdimask_id;   /* id for source grid mask                  */
  int nc_dstgrdimask_id;   /* id for dest grid mask                    */
  int nc_srcgrdcrnrlat_id; /* id for latitude of source grid corners   */
  int nc_srcgrdcrnrlon_id; /* id for longitude of source grid corners  */
  int nc_dstgrdcrnrlat_id; /* id for latitude of dest grid corners     */
  int nc_dstgrdcrnrlon_id; /* id for longitude of dest grid corners    */
  int nc_srcgrdarea_id;    /* id for area of source grid cells         */
  int nc_dstgrdarea_id;    /* id for area of dest grid cells           */
  int nc_srcgrdfrac_id;    /* id for area fraction on source grid      */
  int nc_dstgrdfrac_id;    /* id for area fraction on dest grid        */
  int nc_srcadd_id;        /* id for map source address                */
  int nc_dstadd_id;        /* id for map destination address           */
  int nc_rmpmatrix_id;     /* id for remapping matrix                  */

  char map_name[1024];
  char normalize_opt[64]; /* character string for normalization option */
  char convention[64];    /* character string for output convention    */
  char src_grid_name[64]; /* grid name for source grid                 */
  char tgt_grid_name[64]; /* grid name for dest   grid                 */
  char src_grid_units[64];
  char tgt_grid_units[64];
  size_t attlen, dimlen;

  int gridID1_gme_c = -1;

  // Open file and read some global information

  // nce(nc_open(interp_file, NC_NOWRITE, &nc_file_id));
  int nc_file_id = cdf_openread(interp_file);

  // Map name

  nce(nc_get_att_text(nc_file_id, NC_GLOBAL, "title", map_name));
  nce(nc_inq_attlen(nc_file_id, NC_GLOBAL, "title", &attlen));
  map_name[attlen] = 0;

  if (cdoVerbose)
    {
      cdoPrint("Reading remapping: %s", map_name);
      cdoPrint("From file: %s", interp_file);
    }

  // Map Tyoe
  *mapType = getMapType(nc_file_id, submapType, numNeighbors, remapOrder);

  if (*mapType == RemapMethod::CONSERV) lgridarea = true;

  remapVarsInit(*mapType, rv);

  rv.mapType = *mapType;
  rv.links_per_value = -1;
  rv.sort_add = false;

  // Normalization option
  nce(nc_get_att_text(nc_file_id, NC_GLOBAL, "normalization", normalize_opt));
  nce(nc_inq_attlen(nc_file_id, NC_GLOBAL, "normalization", &attlen));
  normalize_opt[attlen] = 0;

  if (strcmp(normalize_opt, "none") == 0)
    rv.normOpt = NormOpt::NONE;
  else if (strcmp(normalize_opt, "fracarea") == 0)
    rv.normOpt = NormOpt::FRACAREA;
  else if (strcmp(normalize_opt, "destarea") == 0)
    rv.normOpt = NormOpt::DESTAREA;
  else
    {
      cdoPrint("normalize_opt = %s", normalize_opt);
      cdoAbort("Invalid normalization option");
    }

  if (cdoVerbose) cdoPrint("normalize_opt = %s", normalize_opt);

  // File convention
  nce(nc_get_att_text(nc_file_id, NC_GLOBAL, "conventions", convention));
  nce(nc_inq_attlen(nc_file_id, NC_GLOBAL, "conventions", &attlen));
  convention[attlen] = 0;

  if (strcmp(convention, "SCRIP") != 0)
    {
      cdoPrint("convention = %s", convention);
      if (strcmp(convention, "NCAR-CSM") == 0)
        cdoAbort("Unsupported file convention!");
      else
        cdoAbort("Unknown file convention!");
    }

  // Read some additional global attributes

  // Source and destination grid names

  nce(nc_get_att_text(nc_file_id, NC_GLOBAL, "source_grid", src_grid_name));
  nce(nc_inq_attlen(nc_file_id, NC_GLOBAL, "source_grid", &attlen));
  src_grid_name[attlen] = 0;

  nce(nc_get_att_text(nc_file_id, NC_GLOBAL, "dest_grid", tgt_grid_name));
  nce(nc_inq_attlen(nc_file_id, NC_GLOBAL, "dest_grid", &attlen));
  tgt_grid_name[attlen] = 0;

  if (cdoVerbose) cdoPrint("Remapping between: %s and %s", src_grid_name, tgt_grid_name);

  // Initialize remapgrid structure
  remapGridInit(src_grid);
  remapGridInit(tgt_grid);

  // Read dimension information

  nce(nc_inq_dimid(nc_file_id, "src_grid_size", &nc_srcgrdsize_id));
  nce(nc_inq_dimlen(nc_file_id, nc_srcgrdsize_id, &dimlen));
  src_grid.size = dimlen;
  /*
  if (  src_grid.size != gridInqSize(gridID1) )
    cdoAbort("Source grids have different size!");
  */
  nce(nc_inq_dimid(nc_file_id, "dst_grid_size", &nc_dstgrdsize_id));
  nce(nc_inq_dimlen(nc_file_id, nc_dstgrdsize_id, &dimlen));
  tgt_grid.size = dimlen;
  /*
  if ( tgt_grid.size != gridInqSize(gridID2) )
    cdoAbort("Target grids have different size!");
  */
  status = nc_inq_dimid(nc_file_id, "src_grid_corners", &nc_srcgrdcorn_id);
  if (status == NC_NOERR)
    {
      nce(nc_inq_dimlen(nc_file_id, nc_srcgrdcorn_id, &dimlen));
      src_grid.num_cell_corners = dimlen;
      src_grid.luse_cell_corners = true;
      src_grid.lneed_cell_corners = true;
    }

  status = nc_inq_dimid(nc_file_id, "dst_grid_corners", &nc_dstgrdcorn_id);
  if (status == NC_NOERR)
    {
      nce(nc_inq_dimlen(nc_file_id, nc_dstgrdcorn_id, &dimlen));
      tgt_grid.num_cell_corners = dimlen;
      tgt_grid.luse_cell_corners = true;
      tgt_grid.lneed_cell_corners = true;
    }

  nce(nc_inq_dimid(nc_file_id, "src_grid_rank", &nc_srcgrdrank_id));
  nce(nc_inq_dimlen(nc_file_id, nc_srcgrdrank_id, &dimlen));
  src_grid.rank = dimlen;

  nce(nc_inq_dimid(nc_file_id, "dst_grid_rank", &nc_dstgrdrank_id));
  nce(nc_inq_dimlen(nc_file_id, nc_dstgrdrank_id, &dimlen));
  tgt_grid.rank = dimlen;

  nce(nc_inq_dimid(nc_file_id, "num_links", &nc_numlinks_id));
  nce(nc_inq_dimlen(nc_file_id, nc_numlinks_id, &dimlen));
  rv.num_links = dimlen;
  /*
  if ( rv.num_links == 0 )
    cdoAbort("Number of remap links is 0, no remap weights found!");
  */
  nce(nc_inq_dimid(nc_file_id, "num_wgts", &nc_numwgts_id));
  nce(nc_inq_dimlen(nc_file_id, nc_numwgts_id, &dimlen));
  rv.num_wts = dimlen;

  src_grid.gridID = gridID1;
  tgt_grid.gridID = gridID2;

  if (gridInqType(gridID1) == GRID_GME)
    {
      src_grid.nvgp = gridInqSize(gridID1);
      gridID1_gme_c = gridToUnstructured(gridID1, 1);
    }

  remapGridAlloc(rv.mapType, src_grid);
  remapGridAlloc(rv.mapType, tgt_grid);

  if (gridInqType(gridID1) == GRID_GME) gridInqMaskGME(gridID1_gme_c, &src_grid.vgpm[0]);

  // Allocate address and weight arrays for mapping 1
  if (rv.num_links > 0)
    {
      rv.max_links = rv.num_links;
      rv.src_cell_add.resize(rv.num_links);
      rv.tgt_cell_add.resize(rv.num_links);
      rv.wts.resize(rv.num_wts * rv.num_links);
    }

  // Get variable ids

  nce(nc_inq_varid(nc_file_id, "src_grid_dims", &nc_srcgrddims_id));
  nce(nc_inq_varid(nc_file_id, "src_grid_imask", &nc_srcgrdimask_id));
  nce(nc_inq_varid(nc_file_id, "src_grid_center_lat", &nc_srcgrdcntrlat_id));
  nce(nc_inq_varid(nc_file_id, "src_grid_center_lon", &nc_srcgrdcntrlon_id));
  if (src_grid.num_cell_corners)
    {
      nce(nc_inq_varid(nc_file_id, "src_grid_corner_lat", &nc_srcgrdcrnrlat_id));
      nce(nc_inq_varid(nc_file_id, "src_grid_corner_lon", &nc_srcgrdcrnrlon_id));
    }
  if (lgridarea)
    {
      nce(nc_inq_varid(nc_file_id, "src_grid_area", &nc_srcgrdarea_id));
    }
  nce(nc_inq_varid(nc_file_id, "src_grid_frac", &nc_srcgrdfrac_id));

  nce(nc_inq_varid(nc_file_id, "dst_grid_dims", &nc_dstgrddims_id));
  nce(nc_inq_varid(nc_file_id, "dst_grid_imask", &nc_dstgrdimask_id));
  nce(nc_inq_varid(nc_file_id, "dst_grid_center_lat", &nc_dstgrdcntrlat_id));
  nce(nc_inq_varid(nc_file_id, "dst_grid_center_lon", &nc_dstgrdcntrlon_id));
  if (tgt_grid.num_cell_corners)
    {
      nce(nc_inq_varid(nc_file_id, "dst_grid_corner_lat", &nc_dstgrdcrnrlat_id));
      nce(nc_inq_varid(nc_file_id, "dst_grid_corner_lon", &nc_dstgrdcrnrlon_id));
    }
  if (lgridarea)
    {
      nce(nc_inq_varid(nc_file_id, "dst_grid_area", &nc_dstgrdarea_id));
    }
  nce(nc_inq_varid(nc_file_id, "dst_grid_frac", &nc_dstgrdfrac_id));

  nce(nc_inq_varid(nc_file_id, "src_address", &nc_srcadd_id));
  nce(nc_inq_varid(nc_file_id, "dst_address", &nc_dstadd_id));
  nce(nc_inq_varid(nc_file_id, "remap_matrix", &nc_rmpmatrix_id));

  // Read all variables

  int dims[2];
  nce(nc_get_var_int(nc_file_id, nc_srcgrddims_id, dims));
  src_grid.dims[0] = dims[0];
  src_grid.dims[1] = dims[1];

  nce(nc_get_var_int(nc_file_id, nc_srcgrdimask_id, &src_grid.mask[0]));

  nce(nc_get_var_double(nc_file_id, nc_srcgrdcntrlat_id, src_grid.cell_center_lat));
  nce(nc_get_var_double(nc_file_id, nc_srcgrdcntrlon_id, src_grid.cell_center_lon));

  nce(nc_get_att_text(nc_file_id, nc_srcgrdcntrlat_id, "units", src_grid_units));
  nce(nc_inq_attlen(nc_file_id, nc_srcgrdcntrlat_id, "units", &attlen));
  src_grid_units[attlen] = 0;

  grid_to_radian(src_grid_units, src_grid.size, src_grid.cell_center_lon, "source grid center lon");
  grid_to_radian(src_grid_units, src_grid.size, src_grid.cell_center_lat, "source grid center lat");

  if (src_grid.num_cell_corners)
    {
      nce(nc_get_var_double(nc_file_id, nc_srcgrdcrnrlat_id, src_grid.cell_corner_lat));
      nce(nc_get_var_double(nc_file_id, nc_srcgrdcrnrlon_id, src_grid.cell_corner_lon));

      nce(nc_get_att_text(nc_file_id, nc_srcgrdcrnrlat_id, "units", src_grid_units));
      nce(nc_inq_attlen(nc_file_id, nc_srcgrdcrnrlat_id, "units", &attlen));
      src_grid_units[attlen] = 0;

      grid_to_radian(src_grid_units, src_grid.num_cell_corners * src_grid.size, src_grid.cell_corner_lon, "source grid corner lon");
      grid_to_radian(src_grid_units, src_grid.num_cell_corners * src_grid.size, src_grid.cell_corner_lat, "source grid corner lat");
    }

  if (lgridarea) nce(nc_get_var_double(nc_file_id, nc_srcgrdarea_id, &src_grid.cell_area[0]));

  nce(nc_get_var_double(nc_file_id, nc_srcgrdfrac_id, &src_grid.cell_frac[0]));

  nce(nc_get_var_int(nc_file_id, nc_dstgrddims_id, dims));
  tgt_grid.dims[0] = dims[0];
  tgt_grid.dims[1] = dims[1];

  nce(nc_get_var_int(nc_file_id, nc_dstgrdimask_id, &tgt_grid.mask[0]));

  nce(nc_get_var_double(nc_file_id, nc_dstgrdcntrlat_id, tgt_grid.cell_center_lat));
  nce(nc_get_var_double(nc_file_id, nc_dstgrdcntrlon_id, tgt_grid.cell_center_lon));

  nce(nc_get_att_text(nc_file_id, nc_dstgrdcntrlat_id, "units", tgt_grid_units));
  nce(nc_inq_attlen(nc_file_id, nc_dstgrdcntrlat_id, "units", &attlen));
  tgt_grid_units[attlen] = 0;

  grid_to_radian(tgt_grid_units, tgt_grid.size, tgt_grid.cell_center_lon, "target grid center lon");
  grid_to_radian(tgt_grid_units, tgt_grid.size, tgt_grid.cell_center_lat, "target grid center lat");

  if (tgt_grid.num_cell_corners)
    {
      nce(nc_get_var_double(nc_file_id, nc_dstgrdcrnrlat_id, tgt_grid.cell_corner_lat));
      nce(nc_get_var_double(nc_file_id, nc_dstgrdcrnrlon_id, tgt_grid.cell_corner_lon));

      nce(nc_get_att_text(nc_file_id, nc_dstgrdcrnrlat_id, "units", tgt_grid_units));
      nce(nc_inq_attlen(nc_file_id, nc_dstgrdcrnrlat_id, "units", &attlen));
      tgt_grid_units[attlen] = 0;

      grid_to_radian(tgt_grid_units, tgt_grid.num_cell_corners * tgt_grid.size, tgt_grid.cell_corner_lon, "target grid corner lon");
      grid_to_radian(tgt_grid_units, tgt_grid.num_cell_corners * tgt_grid.size, tgt_grid.cell_corner_lat, "target grid corner lat");
    }

  if (lgridarea) nce(nc_get_var_double(nc_file_id, nc_dstgrdarea_id, &tgt_grid.cell_area[0]));

  nce(nc_get_var_double(nc_file_id, nc_dstgrdfrac_id, &tgt_grid.cell_frac[0]));

  if (rv.num_links > 0)
    {
      readLinks(nc_file_id, nc_srcadd_id, rv.num_links, &rv.src_cell_add[0]);
      readLinks(nc_file_id, nc_dstadd_id, rv.num_links, &rv.tgt_cell_add[0]);

      for (size_t i = 0; i < rv.num_links; ++i)
        {
          rv.src_cell_add[i]--;
          rv.tgt_cell_add[i]--;
        }

      nce(nc_get_var_double(nc_file_id, nc_rmpmatrix_id, &rv.wts[0]));
    }

  // Close input file

  nce(nc_close(nc_file_id));

#else
  cdoAbort("NetCDF support not compiled in!");
#endif

}  // remapReadDataScrip
