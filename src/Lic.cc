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

#include <cdi.h>

#include "cdo_int.h"
#include "grid.h"
#include "pstream_int.h"
#include "specspace.h"
#include "listarray.h"
#include "util_string.h"

#include "color.h"

CPT cpt;

#define SQUARE_FLOW_FIELD_SZ 400
#define DISCRETE_FILTER_SIZE 2048
#define LOWPASS_FILTR_LENGTH 10.00000f
#define LINE_SQUARE_CLIP_MAX 100000.0f
#define VECTOR_COMPONENT_MIN 0.050000f

///		normalize the vector field     ///
void
NormalizVectrs(int n_xres, int n_yres, float *pVectr)
{
  for (int j = 0; j < n_yres; j++)
    for (int i = 0; i < n_xres; i++)
      {
        int index = (j * n_xres + i) << 1;
        float vcMag = (float) (sqrt((double) (pVectr[index] * pVectr[index] + pVectr[index + 1] * pVectr[index + 1])));

        float scale = (vcMag == 0.0f) ? 0.0f : 1.0f / vcMag;
        pVectr[index] *= scale;
        pVectr[index + 1] *= scale;
      }
}

///		make white noise as the LIC input texture     ///
void
MakeWhiteNoise(int n_xres, int n_yres, unsigned char *pNoise)
{
  for (int j = 0; j < n_yres; j++)
    for (int i = 0; i < n_xres; i++)
      {
        int r = rand();
        r = ((r & 0xff) + ((r & 0xff00) >> 8)) & 0xff;
        pNoise[j * n_xres + i] = (unsigned char) r;
      }
}

///		generate box filter LUTs     ///
void
GenBoxFiltrLUT(int LUTsiz, float *p_LUT0, float *p_LUT1)
{
  for (int i = 0; i < LUTsiz; i++) p_LUT0[i] = p_LUT1[i] = i;
}

///		write the LIC image to a PPM file     ///
void
WriteImage2PPM(int n_xres, int n_yres, unsigned char *pImage, const char *f_name)
{
  FILE *o_file;
  if ((o_file = fopen(f_name, "w")) == NULL)
    {
      printf("Can't open output file\n");
      return;
    }

  fprintf(o_file, "P6\n%d %d\n255\n", n_xres, n_yres);

  for (int j = 0; j < n_yres; j++)
    for (int i = 0; i < n_xres; i++)
      {
        unsigned char unchar = pImage[j * n_xres + i];
        fprintf(o_file, "%c%c%c", unchar, unchar, unchar);
      }

  fclose(o_file);
}

///		write the LIC image to a PPM file     ///
#include <png.h>

void
setRGBmag(png_byte *ptr, float val, float mag)
{
  int r = 0, g = 0, b = 0, n;

  for (n = 0; n < cpt.ncolors; n++)
    if (mag > cpt.lut[n].z_low && mag <= cpt.lut[n].z_high) break;

  if (n == cpt.ncolors)
    {
      r = cpt.bfn[0].rgb[0];
      g = cpt.bfn[0].rgb[1];
      b = cpt.bfn[0].rgb[2];
    }
  else
    {
      r = cpt.lut[n].rgb_high[0];
      g = cpt.lut[n].rgb_high[1];
      b = cpt.lut[n].rgb_high[2];
    }

  float foregroundAlpha = .8;
  r = (val * foregroundAlpha) + (r * (1.0 - foregroundAlpha));
  g = (val * foregroundAlpha) + (g * (1.0 - foregroundAlpha));
  b = (val * foregroundAlpha) + (b * (1.0 - foregroundAlpha));
  ptr[0] = r;
  ptr[1] = g;
  ptr[2] = b;
}

int
WriteImage2PNG(int width, int height, unsigned char *pImage, float *mag, const char *filename)
{
  int code;
  png_structp png_ptr = NULL;
  png_infop info_ptr = NULL;
  png_bytep row = NULL;

  // Open file for writing (binary mode)
  FILE *fp = fopen(filename, "wb");
  if (fp == NULL)
    {
      fprintf(stderr, "Could not open file %s for writing\n", filename);
      code = 1;
      goto finalise;
    }

  // Initialize write structure
  png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
  if (png_ptr == NULL)
    {
      fprintf(stderr, "Could not allocate write struct\n");
      code = 1;
      goto finalise;
    }

  // Initialize info structure
  info_ptr = png_create_info_struct(png_ptr);
  if (info_ptr == NULL)
    {
      fprintf(stderr, "Could not allocate info struct\n");
      code = 1;
      goto finalise;
    }

  // Setup Exception handling
  if (setjmp(png_jmpbuf(png_ptr)))
    {
      fprintf(stderr, "Error during png creation\n");
      code = 1;
      goto finalise;
    }

  png_init_io(png_ptr, fp);

  // Write header (8 bit colour depth)
  png_set_IHDR(png_ptr, info_ptr, width, height, 8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_BASE,
               PNG_FILTER_TYPE_BASE);

  // Set title
  /*
  if (title != NULL) {
          png_text title_text;
          title_text.compression = PNG_TEXT_COMPRESSION_NONE;
          title_text.key = "Title";
          title_text.text = title;
          png_set_text(png_ptr, info_ptr, &title_text, 1);
  }
  */
  png_write_info(png_ptr, info_ptr);

  // Allocate memory for one row (3 bytes per pixel - RGB)
  row = (png_bytep) malloc(3 * width * sizeof(png_byte));

  // Write image data
  for (int y = 0; y < height; y++)
    {
      for (int x = 0; x < width; x++)
        {
          setRGBmag(&(row[x * 3]), pImage[y * width + x], mag[y * width + x]);
        }
      png_write_row(png_ptr, row);
    }

  // End write
  png_write_end(png_ptr, NULL);

finalise:
  if (fp != NULL) fclose(fp);
  if (info_ptr != NULL) png_free_data(png_ptr, info_ptr, PNG_FREE_ALL, -1);
  if (png_ptr != NULL) png_destroy_write_struct(&png_ptr, (png_infopp) NULL);
  if (row != NULL) free(row);

  return code;
}

///		flow imaging (visualization) through Line Integral Convolution     ///
void
FlowImagingLIC(int n_xres, int n_yres, float *pVectr, unsigned char *pNoise, unsigned char *pImage, float *p_LUT0, float *p_LUT1,
               float krnlen)
{

  /// for each pixel in the 2D output LIC image///

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(pNoise, pVectr, n_xres, n_yres, pImage, p_LUT0, p_LUT1, krnlen)
#endif
  for (int j = 0; j < n_yres; j++)
    {
  int vec_id;                       /// ID in the VECtor buffer (for the input flow field)
  int advDir;                       /// ADVection DIRection (0: positive;  1: negative)
  int advcts;                       /// number of ADVeCTion stepS per direction (a step counter)
  int ADVCTS = (int) (krnlen * 3);  /// MAXIMUM number of advection steps per direction to break dead loops

  float vctr_x;                                        /// x-component  of the VeCToR at the forefront point
  float vctr_y;                                        /// y-component  of the VeCToR at the forefront point
  float clp0_x;                                        /// x-coordinate of CLiP point 0 (current)
  float clp0_y;                                        /// y-coordinate of CLiP point 0	(current)
  float clp1_x;                                        /// x-coordinate of CLiP point 1 (next   )
  float clp1_y;                                        /// y-coordinate of CLiP point 1 (next   )
  float samp_x;                                        /// x-coordinate of the SAMPle in the current pixel
  float samp_y;                                        /// y-coordinate of the SAMPle in the current pixel
  float tmpLen;                                        /// TeMPorary LENgth of a trial clipped-segment
  float segLen;                                        /// SEGment   LENgth
  float curLen;                                        /// CURrent   LENgth of the streamline
  float prvLen;                                        /// PReVious  LENgth of the streamline
  float W_ACUM;                                        /// ACcuMulated Weight from the seed to the current streamline forefront
  float texVal;                                        /// TEXture VALue
  float smpWgt;                                        /// WeiGhT of the current SaMPle
  float t_acum[2];                                     /// two ACcUMulated composite Textures for the two directions, perspectively
  float w_acum[2];                                     /// two ACcUMulated Weighting values   for the two directions, perspectively
  float *wgtLUT = NULL;                                /// WeiGhT Look Up Table pointing to the target filter LUT
  float len2ID = (DISCRETE_FILTER_SIZE - 1) / krnlen;  /// map a curve LENgth TO an ID in the LUT
    for (int i = 0; i < n_xres; i++)
      {
        /// init the composite texture accumulators and the weight accumulators///
        t_acum[0] = t_acum[1] = w_acum[0] = w_acum[1] = 0.0f;

        /// for either advection direction///
        for (advDir = 0; advDir < 2; advDir++)
          {
            /// init the step counter, curve-length measurer, and streamline seed///
            advcts = 0;
            curLen = 0.0f;
            clp0_x = i + 0.5f;
            clp0_y = j + 0.5f;

            /// access the target filter LUT///
            wgtLUT = (advDir == 0) ? p_LUT0 : p_LUT1;

            /// until the streamline is advected long enough or a tightly  spiralling center / focus is encountered///
            while (curLen < krnlen && advcts < ADVCTS)
              {
                /// access the vector at the sample///
                vec_id = ((int) (clp0_y) *n_xres + (int) (clp0_x)) << 1;
                vctr_x = pVectr[vec_id];
                vctr_y = pVectr[vec_id + 1];

                /// in case of a critical point///
                if (vctr_x == 0.0f && vctr_y == 0.0f)
                  {
                    t_acum[advDir] = (advcts == 0) ? 0.0f : t_acum[advDir];  /// this line is indeed unnecessary
                    w_acum[advDir] = (advcts == 0) ? 1.0f : w_acum[advDir];
                    break;
                  }

                /// negate the vector for the backward-advection case///
                vctr_x = (advDir == 0) ? vctr_x : -vctr_x;
                vctr_y = (advDir == 0) ? vctr_y : -vctr_y;

                /// clip the segment against the pixel boundaries --- find the shorter from the two clipped segments///
                /// replace  all  if-statements  whenever  possible  as  they  might  affect the computational speed///
                segLen = LINE_SQUARE_CLIP_MAX;
                segLen = (vctr_x < -VECTOR_COMPONENT_MIN) ? ((int) (clp0_x) -clp0_x) / vctr_x : segLen;
                segLen = (vctr_x > VECTOR_COMPONENT_MIN) ? ((int) ((int) (clp0_x) + 1.5f) - clp0_x) / vctr_x : segLen;
                segLen = (vctr_y < -VECTOR_COMPONENT_MIN)
                             ? (((tmpLen = ((int) (clp0_y) -clp0_y) / vctr_y) < segLen) ? tmpLen : segLen)
                             : segLen;
                segLen = (vctr_y > VECTOR_COMPONENT_MIN)
                             ? (((tmpLen = ((int) ((int) (clp0_y) + 1.5f) - clp0_y) / vctr_y) < segLen) ? tmpLen : segLen)
                             : segLen;

                /// update the curve-length measurers///
                prvLen = curLen;
                curLen += segLen;
                segLen += 0.0004f;

                /// check if the filter has reached either end///
                segLen = (curLen > krnlen) ? ((curLen = krnlen) - prvLen) : segLen;

                /// obtain the next clip point///
                clp1_x = clp0_x + vctr_x * segLen;
                clp1_y = clp0_y + vctr_y * segLen;

                /// obtain the middle point of the segment as the texture-contributing sample///
                samp_x = (clp0_x + clp1_x) * 0.5f;
                samp_y = (clp0_y + clp1_y) * 0.5f;

                /// obtain the texture value of the sample///
                texVal = pNoise[(int) (samp_y) *n_xres + (int) (samp_x)];

                /// update the accumulated weight and the accumulated composite texture (texture x weight)///
                W_ACUM = wgtLUT[(int) (curLen * len2ID)];
                smpWgt = W_ACUM - w_acum[advDir];
                w_acum[advDir] = W_ACUM;
                t_acum[advDir] += texVal * smpWgt;

                /// update the step counter and the "current" clip point///
                advcts++;
                clp0_x = clp1_x;
                clp0_y = clp1_y;

                /// check if the streamline has gone beyond the flow field///
                if (clp0_x < 0.0f || clp0_x >= n_xres || clp0_y < 0.0f || clp0_y >= n_yres) break;
              }
          }

        /// normalize the accumulated composite texture///
        texVal = (t_acum[0] + t_acum[1]) / (w_acum[0] + w_acum[1]);

        /// clamp the texture value against the displayable intensity range [0, 255]
        texVal = (texVal < 0.0f) ? 0.0f : texVal;
        texVal = (texVal > 255.0f) ? 255.0f : texVal;
        pImage[j * n_xres + i] = (unsigned char) texVal;
      }
    }
}

void
lic1(int num, int n_xres, int n_yres, float *pVectr)
{
  float *p_LUT0 = (float *) malloc(sizeof(float) * DISCRETE_FILTER_SIZE);
  float *p_LUT1 = (float *) malloc(sizeof(float) * DISCRETE_FILTER_SIZE);
  unsigned char *pNoise = (unsigned char *) malloc(sizeof(unsigned char) * n_xres * n_yres);
  unsigned char *pImage = (unsigned char *) malloc(sizeof(unsigned char) * n_xres * n_yres);

  float *mag = (float *) malloc(n_xres * n_yres * sizeof(float));
  int n = n_xres * n_yres;

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(mag, pVectr, n)
#endif
  for (int i = 0; i < n; ++i)
    mag[i] = (float)sqrt((double) pVectr[i * 2] * pVectr[i * 2] + (double) pVectr[i * 2 + 1] * pVectr[i * 2 + 1]);

  if (cdoVerbose)
    {
      float vmin = 1.e33;
      float vmax = -1.e33;
      for (int i = 0; i < n_xres * n_yres; ++i)
        {
          if (mag[i] < vmin) vmin = mag[i];
          if (mag[i] > vmax) vmax = mag[i];
        }
      cdoPrint("ts=%d minval=%g  maxval=%g", num+1, vmin, vmax);
    }

  NormalizVectrs(n_xres, n_yres, pVectr);
  MakeWhiteNoise(n_xres, n_yres, pNoise);
  GenBoxFiltrLUT(DISCRETE_FILTER_SIZE, p_LUT0, p_LUT1);
  FlowImagingLIC(n_xres, n_yres, pVectr, pNoise, pImage, p_LUT0, p_LUT1, LOWPASS_FILTR_LENGTH);
  //const char *fname = "LIC.ppm";
  //WriteImage2PPM(n_xres, n_yres, pImage, fname);
  char filename[256];
  sprintf(filename, "image%04d.png", num);
  WriteImage2PNG(n_xres, n_yres, pImage, mag, filename);

  free(p_LUT0);
  free(p_LUT1);
  free(pNoise);
  free(pImage);
  free(mag);
}

void *
Lic(void *process)
{
  int nrecs;
  int varID, levelID;
  int nlev = 0;
  size_t nmiss;
  int varID1 = -1, varID2 = -1;
  size_t offset;

  cdoInitialize(process);

  // clang-format off
  // int LIC  = cdoOperatorAdd("lic",  0, 0, NULL);
  // clang-format on

  // int operatorID = cdoOperatorID();

  operatorCheckArgc(1);
  char *cpt_file = operatorArgv()[0];

  FILE *cpt_fp = fopen(cpt_file, "r");
  if (cpt_fp == NULL) cdoAbort("Open failed on color palette table %s", cpt_file);

  int status = cptRead(cpt_fp, &cpt);
  if (status != 0) cdoAbort("Error during read of color palette table %s", cpt_file);

  int streamID1 = cdoStreamOpenRead(cdoStreamName(0));

  int vlistID1 = cdoStreamInqVlist(streamID1);

  int taxisID1 = vlistInqTaxis(vlistID1);

  /* find variables */
  int nvars = vlistNvars(vlistID1);
  if (nvars != 2) cdoAbort("Need 2 input variable, found %d\n", nvars);
  varID1 = 0;
  varID2 = 1;

  int ngrids = vlistNgrids(vlistID1);
  if (ngrids != 1) cdoAbort("Need 1 grid, found %d\n", ngrids);

  int gridID = vlistInqVarGrid(vlistID1, varID1);
  int zaxisID = vlistInqVarZaxis(vlistID1, varID1);
  /*
  int vlistID2 = vlistCreate();
  int outvarID = vlistDefVar(vlistID2, gridID, zaxisID, TIME_VARYING);

  int streamID2 = cdoStreamOpenWrite(cdoStreamName(1), cdoFiletype());

  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  pstreamDefVlist(streamID2, vlistID2);
  */
  size_t gridsize = vlistGridsizeMax(vlistID1);
  std::vector<double> array1(gridsize);
  std::vector<float> varray(2 * gridsize);

  std::vector<double> ivar1, ivar2, ovar1;
  if (varID1 != -1 && varID2 != -1)
    {
      nlev = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID1));

      gridsize = gridInqSize(gridID);
      ivar1.resize(nlev * gridsize);
      ivar2.resize(nlev * gridsize);

      gridsize = gridInqSize(gridID);
      ovar1.resize(nlev * gridsize);
    }

  int gridtype = gridInqType(gridID);
  if (gridtype != GRID_LONLAT && gridtype != GRID_GAUSSIAN) cdoAbort("Unsupported grid!\n");

  size_t nx = gridInqXsize(gridID);
  size_t ny = gridInqYsize(gridID);

  int tsID = 0;
  while ((nrecs = cdoStreamInqTimestep(streamID1, tsID)))
    {
      /*
      taxisCopyTimestep(taxisID2, taxisID1);
      pstreamDefTimestep(streamID2, tsID);
      */
      for (int recID = 0; recID < nrecs; recID++)
        {
          pstreamInqRecord(streamID1, &varID, &levelID);

          if (varID == varID1 || varID == varID2)
            {
              pstreamReadRecord(streamID1, array1.data(), &nmiss);
              if (nmiss) cdoAbort("Missing values unsupported for spectral data!");

              gridsize = gridInqSize(gridID);
              offset = gridsize * levelID;

              if (varID == varID1)
                for (size_t i = 0; i < gridsize; ++i) varray[i * 2] = (float) array1[i];
              else if (varID == varID2)
                for (size_t i = 0; i < gridsize; ++i) varray[i * 2 + 1] = (float) array1[i];
            }
        }

      lic1(tsID, nx, ny, varray.data());
      /*
      gridsize = gridInqSize(gridID);
      for (levelID = 0; levelID < nlev; levelID++)
        {
          offset = gridsize * levelID;
          pstreamDefRecord(streamID2, outvarID, levelID);
          pstreamWriteRecord(streamID2, &ovar1[offset], 0);
        }
      */

      tsID++;
    }

  // pstreamClose(streamID2);
  pstreamClose(streamID1);

  cdoFinish();

  return 0;
}

// mencoder
// png + alpha canel libpng
// libpng.org
// global_.06
// ls -1v | grep JPG > files.txt
// mencoder -nosound -ovc lavc -lavcopts vcodec=mpeg4:vbitrate=21600000 -o windowsill_flowers_7.avi -mf type=jpeg:fps=24
// mf://@files.txt -vf scale=1920:1080 convert imageA.png imageB.png -alpha off -compose CopyOpacity -composite out.png makecpt
// -Cjet -T-10/30/.2  > cjet.cpt makecpt -Chot -T0/30/.2 -I > chot.cpt ffmpeg -r 1/5 -start_number 0 -i "image%4d.png" -c:v libx264
// -video_size 4k -vf "fps=25,format=yuv420p" movie.mp4
