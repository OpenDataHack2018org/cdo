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
#ifndef _STDNAMETABLE_H
#define _STDNAMETABLE_H

enum stdnameid
{
  air_pressure,
  pressure_thickness,
  surface_geopotential,
  geopotential,
  air_temperature,
  specific_humidity,
  surface_air_pressure,
  air_pressure_at_sea_level,
  geopotential_height
};

int var_echamcode(int varid);
const char *var_name(int varid);
const char *var_stdname(int varid);
const char *var_units(int varid);

int echamcode_from_stdname(const char *stdname);

typedef struct
{
  int geopot;
  int temp;
  int hum;
  int ps;
  int lsp;
  int gheight;
  int wind;
  int uwind;
  int vwind;
} gribcode_t;

enum struct ModelMode
{
  UNDEF,
  ECHAM,
  WMO,
  HIRLAM
};

void echam_gribcodes(gribcode_t *gribcodes);
void wmo_gribcodes(gribcode_t *gribcodes);
void hirlam_harmonie_gribcodes(gribcode_t *gribcodes);

#endif
