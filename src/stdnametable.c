#include <stdio.h>
#include <assert.h>
#include "stdnametable.h"


typedef struct
{
  int   varid;
  int   echamcode;
  char *name;
  char *stdname;     /* Standard name */
  char *units;       /* Units         */
}
stdnametable_t;


const stdnametable_t stdnametable[] = {
  /* varid                       code    name                standard name                 units */
  { surface_geopotential,         129,  "geosp",            "surface_geopotential",       "m2 s-2" },
  { air_temperature,              130,  "airtemperature",   "air_temperature",            "K" },
  { surface_air_pressure,         134,  "aps",              "surface_air_pressure",       "Pa" },
  { air_pressure_at_sea_level,    151,  "sealevelpressure", "air_pressure_at_sea_level",  "Pa" },
  { geopotential_height,          156,  "geopotheight",     "geopotential_height",        "m" },
};


static int stdnametable_idx(int varid)
{
  int idx;
  int num_entries = (int) (sizeof(stdnametable)/sizeof(stdnametable_t));

  for ( idx = 0; idx < num_entries; ++idx )
    if ( stdnametable[idx].varid == varid ) break;

  assert( idx < num_entries );

  return (idx);
}


int var_echamcode(int varid)
{
  return (stdnametable[stdnametable_idx(varid)].echamcode);
}

const char* var_name(int varid)
{
  return (stdnametable[stdnametable_idx(varid)].name);
}

const char* var_stdname(int varid)
{
  return (stdnametable[stdnametable_idx(varid)].stdname);
}

const char* var_units(int varid)
{
  return (stdnametable[stdnametable_idx(varid)].units);
}
