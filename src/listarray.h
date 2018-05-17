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

#ifndef LISTARRAY_H
#define LISTARRAY_H

void split_intstring(const char *intstr, int *first, int *last, int *inc);

template <typename T>
class ListArray
{
private:
  std::vector<T> array;
  int nalloc;
  int allinc;

  void
  ensureArraySize(int num)
  {
    while (nalloc <= num)
      {
        nalloc += allinc;
        array.resize(nalloc);
      }
  }

public:
  ListArray()
  {
    nalloc = 0;
    allinc = 1024;
  }

  T *
  data()
  {
    return array.data();
  }

  void
  setValue(int num, T value)
  {
    ensureArraySize(num);
    array[num] = value;
  }

  T
  getValue(int num)
  {
    return array[num];
  }

  int
  argvToInt(int argc, char **argv)
  {
    int nint = 0;

    for (int iarg = 0; iarg < argc; iarg++)
      {
        int first, last, inc;
        split_intstring(argv[iarg], &first, &last, &inc);

        if (inc >= 0)
          {
            for (int ival = first; ival <= last; ival += inc) setValue(nint++, ival);
          }
        else
          {
            for (int ival = first; ival >= last; ival += inc) setValue(nint++, ival);
          }
      }

    return nint;
  }

  int
  argvToFlt(int argc, char **argv)
  {
    int nint = 0;

    for (int iarg = 0; iarg < argc; iarg++)
      {
        int len = (int) strlen(argv[iarg]);
        int i;
        for (i = 0; i < len; i++)
          if (argv[iarg][i] != '/' && argv[iarg][i] != '-' && !isdigit(argv[iarg][i])) break;

        if (i != len)
          {
            /*
              if      ( strcmp(argv[iarg],  "inf") == 0 )
              tmp_val =  DBL_MAX;
              else if ( strcmp(argv[iarg], "-inf") == 0 )
              tmp_val = -DBL_MAX;
              else
            */
            double tmp_val = parameter2double(argv[iarg]);
            setValue(nint++, tmp_val);
          }
        else
          {
            int first, last, inc;
            split_intstring(argv[iarg], &first, &last, &inc);

            if (inc >= 0)
              {
                for (int ival = first; ival <= last; ival += inc) setValue(nint++, (double) ival);
              }
            else
              {
                for (int ival = first; ival >= last; ival += inc) setValue(nint++, (double) ival);
              }
          }
      }

    return nint;
  }
};

#endif /* LISTARRAY_H */
