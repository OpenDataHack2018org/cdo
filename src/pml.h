/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2016 Uwe Schulzweida, <uwe.schulzweida AT mpimet.mpg.de>
  See COPYING file for copying and redistribution conditions.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/

#ifndef _PML_H
#define _PML_H

#define  PML_INT         1
#define  PML_FLT         2
#define  PML_WORD        3

#define  PML_INIT(name, size)                 memset(flag_##name, 0, size * sizeof(bool))
#define  PML_DEF(name, size)                  bool flag_##name[size]
#define  PML_DEF_INT(name, size)              int par_##name[size]; int name = 0; PML_DEF(name, size); PML_INIT(name, size)
#define  PML_DEF_FLT(name, size)              double par_##name[size]; double name = 0; PML_DEF(name, size); PML_INIT(name, size)
#define  PML_DEF_WORD(name, size)             char *par_##name[size]; const char *name = 0; PML_DEF(name, size); PML_INIT(name, size)
#define  PML_ADD_INT(pml, name, size, txt)    PML_DEF_INT(name, size);  int pid_##name = pml_add(pml, txt, #name, PML_INT,  par_##name, sizeof(par_##name)/sizeof(int))
#define  PML_ADD_FLT(pml, name, size, txt)    PML_DEF_FLT(name, size);  int pid_##name = pml_add(pml, txt, #name, PML_FLT,  par_##name, sizeof(par_##name)/sizeof(double))
#define  PML_ADD_WORD(pml, name, size, txt)   PML_DEF_WORD(name, size); int pid_##name = pml_add(pml, txt, #name, PML_WORD, par_##name, sizeof(par_##name)/sizeof(char *))
#define  PML_NOCC(pml, name)                  pml_num_par(pml, pid_##name)
//#define  PML_NUM(pml, name)                   npar_##name = pml_num(pml, #name)

#define  PML_CHECK_INT_FLAG(pml, name)        pml_check_int_flag(pml, pid_##name, par_##name, flag_##name)
#define  PML_CHECK_FLT_FLAG(pml, name)        pml_check_flt_flag(pml, pid_##name, par_##name, flag_##name)
#define  PML_CHECK_WORD_FLAG(pml, name)       pml_check_word_flag(pml, pid_##name, par_##name, flag_##name)
#define  PML_CHECK_INT(pml, name)             pml_check_int(pml, pid_##name, par_##name, flag_##name, name)
#define  PML_CHECK_FLT(pml, name)             pml_check_flt(pml, pid_##name, par_##name, flag_##name, name)
#define  PML_CHECK_WORD(pml, name)            pml_check_word(pml, pid_##name, par_##name, flag_##name, name)
#define  PML_CHECK_DATE(pml, name)            pml_check_date(pml, pid_##name, par_##name, flag_##name, name)
#define  PML_CHECK_SEASON(pml, name, month)   pml_check_season(pml, pid_##name, par_##name, flag_##name, month)

#define  MAX_PML_ENTRY    256

typedef struct
{
  char *txt;
  char *name;
  size_t len;
  void *ptr;
  int type;
  int occ;
  size_t size;
} pml_entry_t;


typedef struct
{
  int size;
  char *name;
  pml_entry_t *entry[MAX_PML_ENTRY];
} pml_t;


pml_t *pml_create(const char *name);
void pml_destroy(pml_t *pml);

void pml_print(pml_t *pml);
int pml_add(pml_t *pml, const char *txt, const char *name, int type, void *ptr, size_t size);
int pml_num_par(pml_t *pml, int pid);
int pml_add_entry(pml_entry_t *entry, char *arg);
int pml_read(pml_t *pml, int argc, char **argv);

bool pml_check_int(pml_t *pml, int pid, int *parlist, bool *flaglist, int par);
bool pml_check_flt(pml_t *pml, int pid, double *parlist, bool *flaglist, double par);
bool pml_check_word(pml_t *pml, int pid, char **parlist, bool *flaglist, const char *par);
bool pml_check_date(pml_t *pml, int pid, char **parlist, bool *flaglist, const char *par);
bool pml_check_season(pml_t *pml, int pid, char **parlist, bool *flaglist, int month);

void pml_check_int_flag(pml_t *pml, int pid, int *parlist, bool *flaglist);
void pml_check_flt_flag(pml_t *pml, int pid, double *parlist, bool *flaglist);
void pml_check_word_flag(pml_t *pml, int pid, char **parlist, bool *flaglist);

#endif  /* _PMLIST_H */
