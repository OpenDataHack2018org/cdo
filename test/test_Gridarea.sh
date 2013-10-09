#!/bin/sh
#
CDODEBUG=0
#
if [ "$CDODEBUG" == 0 ]; then CDO="$CDO -s"; fi
CDOOUT=cout
CDOERR=cerr
REFVAL="  5.10064e+14"
GRIDS="global_5 global_2 global_1 global_.5 n32 n80 n160"
RSTAT=0;
#
for GRID in $GRIDS; do
  GLOBAREA=`$CDO output -fldsum -gridarea -random,$GRID`
#  echo "$GRID: >$GLOBAREA< >$REFVAL<"
  if [ "$GLOBAREA" != "$REFVAL" ]; then RSTAT=`expr $RSTAT + 1`; fi
done
#
rm -f $CDOOUT $CDOERR
#
if [ "$CDODEBUG" == 1 ]; then
  echo "rstat: $RSTAT"
fi
#
exit $RSTAT
