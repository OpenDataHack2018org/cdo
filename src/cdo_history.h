#ifndef CDO_HISTORY_H
#define CDO_HISTORY_H

void cdoInqHistory(int fileID);
void cdoDefHistory(int fileID, char *histstring);
void cdo_def_creation_date(int vlistID);
void cdo_def_tracking_id(int vlistID, const char *uuid_attribute);
#endif
