void processDefVarNum(int nvars);
int processInqVarNum();
int processInqTimesteps(void);
void processDefTimesteps(int streamID);
int operatorArgc(void);
char ** operatorArgv(void);
void operatorCheckArgc(int numargs);
void operatorInputArg(const char *enter);
int cdoOperatorAdd(const char *name, int f1, int f2, const char *enter);
int cdoOperatorID(void);
int cdoOperatorF1(int operID);
int cdoOperatorF2(int operID);
const char *cdoOperatorName(int operID);
const char *cdoOperatorEnter(int operID);
int cdoStreamNumber();
int cdoStreamCnt(void);
int cdoStreamName(int cnt);

