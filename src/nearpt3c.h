#ifndef  NEARPT3_H_
#define  NEARPT3_H_

typedef float Coord_T;

#define NPT3SFACT 32000
#define NPT3SCALE(x) (0.5+(x+1)*NPT3SFACT)
//#define NPT3SFACT 1
//#define NPT3SCALE(x) x

void *nearpt3_preprocess(const int nfixpts, Coord_T **pts);
int nearpt3_query(void *g, const Coord_T *q);
void nearpt3_destroy(void *g);

#endif
