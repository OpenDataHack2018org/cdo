//typedef unsigned short int Coord_T;
typedef float Coord_T;

#if defined(__cplusplus)
extern "C" {
#endif
void *nearpt3_preprocess(const int nfixpts, Coord_T **pts);
int nearpt3_query(void *g, const Coord_T *q);
#if defined(__cplusplus)
}
#endif
