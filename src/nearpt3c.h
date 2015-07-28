typedef unsigned short int Coord_T;

#if defined(__cplusplus)
extern "C" {
#endif
void *nearpt3_preprocess(const int nfixpts, unsigned short int **pts);
int nearpt3_query(void *g, const unsigned short int *q);
#if defined(__cplusplus)
}
#endif
