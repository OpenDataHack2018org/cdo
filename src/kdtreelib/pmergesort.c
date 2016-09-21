#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "kdtree.h"


typedef struct param_t {
    char *a;
    char *b;
    size_t first;
    size_t nmemb;
    size_t size;
#if defined(KDTEST)
    int axis;
#else
    int (*cmp) (const void *, const void *);
#endif
    int max_threads;
} param_t;

extern void qsortR(const void *base0, size_t n, size_t size, int axis);

void pm_buildparams(struct param_t *p, void *a, void *b, size_t first,
                    size_t nmemb, size_t size,
#if defined(KDTEST)
                    int axis,
#else
                    int (*cmp) (const void *, const void *),
#endif                   
                    int max_threads);
int pmergesort(void *base, size_t nmemb, size_t size,
#if defined(KDTEST)
               int axis,
#else
               int (*compar) (const void *, const void *),
#endif
               int max_threads);
void *mergesort_t(void *args);


int
pmergesort(void *base, size_t nmemb, size_t size,
#if defined(KDTEST)
           int axis,
#else
           int (*compar) (const void *, const void *),
#endif
           int max_threads)
{
    void *tmp;
    param_t args;

    if ((tmp = calloc(nmemb, size)) == NULL) {
        perror("malloc");
        return 0;
    }
    args.a = (char *) base;
    args.b = (char *) tmp;
    args.first = 0;
    args.nmemb = nmemb;
    args.size = size;
#if defined(KDTEST)
    args.axis = axis;
#else
    args.cmp = compar;
#endif
    args.max_threads = max_threads;

    mergesort_t(&args);

    free(tmp);
    return 1;
}

void
pm_buildparams(struct param_t *p, void *a, void *b, size_t first,
               size_t nmemb, size_t size,
#if defined(KDTEST)
               int axis,
#else
               int (*cmp) (const void *, const void *),
#endif
               int max_threads)
{

    p->a = (char *)a;
    p->b = (char *)b;
    p->first = first;
    p->nmemb = nmemb;
    p->size = size;
#if defined(KDTEST)
    p->axis = axis;
#else
    p->cmp = cmp;
#endif
    p->max_threads = max_threads;
}


void *
mergesort_t(void *args)
{
    size_t i, li, ri;
    struct param_t larg, rarg;
    pthread_t thr[2];
    param_t *mya = (param_t *) args;

    if (mya->max_threads < 2) {
        /*
         * Reached maximum number of threads allocated to this
         * branch. Proceed with sequential sort of this chunk. 
         */
#if defined(KDTEST)
        qsortR(mya->a + mya->first * mya->size, mya->nmemb, mya->size, mya->axis);
#else
        qsort(mya->a + mya->first * mya->size, mya->nmemb, mya->size, mya->cmp);
#endif
    } else {
        /*
         * Start two new threads, each sorting half of array a 
         */
        pm_buildparams(&larg, mya->a, mya->b, mya->first, mya->nmemb / 2,
                       mya->size,
#if defined(KDTEST)
                       mya->axis,
#else
                       mya->cmp,
#endif
                       mya->max_threads / 2);
        /*
         * Recursively sort the left half 
         */
        if (pthread_create(&thr[0], NULL, mergesort_t, (void *) &larg)) {
            perror("pthread_create");
            return NULL;
        }

        pm_buildparams(&rarg, mya->a, mya->b, mya->first + mya->nmemb / 2,
                       mya->nmemb - mya->nmemb / 2, mya->size,
#if defined(KDTEST)
                       mya->axis,
#else
                       mya->cmp,
#endif
                       mya->max_threads / 2);
        /*
         * Recursively sort the right half 
         */
        if (pthread_create(&thr[1], NULL, mergesort_t, (void *) &rarg)) {
            perror("pthread_create");
            return NULL;
        }

        pthread_join(thr[0], NULL);
        pthread_join(thr[1], NULL);

        /*
         * Merge the two sorted chunks of array a into array b 
         */
        li = larg.first;
        ri = rarg.first;
        for(i = mya->first; i < mya->first + mya->nmemb; i++) {
            if (li >= larg.first + larg.nmemb) {
                /*
                 * We already copied everything from the left chunk,
                 * now copy from the right 
                 */
                memcpy(mya->b + i * mya->size, mya->a + ri * mya->size,
                       mya->size);
                ri++;
            } else if (ri >= rarg.first + rarg.nmemb) {
                /*
                 * We already copied everything from the right chunk,
                 * now copy from the left 
                 */
                memcpy(mya->b + i * mya->size, mya->a + li * mya->size,
                       mya->size);
                li++;
            }
            /*
             * We can still copy from both chunks, copy the smaller
             * element 
             */
            else if (
#if defined(KDTEST)
                     qcmp(mya->a + li * mya->size,
                          mya->a + ri * mya->size,
                          mya->axis)
#else
                     mya->cmp(mya->a + li * mya->size,
                              mya->a + ri * mya->size)
#endif
                               < 1) {
                memcpy(mya->b + i * mya->size, mya->a + li * mya->size,
                       mya->size);
                li++;
            } else {
                memcpy(mya->b + i * mya->size, mya->a + ri * mya->size,
                       mya->size);
                ri++;
            }
        }
        /*
         * Now b is sorted, copy it back to a 
         */
        memcpy(mya->a + mya->size * mya->first,
               mya->b + mya->size * mya->first, mya->size * mya->nmemb);
    }
    return NULL;
}
