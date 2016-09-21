/*!\file kdtree_common.c
 * \brief Routines common to the Cartesian and spherical kd-tree versions
 *
 */

#include "kdtree.h"
#include "pqueue.h"


extern int pmergesort(void *base, size_t nmemb, size_t size,
#if defined(KDTEST)
                      int (*compar) (const void *, const void *, int), int axis,
#else
                      int (*compar) (const void *, const void *),
#endif
                      int max_threads);


/* ********************************************************************

   general utility functions 

   ********************************************************************* */

void *
kd_malloc(size_t size, const char *msg)
{
    void *ptr;
    if ((ptr = malloc(size)) == NULL)
        perror(msg);
    return ptr;
}

kdata_t
kd_sqr(kdata_t x)
{
    return x == 0 ? 0 : x * x;
}

int
kd_isleaf(struct kdNode *n)
{
    return (n->left || n->right) ? 0 : 1;
}

/* end utility functions */

/* *******************************************************************
   
   Helper functions for debugging

   ******************************************************************* */

void
kd_printNode(struct kdNode *node)
{
    if (!node) {
        fprintf(stderr, "Node is empty.\n");
        return;
    }
    printf("Node %p at (%f, %f)\n", (void *) node, node->location[0],
           node->location[1]);

    printf("Split axis: %d\n", node->split);
    printf("Corners: (%f, %f)\t(%f, %f)\n", node->min[0], node->min[1],
           node->max[0], node->max[1]);
    printf("Children: %p\t%p\n", (void *) node->left, (void *) node->right);
    printf("Index: %u\n", node->index);

    printf("\n");

}

void
kd_printTree(struct kdNode *node)
{
    if (node == NULL)
        return;

    kd_printTree(node->left);
    if (kd_isleaf(node) )
        printf("%f\t%f\n", node->location[0], node->location[1]);
    kd_printTree(node->right);
}

/* End helper functions */

/* ******************************************************************
   
   Functions for building and destroying trees 

   ******************************************************************** */

static int
_compPoints(const void *p1, const void *p2, int axis)
{
    struct kd_point *a = (struct kd_point *) p1;
    struct kd_point *b = (struct kd_point *) p2;

    int ret = (a->point[axis] > b->point[axis]) ? 1 : (a->point[axis] < b->point[axis]) ? -1 : 0;
    if ( ret == 0 ) ret = (a->index > b->index) ? 1 : (a->index < b->index) ? -1 : 0;

    return ret;
}

static int
_compPoints0(const void *p1, const void *p2)
{
    struct kd_point *a = (struct kd_point *) p1;
    struct kd_point *b = (struct kd_point *) p2;

    int ret = (a->point[0] > b->point[0]) ? 1 : (a->point[0] < b->point[0]) ? -1 : 0;
    if ( ret == 0 ) ret = (a->index > b->index) ? 1 : (a->index < b->index) ? -1 : 0;

    return ret;
}

static int
_compPoints1(const void *p1, const void *p2)
{
    struct kd_point *a = (struct kd_point *) p1;
    struct kd_point *b = (struct kd_point *) p2;

    int ret = (a->point[1] > b->point[1]) ? 1 : (a->point[1] < b->point[1]) ? -1 : 0;
    if ( ret == 0 ) ret = (a->index > b->index) ? 1 : (a->index < b->index) ? -1 : 0;

    return ret;
}

static int
_compPoints2(const void *p1, const void *p2)
{
    struct kd_point *a = (struct kd_point *) p1;
    struct kd_point *b = (struct kd_point *) p2;

    int ret = (a->point[2] > b->point[2]) ? 1 : (a->point[2] < b->point[2]) ? -1 : 0;
    if ( ret == 0 ) ret = (a->index > b->index) ? 1 : (a->index < b->index) ? -1 : 0;

    return ret;
}

/*
void kd_quick_sort(struct kd_point *a, int n, int idx)
{
    int i, j;
    struct kd_point p, t;
    if ( n < 2 ) return;
    p = a[n / 2];
    for ( i = 0, j = n - 1;; i++, j-- )
      {
        while ( a[i].point[idx] < p.point[idx] ) i++;
        while ( p.point[idx] < a[j].point[idx] ) j--;
        if ( i >= j ) break;
        t = a[i];
        a[i] = a[j];
        a[j] = t;
      }
    kd_quick_sort(a, i, idx);
    kd_quick_sort(a + i, n - i, idx);
}
*/

void *
kd_doBuildTree(void *threadarg)
{
    int sortaxis;
    unsigned long pivot;
    kdata_t tmpMinLeft[KD_MAX_DIM], tmpMaxLeft[KD_MAX_DIM], tmpMinRight[KD_MAX_DIM], tmpMaxRight[KD_MAX_DIM];
    struct kdNode *node;
    struct kd_point *points;
    unsigned long nPoints;
    kdata_t *min, *max;
    int depth, dim, max_threads;
    pthread_t threads[2];
    pthread_attr_t attr;
    struct kd_thread_data *argleft, *argright;
    struct kd_thread_data *my_data = (struct kd_thread_data *) threadarg;

    points = my_data->points;
    nPoints = my_data->nPoints;
    min = my_data->min;
    max = my_data->max;
    depth = my_data->depth;
    max_threads = my_data->max_threads;
    dim = my_data->dim;

    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    sortaxis = depth % dim;

    if (nPoints == 1) {
        if ((node = kd_allocNode(points, 0, min, max, sortaxis, dim)) == NULL)
            return NULL;
        node->index = points[0].index;
        return node;
    }

    /*
     * If this iteration is allowed to start more threads, we first
     * use them to parallelize the sorting 
     */
    int (*qcomp) (const void *, const void *);
    if      ( sortaxis == 0 ) qcomp = _compPoints0;
    else if ( sortaxis == 1 ) qcomp = _compPoints1;
    else                      qcomp = _compPoints2;

#if defined(KDTEST)
    pmergesort(points, nPoints, sizeof(struct kd_point), _compPoints, sortaxis, max_threads);
#else
    pmergesort(points, nPoints, sizeof(struct kd_point), qcomp, max_threads);
#endif

    pivot = nPoints / 2;
    if ((node = kd_allocNode(points, pivot, min, max, sortaxis, dim)) == NULL)
        return NULL;

    memcpy(tmpMinLeft, min, dim * sizeof(kdata_t));
    memcpy(tmpMaxLeft, max, dim * sizeof(kdata_t));
    tmpMaxLeft[sortaxis] = node->location[sortaxis];
    argleft =
        kd_buildArg(points, pivot, tmpMinLeft, tmpMaxLeft,
                    depth + 1, max_threads / 2, dim);
    if (!argleft) {
        kd_destroyTree(node);
        return NULL;
    }
    if (max_threads > 1) {
        pthread_create(&threads[0], &attr, kd_doBuildTree, (void *) argleft);
    } else {
        node->left = (kdNode *)kd_doBuildTree((void *) argleft);
        free(argleft);
        if (!node->left) {
            kd_destroyTree(node);
            return NULL;
        }
    }

    memcpy(tmpMinRight, min, dim * sizeof(kdata_t));
    memcpy(tmpMaxRight, max, dim * sizeof(kdata_t));
    tmpMinRight[sortaxis] = node->location[sortaxis];
    argright = kd_buildArg(&points[pivot], nPoints - pivot,
                           tmpMinRight, tmpMaxRight, depth + 1,
                           max_threads / 2, dim);
    if (!argright) {
        kd_destroyTree(node);
        return NULL;
    }
    if (max_threads > 1) {
        pthread_create(&threads[1], &attr, kd_doBuildTree, (void *) argright);
    } else {
        node->right = (kdNode *)kd_doBuildTree((void *) argright);
        free(argright);
        if (!node->right) {
            kd_destroyTree(node);
            return NULL;
        }
    }
    if (max_threads > 1) {
        pthread_join(threads[0], (void **) (&node->left));
        free(argleft);
        pthread_join(threads[1], (void **) (&node->right));
        free(argright);
        if (!node->left || !node->right) {
            kd_destroyTree(node);
            return NULL;
        }
    }
    return (void *) node;
}


void
kd_freeNode(kdNode *node)
{
    if ( node ) free(node);
    return;
}


struct kd_thread_data *
kd_buildArg(struct kd_point *points,
            unsigned long nPoints,
            kdata_t *min, kdata_t *max,
            int depth, int max_threads, int dim)
{
    struct kd_thread_data *d;

    if ((d = (kd_thread_data *)kd_malloc(sizeof(kd_thread_data), "kd_thread_data")) == NULL)
        return NULL;
    d->points = points;
    d->nPoints = nPoints;
    memcpy(d->min, min, dim*sizeof(kdata_t));
    memcpy(d->max, max, dim*sizeof(kdata_t));
    d->depth = depth;
    d->max_threads = max_threads;
    d->dim = dim;
    return d;
}


struct kdNode *
kd_allocNode(struct kd_point *points, unsigned long pivot,
             kdata_t *min, kdata_t *max, int axis, int dim)
{
    struct kdNode *node;

    if ((node = (kdNode *)kd_malloc(sizeof(kdNode), "kd_allocNode (node): ")) == NULL)
        return NULL;
    node->split = axis;
    memcpy(node->location, points[pivot].point, dim * sizeof(kdata_t));
    memcpy(node->min, min, dim * sizeof(kdata_t));
    memcpy(node->max, max, dim * sizeof(kdata_t));
    node->left = node->right = NULL;
    node->index = 0;
    return node;
}

/*!
 * \brief free the kd-tree data structure, 
 *
 * \param node the root node of the tree to be destroyed
 *
 * \param *destr a pointer to the destructor function for the data container.
 *
 * \return This function does not return a value
 */
void
kd_destroyTree(struct kdNode *node)
{
    if (node == NULL)
        return;

    kd_destroyTree(node->left);
    kd_destroyTree(node->right);
    kd_freeNode(node);
}

/* end of tree construction and destruction */

/* Functions dealing with result heaps */

/* Insert the sub-tree starting at node into the result heap res */
int
kd_insertResTree(struct kdNode *node, struct pqueue *res)
{
    struct resItem *point;

    if (node == NULL)
        return 1;

    if (!kd_insertResTree(node->left, res))
        return 0;
    if (kd_isleaf(node)) {
        if ((point = (struct resItem *)kd_malloc(sizeof(struct resItem), "kd_insertResTree: "))
            == NULL)
            return 0;

        point->node = node;
        point->dist_sq = -1;
        pqinsert(res, point);
    }
    if (!kd_insertResTree(node->right, res))
        return 0;
    return 1;
}
