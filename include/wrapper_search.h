#ifndef WRAPPER_SEARCH_H
#define WRAPPER_SEARCH_H

static void * pb_global;

int
cerca_atomo_wrapper (p8est_t * p4est,
                     p4est_topidx_t which_tree,
                     p8est_quadrant_t * quadrant,
                     p4est_locidx_t local_num,
                     void *point);

#endif