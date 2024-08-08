/*
#######################################################################

Author: Thang V. Pham, t.pham@amsterdamumc.nl

All rights reserved.

Citation for the beta-binomial test:

T.V. Pham, S.R. Piersma, M. Warmoes, C.R. Jimenez (2010)
On the beta binomial model for analysis of spectral count data
in label-free tandem mass spectrometry-based proteomics.
Bioinformatics, 26(3):363-369.

Software version: 1.3

#######################################################################
*/

#ifdef _OPENMP
#include <omp.h>
#endif

#include <R.h>
#include <Rmath.h>

#define SMAX (TYPE)1e6
#define SMIN (TYPE)1e-6

#define MAX_INNER 200

#define ALPHA  (TYPE)0.1

#define BETA (TYPE)0.7

#define SMALL (TYPE)1e-12

#define forint(i, a, b) for (int i=(int)(a); i<(int)(b); i++)

#define TYPE double

Rboolean R_ToplevelExec(void (*fun)(void*), void* data);

typedef struct _bb_t {

    TYPE* a;
    TYPE* ta;

    int N; // total size
    int M; // no. of groups
    int* g_size;
    int* g_ind;

    int theta_equal;

    // for optimization
    TYPE _m1;
    TYPE _size;
    TYPE* _a;
    TYPE* _ta;

    TYPE* _m1_array;

    int comp; // for one-sided test

    TYPE f;
    TYPE f0;
} bb_t;

/* One dimensional optimization */
TYPE opt_find_eta(TYPE* v, TYPE g, TYPE h, TYPE p, TYPE lb, TYPE ub) {

    if (fabs(h) < SMALL) {
        h = (h > (TYPE)0 ? SMALL : -SMALL);
    }

    *v = -g / h;		// Newton step

    // if h < 0, switch direction. Make sure v is a descent direction
    if (g * (*v) > 0) *v = -(*v);

    TYPE eta = (TYPE)1.0;

    if (p + *v >= ub) eta = (ub - SMALL - p) / *v;
    if (p + *v <= lb) eta = (lb + SMALL - p) / *v;

    if (eta < (TYPE)0.0) eta = (TYPE)0.0;

    return eta;
}


TYPE minimize1d(TYPE(*fval)(TYPE*, TYPE*, bb_t*, TYPE),
    bb_t* x,
    TYPE p0,
    TYPE lb,
    TYPE ub) {

    TYPE p = p0;

    TYPE old_p = p;

    int cc = 0;

    forint(i, 0, MAX_INNER) {

        TYPE f, g, h;

        f = fval(&g, &h, x, p);

        TYPE v;

        TYPE eta = opt_find_eta(&v, g, h, p, lb, ub);

        TYPE lambda = g * v;	// Newton decrement

        TYPE fnew = fval(0, 0, x, p + eta * v);

        while ((fnew - f) > ALPHA* lambda* eta) {

            eta *= BETA;

            fnew = fval(0, 0, x, p + eta * v);

            if (eta < (TYPE)1e-15) {
                eta = (TYPE)0.0;
                fnew = f;
                break;
            };
        }

        p += eta * v;

        if (fabs(fnew - f) < (TYPE)1e-8 && fabs(p - old_p) < (TYPE)1e-8) {
            cc++;
        }
        else {
            cc = 0;
        }

        if (cc > 1) {
            break;
        }

        old_p = p;
    }
    return p;
}


/* method of moments */
void bb_simple_estimate_tm(bb_t* x,
    TYPE* alp,
    TYPE* bet) {

    TYPE p = (TYPE)0.0;
    TYPE m2 = (TYPE)0.0;

    forint(i, 0, x->_size) {

        TYPE tmp = x->_a[i] / x->_ta[i];

        p += tmp;

        m2 += tmp * tmp;

    }

    p /= (TYPE)(x->_size);
    m2 /= (TYPE)(x->_size);

    if (p < SMALL) {
        *alp = (TYPE)1.0;
        *bet = (TYPE)1e4;
    }
    else {

        TYPE s;

        if ((m2 - p * p) < SMALL) {
            s = (TYPE)1e4;
        }
        else {
            s = (p - m2) / (m2 - p * p);

            if (s < SMIN) {
                s = SMIN;
            }
            if (s > SMAX) {
                s = SMAX;
            }
        }

        *alp = p * s;
        *bet = ((TYPE)1.0 - p) * s;
    }
}

/* evaluate as a function of s */
TYPE fval_s(TYPE* g, TYPE* h, bb_t* x, TYPE s) {

    TYPE m2 = ((TYPE)1.0 - x->_m1);

    if (g) {

        TYPE f = (TYPE)0.0;
        *g = (TYPE)0.0;
        *h = (TYPE)0.0;

        TYPE lg_a0 = lgammafn(s);
        TYPE di_a0 = digamma(s);
        TYPE tri_a0 = trigamma(s);

        TYPE a1 = s * x->_m1;
        TYPE a2 = s * m2;

        TYPE lg_a1 = lgammafn(a1);
        TYPE di_a1 = x->_m1 * digamma(a1);
        TYPE tri_a1 = x->_m1 * x->_m1 * trigamma(a1);

        TYPE lg_a2 = lgammafn(a2);
        TYPE di_a2 = m2 * digamma(a2);
        TYPE tri_a2 = m2 * m2 * trigamma(a2);

        forint(i, 0, x->_size) {

            TYPE term2 = x->_ta[i] + s;

            TYPE lg_term2 = lgammafn(term2);
            TYPE di_term2 = digamma(term2);
            TYPE tri_term2 = trigamma(term2);

            TYPE k1 = x->_a[i] + a1;
            TYPE lg_k1 = lgammafn(k1);
            TYPE di_k1 = x->_m1 * digamma(k1);
            TYPE tri_k1 = x->_m1 * x->_m1 * digamma(k1);

            TYPE k2 = x->_ta[i] - x->_a[i] + a2;
            TYPE lg_k2 = lgammafn(k2);
            TYPE di_k2 = m2 * digamma(k2);
            TYPE tri_k2 = m2 * m2 * trigamma(k2);

            f -= (lg_a0 - lg_term2 + lg_k1 - lg_a1 + lg_k2 - lg_a2);

            *g -= (di_a0 - di_term2 + di_k1 - di_a1 + di_k2 - di_a2);

            *h -= (tri_a0 - tri_term2 + tri_k1 - tri_a1 + tri_k2 - tri_a2);
        }
        return(f);
    }
    else {

        TYPE a1 = s * x->_m1;
        TYPE a2 = s * m2;

        TYPE f = -(lgammafn(s) - lgammafn(a1) - lgammafn(a2)) * x->_size;

        forint(i, 0, x->_size) {
            f -= (-lgammafn(x->_ta[i] + s)
                + lgammafn(x->_a[i] + a1)
                + lgammafn(x->_ta[i] - x->_a[i] + a2));
        }
        return(f);
    }
}

/* as a function of 1/s */
TYPE fval_s_inv(TYPE* g, TYPE* h, bb_t* x, TYPE s_inv) {

    TYPE m2 = ((TYPE)1.0 - x->_m1);

    TYPE s = (TYPE)1.0 / s_inv;

    if (g) {

        TYPE f = (TYPE)0.0;
        *g = (TYPE)0.0;
        *h = (TYPE)0.0;

        TYPE lg_a0 = lgammafn(s);
        TYPE di_a0 = digamma(s);
        TYPE tri_a0 = trigamma(s);

        TYPE a1 = s * x->_m1;
        TYPE a2 = s * m2;

        TYPE lg_a1 = lgammafn(a1);
        TYPE di_a1 = x->_m1 * digamma(a1);
        TYPE tri_a1 = x->_m1 * x->_m1 * trigamma(a1);

        TYPE lg_a2 = lgammafn(a2);
        TYPE di_a2 = m2 * digamma(a2);
        TYPE tri_a2 = m2 * m2 * trigamma(a2);

        forint(i, 0, x->_size) {

            TYPE term2 = x->_ta[i] + s;

            TYPE lg_term2 = lgammafn(term2);
            TYPE di_term2 = digamma(term2);
            TYPE tri_term2 = trigamma(term2);

            TYPE k1 = x->_a[i] + a1;
            TYPE lg_k1 = lgammafn(k1);
            TYPE di_k1 = x->_m1 * digamma(k1);
            TYPE tri_k1 = x->_m1 * x->_m1 * digamma(k1);

            TYPE k2 = x->_ta[i] - x->_a[i] + a2;
            TYPE lg_k2 = lgammafn(k2);
            TYPE di_k2 = m2 * digamma(k2);
            TYPE tri_k2 = m2 * m2 * trigamma(k2);

            f -= (lg_a0 - lg_term2 + lg_k1 - lg_a1 + lg_k2 - lg_a2);

            TYPE tmp_g = (di_a0 - di_term2 + di_k1 - di_a1 + di_k2 - di_a2);

            TYPE s_inv_sqr = s_inv * s_inv;

            *g += s_inv_sqr * tmp_g;

            *h -= ((TYPE)2.0 * s_inv_sqr * s_inv * tmp_g
                + s_inv_sqr * s_inv_sqr *
                (tri_a0 - tri_term2 + tri_k1 - tri_a1 + tri_k2 - tri_a2));
        }
        return(f);
    }
    else {

        TYPE a1 = s * x->_m1;
        TYPE a2 = s * m2;

        TYPE f = -(lgammafn(s) - lgammafn(a1) - lgammafn(a2)) * x->_size;

        forint(i, 0, x->_size) {
            f -= (-lgammafn(x->_ta[i] + s)
                + lgammafn(x->_a[i] + a1)
                + lgammafn(x->_ta[i] - x->_a[i] + a2));
        }
        return(f);
    }
}

TYPE fval_s_equal_inv(TYPE* dx, TYPE* dxx, bb_t* x, TYPE s_inv) {

    if (dx) {

        TYPE f = (TYPE)0.0;
        *dx = (TYPE)0.0;
        *dxx = (TYPE)0.0;
        TYPE _dx, _dxx;

        forint(g, 0, x->M) {

            x->_size = x->g_size[g];
            x->_a = x->a + x->g_ind[g];
            x->_ta = x->ta + x->g_ind[g];

            x->_m1 = x->_m1_array[g];

            f += fval_s_inv(&_dx, &_dxx, x, s_inv);
            *dx += _dx;
            *dxx += _dxx;
        }

        return(f);

    }

    else {

        TYPE f = (TYPE)0.0;

        forint(g, 0, x->M) {

            x->_size = x->g_size[g];
            x->_a = x->a + x->g_ind[g];
            x->_ta = x->ta + x->g_ind[g];

            x->_m1 = x->_m1_array[g];

            f += fval_s_inv(0, 0, x, s_inv);
        }
        return(f);
    }
}

/* global optimization with DC programming? */
void fit_m(bb_t* x, const TYPE s) {

    forint(i, 0, 200) {

        TYPE m2 = ((TYPE)1.0 - x->_m1);

        TYPE alp = s * x->_m1;
        TYPE bet = s * m2;

        TYPE new_m1 = -digamma(alp) * x->_size;
        TYPE new_m2 = -digamma(bet) * x->_size;

        forint(j, 0, x->_size) {
            new_m1 += digamma(alp + x->_a[j]);
            new_m2 += digamma(bet + x->_ta[j] - x->_a[j]);
        }

        new_m1 *= alp;
        new_m2 *= bet;

        TYPE tmp = new_m1 + new_m2;

        new_m1 /= tmp;

        if (fabs(x->_m1 - new_m1) < (TYPE)1e-8) {

            x->_m1 = new_m1;

            if (x->_m1 < SMALL) {
                x->_m1 = SMALL;
            }

            if ((x->_m1 + SMALL) > (TYPE)1.0) {
                x->_m1 = (TYPE)1.0 - SMALL;
            }

            break;
        }
        else {

            x->_m1 = new_m1;

            if (x->_m1 < SMALL) {
                x->_m1 = SMALL;
            }

            if ((x->_m1 + SMALL) > (TYPE)1.0) {
                x->_m1 = (TYPE)1.0 - SMALL;
            }
        }
    }
}

/* return MLE */
TYPE bbmle(bb_t* x, int g, TYPE alp0, TYPE bet0, TYPE* alp, TYPE* bet) {

    // optimize for group g
    if (g >= 0) {
        x->_size = x->g_size[g];
        x->_a = x->a + x->g_ind[g];
        x->_ta = x->ta + x->g_ind[g];
    }
    else {
        x->_size = x->N;
        x->_a = x->a;
        x->_ta = x->ta;
    }


    TYPE s_inv = (TYPE)1.0 / (alp0 + bet0);

    x->_m1 = alp0 * s_inv;

    forint(i, 0, 5000) {

        TYPE old_m1 = x->_m1;

        fit_m(x, (TYPE)1.0 / s_inv);

        /* update s */
        TYPE old_s_inv = s_inv;
        TYPE old_f = fval_s_inv(0, 0, x, old_s_inv);

        s_inv = minimize1d(fval_s_inv, x, old_s_inv, SMIN, SMAX);

        TYPE new_f = fval_s_inv(0, 0, x, s_inv);


        if (fabs(s_inv - old_s_inv) < 1e-12
            && fabs(x->_m1 - old_m1) < 1e-12
            && fabs(old_f - new_f) < 1e-12) {
            break;
        }
    }

    *alp = x->_m1 / s_inv;
    *bet = ((TYPE)1.0 - x->_m1) / s_inv;

    return(-fval_s_inv(0, 0, x, s_inv));


}


TYPE fval_s_equal(TYPE* dx, TYPE* dxx, bb_t* x, TYPE s) {

    if (dx) {

        TYPE f = (TYPE)0.0;
        *dx = (TYPE)0.0;
        *dxx = (TYPE)0.0;
        TYPE _dx, _dxx;

        forint(g, 0, x->M) {

            x->_size = x->g_size[g];
            x->_a = x->a + x->g_ind[g];
            x->_ta = x->ta + x->g_ind[g];

            x->_m1 = x->_m1_array[g];

            f += fval_s(&_dx, &_dxx, x, s);
            *dx += _dx;
            *dxx += _dxx;
        }
        return(f);
    }

    else {

        TYPE f = (TYPE)0.0;

        forint(g, 0, x->M) {

            x->_size = x->g_size[g];
            x->_a = x->a + x->g_ind[g];
            x->_ta = x->ta + x->g_ind[g];

            x->_m1 = x->_m1_array[g];

            f += fval_s(0, 0, x, s);
        }
        return(f);
    }
}

/* return MLE equal theta = equal s*/
TYPE bbmle_equal(bb_t* x, TYPE alp, TYPE bet) {

    TYPE alp0;
    TYPE bet0;

    TYPE s;

    if (alp < (TYPE)0.0) { // not initialize

        s = (TYPE)0.0;

        forint(g, 0, x->M) {

            x->_size = x->g_size[g];
            x->_a = x->a + x->g_ind[g];
            x->_ta = x->ta + x->g_ind[g];

            bb_simple_estimate_tm(x, &alp0, &bet0);

            TYPE tmp_s = alp0 + bet0;
            x->_m1_array[g] = alp0 / tmp_s;

            s += tmp_s;

        }

        s /= (TYPE)x->M;

    }
    else {

        alp0 = alp;
        bet0 = bet;

        s = alp0 + bet0;

        forint(g, 0, x->M) {
            x->_m1_array[g] = alp0 / s;
        }

    }


    TYPE s_inv = (TYPE)1.0 / s;

    forint(i, 0, 5000) {

        /* update _m1 */
        TYPE max_diff = (TYPE)0.0;

        forint(g, 0, x->M) {

            x->_size = x->g_size[g];
            x->_a = x->a + x->g_ind[g];
            x->_ta = x->ta + x->g_ind[g];

            x->_m1 = x->_m1_array[g];

            TYPE old_m1 = x->_m1;

            fit_m(x, (TYPE)1.0 / s_inv);

            /* store the values */
            x->_m1_array[g] = x->_m1;

            /* check the difference */
            if (max_diff < fabs(x->_m1 - old_m1)) {
                max_diff = fabs(x->_m1 - old_m1);
            }

        }

        /* update s */
        TYPE old_s_inv = s_inv;

        s_inv = minimize1d(fval_s_equal_inv, x, old_s_inv, SMIN, SMAX);


        if (fabs(s_inv - old_s_inv) < (TYPE)1e-12 && max_diff < (TYPE)1e-12) {
            break;
        }
    }

    return(-fval_s_equal_inv(0, 0, x, s_inv));


}

void do_bb_test(bb_t* x) {

    TYPE f0, f;

    TYPE alp;
    TYPE bet;

    x->_size = x->N;
    x->_a = x->a;
    x->_ta = x->ta;

    //bb_moment_estimate(x, &alp, &bet);
    bb_simple_estimate_tm(x, &alp, &bet);

    TYPE alp0 = alp;
    TYPE bet0 = bet;

    f0 = bbmle(x, -1, alp0, bet0, &alp, &bet);

    alp0 = alp;
    bet0 = bet;

    if (x->theta_equal > 0) {

        TYPE f_init1 = bbmle_equal(x, alp0, bet0);

        int comp = x->_m1_array[0] > x->_m1_array[1];

        TYPE f_init2 = bbmle_equal(x, (TYPE)-1.0, (TYPE)-1.0);

        if (f_init1 > f_init2) {
            f = f_init1;
            x->comp = comp;
        }
        else {
            f = f_init2;
            x->comp = x->_m1_array[0] > x->_m1_array[1];
        }
    }

    else {
        f = (TYPE)0.0;

        TYPE v0 = 0;
        TYPE v1 = 0;

        forint(g, 0, x->M) {

            TYPE tmp = bbmle(x, g, alp0, bet0, &alp, &bet);

            TYPE vv = x->_m1;

            // try  a different initialization
            x->_size = x->g_size[g];
            x->_a = x->a + x->g_ind[g];
            x->_ta = x->ta + x->g_ind[g];

            bb_simple_estimate_tm(x, &alp, &bet);

            TYPE alp_tmp, bet_tmp;

            TYPE tmp2 = bbmle(x, g, alp, bet, &alp_tmp, &bet_tmp);

            if (tmp > tmp2) {
                f += tmp;
            }
            else {
                f += tmp2;
                vv = x->_m1;
            }
            if (g == 0) {
                v0 = vv;
            }
            if (g == 1) {
                v1 = vv;
            }
        }

        x->comp = v0 > v1;
    }

    x->f0 = f0;
    x->f = f;
}

void bbCores(int* n_procs) {
    #ifdef _OPENMP
    * n_procs = omp_get_num_procs();
    #endif
}

static void tp_user(void* dummy) {
    R_CheckUserInterrupt();
}

// won't longjmp-out of your context
int tp_check(void) {
    return (R_ToplevelExec(tp_user, NULL) == FALSE);
}

void bb(int* lK,
    double* a,
    double* ta,
    int* lM,
    int* g_size,
    int* g_ind,
    double* mem,
    int* no_threads,
    double* pval) {

    int verbose = *no_threads;

    if (*no_threads < 0) {
        *no_threads = -*no_threads;
    }

    int theta_equal = ((TYPE)mem[0] > 0 ? 1 : 0);

    TYPE tail = (TYPE)(mem[1]);

    int N = 0;
    forint(i, 0, *lM) {
        N += g_size[i];
    }

    int init_block = 0; // global variables

    int block = 2 * N + (*lM);

    TYPE* works = (TYPE*)mem;

    int thres_display = 0;

    // initial block

    // openmp
    #ifdef _OPENMP
    omp_set_dynamic(0);

    omp_set_num_threads(*no_threads);

    if (*no_threads > 1) {
        if (verbose > 0) Rprintf("Using %d threads ...\nNo. of data rows = %d, no. of groups = %d, no. of samples = %d.\n", *no_threads, *lK, *lM, N);
    }
    else {
        if (verbose > 0) Rprintf("Using a single CPU core ...\nNo. of data rows = %d, no. of groups = %d, no. of samples = %d.\n", *lK, *lM, N);
    }
    #else
    if (verbose > 0) Rprintf("Using a single CPU core ...\nNo. of data rows = %d, no. of groups = %d, no. of samples = %d.\n", *lK, *lM, N);
    #endif

    int k = 0;

    int stop_sig = 0;

    #pragma omp parallel for schedule(dynamic)
    forint(i, 0, *lK) {

        if (stop_sig) continue;

        int thread_id = 0;

        #ifdef _OPENMP
        thread_id = omp_get_thread_num();
        #endif

        bb_t x;

        x.N = N;
        x.theta_equal = theta_equal;

        x.M = *lM;
        x.g_size = g_size;
        x.g_ind = g_ind;

        int mem_start = init_block + thread_id * block;

        x.a = works + mem_start;
        x.ta = works + mem_start + 1 * N;

        x._m1_array = works + mem_start + 2 * N;

        /*
          end pointer = works + mem_start + 2 * N + x.M;
        */

        int ind = i * N;

        forint(j, 0, N) {
            x.a[j] = (TYPE)a[ind];
            x.ta[j] = (TYPE)ta[ind];
            ind++;
        }

        do_bb_test(&x);


        double g = 2.0 * ((double)x.f - (double)x.f0);

        if (tail > 0.5) {
            if (!x.comp) {
                pval[i] = pnorm(sqrt(g), 0.0, 1.0, 0, 0);
            }
            else {
                pval[i] = pnorm(-sqrt(g), 0.0, 1.0, 0, 0);
            }
        }
        else {
            if (tail < -0.5) {
                if (!x.comp) {
                    pval[i] = pnorm(sqrt(g), 0.0, 1.0, 1, 0);
                }
                else {
                    pval[i] = pnorm(-sqrt(g), 0.0, 1.0, 1, 0);
                }
            }
            else {
                pval[i] = pchisq(g, (double)*lM - 1.0, 0, 0);
            }
        }

        #pragma omp atomic
        k++;

        if (thread_id == 0) {
            if (verbose > 0 && k > thres_display) {

                Rprintf("%d%%\n", k * 100 / (*lK));
                R_FlushConsole();
                thres_display = k + (*lK) / 20;
            }

            if (tp_check()) { // user interrupted ...
                stop_sig = 1;
                #pragma omp flush(stop_sig)
            }
        }
    }

    if (!stop_sig) {
        if (verbose > 0) Rprintf("Done.\n");
    }

}
