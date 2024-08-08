/*
#######################################################################

Author: Thang V. Pham, t.pham@amsterdamumc.nl

All rights reserved.

Citation for the beta-binomial test:

T.V. Pham, C.R. Jimenez (2012)
An accurate paired sample test for count data.
Bioinformatics, 28(18):i596-i602.

Software version: 1.3

#######################################################################
*/

#ifdef _OPENMP
#include <omp.h>
#endif

#include <R.h>
#include <Rmath.h>

/* ln(2) */
#define T_LN2 0.693147180559945309417232121458


#define ZERO_REG (TYPE)0.5

#define MAX_INNER 2000

#define SMALL 1e-12

#define ALPHA  0.1

#define BETA 0.7

#define SMAX (TYPE)1e6
#define SMIN (TYPE)1e-6

#define forint(i, a, b) for (int i=(int)(a); i<(int)(b); i++)

#define TYPE double

LibExtern int R_interrupts_pending;

#include "q15.h"

typedef struct _data_t {

    TYPE* a;
    TYPE* b;
    TYPE* ta;
    TYPE* tb;
    int N;

    const TYPE* qz;
    const TYPE* qlw;
    int Z;

    TYPE* ll;
    TYPE* log_one_plus_z;
    TYPE* log_one_minus_z;

    TYPE* works;

    TYPE f0;
    TYPE f;
    TYPE alp;
    TYPE bet;

} data_t;



int tp_check(void);

/* fval2 */
TYPE fval2(TYPE* ga, TYPE* gb, TYPE* haa, TYPE* hab, TYPE* hbb, data_t* x, TYPE alp, TYPE bet) {

    TYPE* y = x->works;

    TYPE* w = x->works + x->Z; // initial = tmp

    TYPE* da_y = x->works + 2 * x->Z;

    TYPE* db_y = x->works + 3 * x->Z;

    TYPE* tmp = x->works + 4 * x->Z;

    TYPE a1 = alp - (TYPE)1.0;
    TYPE b1 = bet - (TYPE)1.0;
    //TYPE c = -(alp + bet -(TYPE)1.0) * log((TYPE)2.0) + lgammafn(alp+bet) - lgammafn(alp) - lgammafn(bet);
    TYPE c = -(alp + bet - (TYPE)1.0) * T_LN2 + lgammafn(alp + bet) - lgammafn(alp) - lgammafn(bet);

    forint(i, 0, x->Z) {
        tmp[i] = x->qlw[i] + a1 * x->log_one_plus_z[i] + b1 * x->log_one_minus_z[i] + c;
    }

    TYPE digamma_ab;
    TYPE digamma_a;
    TYPE digamma_b;

    TYPE trigamma_ab;
    TYPE trigamma_a;
    TYPE trigamma_b;

    TYPE f = (TYPE)0.0;

    if (ga) {

        *ga = (TYPE)0.0;
        *gb = (TYPE)0.0;

        *haa = (TYPE)0.0;
        *hab = (TYPE)0.0;
        *hbb = (TYPE)0.0;

        digamma_ab = digamma(alp + bet);
        digamma_a = digamma(alp);
        digamma_b = digamma(bet);

        trigamma_ab = trigamma(alp + bet);
        trigamma_a = trigamma(alp);
        trigamma_b = trigamma(bet);
    }


    TYPE* ptr = x->ll;

    forint(n, 0, x->N) {

        TYPE max_y = tmp[0] + ptr[0];

        forint(i, 0, x->Z) {

            y[i] = tmp[i] + *ptr++;

            if (max_y < y[i]) max_y = y[i];

        }

        TYPE se = (TYPE)0.0;

        forint(i, 0, x->Z) {
            se += exp(y[i] - max_y);
        }

        TYPE fx = log(se) + max_y;

        f = f - fx;

        if (ga) {

            TYPE da_f = (TYPE)0.0;
            TYPE db_f = (TYPE)0.0;

            forint(i, 0, x->Z) {

                w[i] = exp(y[i] - fx);

                //da_y[i] = x->log_one_plus_z[i] - log(2.0) + digamma_ab - digamma_a;
                da_y[i] = x->log_one_plus_z[i] - T_LN2 + digamma_ab - digamma_a;

                //db_y[i] = x->log_one_minus_z[i] - log(2.0) + digamma_ab - digamma_b;
                db_y[i] = x->log_one_minus_z[i] - T_LN2 + digamma_ab - digamma_b;

                da_f += w[i] * da_y[i];

                db_f += w[i] * db_y[i];

            }

            TYPE daa_f = trigamma_ab - trigamma_a;
            TYPE dab_f = trigamma_ab;
            TYPE dbb_f = trigamma_ab - trigamma_b;

            forint(i, 0, x->Z) {

                TYPE da_w = w[i] * (da_y[i] - da_f);

                TYPE db_w = w[i] * (db_y[i] - db_f);

                daa_f += da_w * da_y[i];

                dab_f += da_w * db_y[i];

                dbb_f += db_w * db_y[i];

            }

            *ga -= da_f;
            *gb -= db_f;

            *haa -= daa_f;
            *hab -= dab_f;
            *hbb -= dbb_f;

        }

    }

    return f;

}



/* fval */

TYPE fval(TYPE* g, TYPE* h, data_t* x, TYPE alp, TYPE bet, int derivative) {


    TYPE* y = x->works;

    TYPE* w = x->works + x->Z; // initial = tmp

    TYPE* d_y = x->works + 2 * x->Z;

    TYPE* tmp = x->works + 3 * x->Z;


    TYPE a1 = alp - (TYPE)1.0;
    TYPE b1 = bet - (TYPE)1.0;
    //TYPE c = -(alp + bet -(TYPE)1.0) * log((TYPE)2.0) + lgammafn(alp+bet) - lgammafn(alp) - lgammafn(bet);
    TYPE c = -(alp + bet - (TYPE)1.0) * T_LN2 + lgammafn(alp + bet) - lgammafn(alp) - lgammafn(bet);

    forint(i, 0, x->Z) {
        tmp[i] = x->qlw[i] + a1 * x->log_one_plus_z[i] + b1 * x->log_one_minus_z[i] + c;
    }

    TYPE digamma_ab;
    TYPE digamma_a;
    TYPE digamma_b;

    TYPE trigamma_ab;
    TYPE trigamma_a;
    TYPE trigamma_b;


    TYPE f = (TYPE)0.0;


    if (g) {

        *g = (TYPE)0.0;
        *h = (TYPE)0.0;

        digamma_ab = digamma(alp + bet);
        digamma_a = digamma(alp);
        digamma_b = digamma(bet);

        trigamma_ab = trigamma(alp + bet);
        trigamma_a = trigamma(alp);
        trigamma_b = trigamma(bet);
    }


    TYPE* ptr = x->ll;

    forint(n, 0, x->N) {

        TYPE max_y = tmp[0] + ptr[0];

        forint(i, 0, x->Z) {

            y[i] = tmp[i] + *ptr++;

            if (max_y < y[i]) max_y = y[i];

        }

        TYPE se = (TYPE)0.0;

        forint(i, 0, x->Z) {
            se += exp(y[i] - max_y);
        }

        TYPE fx = log(se) + max_y;

        f = f - fx;

        if (g) {

            if (derivative == 0) {

                TYPE df = (TYPE)0.0;

                for (int i = 0; i < x->Z; i++) {

                    w[i] = exp(y[i] - fx);

                    //d_y[i] = x->log_one_plus_z[i] - log((TYPE)2.0) + digamma_ab - digamma_a;
                    d_y[i] = x->log_one_plus_z[i] - T_LN2 + digamma_ab - digamma_a;

                    df += w[i] * d_y[i];
                }

                TYPE ddf = trigamma_ab - trigamma_a;

                for (int i = 0; i < x->Z; i++) {

                    TYPE d_w = w[i] * (d_y[i] - df);

                    ddf += d_w * d_y[i];
                }

                *g = *g - df;
                *h = *h - ddf;
            }

            else if (derivative == 1) {

                TYPE df = (TYPE)0.0;

                for (int i = 0; i < x->Z; i++) {

                    w[i] = exp(y[i] - fx);

                    //d_y[i] = x->log_one_minus_z[i] - log((TYPE)2.0) + digamma_ab - digamma_b;
                    d_y[i] = x->log_one_minus_z[i] - T_LN2 + digamma_ab - digamma_b;

                    df += w[i] * d_y[i];
                }

                TYPE ddf = trigamma_ab - trigamma_b;

                for (int i = 0; i < x->Z; i++) {

                    TYPE d_w = w[i] * (d_y[i] - df);

                    ddf += d_w * d_y[i];
                }

                *g = *g - df;
                *h = *h - ddf;

            }
            else if (derivative == 2) {

                TYPE df = (TYPE)0.0;

                for (int i = 0; i < x->Z; i++) {

                    w[i] = exp(y[i] - fx);

                    d_y[i] = x->log_one_plus_z[i]
                        + x->log_one_minus_z[i]
                        //- (TYPE)2.0 * log((TYPE)2.0)
                        - (TYPE)2.0 * T_LN2
                        + (TYPE)2.0 * digamma_ab
                        - digamma_a - digamma_b;

                    df += w[i] * d_y[i];
                }

                TYPE ddf = 4 * trigamma_ab - trigamma_a - trigamma_b;

                for (int i = 0; i < x->Z; i++) {

                    TYPE d_w = w[i] * (d_y[i] - df);

                    ddf += d_w * d_y[i];
                }

                *g = *g - df;
                *h = *h - ddf;
            }
        }
    }

    return f;

}


TYPE fval_ab(TYPE* g, TYPE* h, data_t* x, TYPE p) {

    if (g) return fval(g, h, x, p, p, 2);
    else return fval(g, h, x, p, p, -1);

}


TYPE find_eta(TYPE* v, TYPE g, TYPE h, TYPE p, TYPE lb, TYPE ub) {

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


void update_ab(TYPE* new_f,
    data_t* x,
    TYPE* a, TYPE* b,
    TYPE f, TYPE fcond,
    TYPE ga, TYPE gb,
    TYPE eta,
    TYPE va, TYPE vb,
    TYPE ua, TYPE la,
    TYPE ub, TYPE lb,
    TYPE ab_bound) {

    TYPE lambda = ga * va + gb * vb;	// Newton decrement, or gtp

    *new_f = fval2(0, 0, 0, 0, 0, x, *a + eta * va, *b + eta * vb) * fcond;

    while ((*new_f - f) > ALPHA* lambda* eta) {

        eta *= BETA;

        *new_f = fval2(0, 0, 0, 0, 0, x, *a + eta * va, *b + eta * vb) * fcond;

        if (eta < (TYPE)1e-15) {
            eta = (TYPE)0.0;
            *new_f = f;
            break;
        };
    }

    *a += eta * va;
    *b += eta * vb;

}

// optimization with linear bound
void nr2b_projection(TYPE(*fval2)(TYPE*, TYPE*, TYPE*, TYPE*, TYPE*, data_t*, TYPE, TYPE),
    data_t* x,
    TYPE* a,
    TYPE la,
    TYPE ua,
    TYPE* b,
    TYPE lb,
    TYPE ub,
    TYPE ab_bound, TYPE fcond) {

    TYPE old_a = *a;
    TYPE old_b = *b;

    int cc = 0;

    int slide_ok = 1;
    int newton_ok = 1;

    forint(i, 0, MAX_INNER) {

        TYPE ga;
        TYPE gb;
        TYPE haa;
        TYPE hab;
        TYPE hbb;

        TYPE f = fval2(&ga, &gb, &haa, &hab, &hbb, x, *a, *b);


        f *= fcond;
        ga *= fcond;
        gb *= fcond;
        haa *= fcond;
        hab *= fcond;
        hbb *= fcond;

        if (fabs(ga) < 1e-20 && fabs(gb) < 1e-20) {
            break;
        };


        TYPE va = (TYPE)0.0;
        TYPE vb = (TYPE)0.0;

        TYPE new_f = f;

        // Newton method for positive definite H with a good condition number
        TYPE det = haa * hbb - hab * hab;


        TYPE eta = (TYPE)1.0;

        // try sliding first
        TYPE da = gb - ga;
        TYPE db = ga - gb;

        TYPE d_ddf_d = da * (haa * da + hab * db) + db * (hab * da + hbb * db);
        TYPE df_d = ga * da + gb * db;



        if ((*a + *b - ab_bound) < 1e-9 && d_ddf_d > 1e-8 && slide_ok) {

            TYPE delta = -df_d / d_ddf_d;

            if (delta < 0) { // make sure delta is positive
                delta = -delta;
            }

            va = da * delta;
            vb = db * delta;

            if (*a + eta * va >= ua) eta = (ua - SMALL - *a) / va;
            if (*a + eta * va <= la) eta = (la + SMALL - *a) / va;

            if (*b + eta * vb >= ub) eta = (ub - SMALL - *b) / vb;
            if (*b + eta * vb <= lb) eta = (lb + SMALL - *b) / vb;

            update_ab(&new_f, x, a, b, f, fcond, ga, gb, eta, va, vb, ua, la, ub, lb, ab_bound);

            if (fabs(new_f - f) < (TYPE)1e-10 && fabs(*a - old_a) < (TYPE)1e-10 && fabs(*b - old_b) < (TYPE)1e-10) {
                slide_ok = 0;
            }
            else {
                slide_ok = 1;
                newton_ok = 1;
            }

            old_a = *a;
            old_b = *b;
            continue;

        }

        // try newton
        if (haa > 1e-30 && det > 1e-30 && newton_ok) { // positive definite and still improving

            va = -(hbb * ga - hab * gb) / det;
            vb = -(-hab * ga + haa * gb) / det;

            if (*a + eta * va >= ua) eta = (ua - SMALL - *a) / va;
            if (*a + eta * va <= la) eta = (la + SMALL - *a) / va;

            if (*b + eta * vb >= ub) eta = (ub - SMALL - *b) / vb;
            if (*b + eta * vb <= lb) eta = (lb + SMALL - *b) / vb;

            if (*a + eta * va + *b + eta * vb <= ab_bound) eta = (ab_bound + SMALL - *a - *b) / (va + vb);


            update_ab(&new_f, x, a, b, f, fcond, ga, gb, eta, va, vb, ua, la, ub, lb, ab_bound);

            if (fabs(new_f - f) < (TYPE)1e-10 && fabs(*a - old_a) < (TYPE)1e-10 && fabs(*b - old_b) < (TYPE)1e-10) {
                newton_ok = 0;
            }
            else {
                slide_ok = 1;
                newton_ok = 1;
            }

            old_a = *a;
            old_b = *b;
            continue;


        }


        // try gradient
        TYPE eta_a = find_eta(&va, ga, haa, *a, la > (ab_bound - *b) ? la : (ab_bound - *b), ua);

        TYPE eta_b = find_eta(&vb, gb, hbb, *b, lb > (ab_bound - *a) ? lb : (ab_bound - *a), ub);

        if (fabs(eta_a * va) > fabs(eta_b * vb)) {
            eta = eta_a;
            vb = (TYPE)0.0;
        }
        else {
            eta = eta_b;
            va = (TYPE)0.0;
        }

        update_ab(&new_f, x, a, b, f, fcond, ga, gb, eta, va, vb, ua, la, ub, lb, ab_bound);

        if (fabs(new_f - f) < (TYPE)1e-10 && fabs(*a - old_a) < (TYPE)1e-10 && fabs(*b - old_b) < (TYPE)1e-10) {
            cc++;
        }
        else {
            cc = 0;

            slide_ok = 1;
            newton_ok = 1;

        }

        if (cc > 1) {
            break;
        }

        old_a = *a;
        old_b = *b;

    }
}


/* 1-dimensional optimization */

TYPE nr(TYPE(*fval)(TYPE*, TYPE*, data_t*, TYPE),
    data_t* x,
    TYPE p0,
    TYPE lb,
    TYPE ub) {

    TYPE p = p0;

    TYPE old_p = p;

    int cc = 0;

    for (int i = 0; i < MAX_INNER; i++) {

        TYPE f, g, h;

        f = fval(&g, &h, x, p);

        TYPE v;

        TYPE eta = find_eta(&v, g, h, p, lb, ub);


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


        if (fabs(fnew - f) < (TYPE)1e-10 && fabs(p - old_p) < (TYPE)1e-10) {
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



TYPE max_val3(TYPE a, TYPE b, TYPE c) {
    if (c > a&& c > b) {
        return c;
    }
    else {
        return (a > b ? a : b);
    }
}

TYPE min_val3(TYPE a, TYPE b, TYPE c) {
    if (c < a && c < b) {
        return c;
    }
    else {
        return (a > b ? b : a);
    }
}




void do_ibb_test(data_t* x) {

    TYPE ab_bound = x->works[0];

    TYPE lower_bound = (TYPE)1.0;
    TYPE upper_bound = (TYPE)1e5;

    // real lower bound
    TYPE tmp = (lower_bound > (ab_bound / (TYPE)2.0) ? lower_bound : (ab_bound / (TYPE)2.0));

    /****/
    // try 3 initializations

    TYPE ab0 = tmp + (TYPE)1.0;
    ab0 = nr(fval_ab, x, ab0, tmp, upper_bound);
    TYPE fab0 = fval(0, 0, x, ab0, ab0, -1);

    TYPE ab1 = upper_bound - (TYPE)1.0;
    ab1 = nr(fval_ab, x, ab1, tmp, upper_bound);
    TYPE fab1 = fval(0, 0, x, ab1, ab1, -1);

    TYPE ab2 = (upper_bound + tmp) / (TYPE)2.0;
    ab2 = nr(fval_ab, x, ab2, tmp, upper_bound);
    TYPE fab2 = fval(0, 0, x, ab2, ab2, -1);


    TYPE f0 = fab0;
    TYPE ab = ab0;

    if (f0 > fab1) {
        f0 = fab1;
        ab = ab1;
    }

    if (f0 > fab2) {
        f0 = fab2;
        ab = ab2;
    }

    TYPE alp = ab;
    TYPE bet = ab;

    nr2b_projection(fval2, x, &alp, lower_bound, upper_bound, &bet, lower_bound, upper_bound, ab_bound, (TYPE)1.0);

    TYPE f = fval2(0, 0, 0, 0, 0, x, alp, bet);

    TYPE fmin = f;
    TYPE amin = alp;
    TYPE bmin = bet;


    alp = tmp + (TYPE)1.0;
    bet = upper_bound - (TYPE)1.0;

    nr2b_projection(fval2, x, &alp, lower_bound, upper_bound, &bet, lower_bound, upper_bound, ab_bound, (TYPE)1.0);

    f = fval2(0, 0, 0, 0, 0, x, alp, bet);

    if (f < fmin) {
        fmin = f;
        amin = alp;
        bmin = bet;
    }


    alp = upper_bound - (TYPE)1.0;
    bet = tmp + (TYPE)1.0;

    nr2b_projection(fval2, x, &alp, lower_bound, upper_bound, &bet, lower_bound, upper_bound, ab_bound, (TYPE)1.0);

    f = fval2(0, 0, 0, 0, 0, x, alp, bet);

    if (f < fmin) {
        fmin = f;
        amin = alp;
        bmin = bet;
    }

    x->f0 = -f0;

    x->f = -fmin;

    x->alp = amin;

    x->bet = bmin;

}


void ibb(int* lK,
    double* aa,
    double* bb,
    double* taa,
    double* tbb,
    int* lN,
    double* mem,
    int* no_threads,
    double* pval,
    double* fc) {

    int verbose = *no_threads;

    if (*no_threads < 0) {
        *no_threads = -*no_threads;
    }

    int stop_sig = 0;

    TYPE lower_bound = (TYPE)(mem[0]);
    TYPE tail = (TYPE)(mem[1]);

    int init_block = 3 * _Z; // qz + qlw +
    int block = 4 * (*lN) + (*lN) * _Z + 5 * _Z;

    TYPE* works = (TYPE*)mem;

    // initial block
    TYPE* log_one_minus_z = works;
    TYPE* log_one_plus_z = works + _Z;
    TYPE* phi = works + 2 * _Z;

    int thres_display = 0;

    forint(i, 0, _Z) {
        log_one_minus_z[i] = log1p(-_qz[i]);
        log_one_plus_z[i] = log1p(_qz[i]);
        phi[i] = ((TYPE)1.0 + _qz[i]) / ((TYPE)1.0 - _qz[i]);
    }

    #ifdef _OPENMP
    omp_set_dynamic(0);

    omp_set_num_threads(*no_threads);

    if (*no_threads > 1) {
        if (verbose > 0) Rprintf("Using %d threads ...\nNo. of data rows = %d, no. of pair(s) = %d.\n", *no_threads, *lK, *lN);
    }
    else {
        if (verbose > 0) Rprintf("Using a single CPU core ...\nNo. of data rows = %d, no. of pair(s) = %d.\n", *lK, *lN);
    }
    #else
    if (verbose > 0) Rprintf("Using a single CPU core ...\nNo. of data rows = %d, no. of pair(s) = %d.\n", *lK, *lN);
    #endif

    int k = 0;

    //#pragma omp parallel for schedule(static)
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < *lK; i++) {

        if (stop_sig) continue;

        int thread_id = 0;

        #ifdef _OPENMP
        thread_id = omp_get_thread_num();
        #endif

        data_t x;

        x.N = *lN;
        x.qz = _qz;
        x.qlw = _qw;
        x.Z = _Z;

        x.log_one_minus_z = log_one_minus_z;
        x.log_one_plus_z = log_one_plus_z;

        int mem_start = init_block + thread_id * block;

        x.a = works + mem_start;
        x.b = works + mem_start + (*lN);
        x.ta = works + mem_start + 2 * (*lN);
        x.tb = works + mem_start + 3 * (*lN);

        x.ll = works + mem_start + 4 * (*lN);

        x.works = works + mem_start + 4 * (*lN) + (*lN) * _Z;

        int ind = i * (*lN);

        TYPE* ptr = x.ll;

        forint(j, 0, (*lN)) {

            x.a[j] = (TYPE)aa[ind];
            x.b[j] = (TYPE)bb[ind];
            x.ta[j] = (TYPE)taa[ind];
            x.tb[j] = (TYPE)tbb[ind];

            if ((x.a[j] + x.b[j]) > 0) {

                TYPE factor = x.ta[j] / x.tb[j];

                forint(k, 0, _Z) {

                    TYPE p = phi[k] / (phi[k] + factor);
                    *ptr++ = x.a[j] * log1p(-p) + x.b[j] * log(p);

                }
            }
            else {
                forint(k, 0, _Z) {

                    TYPE p = phi[k] / (phi[k] + (TYPE)1.0);
                    *ptr++ = ZERO_REG * (log1p(-p) + log(p));
                }
            }

            ind++;
        }

        // initial values
        x.works[0] = lower_bound;

        do_ibb_test(&x);

        fc[i] = (double)(x.alp / x.bet);

        double g = 2.0 * ((double)x.f - (double)x.f0);

        if (tail > 0.5) {
            if (x.alp > x.bet) {
                pval[i] = pnorm(sqrt(g), 0.0, 1.0, 0, 0);
            }
            else {
                pval[i] = pnorm(-sqrt(g), 0.0, 1.0, 0, 0);
            }
        }
        else {
            if (tail < -0.5) {
                if (x.alp > x.bet) {
                    pval[i] = pnorm(sqrt(g), 0.0, 1.0, 1, 0);
                }
                else {
                    pval[i] = pnorm(-sqrt(g), 0.0, 1.0, 1, 0);
                }
            }
            else {
                pval[i] = pchisq(g, 1.0, 0, 0);
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
