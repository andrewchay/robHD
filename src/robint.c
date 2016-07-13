#include <R.h>
#include <Rinternals.h>
#include <string.h>
#include <Rmath.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Applic.h>
#define EPSILON (1e-9)
static double *vector(int n)
{
  double *v;
  v = Calloc(n, double);
  return(v);
}

static void free_vector(double *v)
{
  Free(v);
}

static void fww(double *b, double theta, double *x, double *y, double *omega,
                int n, int p, double *val)
{
  int i, k;
  double *s;
  s = vector(n);
  for (i = 0; i < n; i++) s[i] = 0;
  for (i = 0; i < n; i++)
  {
    for (k = 0; k < p; k++)
      s[i] += x[i + k * n] * b[k];
  }
  for (i = 0; i < n; i++)
  {
    val[i] = exp(-pow((y[i] - s[i]), 2) / theta + log(omega[i])) * 2 / theta;
  }
  free_vector(s);
}


static int sgn(double z)
{
    if (z > 0) return(1);
    else if (z < 0) return(-1);
    else return(0);
}

static double soft(double lambda, double z)
{
    if (fabs(z) > lambda) return(sgn(z) * (fabs(z) - lambda));
    else return(0);
}

static double flasso(double lambda, double b, double w, double r)
{
  return(soft(lambda / w, b + r / w));
}

static double fmcp(double lambda, double kappa, double b, double w, double r)
{
  if (fabs(b + r / w) > lambda / kappa) return(b + r / w); else
  {
    if (w > kappa) return(soft(lambda / w, b + r / w) / (1 - kappa / w));
    else return(0);
  }
}


static int equalZero(double num)
{
  if (fabs(num) < EPSILON) return(1); else return(0);
}

static int checkConvergence(double *beta, double *beta_old, double eps,
                            int len)
{
  int i, converged = 1;
  for (i = 0; i <= len; i++)
  {
    if (!equalZero(beta[i]) && !equalZero(beta_old[i]))
    {
      if (fabs((beta[i]-beta_old[i]) / beta_old[i]) > eps)
      {
        converged = 0;
        break;
      }
    }
    else if (!equalZero(beta[i]) && equalZero(beta_old[i]))
    {converged = 0;break;}
    else if (equalZero(beta[i]) && !equalZero(beta_old[i]))
    {converged = 0;break;}
  }
  return(converged);
}

static void cdfit(double *x, double *y, double *omega, double *init, double lambda,
                  double kappa, double theta, char *penalty, double eps,
                  int max, int inmax, int n, int p, int *iter)
{
  double *beta_old, *beta, *w, *ww, *re;
  double  z, r = 0, s = 0;
  int i, j, itetime = 0, flag = 1, itetime2 = 0, flag2 = 1;
  beta = vector(p);
  beta_old = vector(p);
  w = vector(p); /* p dimensional vector for X'WX. */
  re = vector(n);
  ww = vector(n); /* n dimensional weight vector for diag(W). */
  for (i = 0; i < p; i++)
  {
    beta_old[i] = init[i];
    beta[i] = init[i];
  }
  for (i = 0; i < n; i++)
  {
    s = 0;
    for (j = 0; j < p; j++)
      s += x[i + j * n] * beta[j];
    re[i] = y[i] - s; /* re is the r vector. */
  }

  /* This is the MM layer. */
  while (flag && itetime <= max)
  {
    if (itetime > 0) for (i = 0; i < p; i++) beta_old[i] = beta[i];
    fww(beta, theta, x, y, omega, n, p, ww);
    for (j = 0; j < p; j++)
    {
      s = 0;
      for (i = 0; i < n; i++)
        s += x[i + j * n] * ww[i] * x[i + j * n];
      w[j] = s;
    }
    /* This is one coordinate descent layer. */
    itetime2 = 0;
    flag2 = 1;
    while (flag2 && itetime2 < inmax)
    {
      for (j = 0; j < p; j++)
      {
        r = 0;
        for (i = 0; i < n; i++) r += x[i + j * n] * ww[i] * re[i];
        if (strcmp(penalty, "MCP") == 0) z = fmcp(lambda, kappa, beta[j], w[j], r);
        if (strcmp(penalty, "LASSO") == 0) z = flasso(lambda, beta[j], w[j], r);
        if (fabs(beta[j] - z) > eps) {
          for (i = 0; i < n; i++) re[i] += x[i + j * n] * (beta[j] - z);}
        beta[j] = z;
      }
      itetime2++;
      flag2 = !checkConvergence(beta, beta_old, eps, p);
    }
    itetime++;
    flag = !checkConvergence(beta, beta_old, eps, p);
  }
  for (i = 0; i < p; i++) init[i] = beta[i];
  iter[0] = itetime;
  free_vector(beta);
  free_vector(beta_old);
  free_vector(w);
  free_vector(ww);
  free_vector(re);
}

SEXP pathSearch(SEXP Rx, SEXP Ry, SEXP Romega, SEXP Rlambda, SEXP Rkappa,
                SEXP Rtheta, SEXP Rn_lambda, SEXP Rn_theta,
                SEXP Reps, SEXP Rmax, SEXP Rinmax, SEXP Rpenalty, SEXP Rinit)
{
  const char *penalty;
  int  i, j, l, n, p, max, iter[1], n_lambda, n_theta, inmax;
  double eps;
  double *lambda, *y, *x, *init, *in, kappa, *theta, *omega;
  void cdfit(), free_vec();
  SEXP betahat,  betatild, Riter, return_list;
  p = Rf_ncols(Rx);
  n = Rf_nrows(Rx);
  Rx = coerceVector(Rx, REALSXP);
  Ry = coerceVector(Ry, REALSXP);
  Romega = coerceVector(Romega, REALSXP);
  Rlambda = coerceVector(Rlambda, REALSXP);
  Rkappa = coerceVector(Rkappa, REALSXP);
  Rtheta = coerceVector(Rtheta, REALSXP);
  Rn_lambda = coerceVector(Rn_lambda, INTSXP);
  Rn_theta = coerceVector(Rn_theta, INTSXP);
  Reps = coerceVector(Reps, REALSXP);
  Rinit = coerceVector(Rinit, REALSXP);
  Rmax = coerceVector(Rmax, INTSXP);
  Rinmax = coerceVector(Rinmax, INTSXP);


  lambda = REAL(Rlambda);
  kappa = REAL(Rkappa)[0];
  theta = REAL(Rtheta);
  penalty = CHAR(STRING_ELT(Rpenalty, 0));
  n_lambda = INTEGER(Rn_lambda)[0];
  n_theta = INTEGER(Rn_theta)[0];
  eps = REAL(Reps)[0];
  max = INTEGER(Rmax)[0];
  inmax = INTEGER(Rinmax)[0];
  x = REAL(Rx);
  y = REAL(Ry);
  omega = REAL(Romega);
  in = REAL(Rinit);
  init = vector(p);
  for (j = 0; j < p; j++) init[j] = in[j];
  PROTECT(betahat = Rf_allocMatrix(REALSXP, p, n_lambda * n_theta));
  PROTECT(betatild = Rf_allocMatrix(REALSXP, p, n_lambda * n_theta));
  PROTECT(Riter = Rf_allocMatrix(INTSXP, n_lambda, n_theta));
  PROTECT(return_list = Rf_allocVector(VECSXP, 6));
  /*LASSO path along lambda*/
  for (i = 0; i < n_lambda; i++){
    for (l = 0; l < n_theta; l++)
    {
      cdfit(x, y, omega, init, lambda[i], 0.0, theta[l], "LASSO", eps, max,
      inmax, n, p, iter);
      for (j = 0; j < p; j++)
        REAL(betahat)[(i + l * n_lambda) * p + j] = init[j];
      INTEGER(Riter)[i + l * n_lambda] = iter[0];
    }}
  SET_VECTOR_ELT(return_list, 0, betahat);
  /*MCP path along kappa with n_lambda interim points*/
  if (strcmp(penalty, "MCP") == 0)
  {
    for (i = 0; i < n_lambda; i++)
    {
      for (l = 0; l < n_theta; l++)
      {
        for (j = 0; j < p; j++) init[j] = REAL(betahat)[(i + l * n_lambda) * p + j];
        cdfit(x, y, omega, init, lambda[i], kappa, theta[l], "MCP", eps, max,
        inmax, n, p, iter);
        for (j = 0; j < p; j++)
          REAL(betatild)[(i + l * n_lambda) * p + j] = init[j];
        INTEGER(Riter)[i + l * n_lambda] = iter[0];
      }
    }
  }
  SET_VECTOR_ELT(return_list, 1, betatild);
  SET_VECTOR_ELT(return_list, 2, Riter);
  SET_VECTOR_ELT(return_list, 3, Rlambda);
  SET_VECTOR_ELT(return_list, 4, Rkappa);
  SET_VECTOR_ELT(return_list, 5, Rtheta);
  UNPROTECT(4);
  free_vector(init);
  return(return_list);
}

static R_CMethodDef DotCEntries[] = {
  {"pathSearch", (DL_FUNC) &pathSearch, 13},
  {NULL}
};

void R_init_pareto(DllInfo *info)
{
  R_registerRoutines(info, DotCEntries, NULL, NULL, NULL);
}
