/* funcionals.c */
void cac_compute_mean(struct resources_common *res_common);
void cac_compute_gradient_mean(struct resources_common *res_common, double *deriv_int, double *deriv_ext);
void cac_compute_variance2(struct resources_common *res_common, float *s2int, float *s2ext);
double cac_compute_dist_constraints(struct resources_common *res_common, double *xin);
int cac_mean_evaluate(double *xin, double *f, double *deriv, void *instance);
int cac_gaussian_evaluate(double *xin, double *f, double *deriv, void *instance);
int cac_hist_evaluate(double *xin, double *f, double *deriv, void *instance);
double cac_gaussfunc(double x, double sigma2);
double cac_gaussfunc_deriv(double x, double sigma2);
double cac_histogram_parzen_window(int type, double k, struct resources_common *res_common);
void cac_histogram_parzen_window_deriv(int type, double k, int j, struct resources_common *res_common, double *result_x, double *result_y);
