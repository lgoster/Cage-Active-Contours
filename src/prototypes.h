/* affcoord.c */
void cac_get_affcoord_pixel(float v0x, float v0y, float *cage, int cage_size, float *affcoord);
float **cac_get_affcoord(struct vector *points, struct vector *cagepoints);
void cac_setup_omega1_omega2(double *xin, int init, struct resources_common *res_common);
void cac_update_points_non_uniform_sampling(struct resources_common *res_common, double *xin);
void cac_interpolate_points_non_uniform_sampling(struct resources_common *res_common);
void cac_interpolate_points_uniform_sampling(struct resources_common *res_common);
/* energy_functions.c */
void cac_compute_mean(struct resources_common *res_common);
void cac_compute_gradient_mean(struct resources_common *res_common, double *deriv_int, double *deriv_ext);
void cac_compute_variance2(struct resources_common *res_common, float *s2int, float *s2ext);
double cac_compute_dist_constraints(struct resources_common *res_common, double *xin);
int cac_mean_evaluate(double *xin, double *f, double *deriv, void *instance);
int cac_gaussian_evaluate(double *xin, double *f, double *deriv, void *instance);
/* image.c */
struct image *cac_image_new(void);
struct image *cac_image_alloc(int nrow, int ncol);
void cac_image_delete(struct image *image);
void cac_image_copy(struct image *in, struct image *out);
void cac_image_clear(struct image *im, double value);
void cac_image_line_draw(struct image *image, int a0, int b0, int a1, int b1, float c);
/* image_gradient.c */
void cac_get_gradient_images(struct resources_common *res);
struct image *cac_gradient_images_finite_differences_convolution(struct image *in, struct vector *filter_x, struct vector *filter_y);
struct vector *cac_get_filter_gradient(int filter_type, int derivative_order, int length, int index);
/* image_io.c */
struct image *cac_pgm_read_image(char *name);
void cac_pgm_write_image(struct image *image, char *name);
void cac_write_raw_image(struct image *image, char *name);
/* main.c */
int main(int argc, char *argv[]);
/* models.c */
void cac_mean(struct image *img, struct image *mask_in, struct vector *cage_init, struct vector *cage_curr, struct vector *cage_out);
void cac_gaussian(struct image *img, struct image *mask_in, struct vector *cage_init, struct vector *cage_curr, struct vector *cage_out);
void cac_common(int algorithm, struct image *img, struct image *mask_in, struct vector *cage_init, struct vector *cage_curr, struct vector *cage_out, struct resources_common *res_common);
/* newton_algorithm.c */
void cac_newton_algorithm(int n, double *x, tnc_function *proc_evaluate, struct resources_common *res_common);
/* queue.c */
struct queue *cac_new_queue(int size_elem, int expand_size);
void cac_delete_queue(struct queue *q);
int cac_is_queue_empty(struct queue *q);
void cac_put_elem_queue(char *elem, struct queue *q);
void cac_get_elem_queue(char *elem, struct queue *q);
int cac_get_queue_nb_elements(struct queue *q);
int cac_is_queue_full(struct queue *q);
void cac_expand_queue(struct queue *q);
/* tnc.c */
int tnc(int n, double x[], double *f, double g[], tnc_function *function, void *state, double low[], double up[], double scale[], double offset[], int messages, int maxiter, int maxCGit, int maxnfeval, double eta, double stepmx, double accuracy, double fmin, double ftol, double xtol, double pgtol, double rescale, int *nfeval);
void tnc_coercex(int n, double x[], double low[], double up[]);
void tnc_unscalex(int n, double x[], double xscale[], double xoffset[]);
void tnc_scalex(int n, double x[], double xscale[], double xoffset[]);
void tnc_scaleg(int n, double g[], double xscale[], double fscale);
void tnc_setConstraints(int n, double x[], int pivot[], double xscale[], double xoffset[], double low[], double up[]);
tnc_rc tnc_minimize(int n, double x[], double *f, double gfull[], tnc_function *function, void *state, double xscale[], double xoffset[], double *fscale, double low[], double up[], tnc_message messages, int maxCGit, int maxiter, int maxnfeval, int *nfeval, double eta, double stepmx, double accuracy, double fmin, double ftol, double xtol, double pgtol, double rescale);
void tnc_project(int n, double x[], int pivot[]);
void tnc_projectConstants(int n, double x[], double xscale[]);
double tnc_stepMax(double step, int n, double x[], double dir[], int pivot[], double low[], double up[], double xscale[], double xoffset[]);
tnc_logical tnc_addConstraint(int n, double x[], double p[], int pivot[], double low[], double up[], double xscale[], double xoffset[]);
tnc_logical tnc_removeConstraint(double gtpnew, double gnorm, double pgtolfs, double f, double fLastConstraint, double g[], int pivot[], int n);
int tnc_direction(double *zsol, double *diagb, double *x, double g[], int n, int maxCGit, int *niter, int maxnfeval, int *nfeval, tnc_logical upd1, double yksk, double yrsr, double *sk, double *yk, double *sr, double *yr, tnc_logical lreset, tnc_function *function, void *state, double xscale[], double xoffset[], double fscale, int *pivot, double accuracy, double gnorm, double xnorm, double low[], double up[]);
void tnc_diagonalScaling(int n, double e[], double v[], double gv[], double r[]);
double tnc_initialStep(double fnew, double fmin, double gtp, double smax);
int tnc_hessianTimesVector(double v[], double gv[], int n, double x[], double g[], tnc_function *function, void *state, double xscale[], double xoffset[], double fscale, double accuracy, double xnorm, double low[], double up[]);
int tnc_msolve(double g[], double y[], int n, double sk[], double yk[], double diagb[], double sr[], double yr[], tnc_logical upd1, double yksk, double yrsr, tnc_logical lreset);
void tnc_ssbfgs(int n, double gamma, double sj[], double hjv[], double hjyj[], double yjsj, double yjhyj, double vsj, double vhyj, double hjp1v[]);
int tnc_initPreconditioner(double diagb[], double emat[], int n, tnc_logical lreset, double yksk, double yrsr, double sk[], double yk[], double sr[], double yr[], tnc_logical upd1);
tnc_ls_rc tnc_linearSearch(int n, double *low, double *up, double *xscale, double *xoffset, double fscale, int *pivot, double *x, double *f, double *gfull, double *s, double *stp, double max_step, double dginit, double xtol, int maxnfeval, int *count, tnc_function *function, void *state, int only_line_search);
int tnc_update_trial_interval(double *x, double *fx, double *dx, double *y, double *fy, double *dy, double *t, double *ft, double *dt, const double tmin, const double tmax, int *brackt);
double tnc_mchpr1(void);
void tnc_dxpy1(int n, double dx[], double dy[]);
void tnc_daxpy1(int n, double da, double dx[], double dy[]);
void tnc_dcopy1(int n, double dx[], double dy[]);
void tnc_dneg1(int n, double v[]);
double tnc_ddot1(int n, double dx[], double dy[]);
double tnc_dnrm21(int n, double dx[]);
/* tnc_utils.c */
void cac_tnc_print_msg(char *msg, void *instance);
void cac_tnc_setup_omega1_omega2(double *xin, void *instance);
void cac_tnc_set_uniform_sampling(int value, void *instance);
void cac_tnc_set_project_gradient(int value, void *instance);
int cac_tnc_get_project_gradient(void *instance);
void cac_tnc_project_gradient_line_vertex_to_center(double *xin, double *g, void *instance);
double cac_tnc_get_maxstep_linsrch(int n, double *x, double *pk, void *instance);
void cac_tnc_normalize_gradient(int n, double *g);
int cac_progress(double *x, double *pk, double fx, double gnorm, double alpha, int k, void *instance);
/* utils.c */
void *cac_xmalloc(size_t size);
void cac_print_msg(FILE *fp, const char *pFormat, ...);
void cac_error(char *str);
int cac_inside_image_support(struct image *u, double x, double y);
void cac_project_to_imsupport(struct image *u, double *x, double *y);
double cac_get_value_image(struct image *u, double x, double y);
double cac_get_bilinear_interpolation_image(struct image *in, double xin, double yin);
void cac_contour_get_interior_contour(struct vector *contour, struct image *img, int conn);
void cac_gaussian_filter(float sigma, float ratio, int max_samples, float **filter, int *nsamples);
void cac_circular_convolution(float *input, int npoints, float *filter, int length_fil, float *output);
void cac_contour_blur(struct vector *contour);
int cac_imagetovector_count_non_zero(struct image *mask);
struct vector *cac_imagetovector(struct image *mask);
void cac_filloutside_queue(struct image *contour, struct image *mask, struct queue *q);
void cac_filloutside(struct image *contour, struct image *mask, struct resources_common *res_common);
void cac_mask_invert(struct image *in, struct image *out);
/* vector.c */
struct vector *cac_vector_new(void);
struct vector *cac_vector_alloc(int N, int dim);
void cac_vector_delete(struct vector *vector);
void cac_vector_clear(struct vector *vector, float value);
/* vector_io.c */
struct vector *cac_read_points2d(char *fname);
void cac_write_points2d(struct vector *vector, int flag, char *fname);
