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
