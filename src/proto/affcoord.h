/* affcoord.c */
void cac_get_affcoord_pixel(float v0x, float v0y, float *cage, int cage_size, float *affcoord);
float **cac_get_affcoord(struct vector *points, struct vector *cagepoints);
void cac_setup_omega1_omega2(double *xin, int init, struct resources_common *res_common);
void cac_update_points_non_uniform_sampling(struct resources_common *res_common, double *xin);
void cac_interpolate_points_non_uniform_sampling(struct resources_common *res_common);
void cac_interpolate_points_uniform_sampling(struct resources_common *res_common);
