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
