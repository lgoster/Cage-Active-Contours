/* image_gradient.c */
void cac_get_gradient_images(struct resources_common *res);
struct image *cac_gradient_images_finite_differences_convolution(struct image *in, struct vector *filter_x, struct vector *filter_y);
struct vector *cac_get_filter_gradient(int filter_type, int derivative_order, int length, int index);
