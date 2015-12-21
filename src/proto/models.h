/* models.c */
void cac_mean(struct image *img, struct image *mask_in, struct vector *cage_init, struct vector *cage_curr, struct vector *cage_out);
void cac_gaussian(struct image *img, struct image *mask_in, struct vector *cage_init, struct vector *cage_curr, struct vector *cage_out);
void cac_common(int algorithm, struct image *img, struct image *mask_in, struct vector *cage_init, struct vector *cage_curr, struct vector *cage_out, struct resources_common *res_common);
