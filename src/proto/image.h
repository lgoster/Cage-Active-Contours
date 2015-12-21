/* image.c */
struct image *cac_image_new(void);
struct image *cac_image_alloc(int nrow, int ncol);
void cac_image_delete(struct image *image);
void cac_image_copy(struct image *in, struct image *out);
void cac_image_clear(struct image *im, double value);
void cac_image_line_draw(struct image *image, int a0, int b0, int a1, int b1, float c);
