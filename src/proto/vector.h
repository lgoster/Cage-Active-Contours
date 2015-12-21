/* vector.c */
struct vector *cac_vector_new(void);
struct vector *cac_vector_alloc(int N, int dim);
void cac_vector_delete(struct vector *vector);
void cac_vector_clear(struct vector *vector, float value);
