/* queue.c */
struct queue *cac_new_queue(int size_elem, int expand_size);
void cac_delete_queue(struct queue *q);
int cac_is_queue_empty(struct queue *q);
void cac_put_elem_queue(char *elem, struct queue *q);
void cac_get_elem_queue(char *elem, struct queue *q);
int cac_get_queue_nb_elements(struct queue *q);
int cac_is_queue_full(struct queue *q);
void cac_expand_queue(struct queue *q);
