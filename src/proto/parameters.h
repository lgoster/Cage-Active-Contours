/* parameters.c */
void cac_parameters_init(struct parameters *param);
void cac_parameters_close(struct parameters *param);
void cac_set_property_number(struct parameters *param, int property, double value);
void cac_set_property_string(struct parameters *param, int property, const char *value);
