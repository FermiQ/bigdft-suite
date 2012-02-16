#ifndef BINDINGS_API_H
#define BINDINGS_API_H

void FC_FUNC_(atoms_new, ATOMS_NEW)(void *atoms, void *sym);
void FC_FUNC_(atoms_free, ATOMS_FREE)(void *atoms);
void FC_FUNC_(atoms_new_from_file, ATOMS_NEW_FROM_FILE)(int *lstat, void *atoms,
                                                        f90_pointer_double_2D *rxyz,
                                                        const gchar *filename, int *ln);
void FC_FUNC_(atoms_set_n_atoms, ATOMS_SET_N_ATOMS)(void *atoms,
                                                    f90_pointer_double_2D *rxyz, int *nat);
void FC_FUNC_(atoms_set_n_types, ATOMS_SET_N_TYPES)(void *atoms, int *ntypes);
void FC_FUNC_(atoms_set_symmetries, ATOMS_SET_SYMMETRIES)(void *atoms, double *rxyz,
                                                          int *disable, double *tol,
                                                          double *elecfield);
void FC_FUNC_(atoms_set_displacement, ATOMS_SET_DISPLACEMENT)(void *atoms, double *rxyz,
                                                              double *randdis);
void FC_FUNC_(atoms_set_name, ATOMS_SET_NAME)(void *atoms, int *ityp, gchar *name);
void FC_FUNC_(atoms_sync, ATOMS_SYNC)(void *atoms, double *alat1, double *alat2,
                                      double *alat3, gchar *geocode, gchar *format,
                                      gchar *units);
void FC_FUNC_(atoms_get_geocode, ATOMS_GET_GEOCODE)(void *atoms, char *geocode);
void FC_FUNC_(atoms_get_iatype, ATOMS_GET_IATYPE)(void *atoms, f90_pointer_int *iatype);
void FC_FUNC_(atoms_get_iasctype, ATOMS_GET_IASCTYPE)(void *atoms, f90_pointer_int *iasctype);
void FC_FUNC_(atoms_get_natpol, ATOMS_GET_NATPOL)(void *atoms, f90_pointer_int *natpol);
void FC_FUNC_(atoms_get_ifrztyp, ATOMS_GET_IFRZTYP)(void *atoms, f90_pointer_int *ifrztyp);
void FC_FUNC_(atoms_get_amu, ATOMS_GET_AMU)(void *atoms, f90_pointer_double *amu);
void FC_FUNC_(atoms_get_aocc, ATOMS_GET_AOCC)(void *atoms, f90_pointer_double_2D *aocc);
void FC_FUNC_(atoms_get_nelpsp, ATOMS_GET_NELPSP)(void *atoms, f90_pointer_int *nelpsp);
void FC_FUNC_(atoms_get_npspcode, ATOMS_GET_NPSPCODE)(void *atoms, f90_pointer_int *npspcode);
void FC_FUNC_(atoms_get_nzatom, ATOMS_GET_NZATOM)(void *atoms, f90_pointer_int *nzatom);
void FC_FUNC_(atoms_get_nlcc_ngv, ATOMS_GET_NLCC_NGV)(void *atoms, f90_pointer_int *nlcc_ngv);
void FC_FUNC_(atoms_get_nlcc_ngc, ATOMS_GET_NLCC_NGC)(void *atoms, f90_pointer_int *nlcc_ngc);
void FC_FUNC_(atoms_get_ixcpsp, ATOMS_GET_IXCPSP)(void *atoms, f90_pointer_int *ixcpsp);
void FC_FUNC_(atoms_get_radii_cf, ATOMS_GET_RADII_CF)(void *atoms, f90_pointer_double_2D *radii_cf);
void FC_FUNC_(atoms_get_psppar, ATOMS_GET_PSPPAR)(void *atoms, f90_pointer_double_3D *psppar);
void FC_FUNC_(atoms_get_nlccpar, ATOMS_GET_NLCCPAR)(void *atoms, f90_pointer_double_2D *nlccpar);
void FC_FUNC_(atoms_get_ig_nlccpar, ATOMS_GET_IG_NLCCPAR)(void *atoms, f90_pointer_double_2D *ig_nlccpar);
void FC_FUNC_(atoms_copy_nat, ATOMS_COPY_NAT)(void *atoms, int *nat);
void FC_FUNC_(atoms_copy_ntypes, ATOMS_COPY_NTYPES)(void *atoms, int *ntypes);
void FC_FUNC_(atoms_copy_geometry_data, ATOMS_COPY_GEOMETRY_DATA)
     (void *atoms, gchar *geocode, gchar *format, gchar *units);
void FC_FUNC_(atoms_copy_alat, ATOMS_COPY_ALAT)(void *atoms, double *alat1,
                                                double *alat2, double *alat3);
void FC_FUNC_(atoms_copy_psp_data, ATOMS_COPY_PSP_DATA)(void *atoms, int *natsc, int *donlcc);
void FC_FUNC_(atoms_copy_name, ATOMS_COPY_NAME)(void *atoms, int *ityp, gchar *name, int *ln);
void FC_FUNC_(init_atomic_values, INIT_ATOMIC_VALUES)(int *verb, void *atoms, int *ixc);
void FC_FUNC_(read_radii_variables, READ_RADII_VARIABLES)(void *atoms, double *radii_cf);
void FC_FUNC_(atoms_write, ATOMS_WRITE)(void *atoms, const gchar *filename, int *ln2,
                                        double *rxyz, f90_pointer_double *forces,
                                        const double *energy, const gchar *comment, int *ln);


void FC_FUNC_(localfields_new, LOCALFIELDS_NEW)(void *denspotd,
                                                void *rhod, void *dpcom);
void FC_FUNC_(localfields_free, LOCALFIELDS_FREE)(void *denspotd);
void FC_FUNC(allocaterhopot, ALLOCATERHOPOT)(const guint *iproc,
                                             const void *glr, const double *hxh,
                                             const double *hyh, const double *hzh,
                                             const void *in, const void *atoms,
                                             const double *rxyz,
                                             void *denspotd);
void FC_FUNC(system_createkernels, SYSTEM_CREATEKERNELS)
     (const guint *iproc, const guint *nproc, const guint *verb,
      const gchar *geocode, const void *d, const double *hh, 
      const void *in, void *denspot);

#endif
