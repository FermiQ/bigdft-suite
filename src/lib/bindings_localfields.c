/**************************************/
/* BigDFT_LocalFields data structure. */
/**************************************/

#include <config.h>

#ifdef HAVE_GLIB
#include <glib-object.h>
#endif

#include "bigdft.h"
#include "bindings.h"
#include "bindings_api.h"

#include <string.h>
#include <stdlib.h>
#include <stdio.h>

static void bigdft_localfields_dispose(GObject *localfields);
static void bigdft_localfields_finalize(GObject *localfields);

#ifdef HAVE_GLIB
enum {
  V_EXT_READY_SIGNAL,
  LAST_SIGNAL
};

G_DEFINE_TYPE(BigDFT_LocalFields, bigdft_localfields, G_TYPE_OBJECT)

static guint bigdft_localfields_signals[LAST_SIGNAL] = { 0 };

static void bigdft_localfields_class_init(BigDFT_LocalFieldsClass *klass)
{
  /* Connect the overloading methods. */
  G_OBJECT_CLASS(klass)->dispose      = bigdft_localfields_dispose;
  G_OBJECT_CLASS(klass)->finalize     = bigdft_localfields_finalize;
  /* G_OBJECT_CLASS(klass)->set_property = visu_data_set_property; */
  /* G_OBJECT_CLASS(klass)->get_property = visu_data_get_property; */

  bigdft_localfields_signals[V_EXT_READY_SIGNAL] =
    g_signal_new("v-ext-ready", G_TYPE_FROM_CLASS(klass),
                 G_SIGNAL_RUN_LAST | G_SIGNAL_NO_RECURSE | G_SIGNAL_NO_HOOKS,
		 0, NULL, NULL, g_cclosure_marshal_VOID__VOID,
                 G_TYPE_NONE, 0, NULL);
}
#endif

static void bigdft_localfields_init(BigDFT_LocalFields *obj)
{
#ifdef HAVE_GLIB
  memset(obj + sizeof(GObject), 0, sizeof(BigDFT_LocalFields) - sizeof(GObject));
#else
  memset(obj, 0, sizeof(BigDFT_LocalFields));
#endif
  FC_FUNC_(localfields_new, LOCALFIELDS_NEW)(&obj->data, &obj->rhod, &obj->dpcom);
  FC_FUNC_(initialize_dft_local_fields, INITIALIZE_DFT_LOCAL_FIELDS)(obj->data);
}
static void bigdft_localfields_dispose(GObject *obj)
{
#ifdef HAVE_GLIB
  BigDFT_LocalFields *denspot = BIGDFT_LOCALFIELDS(obj);

  if (denspot->dispose_has_run)
    return;
  denspot->dispose_has_run = TRUE;

  g_object_unref(G_OBJECT(denspot->glr));
  /* Chain up to the parent class */
  G_OBJECT_CLASS(bigdft_localfields_parent_class)->dispose(obj);
#endif
}
static void bigdft_localfields_finalize(GObject *obj)
{
  BigDFT_LocalFields *denspotd = BIGDFT_LOCALFIELDS(obj);

  FC_FUNC_(localfields_free, LOCALFIELDS_FREE)(&denspotd->data);

#ifdef HAVE_GLIB
  G_OBJECT_CLASS(bigdft_localfields_parent_class)->finalize(obj);
#endif
}

BigDFT_LocalFields* bigdft_localfields_new (const BigDFT_LocReg *glr,
                                            const BigDFT_Inputs *in,
                                            guint iproc, guint nproc)
{
  BigDFT_LocalFields *localfields;
  double hh[3];
  guint ndegree_ip = 16, verb = 0;

#ifdef HAVE_GLIB
  localfields = BIGDFT_LOCALFIELDS(g_object_new(BIGDFT_LOCALFIELDS_TYPE, NULL));
  g_object_ref(G_OBJECT(glr));
#else
  localfields = g_malloc(sizeof(BigDFT_LocalFields));
  memset(localfields, 0, sizeof(BigDFT_LocalFields));
  bigdft_localfields_init(localfields);
#endif
  localfields->glr   = glr;

  hh[0] = glr->h[0] * 0.5;
  hh[1] = glr->h[1] * 0.5;
  hh[2] = glr->h[2] * 0.5;
  FC_FUNC_(denspot_communications, DENSPOT_COMMUNICATIONS)
    (&iproc, &nproc, glr->d, hh, hh + 1, hh + 2, in->data,
     glr->atoms->data, glr->atoms->rxyz.data, glr->radii, localfields->dpcom, localfields->rhod);
  FC_FUNC(allocaterhopot, ALLOCATERHOPOT)(&iproc, glr->data,
                                          hh, hh + 1, hh + 2, in->data,
                                          glr->atoms->data, glr->atoms->rxyz.data,
                                          localfields->data);
  FC_FUNC_(localfields_copy_metadata, LOCALFIELDS_COPY_METADATA)
    (localfields->data, &localfields->rhov_is, localfields->h,
     &localfields->psoffset);
  GET_ATTR_DBL   (localfields, LOCALFIELDS, rhov,  RHOV);
  GET_ATTR_DBL_4D(localfields, LOCALFIELDS, v_ext, V_EXT);
  GET_ATTR_DBL_4D(localfields, LOCALFIELDS, v_xc,  V_XC);
  
  FC_FUNC_(system_createkernels, SYSTEM_CREATEKERNELS)
    (&iproc, &nproc, &verb, &glr->geocode, glr->d, hh,
     in->data, localfields->data);
  GET_ATTR_DBL   (localfields, LOCALFIELDS, pkernel,    PKERNEL);
  GET_ATTR_DBL   (localfields, LOCALFIELDS, pkernelseq, PKERNELSEQ);

  return localfields;
}
void bigdft_localfields_free(BigDFT_LocalFields *denspot)
{
#ifdef HAVE_GLIB
  g_object_unref(G_OBJECT(denspot));
#else
  bigdft_localfields_finalize(denspot);
  g_free(denspot);
#endif
}
void bigdft_localfields_create_effective_ionic_pot(BigDFT_LocalFields *denspot,
                                                   const BigDFT_Inputs *in,
                                                   guint iproc, guint nproc)
{
  double hh[3];
  
  hh[0] = denspot->glr->h[0] * 0.5;
  hh[1] = denspot->glr->h[1] * 0.5;
  hh[2] = denspot->glr->h[2] * 0.5;

  FC_FUNC(createeffectiveionicpotential, CREATEEFFECTIVEIONICPOTENTIAL)
    (&iproc, &nproc, in->data, denspot->glr->atoms->data, denspot->glr->atoms->rxyz.data,
     denspot->glr->atoms->shift, denspot->glr->data, hh, hh + 1, hh + 2,
     denspot->dpcom, denspot->pkernel, denspot->v_ext, in->elecfield,
     &denspot->psoffset);

#ifdef HAVE_GLIB
  g_signal_emit(G_OBJECT(denspot), bigdft_localfields_signals[V_EXT_READY_SIGNAL],
		0 /* details */, NULL);
#endif
}
