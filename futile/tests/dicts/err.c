#include "futile.h"

#include <string.h>
#include <stdio.h>

#define ASSERT(T) {if (!(T)) {fprintf(stdout, "%s: failed\n", #T); return FALSE;}}
#define ALLEQ(T,U,N) {int it; for(it=0;it < N;it=it+1) { ASSERT(T[it]==U[it]);}}

gboolean test_f90_err()
{
  ErrId id;

  id = err_define("ERROR_TEST", "an error for test.", NULL);
  ASSERT(id != ERR_NOT_DEFINE);

  err_open_try();
  err_throw_by_id("something is wrong!", id);

  ASSERT(err_check(ERR_ANY));
  ASSERT(err_check(id));
  ASSERT(err_check_by_name("ERROR_TEST"));

  err_close_try();
  ASSERT(!err_check(ERR_ANY));

  err_open_try();
  err_throw_by_name("something is wrong!", "ERROR_TEST");

  ASSERT(err_check(ERR_ANY));
  ASSERT(err_check(id));
  ASSERT(err_check_by_name("ERROR_TEST"));

  err_close_try();
  ASSERT(!err_check(ERR_ANY));

  return TRUE;
}

#define RUN(T) {if (!T) return 1; else fprintf(stdout, "%s: OK\n", #T);}

int main(int argc, char **argv)
{
#ifdef GLIB_MAJOR_VERSION
  g_type_init();
#endif
  futile_initialize();
  RUN(test_f90_err());
  futile_finalize();
  return 0;
}
