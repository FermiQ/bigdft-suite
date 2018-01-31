module pspiof_m
  type pspiof_pspdata_t
    integer :: tutu
  end type pspiof_pspdata_t
  
  integer, parameter, public :: PSPIO_STRLEN_ERROR = 1024
  integer, parameter, public :: PSPIO_STRLEN_LINE = 256
  integer, parameter, public :: PSPIO_STRLEN_TITLE = 80
  integer, parameter, public :: PSPIO_STRLEN_DESCRIPTION = 4096
  integer, parameter, public :: PSPIO_SUCCESS = 0
  integer, parameter, public :: PSPIO_ERROR = -1
  integer, parameter, public :: PSPIO_EFILE_CORRUPT = 1
  integer, parameter, public :: PSPIO_EFILE_FORMAT = 2
  integer, parameter, public :: PSPIO_EGSL = 3
  integer, parameter, public :: PSPIO_EIO = 4
  integer, parameter, public :: PSPIO_ENOFILE = 5
  integer, parameter, public :: PSPIO_ENOMEM = 6
  integer, parameter, public :: PSPIO_ENOSUPPORT = 7
  integer, parameter, public :: PSPIO_ETYPE = 8
  integer, parameter, public :: PSPIO_EVALUE = 9
  integer, parameter, public :: PSPIO_FMT_NFORMATS = 19
  integer, parameter, public :: PSPIO_FMT_UNKNOWN = -1
  integer, parameter, public :: PSPIO_FMT_NONE = 0
  integer, parameter, public :: PSPIO_FMT_ABINIT_1 = 1
  integer, parameter, public :: PSPIO_FMT_ABINIT_2 = 2
  integer, parameter, public :: PSPIO_FMT_ABINIT_3 = 3
  integer, parameter, public :: PSPIO_FMT_ABINIT_4 = 4
  integer, parameter, public :: PSPIO_FMT_ABINIT_5 = 5
  integer, parameter, public :: PSPIO_FMT_ABINIT_6 = 6
  integer, parameter, public :: PSPIO_FMT_ABINIT_7 = 7
  integer, parameter, public :: PSPIO_FMT_ABINIT_8 = 8
  integer, parameter, public :: PSPIO_FMT_ABINIT_9 = 9
  integer, parameter, public :: PSPIO_FMT_ABINIT_10 = 10
  integer, parameter, public :: PSPIO_FMT_ABINIT_11 = 11
  integer, parameter, public :: PSPIO_FMT_ABINIT_17 = 12
  integer, parameter, public :: PSPIO_FMT_ATOM = 13
  integer, parameter, public :: PSPIO_FMT_FHI98PP = 14
  integer, parameter, public :: PSPIO_FMT_OCTOPUS_HGH = 15
  integer, parameter, public :: PSPIO_FMT_SIESTA = 16
  integer, parameter, public :: PSPIO_FMT_UPF = 17
  integer, parameter, public :: PSPIO_FMT_XML = 18
  integer, parameter, public :: PSPIO_EQN_DIRAC = 1
  integer, parameter, public :: PSPIO_EQN_SCALAR_REL = 2
  integer, parameter, public :: PSPIO_EQN_SCHRODINGER = 3
  integer, parameter, public :: PSPIO_SCM_UNKNOWN = 0
  integer, parameter, public :: PSPIO_SCM_BHS = 1
  integer, parameter, public :: PSPIO_SCM_GTH = 2
  integer, parameter, public :: PSPIO_SCM_HAMANN = 3
  integer, parameter, public :: PSPIO_SCM_HGH = 4
  integer, parameter, public :: PSPIO_SCM_HSC = 5
  integer, parameter, public :: PSPIO_SCM_KERKER = 6
  integer, parameter, public :: PSPIO_SCM_MRPP = 7
  integer, parameter, public :: PSPIO_SCM_RRKJ = 8
  integer, parameter, public :: PSPIO_SCM_TM = 9
  integer, parameter, public :: PSPIO_SCM_TM2 = 10
  integer, parameter, public :: PSPIO_SCM_RTM = 11
  integer, parameter, public :: PSPIO_SCM_ONCV = 12
  integer, parameter, public :: PSPIO_MESH_UNKNOWN = -1
  integer, parameter, public :: PSPIO_MESH_NONE = 0
  integer, parameter, public :: PSPIO_MESH_LOG1 = 1
  integer, parameter, public :: PSPIO_MESH_LOG2 = 2
  integer, parameter, public :: PSPIO_MESH_LINEAR = 3
  integer, parameter, public :: PSPIO_DIFF = -1
  integer, parameter, public :: PSPIO_EQUAL = -2
  integer, parameter, public :: PSPIO_MTEQUAL = -3
  integer, parameter, public :: PSPIO_INTERP_GSL_CSPLINE = 1
  integer, parameter, public :: PSPIO_INTERP_JB_CSPLINE = 2
  integer, parameter, public :: PSPIO_NLCC_UNKNOWN = -1
  integer, parameter, public :: PSPIO_NLCC_NONE = 0
  integer, parameter, public :: PSPIO_NLCC_FHI = 1
  integer, parameter, public :: PSPIO_NLCC_LOUIE = 2
  integer, parameter, public :: PSPIO_NLCC_TETER1 = 3
  integer, parameter, public :: PSPIO_NLCC_TETER2 = 4
  integer, parameter, public :: PSPIO_NLCC_ATOM = 5
  integer, parameter, public :: PSPIO_NLCC_ONCV = 5

contains

  function pspiof_pspdata_alloc(data)
    type(pspiof_pspdata_t) :: data
    integer :: pspiof_pspdata_alloc

    stop "FAKE PSPIO."
  end function pspiof_pspdata_alloc

  function pspiof_pspdata_free(data)
    type(pspiof_pspdata_t) :: data
    integer :: pspiof_pspdata_free

    stop "FAKE PSPIO."
  end function pspiof_pspdata_free

  function pspiof_pspdata_read(data, kind, filename)
    type(pspiof_pspdata_t) :: data
    integer :: kind
    character(len=1) :: filename(*)
    integer :: pspiof_pspdata_read

    stop "FAKE PSPIO."
  end function pspiof_pspdata_read
end module pspiof_m
