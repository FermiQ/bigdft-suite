#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import __builtin__

#USE_CHECKOUT_SRC = True

#if USE_CHECKOUT_SRC:
#    sys.path.insert(0, '/home/caliste/local/jhbuild')
#    pkgdatadir = None
#    datadir = None
#    import jhbuild
#    srcdir = os.path.abspath(os.path.join(os.path.dirname(jhbuild.__file__), '..'))
#else:
#    pkgdatadir = "@pkgdatadir@"
#    datadir = "@datadir@"
#    srcdir = "@srcdir@"
#    if '@pythondir@' not in sys.path:
#        sys.path.insert(0, '@pythondir@')
#    try:
#        import jhbuild
#    except ImportError:
#        sys.path.insert(0, srcdir)
pkgdatadir = None
datadir = None
import jhbuild
srcdir = os.path.abspath(os.path.join(os.path.dirname(jhbuild.__file__), '../..'))

__builtin__.__dict__['PKGDATADIR'] = pkgdatadir
__builtin__.__dict__['DATADIR'] = datadir
__builtin__.__dict__['SRCDIR'] = srcdir

import jhbuild.main
jhbuild.main.main(sys.argv[1:])
