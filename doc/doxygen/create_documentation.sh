#!/bin/bash

## Copyright (C) 2002-2006 J. Alberdi-Rodriguez
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2, or (at your option)
## any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
## 02110-1301, USA.
##
## $Id$

#save the current directory
cd `dirname $0`
dox_dir=$PWD
cd ../..
oct_base_dir=$PWD

# go to the source directory
cd $oct_base_dir/src

echo "Configuring"
autoreconf -i
./configure
cd src/include
make
# to produce config_F90.h and options.h
cd $oct_base_dir

# FIXME: set all -DHAVE_ via
# sed s'|/\* \#undef HAVE_|#define HAVE_|' ../config.h | grep '^#define' > include/config_F90.h
# FIXME: provide header files (e.g. fcs_config.h)

echo "Creating documentation"

# Doxyfile is for version 1.8.6
cp $dox_dir/Doxyfile .
cp $dox_dir/octopus.png .
#copy to the new location
sed -i "s|doxygen_doc|$dox_dir/$folder|" Doxyfile
#call doxygen
doxygen

cd ..
echo "Doxygen documentation created in $dox_dir/html/index.html"
