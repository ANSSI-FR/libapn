# Copyright 2017 Jean-Pierre Flori, Jérémy Jean
# 
# This file is part of libapn.
# 
# libapn is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
# 
# libapn is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with libapn; if not, see <http://www.gnu.org/licenses/>.

AM_CPPFLAGS = -I$(top_srcdir)/src
AM_CXXFLAGS = $(PTHREAD_CFLAGS)
AM_LDFLAGS = $(PTHREAD_CFLAGS) $(PTHREAD_LIBS)

noinst_PROGRAMS = search_permutation search_function search_involution search_compressing

search_permutation_SOURCES = search_permutation.cpp
search_function_SOURCES = search_function.cpp
search_involution_SOURCES = search_involution.cpp
search_compressing_SOURCES = search_compressing.cpp
