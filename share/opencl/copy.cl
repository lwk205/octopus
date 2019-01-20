/*
 Copyright (C) 2010 X. Andrade

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2, or (at your option)
 any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 02110-1301, USA.

*/

#include <cl_global.h>

__kernel void copy(const int np,
		   const __global double * restrict xx, const int ldxx,
		   __global double * restrict yy, const int ldyy){
  
  int ist = get_global_id(0);
  int ip = get_global_id(1);
  
  if(ip < np) yy[(ip<<ldyy) + ist] = xx[(ip<<ldxx) + ist];

}

/*
 Local Variables:
 mode: c
 coding: utf-8
 End:
*/
