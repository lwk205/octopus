/*
 Copyright (C) 2019 M. Lueders

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
#include <cl_complex.h>

__kernel void batch_mult_real(const int nst,
			   const int np,
			   const __global double * restrict psi1, 
			   const __global double * restrict psi2, 
         const int ldpsi,
			   __global double * restrict mod_sqr){
  
  int ist = get_global_id(0);
  int ip  = get_global_id(1);

  if(ip >= np || ist >= nst) return;

  double ff1 = psi1[(ip<<ldpsi) + ist];
  double ff2 = psi2[(ip<<ldpsi) + ist];
  mod_sqr[(ip<<ldpsi) + ist] = ff1*ff2;

}

__kernel void batch_mult_complex(const int nst,
			      const int np,
			      const __global double2 * restrict psi1, 
			      const __global double2 * restrict psi2, 
            const int ldpsi,
			      __global double * restrict mod_sqr){
  
  int ist = get_global_id(0);
  int ip  = get_global_id(1);

  if(ip >= np || ist >= nst) return;

  double2 ff1 = psi1[(ip<<ldpsi) + ist];
  double2 ff2 = psi2[(ip<<ldpsi) + ist];
  mod_sqr[(ip<<ldpsi) + ist] = (ff1.x*ff2.x + ff1.y*ff2.y);

}

/*
 Local Variables:
 mode: c
 coding: utf-8
 End:
*/
