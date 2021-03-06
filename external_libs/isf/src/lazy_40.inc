!> @file
!!   Filters for interpolating scaling functions (order 40)
!! @author
!!    Copyright (C) 2009-2011 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
integer, parameter :: m=42
real(kind=8), dimension(-m:m) :: ch = (/ &
     0.d0, 0.d0, 0.d0, &
     -2.33895335687458837d-13,0.d0,9.61499474542229436d-12,0.d0,&
     -1.9312403731519637d-10,0.d0,2.52621846791090201d-9,0.d0,&
     -2.42028027409528354d-8,0.d0,1.81103730854716044d-7,0.d0,&
     -1.10227332483179024d-6,0.d0,5.61214589957214346d-6,0.d0,&
     -0.0000244006343459658411d0,0.d0,0.0000920510703104425647d0,0.d0,&
     -0.000305221969976730609d0,0.d0,0.000899343879450152758d0,0.d0,&
     -0.00237826492565707063d0,0.d0,0.00569939227746221067d0,0.d0,&
     -0.0125090557778066702d0,0.d0,0.0254814099177543281d0,0.d0,&
     -0.0491427191270976329d0,0.d0,0.0930820915230908105d0,0.d0,&
     -0.189611667917407207d0,0.d0,0.628712372568244948d0,1.d0,0.628712372568244948d0,&
     0.d0,-0.189611667917407207d0,0.d0,0.0930820915230908105d0,0.d0,&
     -0.0491427191270976329d0,0.d0,0.0254814099177543281d0,0.d0,&
     -0.0125090557778066702d0,0.d0,0.00569939227746221067d0,0.d0,&
     -0.00237826492565707063d0,0.d0,0.000899343879450152758d0,0.d0,&
     -0.000305221969976730609d0,0.d0,0.0000920510703104425647d0,0.d0,&
     -0.0000244006343459658411d0,0.d0,5.61214589957214346d-6,0.d0,&
     -1.10227332483179024d-6,0.d0,1.81103730854716044d-7,0.d0,&
     -2.42028027409528354d-8,0.d0,2.52621846791090201d-9,0.d0,&
     -1.9312403731519637d-10,0.d0,9.61499474542229436d-12,0.d0,&
     -2.33895335687458837d-13,&
     0.d0, 0.d0, 0.d0   /)
real(kind=8), dimension(-m:m) ::  cg,cht,cgt

!******** coefficients for wavelet transform *********************
do i=-m,m
   cht(i)=0.d0
   cg(i)=0.d0
   cgt(i)=0.d0
enddo

! the normalization is chosen such that a constant function remains the same constant 
! on each level of the transform

cht( 0)=1.d0

! g coefficients from h coefficients
do i=-m,m-1
   cg(i+1)=cht(-i)*(-1.d0)**(i+1)
   cgt(i+1)=ch(-i)*(-1.d0)**(i+1)
enddo


