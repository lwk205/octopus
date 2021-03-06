!! Copyright (C) 2006-2011 M. Marques, U. De Giovannini
!!
!! This program is free software; you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation; either version 2, or (at your option)
!! any later version.
!!
!! This program is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with this program; if not, write to the Free Software
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

#include "global.h"

module pes_oct_m
  use geometry_oct_m
  use global_oct_m
  use grid_oct_m
  use hamiltonian_elec_oct_m
  use mesh_oct_m
  use messages_oct_m
  use mpi_oct_m
  use namespace_oct_m
  use output_oct_m
  use parser_oct_m
  use pes_mask_oct_m
  use pes_spm_oct_m
  use pes_flux_oct_m
  use restart_oct_m
  use simul_box_oct_m
  use states_elec_oct_m
  use varinfo_oct_m
    
  implicit none

  private

  public ::                             &
    pes_t,                              &
    pes_end,                            &
    pes_init,                           &
    pes_init_write,                     &
    pes_calc,                           &
    pes_output,                         &
    pes_load,                           &
    pes_dump

  type pes_t
    private
    logical, public :: calc_spm
    type(pes_spm_t) :: spm

    logical, public :: calc_mask
    type(pes_mask_t) :: mask

    logical, public :: calc_flux
    type(pes_flux_t) :: flux
    
  end type pes_t

  integer, parameter ::     &
    PHOTOELECTRON_NONE = 0, &
    PHOTOELECTRON_SPM  = 2, &
    PHOTOELECTRON_MASK = 4, &
    PHOTOELECTRON_FLUX = 8

contains

  ! ---------------------------------------------------------
  !elemental (PUSH/POP)_SUB are not PURE.
  subroutine pes_spm_nullify(this)
    type(pes_spm_t), intent(out) :: this
    !
    PUSH_SUB(pes_spm_nullify)
    this%wf            =>null()
    POP_SUB(pes_spm_nullify)
    return
  end subroutine pes_spm_nullify

  ! ---------------------------------------------------------
  !elemental (PUSH/POP)_SUB are not PURE.
  subroutine pes_mask_nullify(this)
    type(pes_mask_t), intent(out) :: this
    !
    PUSH_SUB(pes_mask_nullify)
    this%k=>null()
    this%mesh=>null()
    this%vec_pot=>null()
    this%Mk=>null()
    this%mask_R=>null()
    this%Lk=>null()
    POP_SUB(pes_mask_nullify)
    return
  end subroutine pes_mask_nullify

  ! ---------------------------------------------------------
  !elemental (PUSH/POP)_SUB are not PURE.
  subroutine pes_nullify(this)
    type(pes_t), intent(out) :: this
    !
    PUSH_SUB(pes_nullify)
    !this%calc_spm=.false.
    call pes_spm_nullify(this%spm)
    !this%calc_mask=.false.
    call pes_mask_nullify(this%mask)
    POP_SUB(pes_nullify)
    return
  end subroutine pes_nullify

  ! ---------------------------------------------------------
  subroutine pes_init(pes, namespace, mesh, sb, st, save_iter, hm, max_iter, dt)
    type(pes_t),         intent(out)   :: pes
    type(namespace_t),   intent(in)    :: namespace
    type(mesh_t),        intent(in)    :: mesh
    type(simul_box_t),   intent(in)    :: sb
    type(states_elec_t), intent(in)    :: st
    integer,             intent(in)    :: save_iter
    type(hamiltonian_elec_t), intent(in)    :: hm
    integer,             intent(in)    :: max_iter
    FLOAT,               intent(in)    :: dt

    character(len=50)    :: str
    integer :: photoelectron_flags

    PUSH_SUB(pes_init)
    
    call pes_nullify(pes)

    !%Variable PhotoElectronSpectrum
    !%Type integer
    !%Default none
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !% This variable controls the method used for the calculation of
    !% the photoelectron spectrum. You can specify more than one value
    !% by giving them as a sum, for example:
    !% <tt>PhotoElectronSpectrum = pes_spm + pes_mask</tt>
    !%Option none 0
    !% The photoelectron spectrum is not calculated. This is the default.
    !%Option pes_spm 2
    !% Store the wavefunctions at specific points in order to 
    !% calculate the photoelectron spectrum at a point far in the box as proposed in 
    !% A. Pohl, P.-G. Reinhard, and E. Suraud, <i>Phys. Rev. Lett.</i> <b>84</b>, 5090 (2000).
    !%Option pes_mask 4
    !% Calculate the photo-electron spectrum using the mask method.
    !% U. De Giovannini, D. Varsano, M. A. L. Marques, H. Appel, E. K. U. Gross, and A. Rubio,
    !% <i>Phys. Rev. A</i> <b>85</b>, 062515 (2012).
    !%Option pes_flux 8
    !% Calculate the photo-electron spectrum using the t-surff technique, <i>i.e.</i>, 
    !% spectra are computed from the electron flux through a surface close to the absorbing 
    !% boundaries of the box. (Experimental.)
    !% L. Tao and A. Scrinzi, <i>New Journal of Physics</i> <b>14</b>, 013021 (2012).
    !%End

    call parse_variable(namespace, 'PhotoElectronSpectrum', PHOTOELECTRON_NONE, photoelectron_flags)
    if(.not.varinfo_valid_option('PhotoElectronSpectrum', photoelectron_flags, is_flag = .true.)) then
      call messages_input_error(namespace, 'PhotoElectronSpectrum')
    end if
    
    pes%calc_spm  = bitand(photoelectron_flags, PHOTOELECTRON_SPM) /= 0
    pes%calc_mask = bitand(photoelectron_flags, PHOTOELECTRON_MASK) /= 0
    pes%calc_flux = bitand(photoelectron_flags, PHOTOELECTRON_FLUX) /= 0

    !Header Photoelectron info
    if(pes%calc_spm .or. pes%calc_mask .or. pes%calc_flux) then 
      write(str, '(a,i5)') 'Photoelectron'
      call messages_print_stress(stdout, trim(str), namespace=namespace)
    end if 

    
    if(pes%calc_spm)  call pes_spm_init(pes%spm, namespace, mesh, st, save_iter)
    if(pes%calc_mask) call pes_mask_init(pes%mask, namespace, mesh, sb, st, hm, max_iter,dt)
    if(pes%calc_flux) call pes_flux_init(pes%flux, namespace, mesh, st, hm, save_iter, max_iter)


    !Footer Photoelectron info
    if(pes%calc_spm .or. pes%calc_mask .or. pes%calc_flux) then 
      call messages_print_stress(stdout, namespace=namespace)
    end if 

    POP_SUB(pes_init)
  end subroutine pes_init


  ! ---------------------------------------------------------
  subroutine pes_end(pes)
    type(pes_t), intent(inout) :: pes

    PUSH_SUB(pes_end)

    if(pes%calc_spm)   call pes_spm_end  (pes%spm)
    if(pes%calc_mask) call pes_mask_end(pes%mask)
    if(pes%calc_flux) call pes_flux_end(pes%flux)

    POP_SUB(pes_end)
  end subroutine pes_end


  ! ---------------------------------------------------------
  subroutine pes_calc(pes, namespace, mesh, st, dt, iter, gr, hm, stopping)
    type(pes_t),              intent(inout) :: pes
    type(namespace_t),        intent(in)    :: namespace
    type(mesh_t),             intent(in)    :: mesh
    type(states_elec_t),      intent(inout) :: st
    type(grid_t),             intent(in)    :: gr
    FLOAT,                    intent(in)    :: dt
    integer,                  intent(in)    :: iter
    type(hamiltonian_elec_t), intent(in)    :: hm
    logical,                  intent(in)    :: stopping

    PUSH_SUB(pes_calc)

    if(pes%calc_spm)  call pes_spm_calc(pes%spm, st, mesh, dt, iter, hm)
    if(pes%calc_mask) call pes_mask_calc(pes%mask, namespace, mesh, st, dt, iter)
    if(pes%calc_flux) call pes_flux_calc(pes%flux, mesh, st, gr, hm, iter, dt, stopping)

    POP_SUB(pes_calc)
  end subroutine pes_calc


  ! ---------------------------------------------------------
  subroutine pes_output(pes, namespace, mesh, st, iter, outp, dt, gr, geo)
    type(pes_t),         intent(inout) :: pes
    type(namespace_t),   intent(in)    :: namespace
    type(mesh_t),        intent(in)    :: mesh
    type(states_elec_t), intent(in)    :: st
    integer,             intent(in)    :: iter
    type(output_t),      intent(in)    :: outp
    FLOAT,               intent(in)    :: dt
    type(grid_t),        intent(in)    :: gr
    type(geometry_t),    intent(in)    :: geo

    PUSH_SUB(pes_output)
    
    if(pes%calc_spm) call pes_spm_output(pes%spm, mesh, st, namespace, iter, dt)

    if(pes%calc_mask) call pes_mask_output (pes%mask, mesh, st, outp, namespace, "td.general/PESM", gr, geo,iter)

    if(pes%calc_flux) call pes_flux_output(pes%flux, mesh, st, namespace, dt)

    POP_SUB(pes_output)
  end subroutine pes_output


  ! ---------------------------------------------------------
  subroutine pes_dump(pes, namespace, restart, st, mesh, ierr)
    type(pes_t),         intent(in)  :: pes
    type(namespace_t),   intent(in)  :: namespace
    type(restart_t),     intent(in)  :: restart
    type(states_elec_t), intent(in)  :: st
    type(mesh_t),        intent(in)  :: mesh
    integer,             intent(out) :: ierr

    PUSH_SUB(pes_dump)

    ierr = 0

    if (restart_skip(restart)) then
      POP_SUB(pes_dump)
      return
    end if

    if (debug%info) then
      message(1) = "Debug: Writing PES restart."
      call messages_info(1)
    end if

    if (pes%calc_mask) then
      call pes_mask_dump(pes%mask, namespace, restart, st, ierr)
    end if

    if (pes%calc_flux) then
      call pes_flux_dump(restart, pes%flux, mesh, st, ierr)
    end if

    if (pes%calc_spm) then
      call pes_spm_dump(restart, pes%spm, st, ierr)
    end if

    if (debug%info) then
      message(1) = "Debug: Writing PES restart done."
      call messages_info(1)
    end if

    POP_SUB(pes_dump)
  end subroutine pes_dump


  ! ---------------------------------------------------------
  subroutine pes_load(pes, namespace, restart, st, ierr)
    type(pes_t),         intent(inout) :: pes
    type(namespace_t),   intent(in)    :: namespace
    type(restart_t),     intent(in)    :: restart
    type(states_elec_t), intent(inout) :: st
    integer,             intent(out)   :: ierr

    PUSH_SUB(pes_load)

    ierr = 0

    if (restart_skip(restart)) then
      ierr = -1
      POP_SUB(pes_load)
      return
    end if

    if (debug%info) then
      message(1) = "Debug: Reading PES restart."
      call messages_info(1)
    end if

    if (pes%calc_mask) then
      call pes_mask_load(pes%mask, namespace, restart, st, ierr)
    end if

    if(pes%calc_flux) then
      call pes_flux_load(restart, pes%flux, st, ierr)
    end if

    if (pes%calc_spm) then
      call pes_spm_load(restart, pes%spm, st, ierr)
    end if

    if (debug%info) then
      message(1) = "Debug: Reading PES restart done."
      call messages_info(1)
    end if

    POP_SUB(pes_load)
  end subroutine pes_load


  ! ---------------------------------------------------------
  subroutine pes_init_write(pes, mesh, st, namespace)
    type(pes_t),         intent(in)  :: pes
    type(mesh_t),        intent(in)  :: mesh
    type(states_elec_t), intent(in)  :: st
    type(namespace_t),   intent(in)  :: namespace


    PUSH_SUB(pes_init_write)

    if(mpi_grp_is_root(mpi_world)) then

      if(pes%calc_spm)   call pes_spm_init_write (pes%spm, mesh, st, namespace)

    end if

    POP_SUB(pes_init_write)
  end subroutine pes_init_write



end module pes_oct_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
