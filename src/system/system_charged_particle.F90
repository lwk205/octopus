!! Copyright (C) 2020 H. Appel, S. Ohlmann, M. Oliveira, N. Tancogne-Dejean
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

module system_charged_particle_oct_m
  use clock_oct_m
  use global_oct_m
  use interaction_abst_oct_m
  use interaction_coulomb_force_oct_m
  use io_oct_m
  use iso_c_binding
  use messages_oct_m
  use mpi_oct_m
  use namespace_oct_m
  use parser_oct_m
  use profiling_oct_m
  use propagator_abst_oct_m
  use quantity_oct_m
  use space_oct_m
  use system_abst_oct_m
  use system_classical_particle_oct_m
  use write_iter_oct_m

  implicit none

  private
  public ::               &
    charged_particle_t

  type, extends(classical_particle_t) :: charged_particle_t
    private

    FLOAT :: charge

  contains
    procedure :: add_interaction_partner => charged_particle_add_interaction_partner
    procedure :: has_interaction => charged_particle_has_interaction
    procedure :: initial_conditions => charged_particle_initial_conditions
    procedure :: do_td_operation => charged_particle_do_td
    procedure :: iteration_info => charged_particle_iteration_info
    procedure :: output_start => charged_particle_output_start
    procedure :: output_write => charged_particle_output_write
    procedure :: output_finish => charged_particle_output_finish
    procedure :: is_tolerance_reached => charged_particle_is_tolerance_reached
    procedure :: store_current_status => charged_particle_store_current_status
    procedure :: update_quantity => charged_particle_update_quantity
    procedure :: update_exposed_quantity => charged_particle_update_exposed_quantity
    procedure :: set_pointers_to_interaction => charged_particle_set_pointers_to_interaction
    procedure :: update_interactions_start => charged_particle_update_interactions_start
    procedure :: update_interactions_finish => charged_particle_update_interactions_finish
    final :: charged_particle_finalize
  end type charged_particle_t

  interface charged_particle_t
    procedure charged_particle_constructor
  end interface charged_particle_t

contains

  ! ---------------------------------------------------------
  !> The factory routine (or constructor) allocates a pointer of the
  !! corresponding type and then calls the init routine which is a type-bound
  !! procedure of the corresponding type. With this design, also derived
  !! classes can use the init routine of the parent class.
  function charged_particle_constructor(namespace) result(sys)
    class(charged_particle_t), pointer  :: sys
    type(namespace_t),       intent(in) :: namespace

    PUSH_SUB(charged_particle_constructor)

    SAFE_ALLOCATE(sys)

    call charged_particle_init(sys, namespace)

    POP_SUB(charged_particle_constructor)
  end function charged_particle_constructor

  ! ---------------------------------------------------------
  !> The init routine is a module level procedure
  !! This has the advantage that different classes can have different
  !! signatures for the initialization routines because they are not
  !! type-bound and thus also not inherited.
  subroutine charged_particle_init(this, namespace)
    class(charged_particle_t), intent(inout) :: this
    type(namespace_t),      intent(in) :: namespace

    PUSH_SUB(charged_particle_init)

    call classical_particle_init(this%classical_particle_t, namespace)

    !%Variable ClassicalParticleCharge
    !%Type float
    !%Section ClassicalDynamics
    !%Description
    !% Charge of classical particle
    !%End
    call parse_variable(namespace, 'ClassicalParticleCharge', M_ONE, this%charge)
    call messages_print_var_value(stdout, 'ClassicalParticleCharge', this%charge)

    POP_SUB(charged_particle_init)
  end subroutine charged_particle_init

  ! ---------------------------------------------------------
  subroutine charged_particle_add_interaction_partner(this, partner)
    class(charged_particle_t), target, intent(inout) :: this
    class(system_abst_t),            intent(inout) :: partner

    class(interaction_coulomb_force_t), pointer :: coulomb_force
    type(interaction_coulomb_force_t) :: coulomb_force_t

    PUSH_SUB(charged_particle_add_interaction_partner)

    call this%classical_particle_t%add_interaction_partner(partner)
    if (partner%has_interaction(coulomb_force_t)) then
      coulomb_force => interaction_coulomb_force_t(this%space%dim, partner)
      this%quantities(POSITION)%required = .true.
      this%quantities(CHARGE)%required = .true.
      coulomb_force%system_charge => this%charge
      ! coulomb_force%system_pos  => this%pos
      call this%interactions%add(coulomb_force)
    end if

    POP_SUB(charged_particle_add_interaction_partner)
  end subroutine charged_particle_add_interaction_partner

  ! ---------------------------------------------------------
  logical function charged_particle_has_interaction(this, interaction) result(particle_has_interaction)
    class(charged_particle_t), intent(in) :: this
    class(interaction_abst_t), intent(in) :: interaction

    logical classical_particle_has
    logical charged_particle_has

    PUSH_SUB(charged_particle_has_interaction)

    classical_particle_has = this%classical_particle_t%has_interaction(interaction)
    select type (interaction)
    type is (interaction_coulomb_force_t)
      charged_particle_has = .true.
    class default
      charged_particle_has = .false.
    end select

    particle_has_interaction = classical_particle_has .or. charged_particle_has

    POP_SUB(charged_particle_has_interaction)
  end function charged_particle_has_interaction

  ! ---------------------------------------------------------
  subroutine charged_particle_initial_conditions(this, from_scratch)
    class(charged_particle_t), intent(inout) :: this
    logical,                 intent(in)    :: from_scratch

    PUSH_SUB(charged_particle_initial_conditions)

    call this%classical_particle_t%initial_conditions(from_scratch)

    POP_SUB(charged_particle_initial_conditions)
  end subroutine charged_particle_initial_conditions

  ! ---------------------------------------------------------
  subroutine charged_particle_do_td(this, operation)
    class(charged_particle_t), intent(inout) :: this
    integer,                 intent(in)    :: operation

    PUSH_SUB(charged_particle_do_td)

    call this%classical_particle_t%do_td_operation(operation)

    POP_SUB(charged_particle_do_td)
  end subroutine charged_particle_do_td

  ! ---------------------------------------------------------
  logical function charged_particle_is_tolerance_reached(this, tol) result(converged)
    class(charged_particle_t),   intent(in)    :: this
    FLOAT,                     intent(in)    :: tol

    PUSH_SUB(charged_particle_is_tolerance_reached)

    converged = this%classical_particle_t%is_tolerance_reached(tol)

    POP_SUB(charged_particle_is_tolerance_reached)
   end function charged_particle_is_tolerance_reached

   ! ---------------------------------------------------------
   subroutine charged_particle_store_current_status(this)
     class(charged_particle_t),   intent(inout)    :: this

     PUSH_SUB(charged_particle_store_current_status)

     call this%classical_particle_t%store_current_status()

     POP_SUB(charged_particle_store_current_status)
   end subroutine charged_particle_store_current_status

  ! ---------------------------------------------------------
  subroutine charged_particle_iteration_info(this)
    class(charged_particle_t), intent(in) :: this

    PUSH_SUB(charged_particle_iteration_info)

    call this%classical_particle_t%iteration_info()

    POP_SUB(charged_particle_iteration_info)
  end subroutine charged_particle_iteration_info

  ! ---------------------------------------------------------
  subroutine charged_particle_output_start(this)
    class(charged_particle_t), intent(inout) :: this

    PUSH_SUB(charged_particle_output_start)

    call this%classical_particle_t%output_start()

    POP_SUB(charged_particle_output_start)
  end subroutine charged_particle_output_start

  ! ---------------------------------------------------------
  subroutine charged_particle_output_finish(this)
    class(charged_particle_t), intent(inout) :: this

    PUSH_SUB(charged_particle_output_finish)

    call this%classical_particle_t%output_finish()

    POP_SUB(charged_particle_output_finish)
  end subroutine charged_particle_output_finish

  ! ---------------------------------------------------------
  subroutine charged_particle_output_write(this, iter)
    class(charged_particle_t), intent(inout) :: this
    integer,                 intent(in)    :: iter

    PUSH_SUB(charged_particle_output_write)

    call this%classical_particle_t%output_write(iter)

    POP_SUB(charged_particle_output_write)
  end subroutine charged_particle_output_write

  ! ---------------------------------------------------------
  subroutine charged_particle_update_quantity(this, iq, clock)
    class(charged_particle_t),   intent(inout) :: this
    integer,                   intent(in)    :: iq
    class(clock_t),            intent(in)    :: clock

    PUSH_SUB(charged_particle_update_quantity)

    call this%classical_particle_t%update_quantity(iq, clock)

    POP_SUB(charged_particle_update_quantity)
  end subroutine charged_particle_update_quantity

 ! ---------------------------------------------------------
 logical function charged_particle_update_exposed_quantity(this, iq, clock) result(updated)
    class(charged_particle_t),   intent(inout) :: this
    integer,                   intent(in)    :: iq
    class(clock_t),            intent(in)    :: clock

    PUSH_SUB(charged_particle_update_exposed_quantity)

    updated = this%classical_particle_t%update_exposed_quantity(iq, clock)

    POP_SUB(charged_particle_update_exposed_quantity)
  end function charged_particle_update_exposed_quantity

  ! ---------------------------------------------------------
  subroutine charged_particle_set_pointers_to_interaction(this, inter)
    class(charged_particle_t), target,  intent(inout) :: this
    class(interaction_abst_t),        intent(inout) :: inter

    PUSH_SUB(charged_particle_set_pointers_to_interaction)

    call this%classical_particle_t%set_pointers_to_interaction(inter)
    select type(inter)
    type is(interaction_coulomb_force_t)
      this%quantities(POSITION)%required = .true.
      this%quantities(CHARGE)%required = .true.
      ! inter%partner_pos => this%pos
      inter%partner_charge => this%charge
    class default
      message(1) = "Unsupported interaction."
      call messages_fatal(1)
    end select

    POP_SUB(charged_particle_set_pointers_to_interaction)
  end subroutine charged_particle_set_pointers_to_interaction

  ! ---------------------------------------------------------
  subroutine charged_particle_update_interactions_start(this)
    class(charged_particle_t), intent(inout) :: this

    PUSH_SUB(charged_particle_update_interactions_start)

    call this%classical_particle_t%update_interactions_start()

    POP_SUB(charged_particle_update_interactions_start)
  end subroutine charged_particle_update_interactions_start

  ! ---------------------------------------------------------
  subroutine charged_particle_update_interactions_finish(this)
    class(charged_particle_t), intent(inout) :: this

    type(interaction_iterator_t) :: iter

    PUSH_SUB(charged_particle_update_interactions_finish)

    call this%classical_particle_t%update_interactions_finish()
    call iter%start(this%interactions)
    do while (iter%has_next())
      select type (interaction => iter%get_next_interaction())
      type is (interaction_coulomb_force_t)
        ! this%tot_force(1:this%space%dim) = this%tot_force(1:this%space%dim) + interaction%force(1:this%space%dim)
      class default
        message(1) = "Unknown interaction by the charged particle " + this%namespace%get()
        call messages_fatal(1)
      end select
    end do

    POP_SUB(charged_particle_update_interactions_finish)
  end subroutine charged_particle_update_interactions_finish

  ! ---------------------------------------------------------
  subroutine charged_particle_finalize(this)
    type(charged_particle_t), intent(inout) :: this

    PUSH_SUB(charged_particle_finalize)
    POP_SUB(charged_particle_finalize)
  end subroutine charged_particle_finalize

end module system_charged_particle_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
