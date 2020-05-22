!! Copyright (C) 2020 Heiko Appel
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

module interaction_lorentz_force_oct_m
  use clock_oct_m
  use global_oct_m
  use interaction_abst_oct_m
  use messages_oct_m
  use profiling_oct_m
  use quantity_oct_m
  use system_abst_oct_m

  implicit none

  private
  public ::                &
    interaction_lorentz_force_t

  type, extends(interaction_abst_t) :: interaction_lorentz_force_t
    integer :: dim

    FLOAT :: force(MAX_DIM)

    class(system_abst_t), pointer :: partner

    FLOAT, pointer :: system_pos(:)
    FLOAT, pointer :: system_vel(:)
    FLOAT, pointer :: system_charge

    FLOAT, pointer :: partner_E_field(:)
    FLOAT, pointer :: partner_B_field(:)

  contains
    procedure :: update => interaction_lorentz_force_update
    procedure :: end => interaction_lorentz_force_end
    final :: interaction_lorentz_force_finalize
  end type interaction_lorentz_force_t

  interface interaction_lorentz_force_t
    module procedure interaction_lorentz_force_init
  end interface interaction_lorentz_force_t

contains

  ! ---------------------------------------------------------
  function interaction_lorentz_force_init(dim, partner) result(this)
    integer,                      intent(in)    :: dim
    class(system_abst_t), target, intent(inout) :: partner
    class(interaction_lorentz_force_t), pointer :: this

    PUSH_SUB(interaction_lorentz_force_init)

    SAFE_ALLOCATE(this)

    this%dim = dim
    this%partner => partner

    ! Gravity interaction needs two quantities from each system: the position and the mass
    this%n_system_quantities = 3
    this%n_partner_quantities = 2
    SAFE_ALLOCATE(this%system_quantities(this%n_system_quantities))
    SAFE_ALLOCATE(this%partner_quantities(this%n_partner_quantities))
    this%system_quantities(1) = POSITION
    this%system_quantities(2) = VELOCITY
    this%system_quantities(3) = CHARGE
    this%partner_quantities(1) = E_FIELD
    this%partner_quantities(2) = B_FIELD

    call partner%set_pointers_to_interaction(this)

    POP_SUB(interaction_lorentz_force_init)
  end function interaction_lorentz_force_init

  ! ---------------------------------------------------------
  logical function interaction_lorentz_force_update(this, clock) result(updated)
    class(interaction_lorentz_force_t), intent(inout) :: this
    class(clock_t),               intent(in)    :: clock

    logical :: allowed_to_update

    PUSH_SUB(interaction_lorentz_force_update)

    ! We should only try to update the interaction if it is not yet at the desired time
    ASSERT(.not. (this%clock == clock))

    allowed_to_update = this%partner%update_exposed_quantities(clock, this%n_partner_quantities, this%partner_quantities)

    if (allowed_to_update) then
      ! We can now compute the interaction from the updated pointers
      ASSERT(associated(this%system_pos))
      ASSERT(associated(this%system_vel))
      ASSERT(associated(this%system_charge))
      ASSERT(associated(this%partner_E_field))
      ASSERT(associated(this%partner_B_field))

      ! curl is defined only in 3D
      ASSERT(this%dim == 3)

      this%force(1) = this%system_charge * (this%partner_E_field(1) + &
        this%system_vel(2) * this%partner_B_field(3) - this%system_vel(3) * this%partner_B_field(2))
      this%force(2) = this%system_charge * (this%partner_E_field(2) + &
        this%system_vel(3) * this%partner_B_field(1) - this%system_vel(1) * this%partner_B_field(3))
      this%force(3) = this%system_charge * (this%partner_E_field(3) + &
        this%system_vel(1) * this%partner_B_field(2) - this%system_vel(2) * this%partner_B_field(1))

      ! Update was successful, so set new interaction time
      updated = .true.
      call this%clock%set_time(clock)

      if (debug%info) then
        write(message(1), '(a,a)') "Debug: -- Updated interaction with ", trim(this%partner%namespace%get())
        write(message(2), '(a,i3,a,i3)') "Debug: ---- clocks are ", clock%get_tick(), " and ", this%partner%clock%get_tick()
        call messages_info(2)
      end if
    else
      if (debug%info) then
        write(message(1), '(a,a)') "Debug: -- Cannot update yet the interaction with ", trim(this%partner%namespace%get())
        write(message(2), '(a,i3,a,i3)') "Debug: ---- clocks are ", clock%get_tick(), " and ", this%partner%clock%get_tick()
        call messages_info(2)
      end if
      updated = .false.
    end if

    POP_SUB(interaction_lorentz_force_update)
  end function interaction_lorentz_force_update

  ! ---------------------------------------------------------
  subroutine interaction_lorentz_force_end(this)
    class(interaction_lorentz_force_t), intent(inout) :: this

    PUSH_SUB(interaction_lorentz_force_end)

    nullify(this%partner)
    this%force = M_ZERO
    nullify(this%system_charge)
    nullify(this%system_pos)
    nullify(this%system_vel)
    nullify(this%partner_E_field)
    nullify(this%partner_B_field)

    call interaction_abst_end(this)

    POP_SUB(interaction_lorentz_force_end)
  end subroutine interaction_lorentz_force_end

  ! ---------------------------------------------------------
  subroutine interaction_lorentz_force_finalize(this)
    type(interaction_lorentz_force_t), intent(inout) :: this

    PUSH_SUB(interaction_lorentz_force_finalize)

    call this%end()

    POP_SUB(interaction_lorentz_force_finalize)
  end subroutine interaction_lorentz_force_finalize

end module interaction_lorentz_force_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
