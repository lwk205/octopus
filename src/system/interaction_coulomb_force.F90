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

module interaction_coulomb_force_oct_m
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
    interaction_coulomb_force_t

  type, extends(interaction_abst_t) :: interaction_coulomb_force_t
    integer :: dim

    FLOAT :: force(MAX_DIM)

    class(system_abst_t), pointer :: partner

    FLOAT, pointer :: system_charge
    FLOAT, pointer :: system_pos(:)

    FLOAT, pointer :: partner_charge
    FLOAT, pointer :: partner_pos(:)

  contains
    procedure :: update => interaction_coulomb_force_update
    procedure :: end => interaction_coulomb_force_end
    final :: interaction_coulomb_force_finalize
  end type interaction_coulomb_force_t

  interface interaction_coulomb_force_t
    module procedure interaction_coulomb_force_init
  end interface interaction_coulomb_force_t

contains

  ! ---------------------------------------------------------
  function interaction_coulomb_force_init(dim, partner) result(this)
    integer,                      intent(in)    :: dim
    class(system_abst_t), target, intent(inout) :: partner
    class(interaction_coulomb_force_t), pointer       :: this

    PUSH_SUB(interaction_coulomb_force_init)

    SAFE_ALLOCATE(this)

    this%dim = dim
    this%partner => partner

    ! Gravity interaction needs two quantities from each system: the position and the mass
    this%n_system_quantities = 2
    this%n_partner_quantities = 2
    SAFE_ALLOCATE(this%system_quantities(this%n_system_quantities))
    SAFE_ALLOCATE(this%partner_quantities(this%n_partner_quantities))
    this%system_quantities(1) = POSITION
    this%system_quantities(2) = CHARGE
    this%partner_quantities(1) = POSITION
    this%partner_quantities(2) = CHARGE

    call partner%set_pointers_to_interaction(this)

    POP_SUB(interaction_coulomb_force_init)
  end function interaction_coulomb_force_init

  ! ---------------------------------------------------------
  logical function interaction_coulomb_force_update(this, clock) result(updated)
    class(interaction_coulomb_force_t), intent(inout) :: this
    class(clock_t),               intent(in)    :: clock

    logical :: allowed_to_update
    FLOAT, parameter :: eps_0 = CNST(625000)/(CNST(22468879468420441) * M_PI) ! vacuum permittivity in SI units
    FLOAT, parameter :: COULCONST = M_ONE/(M_FOUR * M_PI * eps_0) ! Coulomb constant in SI units
    FLOAT :: dist3

    PUSH_SUB(interaction_coulomb_force_update)

    ! We should only try to update the interaction if it is not yet at the desired time
    ASSERT(.not. (this%clock == clock))

    allowed_to_update = this%partner%update_exposed_quantities(clock, this%n_partner_quantities, this%partner_quantities)

    if (allowed_to_update) then
      ! We can now compute the interaction from the updated pointers
      ASSERT(associated(this%partner_pos))
      ASSERT(associated(this%system_pos))
      ASSERT(associated(this%partner_charge))
      ASSERT(associated(this%system_charge))

      dist3 = sum((this%partner_pos(1:this%dim) - this%system_pos(1:this%dim))**2)**(M_THREE/M_TWO)

      this%force(1:this%dim) = (this%partner_pos(1:this%dim) - this%system_pos(1:this%dim)) &
        / dist3 * (COULCONST * this%system_charge * this%partner_charge)

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

    POP_SUB(interaction_coulomb_force_update)
  end function interaction_coulomb_force_update

  ! ---------------------------------------------------------
  subroutine interaction_coulomb_force_end(this)
    class(interaction_coulomb_force_t), intent(inout) :: this

    PUSH_SUB(interaction_coulomb_force_end)

    nullify(this%partner)
    this%force = M_ZERO
    nullify(this%system_charge)
    nullify(this%system_pos)
    nullify(this%partner_charge)
    nullify(this%partner_pos)

    call interaction_abst_end(this)

    POP_SUB(interaction_coulomb_force_end)
  end subroutine interaction_coulomb_force_end

  ! ---------------------------------------------------------
  subroutine interaction_coulomb_force_finalize(this)
    type(interaction_coulomb_force_t), intent(inout) :: this

    PUSH_SUB(interaction_coulomb_force_finalize)

    call this%end()

    POP_SUB(interaction_coulomb_force_finalize)
  end subroutine interaction_coulomb_force_finalize

end module interaction_coulomb_force_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
