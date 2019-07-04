!! Copyright (C) 2019 M. Oliveira
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

module linked_list_oct_m
  use global_oct_m
  use list_node_oct_m
  use messages_oct_m
  use profiling_oct_m
  implicit none

  private
  public :: linked_list_t

  type :: linked_list_t
    private
    class(list_node_t), pointer :: first_node => null()
    class(list_node_t), pointer :: last_node => null()
    class(list_node_t), pointer :: current_node => null()
  contains
    procedure :: add_node
    procedure :: first
    procedure :: next
    procedure :: current
    procedure :: rewind
    procedure :: has_more_values
    generic   :: add => add_node
    procedure :: finalizer
  end type linked_list_t

contains

  subroutine add_node(this, value)
    class(linked_list_t) :: this
    class(*) :: value
    class(list_node_t), pointer :: new_node

    PUSH_SUB(add_node)

    if (.not. associated(this%first_node)) then
      this%first_node => list_node(value, this%first_node)
      this%last_node => this%first_node
    else
      new_node => list_node(value, this%last_node%next())
      call this%last_node%set_next(new_node)
      this%last_node => new_node
    end if

    POP_SUB(add_node)
  end subroutine add_node

  function first(this)
    class(linked_list_t) :: this
    class(*), pointer    :: first

    PUSH_SUB(first)

    first => this%first_node%get()

    POP_SUB(first)
  end function first

  function current(this)
    class(linked_list_t) :: this
    class(*), pointer    :: current

    PUSH_SUB(current)

    current => this%current_node%get()

    POP_SUB(current)
  end function current

  subroutine next(this)
    class(linked_list_t) :: this

    PUSH_SUB(next)

    this%current_node => this%current_node%next()

    POP_SUB(next)
  end subroutine next

  logical function has_more_values(this)
    class(linked_list_t) :: this

    PUSH_SUB(has_more_values)

    has_more_values = associated(this%current_node)

    POP_SUB(has_more_values)
  end function has_more_values

  subroutine rewind(this)
    class(linked_list_t) :: this

    PUSH_SUB(rewind)

    this%current_node => this%first_node

    POP_SUB(rewind)
  end subroutine rewind

  subroutine finalizer(this)
    class(linked_list_t), intent(inout) :: this

    class(list_node_t), pointer :: next

    PUSH_SUB(finalizer)

    call this%rewind()
    do while (associated(this%current_node))
      next => this%current_node%next()
      ! No safe_deallocate macro here, as the call to sizeof() causes an
      ! internal compiler error with GCC 6.4.0
      deallocate(this%current_node)
      this%current_node => next
    end do

    POP_SUB(finalizer)
  end subroutine finalizer

end module linked_list_oct_m