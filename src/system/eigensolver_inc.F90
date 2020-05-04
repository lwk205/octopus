!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
!! Copyright (C) 2020 N. Tancogne-Dejean
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
 ! ---------------------------------------------------------
  subroutine X(eigensolver_run)(eigens, namespace, gr, st, hm, iter, conv, nstconv)
    type(eigensolver_t),      intent(inout) :: eigens
    type(namespace_t),        intent(in)    :: namespace
    type(grid_t),             intent(in)    :: gr
    type(states_elec_t),      intent(inout) :: st
    type(hamiltonian_elec_t), intent(inout) :: hm
    integer,                  intent(in)    :: iter
    logical,        optional, intent(out)   :: conv
    integer,        optional, intent(in)    :: nstconv !< Number of states considered for 
                                                   !< the convergence criteria

    integer :: maxiter, ik, ist, nstconv_
#ifdef HAVE_MPI
    logical :: conv_reduced
    integer :: outcount, lmatvec
    FLOAT, allocatable :: ldiff(:), leigenval(:)
    integer, allocatable :: lconv(:)
#endif

    PUSH_SUB(X(eigensolver_run))

    if(present(conv)) conv = .false.
    if(present(nstconv)) then 
      nstconv_ = nstconv
    else
      nstconv_ = st%nst
    end if

    eigens%matvec = 0

    if(mpi_grp_is_root(mpi_world) .and. eigensolver_has_progress_bar(eigens) .and. .not. debug%info) then
      call loct_progress_bar(-1, st%lnst*st%d%kpt%nlocal)
    end if

    ik_loop: do ik = st%d%kpt%start, st%d%kpt%end
     if(eigens%skip_finite_weight_kpoints.and. st%d%kweights(ik) > M_ZERO) cycle
      maxiter = eigens%es_maxiter

      !Following the paper of Kresse and Furthmueller, we perform the susbpace rotation before 
      !the CG or RMMDIIS steps.
      if(st%calc_eigenval) then
        if(eigens%converged(ik) == 0 .and. hm%theory_level /= INDEPENDENT_PARTICLES) then
          call X(subspace_diag)(eigens%sdiag, namespace, gr%mesh, st, hm, ik, st%eigenval(:, ik), eigens%diff(:, ik))
        end if
      end if

      select case(eigens%es_type)
      case(RS_CG_NEW)
        call X(eigensolver_cg2_new)(namespace, gr, st, hm, eigens%tolerance, maxiter, eigens%converged(ik), ik, eigens%diff(:, ik))
      case(RS_CG)
        call X(eigensolver_cg2)(namespace, gr, st, hm, hm%xc, eigens%pre, eigens%tolerance, maxiter, &
          eigens%converged(ik), ik, eigens%diff(:, ik), eigens%orthogonalize_to_all, &
          eigens%conjugate_direction, eigens%additional_terms, eigens%energy_change_threshold)
      case(RS_PLAN)
        call X(eigensolver_plan)(namespace, gr, st, hm, eigens%pre, eigens%tolerance, maxiter, eigens%converged(ik), ik, &
          eigens%diff(:, ik))
      case(RS_EVO)
        call X(eigensolver_evolution)(namespace, gr%mesh, st, hm, eigens%exponential_operator, eigens%tolerance, maxiter, &
          eigens%converged(ik), ik, eigens%diff(:, ik), tau = eigens%imag_time)
      case(RS_LOBPCG)
        call X(eigensolver_lobpcg)(namespace, gr, st, hm, eigens%pre, eigens%tolerance, maxiter, &
          eigens%converged(ik), ik, eigens%diff(:, ik), hm%d%block_size)
      case(RS_RMMDIIS)
        if(iter <= eigens%rmmdiis_minimization_iter) then
          maxiter = 2
          call X(eigensolver_rmmdiis_min)(namespace, gr, st, hm, eigens%pre, maxiter, eigens%converged(ik), ik)
        else
          call X(eigensolver_rmmdiis)(namespace, gr, st, hm, eigens%pre, eigens%tolerance, maxiter, &
            eigens%converged(ik), ik, eigens%diff(:, ik))
        end if
      case(RS_PSD)
        call X(eigensolver_rmmdiis_min)(namespace, gr, st, hm, eigens%pre, maxiter, eigens%converged(ik), ik)
      end select
        
      if(st%calc_eigenval .and. .not. eigens%folded_spectrum) then
        ! recheck convergence after subspace diagonalization, since states may have reordered
        eigens%converged(ik) = 0
        do ist = 1, st%nst
          if(eigens%diff(ist, ik) < eigens%tolerance) then
            eigens%converged(ik) = ist
          else
            exit
          end if
        end do
      end if
      
      eigens%matvec = eigens%matvec + maxiter
    end do ik_loop

    if(mpi_grp_is_root(mpi_world) .and. eigensolver_has_progress_bar(eigens) .and. .not. debug%info) then
      write(stdout, '(1x)')
    end if

    if(present(conv)) conv = all(eigens%converged(st%d%kpt%start:st%d%kpt%end) >= nstconv_)

#ifdef HAVE_MPI
    if(st%d%kpt%parallel) then
      if(present(conv)) then
        call MPI_Allreduce(conv, conv_reduced, 1, MPI_LOGICAL, MPI_LAND, st%d%kpt%mpi_grp%comm, mpi_err)
        conv = conv_reduced
      end if

      lmatvec = eigens%matvec
      call MPI_Allreduce(lmatvec, eigens%matvec, 1, MPI_INTEGER, MPI_SUM, st%d%kpt%mpi_grp%comm, mpi_err)

      SAFE_ALLOCATE(lconv(1:st%d%kpt%nlocal))
      lconv(1:st%d%kpt%nlocal) = eigens%converged(st%d%kpt%start:st%d%kpt%end)
      call lmpi_gen_allgatherv(st%d%kpt%nlocal, lconv, outcount, eigens%converged, st%d%kpt%mpi_grp)
      ASSERT(outcount == st%d%nik)
      SAFE_DEALLOCATE_A(lconv)

      ! every node needs to know all eigenvalues (and diff)
      SAFE_ALLOCATE(ldiff(1:st%d%kpt%nlocal))
      SAFE_ALLOCATE(leigenval(1:st%d%kpt%nlocal))
      do ist = st%st_start, st%st_end
        ldiff(1:st%d%kpt%nlocal) = eigens%diff(ist, st%d%kpt%start:st%d%kpt%end)
        leigenval(1:st%d%kpt%nlocal) = st%eigenval(ist, st%d%kpt%start:st%d%kpt%end)
        call lmpi_gen_allgatherv(st%d%kpt%nlocal, ldiff, outcount, &
          eigens%diff(ist, :), st%d%kpt%mpi_grp)
        ASSERT(outcount == st%d%nik)
        call lmpi_gen_allgatherv(st%d%kpt%nlocal, leigenval, outcount, &
          st%eigenval(ist, :), st%d%kpt%mpi_grp)
        ASSERT(outcount == st%d%nik)
      end do
      SAFE_DEALLOCATE_A(ldiff)
      SAFE_DEALLOCATE_A(leigenval)
    end if
#endif

    POP_SUB(X(eigensolver_run))
  end subroutine X(eigensolver_run)
