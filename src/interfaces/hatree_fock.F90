! Hatree-Fock module
! This module computes the final Hartree-Fock energy by adding the electronic energy to the nuclear repulsion energy.
! Written by Eline van de Voort, Vrije Universiteit Amsterdam, march 2025

module Hatree_Fock
  implicit none
  private
  public Fock_Matrix

  contains
    subroutine Fock_Matrix(molecule, ao_basis, E_HF, E_repulsion, F, V, T, S, C, eps, D, ao_integrals, n_AO, n_occ, n_atoms)
      use molecular_structure
      use ao_basis
      use compute_integrals
      use diagonalization
      implicit none

      type(molecular_structure_t), intent(in) :: molecule
      type(basis_set_info_t), intent(in)      :: ao_basis
      real(8), intent(out)                    :: E_HF, E_repulsion
      integer, intent(in)                     :: n_AO, n_occ, n_atoms
      integer                                 :: lambda, kappa, i, j, mu, nu
      real(8)                                 :: E_HF_new, delta_E, r
      real(8), allocatable, intent(inout)     :: ao_integrals(:,:,:,:)
      real(8), allocatable, intent(inout)     :: F(:,:), V(:,:), T(:,:), S(:,:), C(:,:), eps(:), D(:,:)

    ! Initialize the density matrix and the Hartree-Fock energy
      D = 0.D0
      E_HF = 0.D0

    ! SCF iteration loop
      do
        call solve_genev(F, S, C, eps)

      ! Construct the density matrix
        do lambda = 1, n_AO
          do kappa = 1, n_AO
            D(kappa, lambda) = sum(C(kappa, 1:n_occ) * C(lambda, 1:n_occ))
          end do
        end do

      ! Generate two-electron integrals
        call generate_2int(ao_basis, ao_integrals)

      ! Update the Fock matrix 
        do kappa = 1, n_AO
          do lambda = 1, n_AO
            F(kappa, lambda) = T(kappa, lambda) - V(kappa, lambda)
            do mu = 1, n_AO
              do nu = 1, n_AO
                F(kappa, lambda) = F(kappa, lambda) + D(mu, nu) * (2.0D0 * ao_integrals(kappa, lambda, mu, nu) - ao_integrals(kappa, nu, mu, lambda))
              end do
            end do
          end do
        end do

      ! Compute the Hartree-Fock energy before checking convergence
        E_HF_new = 0.D0
        do kappa = 1, n_AO
          do lambda = 1, n_AO
            E_HF_new = E_HF_new + (D(kappa, lambda) * (T(kappa, lambda) - V(kappa, lambda) + F(kappa, lambda)))
          end do
        end do
      
      ! Check convergence
        delta_E = abs(E_HF_new - E_HF)
        if (delta_E < 1.0e-9) then
          exit
        end if
        E_HF = E_HF_new
      end do

      deallocate(eps)

    ! Compute nuclear repulsion energy
      E_repulsion = 0.D0
      do i = 1, n_atoms
        do j = i+1, n_atoms
          r = sqrt(sum((molecule%coord(:,i) - molecule%coord(:,j))**2))
          E_repulsion = E_repulsion + (molecule%charge(i) * molecule%charge(j) / r )
        end do
      end do

    ! Final Hartree-Fock energy
      E_HF = E_HF + E_repulsion

    end subroutine Fock_Matrix

end module Hatree_Fock