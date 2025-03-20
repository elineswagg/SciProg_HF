! Program made by Eline van de Voort, from a forked skeleton of Luuk Vischer
! This program calculates the Hartree-Fock energy for a given molecule
! The molecule is defined in the define_molecule subroutine
! The atomic orbital basis is defined in the define_basis subroutine
! The Hartree-Fock energy is calculated in the Fock_Matrix subroutine
! The results are written to a file in the Output_txt subroutine

program HartreeFock
  use ao_basis
  use compute_integrals
  use diagonalization
  use Input_Output
  use Hatree_Fock
  use molecular_structure

  implicit none

! Variable containing the molecular structure
  type(molecular_structure_t) :: molecule

! Variable containing the atomic orbital basis
  type(basis_set_info_t)      :: ao_basis

! Variables for molecular data
  integer                      :: number_electrons, n_AO, n_occ,  n_atoms
  integer                      :: mu, nu, i, kappa, lambda, j
  integer, allocatable         :: atom(:) 
  real(8)                      :: E_HF, E_repulsion
  real(8), allocatable         :: coord(:,:), ao_integrals(:,:,:,:)
  real(8), allocatable         :: F(:,:), V(:,:), T(:,:), S(:,:), C(:,:), eps(:), D(:,:)

! Definition of the molecule
  call define_molecule(molecule, atom, coord, n_atoms, number_electrons)
  
  print*, 'Coordinate matrix:'
  call PrintMatrix(coord)

! Definition of the GTOs
  call define_basis(ao_basis, atom, coord, n_atoms)
  n_AO = ao_basis%nao
  n_occ = number_electrons / 2

! Allocate and compute the integrals for kinetic energy, potential energy and overlap matrix
  allocate(S(n_AO, n_AO))
  call compute_1e_integrals("OVL", ao_basis, ao_basis, S)

  allocate(T(n_AO, n_AO))
  call compute_1e_integrals("KIN", ao_basis, ao_basis, T)

  allocate(V(n_AO, n_AO))
  call compute_1e_integrals("POT", ao_basis, ao_basis, V, molecule)

! Initialize the Fock matrix (core Hamiltonian)
  allocate(F(n_AO, n_AO))
  F = T - V

! Allocate space for orbital coefficients, density matrix, atomic orbital integrals, orbital energies
  allocate(C(n_AO, n_AO))
  allocate(D(n_AO, n_AO))
  allocate(ao_integrals(n_AO, n_AO, n_AO, n_AO))
  allocate (eps(n_AO))

! Call Fock_Matrix subroutine from Hatree_Fock module
  call Fock_Matrix(molecule, ao_basis, E_HF, E_repulsion, F, V, T, S, C, eps, D, ao_integrals, n_AO, n_occ, n_atoms)

! Output the results to a file
  call Output_txt(E_HF, E_repulsion, n_atoms, atom, coord)

  print*, "Number of atomic orbitals: ", n_AO
  Print*, "Nuclear repulsion Energy: ", E_repulsion, " Hartree"
  print*, "Final Hartree-Fock Energy: ", E_HF, " Hartree"

end program HartreeFock