! This module contains the subroutines for reading input and writing output
! Written by Eline van de Voort, Vrije Universiteit Amsterdam, march 2025

module Input_Output
   implicit none
   private
   public define_basis, define_molecule, PrintMatrix, Output_txt

   contains

      subroutine define_molecule(molecule, atom, coord, number_molecules, number_electrons)
         use molecular_structure
         type(molecular_structure_t), intent(inout) :: molecule
         integer, intent(out)                       :: number_molecules, number_electrons
         integer, allocatable, intent(out)          :: atom(:)
         real(8), allocatable, intent(out)          :: coord(:,:)
         character (len=100)                        :: filename, dummy
         integer                                    :: i, info
         real(8), allocatable                       :: charge(:)
         real(8), dimension(3)                      :: temp_coord
         character(len=2)                           :: atom_str

      ! Open the file containing the molecule specifications
         filename = "molecule.txt" 
         open(unit=10,file= filename,status='old')
  
      ! Read the number of atoms in the molecule
         read(10,*) dummy
         read(10,*) dummy  
         read(10,*) dummy, number_molecules
         read(10,*) dummy, number_electrons
         read(10,*) dummy
      
      ! Allocate the arrays
         allocate(charge(number_molecules))
         allocate(coord(3, number_molecules))
         allocate(atom(number_molecules))
         
     ! Read the atom names, charges and coordinates
         do i = 1, number_molecules
            read(10, *, iostat=info) atom_str, charge(i), temp_coord
   
            if (info /= 0) then
               print*, "Error reading data at row:", i
               stop "Error reading molecule data"
           end if
      
      ! Convert atom name to integer
            if (atom_str == "H") then
               atom(i) = 1
            else
               atom(i) = 0
            end if
         
      ! Should be 3, number_molecules matrix to ensure that add_atoms_to_molecule works
            coord(:, i) = temp_coord
         end do
         close(10)
         print*, "Molecule read successfully"

      ! Add the atoms to the molecule
         call add_atoms_to_molecule(molecule, charge, coord)
      end subroutine
   
      subroutine define_basis(ao_basis, atom, coord, number_molecules)
         use ao_basis
         implicit none
         integer, intent(in)                        :: number_molecules
         real(8), intent(in)                        :: coord(:,:)
         integer, intent(in)                        :: atom(:)
         type(basis_set_info_t), intent(inout)      :: ao_basis
         integer                                    :: i

         do i = 1, number_molecules
            if (atom(i) == 1) then
               ! H:  3 uncontracted s-fun:       l  coord        exp
               call add_shell_to_basis(ao_basis, 0, coord(:,i), 0.1D0)
               call add_shell_to_basis(ao_basis, 0, coord(:,i), 1.D0)
               call add_shell_to_basis(ao_basis, 0, coord(:,i), 3.D0)
            else
               ! Non-Hydrogen elements:
               ! 5 uncontracted s-fun:           l  coord        exp
               call add_shell_to_basis(ao_basis, 0, coord(:,i), 0.1D0)
               call add_shell_to_basis(ao_basis, 0, coord(:,i), 0.35D0)
               call add_shell_to_basis(ao_basis, 0, coord(:,i), 1.D0)
               call add_shell_to_basis(ao_basis, 0, coord(:,i), 3.D0)
               call add_shell_to_basis(ao_basis, 0, coord(:,i), 10.D0)

               ! 3 uncontracted p-fun:           l  coord        exp
               call add_shell_to_basis(ao_basis, 1, coord(:,i), 0.2D0)
               call add_shell_to_basis(ao_basis, 1, coord(:,i), 1.D0)
               call add_shell_to_basis(ao_basis, 1, coord(:,i), 5.D0)

               ! 1 uncontracted d-fun:           l  coord        exp
               call add_shell_to_basis(ao_basis, 2, coord(:,i), 1.D0)
            end if
         end do
      end subroutine define_basis

      subroutine PrintMatrix(array)
         real(8) :: array(:,:)
         integer :: i
     
      ! Loop through each row of the matrix
         do i = 1, size(array, 1)
             print *, array(i,:)
         end do
     end subroutine PrintMatrix

     subroutine Output_txt(E_HF, E_repulsion, n_atoms, coord, cycles, atom, n_AO, n_occ, eps)
      use ao_basis
      implicit none
      real(8), intent(in)                     :: E_HF, E_repulsion, coord(3, n_atoms), eps(n_AO)
      type(basis_set_info_t)                  :: ao_basis
      type(basis_func_info_t)                 :: gtos
      integer, intent(in)                     :: n_atoms, cycles, atom(n_atoms), n_AO, n_occ
      integer                                 :: i, info, j
      real(8)                                 :: charge(n_atoms), temp_coord(3), tolerance
      logical                                 :: exists
      character(len=100)                      :: dummy, atom_str
      real(8)                                 :: homo_energy, homo_minus1_energy, lumo_energy, homo_lumo_gap, homo1_lumo_gap

      ! Delete the file if it exists
      inquire(file="hartree_fock_energy.txt", exist=exists)
      if (exists) then
         call execute_command_line("rm hartree_fock_energy.txt")
      end if

      ! Open a file to write the Hartree-Fock energy, overwriting if it exists
      open(unit=10, file="hartree_fock_energy.txt", status="replace", action="write")
      ! Write a first line to the file
      write(10, '(A)') "! This file contains the final Hartree-Fock energy"

      write(10, '(A)') " "

      ! Write the coordinates of the molecule to the file
      write(10, '(A)') "Coordinates of molecule:"
      open(unit=11, file="molecule.txt", action="read")
      do i = 1, 5
         read(11, *) dummy
      end do

      do i = 1, n_atoms
         read(11, *, iostat=info) atom_str, charge(i), temp_coord
         if (info /= 0) then
         print*, "Error reading data at row:", i
         stop "Error reading molecule data"
         end if
         write(10, '(A5, 3F12.6)') trim(atom_str), coord(1, i), coord(2, i), coord(3, i)
      end do
      close(11)

      write(10, '(A)') " "

      ! Write the number of atomic orbitals, occupied orbitals, cycles and energies to the file
      write(10, '(A, I5, A)') "Number of atoms:               ", n_atoms, " atoms"
      write(10, '(A, I5, A)') "Number of atomic orbitals:      ", n_AO, " orbitals"
      write(10, '(A, I5, A)') "Number of occupied orbitals:   ", n_occ, " orbitals"
      write(10, "(A, I5, A)") "Number of cycles:               ", cycles, " cycles"
      write(10, '(A, F6.2, A)') "Nuclear repulsion Energy:         ", E_repulsion, " Hartree"
      write(10, '(A, F6.2, A)') "Final Hartree-Fock Energy:        ", E_HF, " Hartree"

      ! Write the HOMO, LUMO and HOMO-1 energies to file
      write(10, '(A)') " "
      write(10, '(A)') "HOMO, LUMO and HOMO-1 energies:  "

      homo_energy = eps(n_occ)

      ! take into accound degeracy
      tolerance = 0.01
      if (abs(homo_energy - eps(n_occ - 1)) < tolerance) then
         if (abs(homo_energy - eps(n_occ - 2)) < tolerance) then
         homo_minus1_energy = eps(n_occ - 3)
         else
         homo_minus1_energy = eps(n_occ - 2)
         end if
      else
         homo_minus1_energy = eps(n_occ - 1)
      end if

      if (abs(homo_energy - eps(n_occ + 1)) < tolerance) then
         if (abs(homo_energy - eps(n_occ + 2)) < tolerance) then
         lumo_energy = eps(n_occ + 3)
         else
         lumo_energy = eps(n_occ + 2)
         end if
      else
         lumo_energy = eps(n_occ + 1)
      end if

      write(10, "(A, F6.2, A)") 'HOMO energy:              ', homo_energy, " Hatree"
      write(10, "(A, F6.2, A)") 'HOMO-1 energy:            ', homo_minus1_energy, " Hatree"
      write(10, "(A, F6.2, A)") 'LUMO energy:              ', lumo_energy, " Hatree"
      
      homo_lumo_gap = lumo_energy - homo_energy
      homo1_lumo_gap = lumo_energy - homo_minus1_energy
  
      write(10, "(A, F6.2, A)") "Ionization Potential (IP):", -homo_energy, " Hatree"
      write(10, "(A, F6.2, A)") "Electron Affinity (EA):   ", lumo_energy, " Hatree"
      write(10, "(A, F6.2, A)") "HOMO-LUMO Gap:            ",  homo_lumo_gap, " Hatree"
      write(10, "(A, F6.2, A)") "HOMO-1 to LUMO Gap:       ", homo1_lumo_gap, " Hatree"

      ! Loop through each basis function and print its details
      write(10, '(A)') " "
      write(10, '(A)') "  ----------------------------------------------------"
      write(10, '(A)') "Basis set details:"
      do i = 1, n_atoms
         if (atom(i) == 1) then
         write(10, '(A)') "H:  3 uncontracted s-functions:"
         write(10, '(A, 3F12.6)') "  l=0, coordinates=", coord(1, i), coord(2, i), coord(3, i)
         write(10, '(A, F12.6)') "    exp=", 0.1D0
         write(10, '(A, F12.6)') "    exp=", 1.D0
         write(10, '(A, F12.6)') "    exp=", 3.D0
         else
         write(10, '(A)') " "
         write(10, '(A)') "Non-Hydrogen element:"
         write(10, '(A, 3F12.6)') "  l=0, coordinates=", coord(1, i), coord(2, i), coord(3, i)
         write(10, '(A, F12.6)') "    exp=", 0.1D0
         write(10, '(A, F12.6)') "    exp=", 0.35D0
         write(10, '(A, F12.6)') "    exp=", 1.D0
         write(10, '(A, F12.6)') "    exp=", 3.D0
         write(10, '(A, F12.6)') "    exp=", 10.D0

         write(10, '(A)') "  l=1, "
         write(10, '(A, F12.6)') "    exp=", 0.2D0
         write(10, '(A, F12.6)') "    exp=", 1.D0
         write(10, '(A, F12.6)') "    exp=", 5.D0

         write(10, '(A)') "  l=2, "
         write(10, '(A, F12.6)') "    exp=", 1.D0
         end if
      end do
      write(10, '(A)') "  ----------------------------------------------------"
      write(10, '(A)') " "


      ! Close the file
      close(10)
     end subroutine Output_txt

end module Input_Output
