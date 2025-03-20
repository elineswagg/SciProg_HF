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

     subroutine Output_txt(E_HF, E_repulsion, n_atoms, atom, coord)
      implicit none
      real(8), intent(in) :: E_HF, E_repulsion, coord(3, n_atoms)
      integer, intent(in) :: n_atoms, atom(n_atoms)
      integer             :: i
      logical             :: exists

      ! Delete the file if it exists
      inquire(file="hartree_fock_energy.txt", exist=exists)
      if (exists) then
         call execute_command_line("rm hartree_fock_energy.txt")
      end if

      ! Open a file to write the Hartree-Fock energy, overwriting if it exists
      open(unit=10, file="hartree_fock_energy.txt", status="replace", action="write")

      ! Write the Hartree-Fock energy to the file
      write(10, '(A)') "! This file contains the final Hartree-Fock energy"
      write(10, '(A)') "Coordinates of molecule:"
      do i = 1, n_atoms
         write(10, '(I5, 3F12.6)') atom(i), coord(1, i), coord(2, i), coord(3, i)
      end do
  
      write(10, '(A, F12.6)') "Nuclear repulsion Energy: ", E_repulsion
      write(10, '(A, F12.6)') "Final Hartree-Fock Energy: ", E_HF

      ! Close the file
      close(10)
     end subroutine Output_txt

end module Input_Output
