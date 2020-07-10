!
! Copyright (C) 2001-2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
FUNCTION set_hubbard_n( psd ) RESULT( hubbard_n )
  !---------------------------------------------------------------------------
  !
  USE io_global, ONLY : stdout
  !
  IMPLICIT NONE
  !
  INTEGER                      :: hubbard_n
  CHARACTER(LEN=2), INTENT(IN) :: psd
  !
  !
  SELECT CASE( TRIM(ADJUSTL(psd)) )
     !
     ! ... transition metals, 4-th row 
     !
     CASE( 'Sc', 'Ti', 'V',  'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn')  
         hubbard_n = 3
     !
     !  ... transition metals, 5-th  row
     !
     CASE( 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd') 
         hubbard_n = 4
     ! 
     ! ... transition metals, 6-th  row 
     !
     CASE( 'Lu', 'Hf', 'Ta', 'W',  'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg')
     !
         hubbard_n = 5  
     !
     !
     ! ... rare earths (lanthanoid) 
     !
     CASE('La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb' ) 
     !  
         hubbard_n = 4 
     !  ... rare earths (actinoids )
     CASE('Ac', 'Th','Pa','U', 'Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No')
        !
         hubbard_n = 5
        !
     !
     ! ... other elements
     !
     CASE( 'H' )
        !
        hubbard_n =  1
        !
     CASE( 'Li', 'Be', 'B', 'C', 'N', 'O','F')
        !
        hubbard_n =  2
        !
     CASE( 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ga') 
        !
        hubbard_n =  3
        !
     CASE ( 'K', 'Ca', 'Ge','As','Se','Br','In' )

        ! 
        hubbard_n = 4
        !
     CASE ( 'Rb','Sr', 'Sn', 'Sb', 'Te', 'I')

        !
        hubbard_n = 5   
        !
     CASE ( 'Cs','Ba', 'Tl', 'Pb', 'Bi', 'Po', 'At')
        !
        hubbard_n = 6
        !
     CASE DEFAULT
        !
        hubbard_n = -1
        !
        WRITE( stdout, '(/,"psd = ",A,/)' ) psd
        !
        CALL errore( 'set_hubbard_l', 'pseudopotential not yet inserted', 1 )
        !
  END SELECT
  !
  RETURN  
  !
END FUNCTION set_Hubbard_n
