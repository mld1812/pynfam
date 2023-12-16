!------------------------------------------------------------------------------
!> This module contains a derived type for a 2x2 block matrix, where each quadrant
!> is itself a block matrix. It contains the methods for manipulating these
!> structures, in particular, multiplying three of them.
!>
!> To avoid explicit transposes and conjugations, we store the vectorized
!> matrix elements, the intended tranpose, and the intended sign as separate attributes
!> (in the same way dgemm does not require explicit transposes, just a string indicator)
!>
!> 6/17/21 - ISSUES:
!>     1. Methods that combine bbms should have more precise checks on
!>        block structure (e.g. add_bbm should not work unless strucutres
!>        are the same. If I implemented complete general block structure
!>        then the methods can just change the structure, but while
!>        1 block per row/col is enforced we cannot do this.
!>     2. I should use BLAS throughout to speedup these methods
!>     3. Should contract BBM be incorporating the signs of the subblocks?
!>     4. Checks on imag?
!>
!> @authors E.M. Ney, UNC Chapel Hill, 2019
!------------------------------------------------------------------------------
module type_bigblockmatrix
   use pnfam_logger
   use type_blockmatrix

   implicit none

   integer, parameter, private :: dp = kind(1.0d0)

   type bigblockmatrix
       type(blockmatrix) :: m11, m12, m21, m22
       character(1)      :: t11='n', t12='n', t21='n', t22='n'
       real(dp)          :: s11=+1d0,s12=+1d0,s21=+1d0,s22=+1d0
       logical           :: imag=.false.
   end type bigblockmatrix

contains

   !----------------------------------------------------------------------------
   ! Allocate any number of quadrants at once
   !----------------------------------------------------------------------------
   subroutine allocate_bbm(mat, a11,a12,a21,a22, n)
      implicit none
      type(bigblockmatrix), intent(inout) :: mat
      logical, intent(in) :: a11,a12,a21,a22
      integer, intent(in) :: n

      if (a11) call allocate_blockmatrix(mat%m11, n)
      if (a12) call allocate_blockmatrix(mat%m12, n)
      if (a21) call allocate_blockmatrix(mat%m21, n)
      if (a22) call allocate_blockmatrix(mat%m22, n)

   end subroutine

   !----------------------------------------------------------------------------
   ! If allocated block matrices, deallocate
   !----------------------------------------------------------------------------
   subroutine deallocate_bbm(mat)
      implicit none
      type(bigblockmatrix), intent(inout) :: mat

      if (allocated(mat%m11%elem)) call deallocate_blockmatrix(mat%m11)
      if (allocated(mat%m12%elem)) call deallocate_blockmatrix(mat%m12)
      if (allocated(mat%m21%elem)) call deallocate_blockmatrix(mat%m21)
      if (allocated(mat%m22%elem)) call deallocate_blockmatrix(mat%m22)

   end subroutine

   !----------------------------------------------------------------------------
   ! Set all (allocated) matrix elements to a constant i
   !----------------------------------------------------------------------------
   subroutine set_val_bbm(a,i)
      implicit none
      type(bigblockmatrix), intent(inout) :: a
      real(dp), intent(in) :: i

      logical :: a11,a12,a21,a22

      ! Check which quadrants are active
      a11 = allocated(a%m11%elem); a12 = allocated(a%m12%elem)
      a21 = allocated(a%m21%elem); a22 = allocated(a%m22%elem)

      if (a11) a%m11%elem = i
      if (a12) a%m12%elem = i
      if (a21) a%m21%elem = i
      if (a22) a%m22%elem = i

   end subroutine

   !----------------------------------------------------------------------------
   ! Indicated the intended sign of the quadrant
   !----------------------------------------------------------------------------
   subroutine set_sign_bbm(a,s11,s12,s21,s22)
      implicit none
      type(bigblockmatrix), intent(inout) :: a
      real(dp),intent(in) :: s11,s12,s21,s22

      a%s11 = sign(a%s11,s11)
      a%s12 = sign(a%s12,s12)
      a%s21 = sign(a%s21,s21)
      a%s22 = sign(a%s22,s22)

   end subroutine

   !----------------------------------------------------------------------------
   ! Indicate the intended tranpose of the quadrant
   !----------------------------------------------------------------------------
   subroutine set_trans_bbm(a,t11,t12,t21,t22)
      implicit none
      type(bigblockmatrix), intent(inout) :: a
      character(1),intent(in) :: t11,t12,t21,t22

      a%t11 = t11
      a%t12 = t12
      a%t21 = t21
      a%t22 = t22

   end subroutine

   !----------------------------------------------------------------------------
   ! Return the transpose of a set of transpose indicators
   !----------------------------------------------------------------------------
   subroutine switch_trans(t11,t12,t21,t22)
      implicit none
      character(1),intent(inout) :: t11,t12,t21,t22
      character(1), dimension(4) :: all_t
      integer :: i

      all_t = (/ t11,t12,t21,t22 /)
      do i=1,4
         if (all_t(i)=='n') then
             all_t(i)='t'
         else
             all_t(i)='n'
         end if
     end do
     t11 = all_t(1); t12 = all_t(2)
     t21 = all_t(3); t22 = all_t(4)

   end subroutine

   !----------------------------------------------------------------------------
   ! Transpose a bigblockmatrix
   !----------------------------------------------------------------------------
   subroutine transpose_bbm(a)
      implicit none
      type(bigblockmatrix), intent(inout) :: a

      type(blockmatrix) :: auxb
      real(dp) :: auxs

      auxb = a%m12; a%m12 = a%m21; a%m21 = auxb
      auxs = a%s12; a%s12 = a%s21; a%s21 = auxs
      call switch_trans(a%t11,a%t12,a%t21,a%t22)

   end subroutine

   !----------------------------------------------------------------------------
   ! Copy block structure of A to B
   !----------------------------------------------------------------------------
   subroutine copy_block_structure_bbm(a, b)
      implicit none
      type(bigblockmatrix), intent(in) :: a
      type(bigblockmatrix), intent(inout) :: b

      logical :: a11,a12,a21,a22
      logical :: b11,b12,b21,b22

      ! Check which quadrants are active
      a11 = allocated(a%m11%elem); a12 = allocated(a%m12%elem)
      a21 = allocated(a%m21%elem); a22 = allocated(a%m22%elem)
      b11 = allocated(b%m11%elem); b12 = allocated(b%m12%elem)
      b21 = allocated(b%m21%elem); b22 = allocated(b%m22%elem)

      ! WARN if a and b have mis-matched active quadrants?
      if (a11 .and. b11) then
         call copy_block_structure(a%m11, b%m11)
      end if
      if (a12 .and. b12) then
         call copy_block_structure(a%m12, b%m12)
      end if
      if (a21 .and. b21) then
         call copy_block_structure(a%m21, b%m21)
      end if
      if (a22 .and. b22) then
         call copy_block_structure(a%m22, b%m22)
      end if

   end subroutine

   !----------------------------------------------------------------------------
   ! Copy block structure of block matrix a to all quadrants of bigblock matrix B
   !----------------------------------------------------------------------------
   subroutine copyall_block_structure_bbm(a, b)
      implicit none
      type(blockmatrix), intent(in) :: a
      type(bigblockmatrix), intent(inout) :: b

      logical :: b11,b12,b21,b22

      ! Check which quadrants are active
      b11 = allocated(b%m11%elem); b12 = allocated(b%m12%elem)
      b21 = allocated(b%m21%elem); b22 = allocated(b%m22%elem)

      ! WARN if a and b have mis-matched active quadrants?
      if (b11) call copy_block_structure(a, b%m11)
      if (b12) call copy_block_structure(a, b%m12)
      if (b21) call copy_block_structure(a, b%m21)
      if (b22) call copy_block_structure(a, b%m22)

   end subroutine


   !----------------------------------------------------------------------------
   ! Multiply bigblock matrix A by a real scalar
   !----------------------------------------------------------------------------
   subroutine scalar_mult_bbm(a,s)
      implicit none
      type(bigblockmatrix), intent(inout) :: a
      real(dp), intent(in) :: s

      logical :: a11,a12,a21,a22

      ! Check which quadrants are active
      a11 = allocated(a%m11%elem); a12 = allocated(a%m12%elem)
      a21 = allocated(a%m21%elem); a22 = allocated(a%m22%elem)

      ! USE BLAS?
      if (a11) a%m11%elem = a%m11%elem*s
      if (a12) a%m12%elem = a%m12%elem*s
      if (a21) a%m21%elem = a%m21%elem*s
      if (a22) a%m22%elem = a%m22%elem*s

   end subroutine

   !----------------------------------------------------------------------------
   ! Elemental multiplication of two block matrices
   !----------------------------------------------------------------------------
   subroutine emult_bbm(a, b, ab)
      implicit none
      type(bigblockmatrix), intent(in) :: a, b
      type(bigblockmatrix), intent(inout) :: ab

      logical :: a11,a12,a21,a22
      logical :: b11,b12,b21,b22
      integer :: im

      ! Check which quadrants are non-zero
      a11 = allocated(a%m11%elem); a12 = allocated(a%m12%elem)
      a21 = allocated(a%m21%elem); a22 = allocated(a%m22%elem)
      b11 = allocated(b%m11%elem); b12 = allocated(b%m12%elem)
      b21 = allocated(b%m21%elem); b22 = allocated(b%m22%elem)

      ! Handle imaginary
      im = 0
      if (a%imag) im = im + 1
      if (b%imag) im = im + 1
      if (mod(im,2)==0) then
         ab%imag = .false.
      else
         ab%imag = .true.
      end if

      ! USE BLAS?
      if (a11.and.b11) then
         ab%m11%elem(:) = a%m11%elem*b%m11%elem
         if (im==2) ab%m11%elem = - ab%m11%elem
      end if
      if (a12.and.b12) then
         ab%m12%elem(:) = a%m12%elem*b%m12%elem
         if (im==2) ab%m12%elem = - ab%m12%elem
      end if
      if (a21.and.b21) then
         ab%m21%elem(:) = a%m21%elem*b%m21%elem
         if (im==2) ab%m21%elem = - ab%m21%elem
      end if
      if (a22.and.b22) then
         ab%m22%elem(:) = a%m22%elem*b%m22%elem
         if (im==2) ab%m22%elem = - ab%m22%elem
      end if

   end subroutine

   !----------------------------------------------------------------------------
   ! Elemental multiplication of two complex block matrices split to real and
   ! imaginary parts: (reab,imab) = (rea,ima)*(reb,imb) = (rea*reb-ima*imb, 
   !----------------------------------------------------------------------------
   subroutine complex_emult_bbm(rea, ima, reb, imb, reab, imab)
      implicit none
      type(bigblockmatrix), intent(in) :: rea, ima, reb, imb
      type(bigblockmatrix), intent(inout) :: reab, imab

      logical :: ra11,ra12,ra21,ra22, ia11,ia12,ia21,ia22
      logical :: rb11,rb12,rb21,rb22, ib11,ib12,ib21,ib22

      ! Check which quadrants are non-zero
      ra11 = allocated(rea%m11%elem); ra12 = allocated(rea%m12%elem)
      ra21 = allocated(rea%m21%elem); ra22 = allocated(rea%m22%elem)
      ia11 = allocated(ima%m11%elem); ia12 = allocated(ima%m12%elem)
      ia21 = allocated(ima%m21%elem); ia22 = allocated(ima%m22%elem)
      rb11 = allocated(reb%m11%elem); rb12 = allocated(reb%m12%elem)
      rb21 = allocated(reb%m21%elem); rb22 = allocated(reb%m22%elem)
      ib11 = allocated(imb%m11%elem); ib12 = allocated(imb%m12%elem)
      ib21 = allocated(imb%m21%elem); ib22 = allocated(imb%m22%elem)

      ! USE BLAS?
      if (ra11.and.ia11.and.rb11.and.ib11) then
         reab%m11%elem(:) = rea%m11%elem*reb%m11%elem - ima%m11%elem*imb%m11%elem
         imab%m11%elem(:) = ima%m11%elem*reb%m11%elem + imb%m11%elem*rea%m11%elem
      end if
      if (ra12.and.ia12.and.rb12.and.ib12) then
         reab%m12%elem(:) = rea%m12%elem*reb%m12%elem - ima%m12%elem*imb%m12%elem
         imab%m12%elem(:) = ima%m12%elem*reb%m12%elem + imb%m12%elem*rea%m12%elem
      end if
      if (ra21.and.ia21.and.rb21.and.ib21) then
         reab%m21%elem(:) = rea%m21%elem*reb%m21%elem - ima%m21%elem*imb%m21%elem
         imab%m21%elem(:) = ima%m21%elem*reb%m21%elem + imb%m21%elem*rea%m21%elem
      end if
      if (ra22.and.ia22.and.rb22.and.ib22) then
         reab%m22%elem(:) = rea%m22%elem*reb%m22%elem - ima%m22%elem*imb%m22%elem
         imab%m22%elem(:) = ima%m22%elem*reb%m22%elem + imb%m22%elem*rea%m22%elem
      end if

      reab%imag = .false.
      imab%imag = .true.

   end subroutine

   !----------------------------------------------------------------------------
   ! Add bigblock matrix A to B
   !----------------------------------------------------------------------------
   subroutine add_bbm(a,b)
      implicit none
      type(bigblockmatrix), intent(in) :: a
      type(bigblockmatrix), intent(inout) :: b

      logical :: a11,a12,a21,a22
      logical :: b11,b12,b21,b22
      integer :: ierr

      ierr = 0

      ! Check which quadrants are non-zero
      a11 = allocated(a%m11%elem); a12 = allocated(a%m12%elem)
      a21 = allocated(a%m21%elem); a22 = allocated(a%m22%elem)
      b11 = allocated(b%m11%elem); b12 = allocated(b%m12%elem)
      b21 = allocated(b%m21%elem); b22 = allocated(b%m22%elem)

      ! USE BLAS? WARN IF MISMATCHED ACTIVE BLOCKS?
      if (a11) then
         if (.not.b11) then
            b%m11 = a%m11
         else
            if (.not.equal_structure(a%m11,b%m11)) ierr = 1
            b%m11%elem = b%m11%elem + a%m11%elem
         end if
      end if
      if (a12) then
         if (.not.b12) then
            b%m12 = a%m12
         else
            if (.not.equal_structure(a%m12,b%m12)) ierr = 1
            b%m12%elem = b%m12%elem + a%m12%elem
         end if
      end if
      if (a21) then
         if (.not.b21) then
            b%m21 = a%m21
         else
            if (.not.equal_structure(a%m21,b%m21)) ierr = 1
            b%m21%elem = b%m21%elem + a%m21%elem
         end if
      end if
      if (a22) then
         if (.not.b22) then
            b%m22 = a%m22
         else
            if (.not.equal_structure(a%m22,b%m22)) ierr = 1
            b%m22%elem = b%m22%elem + a%m22%elem
         end if
      end if

      if (ierr /= 0 ) then
         call abort(" Error: Attempt to add bigblock matrices with different block structures")
      end if

   end subroutine

   !----------------------------------------------------------------------------
   ! Contract two bigblock matrices into a scalar
   ! s = a11.b11 + a12.b12 + a21.b21 + a22.b22
   !----------------------------------------------------------------------------
   subroutine contract_bbm(a,b,s)
      implicit none
      type(bigblockmatrix), intent(in) :: a,b
      real(dp), intent(out) :: s

      logical :: a11,a12,a21,a22
      logical :: b11,b12,b21,b22
      integer :: im

      ! Handle imaginary
      im = 0
      if (a%imag) im = im + 1
      if (b%imag) im = im + 1

      ! Check which quadrants are non-zero
      a11 = allocated(a%m11%elem); a12 = allocated(a%m12%elem)
      a21 = allocated(a%m21%elem); a22 = allocated(a%m22%elem)
      b11 = allocated(b%m11%elem); b12 = allocated(b%m12%elem)
      b21 = allocated(b%m21%elem); b22 = allocated(b%m22%elem)

      ! USE BLAS?
      s = 0
      if (a11.and.b11) s = s + dot_product(a%m11%elem, b%m11%elem)
      if (a12.and.b12) s = s + dot_product(a%m12%elem, b%m12%elem)
      if (a21.and.b21) s = s + dot_product(a%m21%elem, b%m21%elem)
      if (a22.and.b22) s = s + dot_product(a%m22%elem, b%m22%elem)
      if (im==2) s = -s

   end subroutine

   !----------------------------------------------------------------------------
   ! Multiply 3 bigblock matrices and return one quadrant
   ! - Optionally return the result transposed, conjugated, or negated
   ! - Assumes the resulting block matrix is already allocated
   !----------------------------------------------------------------------------
   subroutine triprod_bbm_quad(ta, ain, tb, bin, tc, cin, qabc, sabc, tabc, abcij)
      implicit none
      type(bigblockmatrix), intent(in)  :: ain, bin, cin
      integer, intent(in) :: qabc
      real(dp), intent(in) :: sabc
      character(1), intent(in) :: ta,tb,tc,tabc
      type(blockmatrix), intent(inout) :: abcij

      integer :: im, quadrant
      logical :: imag
      type(bigblockmatrix) :: a, b, c, aux
      logical :: a11,a12,a21,a22
      logical :: b11,b12,b21,b22
      logical :: c11,c12,c21,c22
      real(dp) :: s

      ! Setup variables based on inputs
      a=ain; b=bin; c=cin; quadrant=qabc; s=sabc

      ! Transpose inputs
      if (ta=='t') call transpose_bbm(a)
      if (tb=='t') call transpose_bbm(b)
      if (tc=='t') call transpose_bbm(c)

      ! Transpose output (A B C)^T = C^T B^T A^T
      ! (Note: by transposing the bigblock matrices, quadrant
      ! 12 moves to 21 and vice versa, we must follow it
      ! to return the requested result)
      if (tabc=='t') then
         aux=a; a=c; c=aux
         call transpose_bbm(a)
         call transpose_bbm(b)
         call transpose_bbm(c)
         if (quadrant==12) then
            quadrant=21
         else if (quadrant==21) then
            quadrant=12
         end if
      end if

      ! Handle imaginary unit (imag is not used...)
      im = 0
      if (a%imag) im=im+1
      if (b%imag) im=im+1
      if (c%imag) im=im+1
      if (mod(im,2)==0) then
         imag=.false.
      else
         imag=.true.
      end if
      if (im>=2) s = -1d0*s

      ! Check which quadrants are active
      a11 = allocated(a%m11%elem); a12 = allocated(a%m12%elem)
      a21 = allocated(a%m21%elem); a22 = allocated(a%m22%elem)
      b11 = allocated(b%m11%elem); b12 = allocated(b%m12%elem)
      b21 = allocated(b%m21%elem); b22 = allocated(b%m22%elem)
      c11 = allocated(c%m11%elem); c12 = allocated(c%m12%elem)
      c21 = allocated(c%m21%elem); c22 = allocated(c%m22%elem)

      abcij%elem=0

      ! DEBUG
      !write(*,'(i0)') q
      select case(quadrant)

         case(11)
            !! DEBUG
            !if(a11.and.b11.and.c11) print *, "a11.and.b11.and.c11 ----------- ", a%t11,b%t11,c%t11,s*a%s11*b%s11*c%s11
            !if(a12.and.b21.and.c11) print *, "a12.and.b21.and.c11 ----------- ", a%t12,b%t21,c%t11,s*a%s12*b%s21*c%s11
            !if(a11.and.b12.and.c21) print *, "a11.and.b12.and.c21 ----------- ", a%t11,b%t12,c%t21,s*a%s11*b%s12*c%s21
            !if(a12.and.b22.and.c21) print *, "a12.and.b22.and.c21 ----------- ", a%t12,b%t22,c%t21,s*a%s12*b%s22*c%s21
            !print *, 

            if(a11.and.b11.and.c11) call triprod(a%t11,a%m11,b%t11,b%m11,c%t11,c%m11,s*a%s11*b%s11*c%s11, 1d0, abcij)
            if(a12.and.b21.and.c11) call triprod(a%t12,a%m12,b%t21,b%m21,c%t11,c%m11,s*a%s12*b%s21*c%s11, 1d0, abcij)
            if(a11.and.b12.and.c21) call triprod(a%t11,a%m11,b%t12,b%m12,c%t21,c%m21,s*a%s11*b%s12*c%s21, 1d0, abcij)
            if(a12.and.b22.and.c21) call triprod(a%t12,a%m12,b%t22,b%m22,c%t21,c%m21,s*a%s12*b%s22*c%s21, 1d0, abcij)

         case(12)
            !! DEBUG
            !if(a11.and.b11.and.c12) print *, "a11.and.b11.and.c12 *********** ", a%t11,b%t11,c%t12,s*a%s11*b%s11*c%s12
            !if(a12.and.b21.and.c12) print *, "a12.and.b21.and.c12 *********** ", a%t12,b%t21,c%t12,s*a%s12*b%s21*c%s12
            !if(a11.and.b12.and.c22) print *, "a11.and.b12.and.c22 *********** ", a%t11,b%t12,c%t22,s*a%s11*b%s12*c%s22
            !if(a12.and.b22.and.c22) print *, "a12.and.b22.and.c22 *********** ", a%t12,b%t22,c%t22,s*a%s12*b%s22*c%s22
            !print *, 

            if(a11.and.b11.and.c12) call triprod(a%t11,a%m11,b%t11,b%m11,c%t12,c%m12,s*a%s11*b%s11*c%s12, 1d0, abcij)
            if(a12.and.b21.and.c12) call triprod(a%t12,a%m12,b%t21,b%m21,c%t12,c%m12,s*a%s12*b%s21*c%s12, 1d0, abcij)
            if(a11.and.b12.and.c22) call triprod(a%t11,a%m11,b%t12,b%m12,c%t22,c%m22,s*a%s11*b%s12*c%s22, 1d0, abcij)
            if(a12.and.b22.and.c22) call triprod(a%t12,a%m12,b%t22,b%m22,c%t22,c%m22,s*a%s12*b%s22*c%s22, 1d0, abcij)

         case(21)
            !! DEBUG
            !if(a21.and.b11.and.c11) print *, "a21.and.b11.and.c11 ----------- ", a%t21,b%t11,c%t11,s*a%s21*b%s11*c%s11
            !if(a22.and.b21.and.c11) print *, "a22.and.b21.and.c11 ----------- ", a%t22,b%t21,c%t11,s*a%s22*b%s21*c%s11
            !if(a21.and.b12.and.c21) print *, "a21.and.b12.and.c21 ----------- ", a%t21,b%t12,c%t21,s*a%s21*b%s12*c%s21
            !if(a22.and.b22.and.c21) print *, "a22.and.b22.and.c21 ----------- ", a%t22,b%t22,c%t21,s*a%s22*b%s22*c%s21
            !print *,

            if(a21.and.b11.and.c11) call triprod(a%t21,a%m21,b%t11,b%m11,c%t11,c%m11,s*a%s21*b%s11*c%s11, 1d0, abcij)
            if(a22.and.b21.and.c11) call triprod(a%t22,a%m22,b%t21,b%m21,c%t11,c%m11,s*a%s22*b%s21*c%s11, 1d0, abcij)
            if(a21.and.b12.and.c21) call triprod(a%t21,a%m21,b%t12,b%m12,c%t21,c%m21,s*a%s21*b%s12*c%s21, 1d0, abcij)
            if(a22.and.b22.and.c21) call triprod(a%t22,a%m22,b%t22,b%m22,c%t21,c%m21,s*a%s22*b%s22*c%s21, 1d0, abcij)

         case(22)
            !! DEBUG
            !if(a21.and.b11.and.c12) print *, "a21.and.b11.and.c12 *********** ", a%t21,b%t11,c%t12,s*a%s21*b%s11*c%s12
            !if(a22.and.b21.and.c12) print *, "a22.and.b21.and.c12 *********** ", a%t22,b%t21,c%t12,s*a%s22*b%s21*c%s12
            !if(a21.and.b12.and.c22) print *, "a21.and.b12.and.c22 *********** ", a%t21,b%t12,c%t22,s*a%s21*b%s12*c%s22
            !if(a22.and.b22.and.c22) print *, "a22.and.b22.and.c22 *********** ", a%t22,b%t22,c%t22,s*a%s22*b%s22*c%s22
            !print *, 

            if(a21.and.b11.and.c12) call triprod(a%t21,a%m21,b%t11,b%m11,c%t12,c%m12,s*a%s21*b%s11*c%s12, 1d0, abcij)
            if(a22.and.b21.and.c12) call triprod(a%t22,a%m22,b%t21,b%m21,c%t12,c%m12,s*a%s22*b%s21*c%s12, 1d0, abcij)
            if(a21.and.b12.and.c22) call triprod(a%t21,a%m21,b%t12,b%m12,c%t22,c%m22,s*a%s21*b%s12*c%s22, 1d0, abcij)
            if(a22.and.b22.and.c22) call triprod(a%t22,a%m22,b%t22,b%m22,c%t22,c%m22,s*a%s22*b%s22*c%s22, 1d0, abcij)

         case default
            call abort(" Value error: quadrant takes integer values 11,12,21,22")

      end select


   end subroutine

   !----------------------------------------------------------------------------
   ! Multiply 3 bigblock matrices and populate all non-zero quadrants
   !----------------------------------------------------------------------------
   subroutine triprod_bbm(ta,a, tb,b, tc,c, abc)
      implicit none
      type(bigblockmatrix), intent(inout) :: abc
      type(bigblockmatrix), intent(in)  :: a, b, c
      character(1), intent(in) :: ta,tb,tc

      logical :: abc11,abc12,abc21,abc22
      logical :: b11,b12,b21,b22
      logical :: c11,c12,c21,c22

      ! Check which quadrants of the result are active
      abc11 = allocated(abc%m11%elem); abc12 = allocated(abc%m12%elem)
      abc21 = allocated(abc%m21%elem); abc22 = allocated(abc%m22%elem)

      if(abc11) call triprod_bbm_quad(ta,a,tb,b,tc,c,11,abc%s11,abc%t11,abc%m11)
      if(abc12) call triprod_bbm_quad(ta,a,tb,b,tc,c,12,abc%s12,abc%t12,abc%m12)
      if(abc21) call triprod_bbm_quad(ta,a,tb,b,tc,c,21,abc%s21,abc%t21,abc%m21)
      if(abc22) call triprod_bbm_quad(ta,a,tb,b,tc,c,22,abc%s22,abc%t22,abc%m22)

   end subroutine

end module type_bigblockmatrix
