module m_misc
    implicit none

    contains

    subroutine convert2nm(nrAtoms, cDOF, nDOF, U, M, masses, prog)
        integer, intent(in)            :: nrAtoms, cDOF, nDOF
        real(kind=8), intent(in)       :: M(nDOF,nDOF), masses(nrAtoms)
        real(kind=8), intent(inout)    :: U(cDOF,nDOF)
        character(len=128), intent(in) :: prog
        integer :: i, j, iAtom

        select case(prog)
            case('gaussian')
            do i = 1,cDOF 
              iAtom = (i-1)/3 + 1
              do j = 1,nDOF
                U(i,j) = (U(i,j)*dsqrt(masses(iAtom))/dsqrt(M(j,j)))
              enddo
            enddo
            case('terachem')
            return
        endselect


    end subroutine convert2nm

end module m_misc
