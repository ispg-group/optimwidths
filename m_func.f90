module m_func
    implicit none


    type myfunc 
      integer :: nrAtoms, cDOF, nDOF 
      logical :: constrained
      real(kind=8), allocatable :: W(:,:), M(:,:), U(:,:)
      integer, allocatable :: atomTypes(:)
    contains 
      procedure, public :: mf_init
      procedure, public :: calc
      procedure, private, nopass :: calc_det 
    end type myfunc

    contains

    subroutine mf_init(self, nrAtoms, cDOF, nDOF, W, M, U, &
                       atomTypes, constrained)
      class(myfunc) :: self 
      integer, intent(in) :: nrAtoms, cDOF, nDOF 
      integer, intent(in), optional :: atomTypes(nrAtoms)
      real(kind=8), intent(in) :: W(nDOF,nDOF), M(nDOF,nDOF), &
                                  U(cDOF,nDOF)
      logical, intent(in) :: constrained
      

      self%nrAtoms = nrAtoms
      self%nDOF = nDOF
      self%cDOF = cDOF
      allocate(self%W(nDOF,nDOF), self%M(nDOF,nDOF), &
               self%U(cDOF,nDOF),self%atomTypes(nrAtoms))
      self%W = W
      self%M = M
      self%U = U  
      self%atomTypes = atomTypes
      self%constrained = constrained

    end subroutine mf_init

    subroutine calc_det(ndim, B, det)
        integer, intent(in) :: ndim
        real(kind=8), intent(in) :: B(ndim, ndim)
        real(kind=8), intent(inout) :: det
        integer, allocatable :: pivot(:)
        integer :: info, i
        
        allocate(pivot(ndim))
        call dgetrf(ndim,ndim,B,ndim,pivot,info)
        if (info == 0) then
          det = 1.d0
          do i=1,ndim
            ! since this a PLU decomposition, we need
            ! take care if the row permutations via
            ! the pivot array
            if (pivot(i) == i) det = det * B(i,i)
            if (pivot(i) /= i) det = det * B(i,i) * (-1.d0)
          enddo
        else
          write(*,*) 'PLU failed'
        endif
    end subroutine calc_det

    real(kind=8) function calc(self, ndim, x)
        class(myfunc) :: self 
        integer, intent(in) :: ndim
        real(kind=8), intent(in) :: x(ndim)
        real(kind=8), allocatable :: mat1(:,:), mat2(:,:),    &
                                     mat3(:,:), A(:,:)
        real(kind=8) :: lnDet1, lnDet2, lnDet3, c1, c2, c3
        integer      :: i, j, k, tmpNdim
        
        tmpNdim = self%nDOF 

        allocate(mat1(tmpNdim,tmpNdim),mat2(tmpNdim,tmpNdim),&
                 mat3(tmpNdim,tmpNdim),A(self%cDOF,self%cDOF))
        
        A = 0.d0
        if (self%constrained) then
          do i=0,self%cDOF-1 
            A(i+1,i+1) = x(self%atomTypes(i/3 + 1)) 
          enddo
        else
          do i=0,self%nrAtoms-1 
            do j=1,3 
                k = 3 * i + j
                A(k,k) = x(i+1) 
            enddo
          enddo
        endif

        mat1 = matmul(transpose(self%U),matmul(A,self%U)) 
        mat2 = matmul(self%W,self%M)
        mat3 = mat2/2.d0 + mat1

        call self%calc_det(tmpNdim,2.d0*mat1,lnDet1)
        call self%calc_det(tmpNdim,mat2,lnDet2)
        call self%calc_det(tmpNdim,mat3,lnDet3)

        c1 = (-1.d0)/(4.d0)
        c2 = c1
        c3 = (1.d0)/(2.d0)
        
        lnDet1 = c1*dlog(lnDet1)
        lnDet2 = c2*dlog(lnDet2)
        lnDet3 = c3*dlog(lnDet3)

        calc = lnDet1 + lnDet2 + lnDet3
        
        deallocate(mat1,mat2,mat3,A)
    end function calc 

end module m_func
