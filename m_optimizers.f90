module m_optimizers
    use m_func
    implicit none

    type nedler_mead  
      private
      ! standard parameters of method 
      real(kind=8) :: alpha, beta, gama, delta
      ! objective function to be optimized
      class(myfunc), pointer :: ufunc 
      ! number of degrees of freedom
      integer :: nDOF 
      ! intial guess for construction of simplex
      real(kind=8), dimension(:), pointer :: x0
      ! size of perturbation of x0
      real(kind=8) :: h 
      ! maimum number of iteration
      integer :: maxiter
      ! convergence criteria:
      !   relative change in simplex
      real(kind=8) :: xtol
      !   relative change in function  
      real(kind=8) :: ftol
    contains
      ! init method
      procedure, public  :: nm_init             
      ! cleanup method
      procedure, public  :: nm_cleanup             
      ! driver 
      procedure, public  :: driveOptimisation   
      ! method to calculate function value for each vertex  
      procedure, private :: simplexFunc         
      ! method to initialise the simplex  
      procedure, private :: initSimplex
      ! method to sort the simplex vertices by function value 
      procedure, private :: orderSimplex
      ! method to perform one iteration of the nm optimisation  
      procedure, private :: propagateSimplex
    end type nedler_mead

    contains
    
    subroutine nm_init(self, nDOF, x0, ufunc, maxiter, h, xtol, ftol, &
                       alpha, beta, gama, delta)
    ! 
    ! Initialisation of the optimisation method 
    !
      use m_func
      implicit none
      integer, intent(in)                :: nDOF
      class(nedler_mead)                 :: self
      class(myfunc), pointer, intent(in) :: ufunc
      real(kind=8), dimension(:), intent(in) :: x0
      real(kind=8), intent(in), optional :: h, xtol, ftol 
      real(kind=8), intent(in), optional :: alpha, beta, &
                                            gama, delta 
      integer, intent(in), optional      :: maxiter 
      
      ! set parameters for simplex propagation steps:
      !  reflection - alpha 
      if (present(alpha)) then 
        self%alpha = alpha 
      else
        self%alpha = 1.d0 
      endif
      !  shrinkage a) - beta
      if (present(beta)) then 
        self%beta = beta 
      else
        self%beta = 5.d-1 
      endif
      !   expansion - gamma
      if (present(gama)) then 
        self%gama = gama 
      else
        self%gama = 2.d0 
      endif
      !   shrinkage b) - delta
      if (present(delta)) then 
        self%delta = delta 
      else
        self%delta = 5.d-1 
      endif
      ! end set parameters

      ! set initial guess
      allocate(self%x0(nDOF))
      self%x0 = x0

      ! set number of degrees of freedom
      self%nDOF = nDOF

      ! set perturbation of initial vertex
      ! for constructing simplex
      if (present(h)) then
        self%h = h
      else
        self%h = 5.d0
      endif

      ! set maximum number of iterations
      if (present(maxiter)) then
        self%maxiter = maxiter
      else
        self%maxiter = 100
      endif

      ! set convergence tolerances:
      !   in simplex relative change
      if (present(xtol)) then
        self%xtol = xtol
      else
        self%xtol = 1.d-7
      endif
      !   in function relative change
      if (present(ftol)) then
        self%ftol = ftol
      else
        self%ftol = 1.d-7
      endif
      ! end set convergence tolerances

      self%ufunc => ufunc 
    end subroutine nm_init

    subroutine nm_cleanup(self)
      implicit none
      class(nedler_mead) :: self
      deallocate(self%x0)
      nullify(self%ufunc)
    end subroutine nm_cleanup

    subroutine simplexFunc(self, simplex, fSimplex)
    ! 
    !  Calculates function values to all vertices; 
    !  callable only from driver
    ! 
      implicit none
      class(nedler_mead) :: self
      !class(myfunc), pointer :: mfunc 
      real(kind=8), dimension(:,:), intent(in)   :: simplex
      real(kind=8), dimension(:), intent(inout)  :: fSimplex
      integer :: i
      
      do i=1,(self%nDOF+1)
        fSimplex(i) = self%ufunc%calc(self%nDOF,simplex(:,i)) 
      enddo

    end subroutine simplexFunc

    subroutine initSimplex(self, simplex, fSimplex)
    !
    !  Initialise the simplex; callable only from driver
    !
      implicit none
      class(nedler_mead) :: self 
      real(kind=8), dimension(:,:), intent(inout) :: simplex
      real(kind=8), dimension(:), intent(inout) :: fSimplex
      real(kind=8), allocatable :: perturb(:,:), cpSimplex(:,:), &
                                   cpFSimplex(:)
      integer :: i

      allocate(perturb(self%nDOF,self%nDOF),     &
               cpSimplex(self%nDOF,self%nDOF+1), &
               cpFSimplex(self%nDOF+1))

      perturb = 0.d0
      do i=1,self%nDOF
        perturb(i,i) = 1.d0
      enddo
      
      simplex(:,1) = self%x0
      do i=2,(self%nDOF+1) 
        simplex(:,i) = self%x0 + self%h * perturb(:,i-1)  
      enddo

      cpSimplex = simplex 
      call self%simplexFunc(cpSimplex, fSimplex)
          
      deallocate(perturb, cpSimplex, cpFSimplex)
    end subroutine initSimplex

    subroutine orderSimplex(self, simplex, fSimplex) 
      ! 
      ! Uses a linked list for insertion sort of vertices
      ! accoring to function value from lowest to highest 
      ! 
      implicit none

      class(nedler_mead) :: self
      real(kind=8), dimension(:,:), intent(inout) :: simplex
      real(kind=8), dimension(:), intent(inout) :: fSimplex
      type :: llist 
        real(kind=8) :: funcValue
        real(kind=8), allocatable :: vertex(:)
        type(llist), pointer :: next => NULL() 
      end type llist
      type(llist), pointer :: tail => NULL(), head => NULL(), &
                              ptr => NULL(), ptr1 => NULL(), &
                              ptr2 => NULL(), nptr => NULL(), &
                              nptr2 => NULL(), current => NULL(), &
                              next => NULL()
      real(kind=8) :: tmpFuncValue
      integer :: i


      allocate(ptr)
      allocate(ptr%vertex(self%nDOF))
      do i=1,(self%nDOF+1)
        ptr%funcValue = fSimplex(i) 
        ptr%vertex = simplex(:,i)
        if (.not. associated(head)) then
          allocate(head)
          allocate(head%vertex(self%nDOF))
          tail => head
          nullify(ptr%next)
          tail%funcValue = ptr%funcValue
          tail%vertex = ptr%vertex
        else
          if (ptr%funcValue < head%funcValue) then  
            ! add to front
            allocate(nptr)
            allocate(nptr%vertex(self%nDOF))
            nptr%funcValue = ptr%funcValue
            nptr%vertex = ptr%vertex
            nptr%next => head
            head => nptr
          elseif (ptr%funcValue >= tail%funcValue) then
            ! add to back
            allocate(tail%next)
            allocate(tail%next%vertex(self%nDOF))
            tail => tail%next
            tail%funcValue = ptr%funcValue
            tail%vertex = ptr%vertex
            nullify(tail%next)
          else
            ptr1 => head
            ptr2 => ptr1%next
            do
              if ((ptr%funcValue >= ptr1%funcValue) .and. &
                  (ptr%funcValue < ptr2%funcValue)) then
                allocate(nptr2)
                allocate(nptr2%vertex(self%nDOF))
                nptr2%funcValue = ptr%funcValue
                nptr2%vertex = ptr%vertex
                nptr2%next => ptr2
                ptr1%next => nptr2
                exit
              endif
              ptr1 => ptr2
              ptr2 => ptr1%next
            enddo
          endif
        endif
      enddo

      deallocate(ptr%vertex)
      deallocate(ptr)

      ptr => head
      i = 1
      do
        fSimplex(i)  = ptr%funcValue
        simplex(:,i) = ptr%vertex
        deallocate(ptr%vertex)
        ptr => ptr%next
        if (.not. associated(ptr)) exit
        i = i + 1
      enddo

      current => head
      next => current%next 
      i = 1
      do 
        deallocate(current)
        nullify(current)
        if(.not. associated(next)) exit
        current => next
        next => current%next
        i = i + 1
      enddo

      nullify(head)
      nullify(tail)
       
    end subroutine

    subroutine propagateSimplex(self, stepType, simplex, fSimplex)
      implicit none
      class(nedler_mead) :: self
      character(len=128), intent(inout) :: stepType 
      real(kind=8), dimension(:,:), intent(inout) :: simplex
      real(kind=8), dimension(:), intent(inout) :: fSimplex
      real(kind=8), allocatable :: cpSimplex(:,:), cpFSimplex(:)
      real(kind=8), allocatable :: c(:), xr(:), xe(:), xc(:)
      real(kind=8) :: fr, fe, fc
      integer :: i
      logical :: shrink

      allocate(cpSimplex(self%nDOF,self%nDOF+1), &
               cpFSimplex(self%nDOF+1), c(self%nDOF), &
               xr(self%nDOF))

      cpSimplex = simplex
      cpFSimplex = fSimplex

      call self%orderSimplex(cpSimplex, cpFSimplex)

      simplex = cpSimplex
      fSimplex = cpFSimplex

      shrink = .false.
      ! calculate simplex centroid
      c = 0.d0
      do i=1,self%nDOF
        c = c + simplex(:,i)
      enddo
      c = c / self%nDOF

      ! reflect worst vertex about centroid  
      xr = c + self%alpha * (c - simplex(:,self%nDOF+1)) 
      fr = self%ufunc%calc(self%nDOF,xr) 
      ! check if reflected point xr is in betweeen second
      ! worst and best vertices in simplex 
      if ((fSimplex(1) <= fr) .and.  &
          (fr < fSimplex(self%nDOF))) then
        simplex(:,self%nDOF+1) = xr
        fSimplex(self%nDOF+1) = fr
        write(steptype,*) 'Simple reflection'
      ! check if xr is better than best vertex in simplex
      elseif (fr < fSimplex(1)) then
        ! expand simplex
        allocate(xe(self%nDOF))
        xe = c + self%gama * (xr - c)
        fe = self%ufunc%calc(self%nDOF,xe)
        ! if expanded point xe better than reflected one
        if (fe < fr) then
          simplex(:,self%nDOF+1) = xe
          fSimplex(self%nDOF+1) = fe
          write(steptype,*) 'Expanded reflection'
        else
          simplex(:,self%nDOF+1) = xr
          fSimplex(self%nDOF+1) = fr
          write(steptype,*) 'Simple reflection'
        endif
        deallocate(xe)
      ! check if it is worse than the second worst vertex 
      elseif (fr >= fSimplex(self%nDOF)) then
        ! check if it is better than the worst
        if (fr < fSimplex(self%nDOF+1)) then
          ! contract simplex
          allocate(xc(self%nDOF))
          xc = c + self%beta * (xr - c)  
          fc = self%ufunc%calc(self%nDOF,xc)
          if (fc <= fr) then
            simplex(:,self%nDOF+1) = xc
            fSimplex(self%nDOF+1) = fc
            write(steptype,*) 'Outside contraction'
          else
            shrink = .true.
          endif
          deallocate(xc)
        ! check if it is worse than the worst
        elseif (fr >= fSimplex(self%nDOF+1)) then 
          xc = c + self%beta * (xr - c) 
          fc = self%ufunc%calc(self%nDOF,xc)
          if (fc < fr) then
            simplex(:,self%nDOF+1) = xc
            fSimplex(self%nDOF+1) = fc
            write(steptype,*) 'Inside contraction'
          else
            shrink = .true.
          endif
        endif
      endif

      if (shrink) then
        do i=2,(self%nDOF+1)
          simplex(:,i) = (simplex(:,1) + self%delta &
                          * (simplex(:,i) - simplex(:,1)))
          fSimplex(i) = self%ufunc%calc(self%nDOF,simplex(:,i))
        enddo
        write(steptype,*) 'Shrink'
      endif

      deallocate(cpSimplex, cpFSimplex, c, xr)
      
    end subroutine propagateSimplex

    subroutine driveOptimisation(self)
      implicit none
      class(nedler_mead) :: self
      character(len=128) :: stepType 
      real(kind=8), allocatable :: simplex(:,:), fSimplex(:)
      real(kind=8) :: xdev, fdev, tmpXdev, tmpFdev, &
                      maxXdev, maxFdev!, xabs
      integer :: i, j, k, ioerr
      
      allocate(simplex(self%nDOF,self%nDOF+1), &
               fSimplex(self%nDOF+1))
      i = 1
      call self%initSimplex(simplex, fSimplex)
      call self%propagateSimplex(stepType, simplex, fSimplex)
      open(11, file='optim.out', action='write', status='REPLACE',&
&          iostat=ioerr)

100   format(A1, A10, 3A30, 5X, A20)
200   format(' ', i10, 3ES30.18, 5X, A20)

      write(11,100) "#", "iteration", "min(f)", "max(dF)",  &
&                   "max(dX)", "step type"
      do 
        i = i + 1
        call self%propagateSimplex(stepType, simplex, fSimplex)
        maxXdev = 0.d0
        maxFdev = 0.d0
        do j=2,(self%nDOF+1)
          tmpXdev = 0.d0  
          do k=1,self%nDOF
            tmpXdev = tmpXdev + (simplex(k,1) - simplex(k,j))**2
          enddo
          tmpFdev = (fSimplex(1) - fSimplex(j))**2
          tmpFdev = dsqrt(tmpFdev)
          tmpXdev = dsqrt(tmpXdev) 
          if (tmpFdev > maxFdev) maxFdev = tmpFdev
          if (tmpXdev > maxXdev) maxXdev = tmpXdev
        enddo 
        write(11,200) i, fSimplex(1), maxFdev, maxXdev, stepType
        if ((maxFdev < self%ftol) .and. &
            (maxXdev < self%xtol))  then
          open(12, file='final.out', action='write', &
               status='REPLACE', iostat=ioerr)
          write(12,*) "Nedler-Mead converged after ", &
                     i, " iterations."
          write(12,*) "The solution is:" 
          write(12,*) ""
          write(12,*) "x = "
          do j=1,self%nDOF
            write(12,'(f20.13)') simplex(j,1)
          enddo
          write(12,*) ""
          write(12,'(A,f20.13)') "f = ", fSimplex(1)
          close(12)
          exit
        elseif (i > self%maxiter) then
          write(12,*) "Nedler-Mead stopped after ", &
                     i, " iterations, without converging."
          write(12,*) "Largest function deviation is"
          write(12,'(A,f20.13)') "df = ",  maxFdev
          close(12)
          exit
        endif
      enddo
      close(11)
    end subroutine driveOptimisation

end module m_optimizers
