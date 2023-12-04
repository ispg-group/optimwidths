program optimwidths
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!                                                             !!
    !! This program calculates the optimal width of the primitive  !!
    !! Gaussian basis functions employed in AIMS. It is an         !! 
    !! implementation of the method described in the SI of         !!
    !! M. P. Esch et al., J. Phys. Chem. A 2019, 123, 2661−2673,   !!
    !! DOI: 10.1021/acs.jpca.9b00952.                              !!
    !!                                                             !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! load modules
    use m_read
    use m_misc
    use m_func
    use m_optimizers
    ! end load modules
    implicit none

    ! input variables 
    integer :: nrAtoms
    integer :: maxiter
    real(kind=8) :: xtol
    real(kind=8) :: ftol
    character(len=128) :: prog
    character(len=128) :: fileName
    ! end input variables

    ! working variables 
    integer :: cartDOF  ! number of Cartesian degrees of freedom 
    integer :: normDOF  ! number of normal modes 
    integer :: nrType   ! number of distinct elements in molecule 
    integer :: i        ! dummy index
    ! normDOF x normDOF diagonal matrix of normal mode frequencies 
    real(kind=8), allocatable :: W(:,:)  
    ! normDOF x normDOF diagonal matrix of normal mode reduced masses 
    real(kind=8), allocatable :: M(:,:)
    ! cartDOF x cartDOF diagonal matrix of TBF widths 
    real(kind=8), allocatable :: A(:,:)
    ! cartDOF x normDOF matrix of normal modes 
    real(kind=8), allocatable :: U(:,:)
    ! nrAtoms dimensional vector of atomic masses
    real(kind=8), allocatable :: masses(:)
    ! nrType  or nrAtoms dimensional TBF width vector
    ! that is mapped onto the diagonal of A
    real(kind=8), allocatable :: x0(:)
    ! nrAtoms dimensional vector of element types
    integer, allocatable :: atomTypes_a(:)
    !real(kind=8) :: d
    logical :: constrained
    ! end working variables
    
    ! classes 
    ! Pointer to objective function class to be passed to nm optimizer class
    class(myfunc), pointer :: pfunc     
    ! Nedler Mead optimizer class
    class(nedler_mead), allocatable :: nm
    ! end classes

    call readInput(nrAtoms, maxiter, xtol, ftol, constrained, &
                   prog, fileName)

    cartDOF = 3 * nrAtoms
    normDOF = cartDOF - 6

    ! allocate working arrays
    allocate(W(normDOF,normDOF), M(normDOF,normDOF), &
             A(cartDOF,cartDOF), U(cartDOF,normDOF), &
             masses(nrAtoms), atomTypes_a(nrAtoms))
    ! allocate memory for pointer to objective function class
    allocate(pfunc)
    allocate(nm)

    nrType=0
    atomTypes_a = 0

    call readFreqOut(nrAtoms, cartDOF, normDOF, W, M, U, &
                     masses, prog, fileName, atomTypes_a,&
                     nrType)
    call convert2nm(nrAtoms, cartDOF, normDOF, U, M, masses)
    call pfunc%mf_init(nrAtoms, cartDOF, normDOF, W, M,  &
                       U, atomTypes_a, constrained)

    if (constrained) then 
      allocate(x0(nrType))
      do i=1,nrType
        write(*, *) "Input initial width of atom ", i
        read(*, *) x0(i)
      enddo
    else
      allocate(x0(nrAtoms))
      do i=1,nrAtoms
        write(*, *) "Input initial width of atom ", atomTypes_a(i)
        read(*, *) x0(i) 
      enddo
    endif

    write(*,*) "The initial guess is:"
    write(*,'(3f8.3)') x0

    if (constrained) then 
      call nm%nm_init(nrType,x0,pfunc,maxiter=maxiter)
    else
      call nm%nm_init(nrAtoms,x0,pfunc,maxiter=maxiter)
    endif
        
    call nm%driveOptimisation()

    ! clean up
    call nm%nm_cleanup()
    deallocate(W, M, A, U, masses, atomTypes_a, &
&              pfunc, x0, nm)
end program optimwidths
