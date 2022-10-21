program optimwidths
    use m_read
    use m_misc
    use m_func
    use m_optimizers
    implicit none

    ! input variables 
    integer            :: nrAtoms, maxiter
    real(kind=8)       :: xtol, ftol
    character(len=128) :: prog, fileName
    ! end input variables

    integer :: cartDOF, normDOF, info, ntype, i
    !real(kind=8) :: S
    real(kind=8), allocatable :: W(:,:), M(:,:), A(:,:), &
                                 U(:,:), widths(:), masses(:),& 
                                 x0(:) 
    integer, allocatable :: pivot(:), atomTypes_a(:)
    real(kind=8) :: d
    logical :: constrained
    class(myfunc), pointer :: pfunc
    type(nedler_mead) :: nm

    call readInput(nrAtoms, maxiter, xtol, ftol, constrained, &
                   prog, fileName)
    cartDOF = 3 * nrAtoms
    normDOF = cartDOF - 6
    allocate(W(normDOF,normDOF), M(normDOF,normDOF), &
             A(cartDOF,cartDOF), U(cartDOF,normDOF), &
             widths(nrAtoms),masses(nrAtoms),        &
             atomTypes_a(nrAtoms))

    allocate(pfunc)

    ntype=0
    atomTypes_a = 0
    call readFreqOut(nrAtoms, cartDOF, normDOF, W, M, U, &
                     masses, prog, fileName, atomTypes_a,&
                     ntype)
    call convert2nm(nrAtoms, cartDOF, normDOF, U, M, masses)
    call pfunc%mf_init(nrAtoms, cartDOF, normDOF, W, M, U, &
                       atomTypes_a, constrained)
    if (constrained) then 
      allocate(x0(ntype))
      do i=1,ntype
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
      call nm%nm_init(ntype,x0,pfunc,maxiter=maxiter)
    else
      call nm%nm_init(nrAtoms,x0,pfunc,maxiter=maxiter)
    endif
        
    call nm%driveOptimisation()

    deallocate(W, M, A, U, widths, masses, atomTypes_a, &
&              pfunc, x0)
end program optimwidths
