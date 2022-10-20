program optimwidths
    use m_func
    use m_optimizers
    implicit none

    ! input variables 
    integer            :: nrAtoms, maxiter
    real(kind=8)       :: xtol, ftol
    character(len=128) :: prog, fileName, detType
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
                   prog, fileName, detType)
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
                     ntype, constrained)
    call convert2nm(nrAtoms, cartDOF, normDOF, U, M, masses)
    call pfunc%mf_init(nrAtoms, cartDOF, normDOF, W, M, U, &
                       atomTypes_a, constrained, detType)
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

subroutine readInput(nrAtoms, maxiter, xtol, ftol, &
                     constrained, prog, fName, detType)
    integer, intent(inout)            :: maxiter, nrAtoms
    real(kind=8), intent(inout)       :: xtol, ftol
    logical, intent(inout)            :: constrained
    character(len=128), intent(inout) :: prog, fName,&
                                         detType
    namelist / optimparams / nrAtoms, maxiter, xtol, &
                             ftol, constrained, prog,&
                             fName, detType
    
    open(10, file='input.in')
    read(unit=10, nml=optimparams)
    write(*, '(A)') 'Read user input:'
    write(unit=*, nml=optimparams)
    close(10)
end subroutine readInput

subroutine readFreqOut(nrAtoms, cDOF, nDOF, W, M, U, masses, &
                       prog, fName, atomTypes_a, ntype)
    use m_func
    implicit none
    integer, intent(in)             :: nrAtoms, cDOF, nDOF
    integer, intent(inout)          :: ntype 
    real(kind=8), intent(inout)     :: W(nDOF,nDOF), M(nDOF,nDOF), &
                                       U(cDOF,nDOF), masses(nrAtoms), &
                                       atomTypes_a(nrAtoms)
    character(len=128), intent(in)  :: prog, fName
    character(len=256)              :: line
    real(kind=8), allocatable       :: atomCoords(:,:)
    integer, allocatable            :: nrMasses(:) 
    character(len=128), allocatable :: atomNames(:)
    character(len=20), allocatable  ::  tempCLine(:)
    real(kind=8), allocatable :: tempRLine(:) 
    real(kind=8), parameter :: freqConv = 4.5563352529119736d-6  
    real(kind=8), parameter :: proton   = 1822.888515d0
    logical :: rMasses, rAtomNames, rAtomCoords, rNormalModes 
    integer :: ioerr, ind, lineNr, i, j, ismass, ifreq, irmass, &
               inmode, jnmode, itmp, jtmp, tmpNrAtoms, nMass,   &
               imass
               

    allocate(atomCoords(nrAtoms,3), atomNames(nrAtoms))
    nMass = 1
    tmpNrAtoms = nrAtoms
    if (tmpNrAtoms > 10) then
      do 
        tmpNrAtoms = tmpNrAtoms - 10 
        nMass = nMass + 1
        if (tmpNrAtoms < 10) exit
      enddo
    endif
    allocate(nrMasses(nMass))
    do i=1,nMass
      if (i < nMass) nrMasses(i) = 10
      if (i == nMass) nrMasses(i) = tmpNrAtoms
    enddo
    rAtomNames = .false.
    rAtomCoords = .false.
    rNormalModes = .false.
    lineNr     = 0
    ismass      = 1
    imass      = 1
    ifreq      = 1
    irmass     = 1
    jnmode     = 0
    inmode     = 1
    W = 0.d0
    M = 0.d0
    U = 0.d0
    select case(prog)
    case('gaussian')
      open(unit=11, file=fName, status='old', action='read', &
           iostat=ioerr) 
      if (ioerr == 0) then 
        do 
          read(11, '(A)', iostat=ioerr) line
          if (ioerr /= 0) exit
          ind = index(trim(line), 'Symbolic Z-matrix')  
          if (ind /= 0) then 
            rAtomNames = .true. 
            lineNr     = lineNr + 1
            cycle
          endif
          
          if (rAtomNames .and. (lineNr > 1)) then
            read(line, '(x,a2,19x,3f9.5)') atomNames(lineNr-1), &
                                           atomCoords(lineNr-1,:)
            atomCoords = 0.d0
            lineNr     = lineNr + 1
            if (lineNr > (nrAtoms + 1)) then
                lineNr     = 0
                rAtomNames = .false. 
                call fillAtomTypes(nrAtoms,ntype,atomTypes_a,&
                                   atomNames)
            endif 
          elseif (rAtomNames .and. (lineNr == 1)) then
            lineNr = lineNr + 1
          endif

          ind = index(trim(line), 'AtmWgt=')
          if (ind /= 0) then
            if (ismass < nMass) allocate(tempCLine(1),tempRLine(10))
            if (ismass == nMass) allocate(tempCLine(1), &
                                         tempRLine(nrMasses(ismass)))
            if (ismass > nMass) cycle
            read(line,*) tempCLine(:), tempRLine(:)
            if (ismass < nMass) then
            do i=1,10
              masses(imass) = tempRLine(i)*proton
              !write(*,*) masses(imass)
              imass = imass + 1
            enddo
            elseif (ismass == nMass) then
            do i=1,nrMasses(ismass)
              masses(imass) = tempRLine(i)*proton
              !write(*,*) masses(imass)
              imass = imass + 1
            enddo
            endif
            ismass = ismass + 1
            deallocate(tempCLine, tempRLine)
          endif

          ind = index(trim(line), 'Standard orientation:')  
          if (ind /= 0) then 
            rAtomCoords = .true. 
            lineNr      = lineNr + 1
            cycle
          endif

          if (rAtomCoords .and. (lineNr > 4)) then
            allocate(tempRLine(6))
            read(line, *) tempRLine 
            atomCoords(lineNr-4,:) = tempRLine(4:)
            deallocate(tempRLine)
            lineNr = lineNr + 1
            if (lineNr > (nrAtoms + 4)) then
                lineNr      = 0
                rAtomCoords = .false. 
            endif 
          elseif (rAtomCoords .and. (lineNr <= 4)) then
            lineNr = lineNr + 1
          endif

          ind = index(trim(line), 'Frequencies')  
          if (ind /= 0) then
            rNormalModes = .true.
            lineNr = lineNr + 1 
            allocate(tempCLine(2),tempRLine(3))
            read(line,*) tempCLine(:), tempRLine(:) 
            do i=1,3
                W(ifreq,ifreq) = tempRLine(i)*freqConv
                ifreq = ifreq + 1
            enddo
            deallocate(tempCLine,tempRLine)
            cycle
          endif

          if (rNormalModes .and. (lineNr == 1)) then
            allocate(tempCLine(3),tempRLine(3))
            read(line,*) tempCLine(:), tempRLine(:) 
            do i=1,3
                M(irmass,irmass) = tempRLine(i)*proton
                irmass = irmass + 1
            enddo
            deallocate(tempCLine,tempRLine)
            lineNr = lineNr + 1
          elseif (rNormalModes .and. (lineNr <= 4)) then
            lineNr = lineNr + 1
          elseif (rNormalModes .and. (lineNr > 4)) then
            lineNr = lineNr + 1
            allocate(tempRLine(2))
            read(line,*) tempRLine(:), &
                         U(inmode:inmode+2,jnmode+1), &
                         U(inmode:inmode+2,jnmode+2), &
                         U(inmode:inmode+2,jnmode+3)
            inmode = inmode + 3
            deallocate(tempRLine)
            if (lineNr > (nrAtoms + 4)) then
                jnmode = jnmode + 3
                inmode = 1
                lineNr = 0
                rNormalModes = .false. 
            endif 
          endif
        enddo 
      endif 
      close(11)
    endselect
    deallocate(atomCoords, atomNames)
end subroutine readFreqOut

subroutine convert2nm(nrAtoms, cDOF, nDOF, U, M, masses)
    integer, intent(in)            :: nrAtoms, cDOF, nDOF
    real(kind=8), intent(in)       :: M(nDOF,nDOF), masses(nrAtoms)
    real(kind=8), intent(inout)    :: U(cDOF,nDOF)
    integer :: i, j, iAtom

    do i = 1,cDOF 
      iAtom = (i-1)/3 + 1
      do j = 1,nDOF
        U(i,j) = (U(i,j)*dsqrt(masses(iAtom))/dsqrt(M(j,j)))
      enddo
    enddo

end subroutine convert2nm

subroutine fillAtomTypes(nrAtoms, ind,  atomTypes_a, atomNames) 
    use m_func
    implicit none
    integer, intent(in) :: nrAtoms
    integer, intent(inout) :: ind
    character(len=128), intent(in) :: atomNames(nrAtoms)
    integer, intent(inout) :: atomTypes_a(nrAtoms)
    type(mf_llist), pointer :: head => NULL(), tail => NULL(), &
                               ptr => NULL(), ptr1 => NULL(),  &
                               next => NULL(), current => NULL() 
    integer :: i
    logical :: inllist

    allocate(ptr)
    ind = ind + 1
    do i=1,nrAtoms
      inllist = .false.
      ptr%atomName = atomNames(i)
      ptr%atomIndex = ind
      if (not(associated(head))) then
        allocate(head)
        tail => head
        nullify(ptr%mf_next)
        tail%atomName = ptr%atomName
        tail%atomIndex = ptr%atomIndex
        ind = ind + 1
      elseif (ptr%atomName /= tail%atomName) then
        ptr1 => head
        do 
          if (ptr%atomName == ptr1%atomName) then 
            inllist = .true.
            exit
          endif
          ptr1 => ptr1%mf_next
          if (not(associated(ptr1))) exit
        enddo
        if (not(inllist)) then
          allocate(tail%mf_next)
          tail => tail%mf_next
          nullify(tail%mf_next)
          tail%atomName = ptr%atomName
          tail%atomIndex = ptr%atomIndex
          ind = ind + 1
        endif
      endif
    enddo
    deallocate(ptr)
    
    ind = ind - 1
    write(*,*) "There are ", ind, " distinct atom types:"
    write(*,'(A,10x,A)') "Name", "Index"
    ptr => head
    do 
      write(*,*) ptr%atomName, ptr%atomIndex
      ptr => ptr%mf_next
      if(not(associated(ptr))) exit
    enddo

    do i=1,nrAtoms
      ptr => head
      do 
        if (ptr%atomName == atomNames(i)) then
          atomTypes_a(i) = ptr%atomIndex
        endif
        ptr => ptr%mf_next
        if(not(associated(ptr))) exit
      enddo
    enddo

    current => head
    next => current%mf_next
    do 
      deallocate(current)
      nullify(current)
      if(not(associated(next))) exit
      current => next
      next => current%mf_next
    enddo

    nullify(head)
    nullify(tail)

end subroutine fillAtomTypes
