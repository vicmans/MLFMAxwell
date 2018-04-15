module modIterativo
    use modGlobalParam
    use modGlobalMetod

    implicit none
    private
    public :: csrmsr, ilut, msrcsr, csrcoo, coocsr, ilu0, csrmsrC, msrcsrC, csrcooC, coocsrC, ilu0C, lusol0C, CGSold, ilutC, amuxC

contains



    subroutine CGSold(func_MatxVec, b, guess, errmax, n, x, mensaje, prec_alu, prec_jlu, prec_ju)
        !dummy
        complex (kind = dp), intent(in), dimension(n) :: guess, b
        complex (kind = dp), intent(out), dimension(n) :: x
        real (kind = dp),  intent(in) :: errmax
        integer (kind = il),  intent(in) :: n
        integer (kind = il), optional, intent(in) :: mensaje
        complex (kind = dp_prec), optional, intent(in), dimension(:) :: prec_alu
        integer (kind = il), optional, intent(in), dimension(:) :: prec_jlu
        integer (kind = ild), optional, intent(in), dimension(:) :: prec_ju

        !
        !local
        complex (kind = dp), allocatable, dimension(:) :: r, R_, u, p, P_, V_, q, U_, Q_ !_ = gorrito
        complex (kind = dp) :: beta, alfa
        complex (kind = dp), dimension(2) :: rho !rho(1) = i-1; rho(2) = i-2
        real (kind = dp) :: resid, residde
        integer (kind = il) :: i
        logical :: prec
        !

        interface
            function func_MatxVec(vecX) result(res)
                complex (kind = 8), intent(in), dimension(:) :: vecX
                complex (kind = 8), dimension(size(vecX,1)) :: res
            end function

        end interface


        allocate(r(n))
        allocate(R_(n))
        allocate(u(n))
        allocate(p(n))
        allocate(P_(n))
        allocate(V_(n))
        allocate(q(n))
        allocate(U_(n))
        allocate(Q_(n))


        prec = .false.

        if (present(prec_alu) .and. present(prec_jlu) .and. present(prec_ju)) then
          prec = .true.
        end if

        call printconsole(' ',  'b_residue/' // trim(realtostrE2(sqrt(real(dot_product(b,b))),5)) )
        
        !inicial
        r = b - func_MatxVec(guess)!r = b - matmul(Z,guess)
        R_ = r
        x = guess
        i = 0
        do
            i = i+1

            rho(1) = sum(r*R_)

            if (rho(1) == 0) then
                !falla
                exit
            end if
            if (i == 1) then
                u = r
                p = u
            else
                beta = rho(1)/rho(2)
                u = r + (beta*q)
                p = u + (beta*(q+(beta*p)))
            end if

            if (prec .eqv. .false.) then
              P_ = p !FALTA PRECONDICIONADOR (P_ = (M⁻1)*p)
            else
              call lusol0C(n, p, P_, prec_alu, prec_jlu, prec_ju)
            end if

            V_ = func_MatxVec(P_)!V_ = matmul(Z,P_)
            !alfa = rho(1)/Nproducto_escalarC(R_,V_,n)
            alfa = rho(1)/sum(R_*V_)
            q = u -(alfa*V_)


            if (prec .eqv. .false.) then
              U_ = (u+q) !FALTA PRECONDICIONADOR (U_ = (M⁻1)*(u+q))
            else
              call lusol0C(n, u+q, U_, prec_alu, prec_jlu, prec_ju)
            end if


            x = x + (alfa*U_)
            Q_ = func_MatxVec(U_)!Q_ = matmul(Z,U_)
            r = r -(alfa*Q_)

            if (present(mensaje)) then

                resid = sqrt(real(dot_product(r,r)))
                residde = sqrt(real(dot_product(b,b)))
                call printconsole( 'Iterac: ' // inttostr(i) // ' residuo ' //  trim(realtostrE2(resid,5)) // ' de ' &
                & // trim(realtostrE2(residde,5)) // ' - ' // trim(gettiempo()), 'iteration_residue/' &
                & // trim(realtostrE2(resid,5)) // '/' // inttostr(i) ) 
call savetolog('Iteracion: ' // inttostr(i) // ' residuo ' // realtostrE(resid) // ' de ' // realtostrE(residde) // '...')
                call savetolog(gettiempo())
            end if

            if (  sqrt(real(dot_product(r,r))) < errmax*sqrt(real(dot_product(b,b))) ) then
                if (present(mensaje)) then
                    call printconsole('El sistema ha convergido con error: ' // &
                      realtostrE(sqrt(real(dot_product(r,r)/dot_product(b,b)))), 'convergence_error/' // &
                      & realtostrE(sqrt(real(dot_product(r,r)/dot_product(b,b)))))
                end if
                exit !Convergio
            end if

            !actualizacion
            rho(2) = rho(1)

        end do


        deallocate(r)
        deallocate(R_)
        deallocate(u)
        deallocate(p)
        deallocate(P_)
        deallocate(V_)
        deallocate(q)
        deallocate(U_)
        deallocate(Q_)

    end subroutine










subroutine lusol0C ( n, y, x, alu, jlu, ju )

!*****************************************************************************80
!
!! LUSOL0 performs a forward followed by a backward solve
! for LU matrix as produced by  ILUT
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer ( kind = il ) N, the order of the matrix.
!
!    Input, real Y(N), the right hand side of the linear system.
!
!    Output, real X(N), the solution.
!
!    ALU, JLU, JU, ...
!
  implicit none

  integer ( kind = il ) n

  complex ( kind = dp_prec ) alu(*)
  integer ( kind = il ) i
  integer ( kind = il ) jlu(*)
  integer ( kind = ild ) ju(*)
  integer ( kind = ild ) k
  complex ( kind = dp ) x(n)
  complex ( kind = dp ) y(n)
!
!  Forward solve
!
  do i = 1, n
    x(i) = y(i)
    do k = jlu(i), ju(i)-1
      x(i) = x(i) - alu(k) * x(jlu(k))
    end do
  end do
!
!  Backward solve.
!
  do i = n, 1, -1
    do k = ju(i), jlu(i+1)-1
      x(i) = x(i) - alu(k) * x(jlu(k))
    end do
    x(i) = alu(i) * x(i)
  end do

  return
end subroutine




subroutine ilu0 ( n, a, ja, ia, alu, jlu, ju, iw, ierr )

!*****************************************************************************80
!
!! ILU0 is an ILU(0) preconditioner.
!
!  Discussion:
!
!    Note that this has been coded in such a way that it can be used
!    with PGMRES.  Normally, since the data structure of a, ja, ia is
!    the same as that of a, ja, ia, savings can be made. In fact with
!    some definitions (not correct for general sparse matrices) all we
!    need in addition to a, ja, ia is an additional diagonal.
!    Ilu0 is not recommended for serious problems. It is only provided
!    here for comparison purposes.
!
!    It is assumed that the the elements in the input matrix are stored
!    in such a way that in each row the lower part comes first and
!    then the upper part. To get the correct ILU factorization, it is
!    also necessary to have the elements of L sorted by increasing
!    column number. It may therefore be necessary to sort the
!    elements of a, ja, ia prior to calling ilu0. This can be
!    achieved by transposing the matrix twice using csrcsc.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer ( kind = il ) N, the order of the matrix.
!
!    Input, real A(*), integer ( kind = il ) JA(*), IA(N+1), the matrix in CSR
!    Compressed Sparse Row format.
!
! on return:
!
! alu,jlu = matrix stored in Modified Sparse Row (MSR) format containing
!           the L and U factors together. The diagonal (stored in
!           alu(1:n) ) is inverted. Each i-th row of the alu,jlu matrix
!           contains the i-th row of L (excluding the diagonal entry=1)
!           followed by the i-th row of U.
!
! ju        = pointer to the diagonal elements in alu, jlu.
!
! ierr        = integer ( kind = il ) indicating error code on return
!           ierr = 0 --> normal return
!           ierr = k --> code encountered a zero pivot at step k.
! work arrays:
!
! iw          = integer ( kind = il ) work array of length n.
!
!
  implicit none

  integer ( kind = il ) n

  real ( kind = dp ) a(*)
  real ( kind = dp ) alu(*)
  integer ( kind = il ) i
  integer ( kind = il ) ia(n+1)
  integer ( kind = il ) ierr
  integer ( kind = il ) ii
  integer ( kind = il ) iw(n)
  integer ( kind = il ) j
  integer ( kind = il ) ja(*)
  integer ( kind = il ) jcol
  integer ( kind = il ) jf
  integer ( kind = il ) jj
  integer ( kind = il ) jlu(*)
  integer ( kind = il ) jm
  integer ( kind = il ) jrow
  integer ( kind = il ) js
  integer ( kind = il ) ju(*)
  integer ( kind = il ) ju0
  integer ( kind = il ) jw
  real ( kind = dp ) tl

  ju0 = n + 2
  jlu(1) = ju0
!
!  Initialize the work vector.
!
  iw(1:n) = 0
!
!  The main loop.
!
  do ii = 1, n

    js = ju0
!
!  Generating row II of L and U.
!
    do j = ia(ii), ia(ii+1)-1
!
!  Copy row II of A, JA, IA into row II of ALU, JLU (L/U) matrix.
!
      jcol = ja(j)

      if ( jcol == ii ) then
        alu(ii) = a(j)
        iw(jcol) = ii
        ju(ii) = ju0
      else
        alu(ju0) = a(j)
        jlu(ju0) = ja(j)
        iw(jcol) = ju0
        ju0 = ju0 + 1
      end if

    end do

    jlu(ii+1) = ju0
    jf = ju0 - 1
    jm = ju(ii) - 1
!
!  Exit if the diagonal element is reached.
!
    do j = js, jm

      jrow = jlu(j)
      tl = alu(j) * alu(jrow)
      alu(j) = tl
!
!  Perform linear combination.
!
      do jj = ju(jrow), jlu(jrow+1)-1
        jw = iw(jlu(jj))
        if ( jw /= 0 ) then
          alu(jw) = alu(jw) - tl * alu(jj)
        end if
      end do

    end do
!
!  Invert and store the diagonal element.
!
    if ( alu(ii) == 0.0D+00 ) then
      ierr = ii
      return
    end if

    alu(ii) = 1.0D+00 / alu(ii)
!
!  Reset pointer IW to zero.
!
    iw(ii) = 0
    do i = js, jf
      iw(jlu(i)) = 0
    end do

  end do

  ierr = 0
  return
end subroutine



subroutine coocsr ( nrow, nnz, a, ir, jc, ao, jao, iao )

!*****************************************************************************80
!
!! COOCSR converts COO to CSR.
!
!  Discussion:
!
!    This routine converts a matrix that is stored in COO coordinate format
!    a, ir, jc into a CSR row general sparse ao, jao, iao format.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer ( kind = il ) NROW, the row dimension of the matrix.
!
!    Input, integer ( kind = il ) NNZ, the number of nonzero elements.
!
! a,
! ir,
! jc    = matrix in coordinate format. a(k), ir(k), jc(k) store the nnz
!         nonzero elements of the matrix with a(k) = actual real value of
!         the elements, ir(k) = its row number and jc(k) = its column
!        number. The order of the elements is arbitrary.
!
! on return:
!
! ir       is destroyed
!
!    Output, real AO(*), JAO(*), IAO(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
  implicit none

  integer ( kind = il ) nrow

  real ( kind = dp ) a(*)
  real ( kind = dp ) ao(*)
  integer ( kind = il ) i
  integer ( kind = il ) iad
  integer ( kind = il ) iao(nrow+1)
  integer ( kind = il ) ir(*)
  integer ( kind = il ) j
  integer ( kind = il ) jao(*)
  integer ( kind = il ) jc(*)
  integer ( kind = il ) k
  integer ( kind = il ) k0
  integer ( kind = il ) nnz
  real ( kind = dp ) x

  iao(1:nrow+1) = 0
!
!  Determine the row lengths.
!
  do k = 1, nnz
    iao(ir(k)) = iao(ir(k)) + 1
  end do
!
!  The starting position of each row.
!
  k = 1
  do j = 1, nrow+1
     k0 = iao(j)
     iao(j) = k
     k = k + k0
  end do
!
!  Go through the structure once more.  Fill in output matrix.
!
  do k = 1, nnz
     i = ir(k)
     j = jc(k)
     x = a(k)
     iad = iao(i)
     ao(iad) = x
     jao(iad) = j
     iao(i) = iad + 1
  end do
!
!  Shift back IAO.
!
  do j = nrow, 1, -1
    iao(j+1) = iao(j)
  end do
  iao(1) = 1

  return
end subroutine



subroutine csrcoo ( nrow, job, nzmax, a, ja, ia, nnz, ao, ir, jc, ierr )

!*****************************************************************************80
!
!! CSRCOO converts Compressed Sparse Row to Coordinate format.
!
!  Discussion:
!
!   This routine converts a matrix that is stored in row general sparse
!   A, JA, IA format into coordinate format AO, IR, JC.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer ( kind = il ) NROW, the row dimension of the matrix.
!
! job   = integer ( kind = il ) serving as a job indicator.
!         if job = 1 fill in only the array ir, ignore jc, and ao.
!         if job = 2 fill in ir, and jc but not ao
!         if job = 3 fill in everything.
!         The reason why these options are provided is that on return
!         ao and jc are the same as a, ja. So when job = 3, a and ja are
!         simply copied into ao, jc.  When job=2, only jc and ir are
!         returned. With job=1 only the array ir is returned. Moreover,
!         the algorithm is in place:
!           call csrcoo (nrow,1,nzmax,a,ja,ia,nnz,a,ia,ja,ierr)
!         will write the output matrix in coordinate format on a, ja,ia.
!         (Important: note the order in the output arrays a, ja, ia. )
!         i.e., ao can be the same as a, ir can be the same as ia
!         and jc can be the same as ja.
!
!    Input, real A(*), integer ( kind = il ) JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
! nzmax = length of space available in ao, ir, jc.
!         the code will stop immediatly if the number of
!         nonzero elements found in input matrix exceeds nzmax.
!
! on return:
!-
! ao, ir, jc = matrix in coordinate format.
!
! nnz        = number of nonzero elements in matrix.
!
! ierr       = integer ( kind = il ) error indicator.
!         ierr == 0 means normal retur
!         ierr == 1 means that the the code stopped
!         because there was no space in ao, ir, jc
!         (according to the value of  nzmax).
!
  implicit none

  integer ( kind = il ) nrow

  real ( kind = dp ) a(*)
  real ( kind = dp ) ao(*)
  integer ( kind = il ) i
  integer ( kind = il ) ia(nrow+1)
  integer ( kind = il ) ierr
  integer ( kind = il ) ir(*)
  integer ( kind = il ) ja(*)
  integer ( kind = il ) jc(*)
  integer ( kind = il ) job
  integer ( kind = il ) k
  integer ( kind = il ) k1
  integer ( kind = il ) k2
  integer ( kind = il ) nnz
  integer ( kind = il ) nzmax

  ierr = 0
  nnz = ia(nrow+1)-1

  if ( nzmax < nnz ) then
    ierr = 1
    return
  end if

  if ( 3 <= job ) then
    ao(1:nnz) = a(1:nnz)
  end if

  if ( 2 <= job ) then
    jc(1:nnz) = ja(1:nnz)
  end if
!
!  Copy backward.
!
  do i = nrow, 1, -1
    k1 = ia(i+1) - 1
    k2 = ia(i)
    do k = k1, k2, -1
      ir(k) = i
    end do
  end do

  return
end subroutine


subroutine msrcsr ( n, a, ja, ao, jao, iao, wk )

!*****************************************************************************80
!
!! MSRCSR converts Modified Sparse Row to Compressed Sparse Row.
!
!  Discussion:
!
!    This routine converts a compressed matrix using a separated diagonal
!    (modified sparse row format) in the Compressed Sparse Row format.
!
!    does not check for zero elements in the diagonal.
!
!    This is an "in place" algorithm (see a, ja, ia).
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer ( kind = il ) N, the row dimension of the matrix.
!
! ao, jao  = sparse matrix in msr sparse storage format
!           see routine csrmsr for details
!
! on return :
!
! a, ja, ia = matrix in csr format. note that the
!           algorithm is in place: ao, jao can be the same
!            as a, ja, in which case it will be overwritten on it
!            upon return.
!
!             here nnz = number of nonzero elements+1
!
!    Workspace, real WK(N).
!
  implicit none

  integer ( kind = il ) n

  real ( kind = dp_prec ) a(*)
  logical added
  real ( kind = dp ) ao(*)
  integer ( kind = il ) iao(n+1)
  integer ( kind = il ) idiag
  integer ( kind = il ) ii
  integer ( kind = il ) iptr
  integer ( kind = il ) j
  integer ( kind = il ) ja(*)
  integer ( kind = il ) jao(*)
  integer ( kind = il ) k
  real ( kind = dp ) wk(n)

  wk(1:n) = a(1:n)

  iao(1) = 1
  iptr = 1

  do ii = 1, n

    added = .false.
    idiag = iptr + ( ja(ii+1) - ja(ii) )

    do k = ja(ii), ja(ii+1)-1

      j = ja(k)

      if ( j < ii ) then
        ao(iptr) = a(k)
        jao(iptr) = j
        iptr = iptr + 1
      else if ( added ) then
        ao(iptr) = a(k)
        jao(iptr) = j
        iptr = iptr + 1
      else
!
!  Add diagonal element.  Only reserve a position for it.
!
        idiag = iptr
        iptr = iptr + 1
        added = .true.
!
!  Then other elements.
!
        ao(iptr) = a(k)
        jao(iptr) = j
        iptr = iptr + 1
      end if

    end do

    ao(idiag) = wk(ii)
    jao(idiag) = ii
    if ( .not. added ) then
      iptr = iptr + 1
    end if
    iao(ii+1) = iptr

  end do

  return
end subroutine



  subroutine csrmsr ( n, a, ja, ia, ao, jao, wk, iwk )

  !*****************************************************************************80
  !
  !! CSRMSR converts Compressed Sparse Row to Modified Sparse Row.
  !
  !  Discussion:
  !
  !    This routine converts a general sparse matrix a, ja, ia into
  !    a compressed matrix using a separated diagonal (referred to as
  !    the bell-labs format as it is used by bell labs semi conductor
  !    group. We refer to it here as the modified sparse row format.
  !
  !    This has been coded in such a way that one can overwrite
  !    the output matrix onto the input matrix if desired by a call of
  !    the form
  !
  !     call csrmsr (n, a, ja, ia, a, ja, wk,iwk)
  !
  !    In case ao, jao, are different from a, ja, then one can
  !    use ao, jao as the work arrays in the calling sequence:
  !
  !     call csrmsr (n, a, ja, ia, ao, jao, ao,jao)
  !
  !    Algorithm is in place.  i.e. both:
  !
  !          call csrmsr (n, a, ja, ia, ao, jao, ao,jao)
  !          (in which  ao, jao, are different from a, ja)
  !           and
  !          call csrmsr (n, a, ja, ia, a, ja, wk,iwk)
  !          (in which  wk, jwk, are different from a, ja)
  !        are OK.
  !
  !  Modified:
  !
  !    07 January 2004
  !
  !  Author:
  !
  !    Youcef Saad
  !
  !  Parameters:
  !
  !    Input, integer ( kind = il ) N, the order of the matrix.
  !
  !    Input, real A(*), integer ( kind = il ) JA(*), IA(N+1), the matrix in CSR
  !    Compressed Sparse Row format.
  !
  ! on return :
  !
  ! ao, jao  = sparse matrix in modified sparse row storage format:
  !         +  ao(1:n) contains the diagonal of the matrix.
  !         +  ao(n+2:nnz) contains the nondiagonal elements of the
  !             matrix, stored rowwise.
  !         +  jao(n+2:nnz) : their column indices
  !         +  jao(1:n+1) contains the pointer array for the nondiagonal
  !             elements in ao(n+1:nnz) and jao(n+2:nnz).
  !             i.e., for i <= n+1 jao(i) points to beginning of row i
  !            in arrays ao, jao.
  !             here nnz = number of nonzero elements+1
  !
  !    Work array, real WK(N).
  !
  !    Work array, integer ( kind = il ) IWK(N+1).
  !
      implicit none

      integer ( kind = il ) n

      real ( kind = dp ) a(*)
      real ( kind = dp ) ao(*)
      integer ( kind = il ) i
      integer ( kind = il ) ia(n+1)
      integer ( kind = il ) icount
      integer ( kind = il ) ii
      integer ( kind = il ) iptr
      integer ( kind = il ) iwk(n+1)
      integer ( kind = il ) j
      integer ( kind = il ) ja(*)
      integer ( kind = il ) jao(*)
      integer ( kind = il ) k
      real ( kind = dp ) wk(n)

      icount = 0
    !
    !  Store away diagonal elements and count nonzero diagonal elements.
    !
      do i = 1, n
        wk(i) = 0.0D+00
        iwk(i+1) = ia(i+1) - ia(i)
        do k = ia(i), ia(i+1)-1
          if ( ja(k) == i ) then
            wk(i) = a(k)
            icount = icount + 1
            iwk(i+1) = iwk(i+1) - 1
          end if
        end do
      end do
    !
    !  Compute total length.
    !
      iptr = n + ia(n+1) - icount
    !
    !  Copy backwards, to avoid collisions.
    !
      do ii = n, 1, -1
        do k = ia(ii+1)-1, ia(ii), -1
          j = ja(k)
          if ( j /= ii ) then
            ao(iptr) = a(k)
            jao(iptr) = j
            iptr = iptr - 1
          end if
        end do
      end do
    !
    !  Compute the pointer values and copy WK.
    !
      jao(1) = n + 2
      do i = 1, n
        ao(i) = wk(i)
        jao(i+1) = jao(i) + iwk(i+1)
      end do

      return
  end subroutine



  subroutine ilut ( n, a, ja, ia, lfil, tol, alu, jlu, ju, iwk, wu, wl, jr, &
    jwl, jwu, ierr )

  !*****************************************************************************80
  !
  !! ILUT is an ILUT preconditioner.
  !
  !  Discussion:
  !
  !    This routine carries ouot incomplete LU factorization with dual
  !    truncation mechanism.  Sorting is done for both L and U.
  !
  !    The dual drop-off strategy works as follows:
  !
  !    1) Theresholding in L and U as set by TOL.  Any element whose size
  !       is less than some tolerance (relative to the norm of current
  !       row in u) is dropped.
  !
  !    2) Keeping only the largest lenl0+lfil elements in L and the
  !       largest lenu0+lfil elements in U, where lenl0=initial number
  !       of nonzero elements in a given row of lower part of A
  !       and lenlu0 is similarly defined.
  !
  !    Flexibility: one can use tol=0 to get a strategy based on keeping the
  !    largest elements in each row of L and U. Taking tol /= 0 but lfil=n
  !    will give the usual threshold strategy (however, fill-in is then
  !    unpredictible).
  !
  !    A must have all nonzero diagonal elements.
  !
  !  Modified:
  !
  !    07 January 2004
  !
  !  Author:
  !
  !    Youcef Saad
  !
  !  Parameters:
  !
  !  Parameters:
  !
  !    Input, integer ( kind = il ) N, the order of the matrix.
  !
  !    Input, real A(*), integer ( kind = il ) JA(*), IA(N+1), the matrix in CSR
  !    Compressed Sparse Row format.
  !
  ! lfil    = integer ( kind = il ). The fill-in parameter. Each row of L and
  !           each row of U will have a maximum of lfil elements
  !           in addition to the original number of nonzero elements.
  !           Thus storage can be determined beforehand.
  !           lfil must be >= 0.
  !
  ! iwk     = integer ( kind = il ). The minimum length of arrays alu and jlu
  !
  ! On return:
  !
  ! alu,jlu = matrix stored in Modified Sparse Row (MSR) format containing
  !           the L and U factors together. The diagonal (stored in
  !           alu(1:n) ) is inverted. Each i-th row of the alu,jlu matrix
  !           contains the i-th row of L (excluding the diagonal entry=1)
  !           followed by the i-th row of U.
  !
  ! ju      = integer ( kind = il ) array of length n containing the pointers to
  !           the beginning of each row of U in the matrix alu,jlu.
  !
  ! ierr    = integer ( kind = il ). Error message with the following meaning.
  !           ierr  = 0    --> successful return.
  !           ierr > 0  --> zero pivot encountered at step number ierr.
  !           ierr  = -1   --> Error. input matrix may be wrong.
  !                            (The elimination process has generated a
  !                            row in L or U whose length is >  n.)
  !           ierr  = -2   --> The matrix L overflows the array al.
  !           ierr  = -3   --> The matrix U overflows the array alu.
  !           ierr  = -4   --> Illegal value for lfil.
  !           ierr  = -5   --> zero pivot encountered.
  !
  ! work arrays:
  !
  ! jr,jwu,jwl, integer ( kind = il ) work arrays of length n.
  ! wu, wl, real work arrays of length n+1, and n resp.
  !
      implicit none

      integer ( kind = il ) n

      real ( kind = dp ) a(*)
      real ( kind = dp ) alu(*)
      real ( kind = dp ) fact
      integer ( kind = il ) ia(n+1)
      integer ( kind = il ) idiag
      integer ( kind = il ) ierr
      integer ( kind = il ) ii
      integer ( kind = il ) iwk
      integer ( kind = il ) j
      integer ( kind = il ) j1
      integer ( kind = il ) j2
      integer ( kind = il ) ja(*)
      integer ( kind = il ) jj
      integer ( kind = il ) jlu(*)
      integer ( kind = il ) jpos
      integer ( kind = il ) jr(*)
      integer ( kind = il ) jrow
      integer ( kind = il ) ju(*)
      integer ( kind = il ) ju0
      integer ( kind = il ) jwl(n)
      integer ( kind = il ) jwu(n)
      integer ( kind = il ) k
      integer ( kind = il ) len
      integer ( kind = il ) lenl
      integer ( kind = il ) lenl0
      integer ( kind = il ) lenu
      integer ( kind = il ) lenu0
      integer ( kind = il ) lfil
      integer ( kind = il ) nl
      real ( kind = dp ) s
      real ( kind = dp ) t
      real ( kind = dp ) tnorm
      real ( kind = dp ) tol
      real ( kind = dp ) wl(n)
      real ( kind = dp ) wu(n)

      if ( lfil < 0 ) then
        ierr = -4
        return
      end if
    !
    !  Initialize JU0 (points to next element to be added to ALU, JLU)
    !  and pointer.
    !
      ju0 = n + 2 
      jlu(1) = ju0
    !
    !  integer ( kind = il ) double pointer array.
    !
      jr(1:n) = 0
    !
    !  The main loop.
    !
      do ii = 1, n

        j1 = ia(ii)
        j2 = ia(ii+1) - 1
        lenu = 0
        lenl = 0

        tnorm = 0.0D+00
        do k = j1, j2
          tnorm = tnorm + abs ( a(k) )
        end do
        tnorm = tnorm / real ( j2-j1+1, kind = dp )
    !
    !  Unpack L-part and U-part of row of A in arrays WL, WU.
    !
        do j = j1, j2

          k = ja(j)
          t = a(j)

          if ( tol * tnorm <= abs ( t ) ) then

            if ( k < ii ) then
              lenl = lenl + 1
              jwl(lenl) = k
              wl(lenl) = t
              jr(k) = lenl
            else
              lenu = lenu+1
              jwu(lenu) = k
              wu(lenu) = t
              jr(k) = lenu
            end if

          end if

        end do

        lenl0 = lenl
        lenu0 = lenu
        jj = 0
        nl = 0
    !
    !  Eliminate previous rows.
    !
    150 continue

        jj = jj + 1

        if ( lenl < jj ) then
          go to 160
        end if
    !
    !  In order to do the elimination in the correct order we need to
    !  exchange the current row number with the one that has
    !  smallest column number, among JJ, JJ+1, ..., LENL.
    !
        jrow = jwl(jj)
        k = jj
    !
    !  Determine the smallest column index.
    !
        do j = jj+1, lenl
           if ( jwl(j) < jrow ) then
              jrow = jwl(j)
              k = j
           end if
        end do
    !
    !  Exchange in JWL.
    !
        j = jwl(jj)
        jwl(jj) = jrow
        jwl(k) = j
    !
    !  Exchange in JR.
    !
        jr(jrow) = jj
        jr(j) = k
    !
    !  Exchange in WL.
    !
        s = wl(k)
        wl(k) = wl(jj)
        wl(jj) = s

        if ( ii <= jrow ) then
          go to 160
        end if
    !
    !  Get the multiplier for row to be eliminated: JROW.
    !
        fact = wl(jj) * alu(jrow)
        jr(jrow) = 0

        if ( abs ( fact ) * wu(n+2-jrow) <= tol * tnorm ) then
          go to 150
        end if
    !
    !  Combine current row and row JROW.
    !
        do k = ju(jrow), jlu(jrow+1)-1
           s = fact * alu(k)
           j = jlu(k)
           jpos = jr(j)
    !
    !  If fill-in element and small disregard.
    !
           if ( abs ( s ) < tol * tnorm .and. jpos == 0 ) then
             cycle
           end if

           if ( ii <= j ) then
    !
    !  Dealing with upper part.
    !
              if ( jpos == 0 ) then
    !
    !  This is a fill-in element.
    !
                 lenu = lenu + 1

                 if ( n < lenu ) then
                   go to 995
                 end if

                 jwu(lenu) = j
                 jr(j) = lenu
                 wu(lenu) = - s
              else
    !
    !  No fill-in element.
    !
                 wu(jpos) = wu(jpos) - s
              end if
           else
    !
    !  Dealing with lower part.
    !
              if ( jpos == 0 ) then
    !
    !  This is a fill-in element.
    !
                 lenl = lenl + 1

                 if ( n < lenl ) then
                   go to 995
                 end if

                 jwl(lenl) = j
                 jr(j) = lenl
                 wl(lenl) = -s
              else
    !
    !  No fill-in element.
    !
                 wl(jpos) = wl(jpos) - s
              end if
           end if

      end do

        nl = nl + 1
        wl(nl) = fact
        jwl(nl) = jrow
      go to 150
    !
    !  Update the L matrix.
    !
     160 continue

        len = min ( nl, lenl0 + lfil )

        call bsort2 ( wl, jwl, nl, len )

        do k = 1, len

           if ( iwk < ju0 ) then
             ierr = -2
             return
           end if

           alu(ju0) =  wl(k)
           jlu(ju0) =  jwl(k)
           ju0 = ju0 + 1

        end do
    !
    !  Save pointer to beginning of row II of U.
    !
        ju(ii) = ju0
    !
    !  Reset double pointer JR to zero (L-part - except first
    !  JJ-1 elements which have already been reset).
    !
      do k = jj, lenl
        jr(jwl(k)) = 0
      end do
    !
    !  Be sure that the diagonal element is first in W and JW.
    !
        idiag = jr(ii)

        if ( idiag == 0 ) then
          go to 900
        end if

        if ( idiag /= 1 ) then

           s = wu(1)
           wu(j) = wu(idiag)
           wu(idiag) = s

           j = jwu(1)
           jwu(1) = jwu(idiag)
           jwu(idiag) = j

        end if

        len = min ( lenu, lenu0 + lfil )

        call bsort2 ( wu(2), jwu(2), lenu-1, len )
    !
    ! Update the U-matrix.
    !
        t = 0.0D+00

        do k = 2, len

           if ( iwk < ju0 ) then
             ierr = -3
             return
           end if

           jlu(ju0) = jwu(k)
           alu(ju0) = wu(k)
           t = t + abs ( wu(k) )
           ju0 = ju0 + 1

        end do
    !
    !  Save norm in WU (backwards). Norm is in fact average absolute value.
    !
        wu(n+2-ii) = t / real ( len + 1, kind = dp )
    !
    !  Store inverse of diagonal element of U.
    !
        if ( wu(1) == 0.0D+00 ) then
          ierr = -5
          return
        end if

        alu(ii) = 1.0D+00 / wu(1)
    !
    !  Update pointer to beginning of next row of U.
    !
      jlu(ii+1) = ju0
    !
    !  Reset double pointer JR to zero (U-part).
    !
      do k = 1, lenu
        jr(jwu(k)) = 0
      end do

      end do

      ierr = 0

      return
    !
    !  Zero pivot :
    !
     900    ierr = ii
        return
    !
    !  Incomprehensible error. Matrix must be wrong.
    !
     995    ierr = -1
        return
  end subroutine


  subroutine bsort2 ( w, ind, n, ncut )

  !*****************************************************************************80
  !
  !! BSORT2 returns the NCUT largest elements of an array, using bubble sort.
  !
  !  Discussion:
  !
  !    This routine carries out a simple bubble sort for getting the NCUT largest
  !    elements in modulus, in array W.  IND is sorted accordingly.
  !    (Ought to be replaced by a more efficient sort especially
  !    if NCUT is not that small).
  !
  !  Modified:
  !
  !    07 January 2004
  !
  !  Author:
  !
  !    Youcef Saad
  !
  !  Parameters:
  !
      implicit none

      integer ( kind = il ) n

      integer ( kind = il ) i
      integer ( kind = il ) ind(*)
      integer ( kind = il ) iswp
      integer ( kind = il ) j
      integer ( kind = il ) ncut
      logical test
      real ( kind = dp ) w(n)
      real ( kind = dp ) wswp

      i = 1

      do

        test = .false.

        do j = n-1, i, -1

          if ( abs ( w(j) ) < abs ( w(j+1) ) ) then
    !
    !  Swap.
    !
            wswp = w(j)
            w(j) = w(j+1)
            w(j+1) = wswp
    !
    !  Reorder the original ind array accordingly.
    !
            iswp = ind(j)
            ind(j) = ind(j+1)
            ind(j+1) = iswp
    !
    !  Set indicator that sequence is still unsorted.
    !
            test = .true.

          end if

        end do

        i = i + 1

        if ( .not. test .or. ncut < i ) then
          exit
        end if

      end do

      return
  end subroutine




  subroutine bsort2C ( w, ind, n, ncut )

  !*****************************************************************************80
  !
  !! BSORT2 returns the NCUT largest elements of an array, using bubble sort.
  !
  !  Discussion:
  !
  !    This routine carries out a simple bubble sort for getting the NCUT largest
  !    elements in modulus, in array W.  IND is sorted accordingly.
  !    (Ought to be replaced by a more efficient sort especially
  !    if NCUT is not that small).
  !
  !  Modified:
  !
  !    07 January 2004
  !
  !  Author:
  !
  !    Youcef Saad
  !
  !  Parameters:
  !
      implicit none

      integer ( kind = ild ) n

      integer ( kind = ild ) i
      integer ( kind = ild ) ind(*)
      integer ( kind = ild ) iswp
      integer ( kind = ild ) j
      integer ( kind = ild ) ncut
      logical test
      complex ( kind = dp ) w(n)
      complex ( kind = dp ) wswp

      i = 1

      do

        test = .false.

        do j = n-1, i, -1

          if ( abs ( w(j) ) < abs ( w(j+1) ) ) then
    !
    !  Swap.
    !
            wswp = w(j)
            w(j) = w(j+1)
            w(j+1) = wswp
    !
    !  Reorder the original ind array accordingly.
    !
            iswp = ind(j)
            ind(j) = ind(j+1)
            ind(j+1) = iswp
    !
    !  Set indicator that sequence is still unsorted.
    !
            test = .true.

          end if

        end do

        i = i + 1

        if ( .not. test .or. ncut < i ) then
          exit
        end if

      end do

      return
  end subroutine





















































  !COMPLEX VERSIONS
  !********************************************************+

  subroutine ilu0C ( n, a, ja, ia, alu, jlu, ju, iw, ierr )

!*****************************************************************************80
!
!! ILU0 is an ILU(0) preconditioner.
!
!  Discussion:
!
!    Note that this has been coded in such a way that it can be used
!    with PGMRES.  Normally, since the data structure of a, ja, ia is
!    the same as that of a, ja, ia, savings can be made. In fact with
!    some definitions (not correct for general sparse matrices) all we
!    need in addition to a, ja, ia is an additional diagonal.
!    Ilu0 is not recommended for serious problems. It is only provided
!    here for comparison purposes.
!
!    It is assumed that the the elements in the input matrix are stored
!    in such a way that in each row the lower part comes first and
!    then the upper part. To get the correct ILU factorization, it is
!    also necessary to have the elements of L sorted by increasing
!    column number. It may therefore be necessary to sort the
!    elements of a, ja, ia prior to calling ilu0. This can be
!    achieved by transposing the matrix twice using csrcsc.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer ( kind = il ) N, the order of the matrix.
!
!    Input, real A(*), integer ( kind = il ) JA(*), IA(N+1), the matrix in CSR
!    Compressed Sparse Row format.
!
! on return:
!
! alu,jlu = matrix stored in Modified Sparse Row (MSR) format containing
!           the L and U factors together. The diagonal (stored in
!           alu(1:n) ) is inverted. Each i-th row of the alu,jlu matrix
!           contains the i-th row of L (excluding the diagonal entry=1)
!           followed by the i-th row of U.
!
! ju        = pointer to the diagonal elements in alu, jlu.
!
! ierr        = integer ( kind = il ) indicating error code on return
!           ierr = 0 --> normal return
!           ierr = k --> code encountered a zero pivot at step k.
! work arrays:
!
! iw          = integer ( kind = il ) work array of length n.
!
!
!ilu0C ( n,       a,      ja, ia,     alu,        jlu,         ju,      iw, ierr )
!ilu0C(num_e, Zesparcida, ja, ia, precond_alu, precond_jlu, precond_ju, iw, ierr)
  implicit none

  integer ( kind = il ) n

  complex ( kind = dp_near ) a(*)
  complex ( kind = dp_prec ) alu(*)
  integer ( kind = ild ) i
  integer ( kind = ild ) ia(n+1)
  integer ( kind = il ) ierr
  integer ( kind = il ) ii
  integer ( kind = ild ) iw(n)
  integer ( kind = ild ) j
  integer ( kind = il ) ja(*)
  integer ( kind = ild ) jcol
  integer ( kind = ild ) jf
  integer ( kind = ild ) jj
  integer ( kind = il ) jlu(*)
  integer ( kind = ild ) jm
  integer ( kind = ild ) jrow
  integer ( kind = ild ) js
  integer ( kind = ild ) ju(*)
  integer ( kind = ild ) ju0
  integer ( kind = ild ) jw
  complex ( kind = dp ) tl
  !integer (kind = il) :: prc, prca

  !prca = -1
  call setprc(n)
  ju0 = n + 2
  jlu(1) = ju0
!
!  Initialize the work vector.
!
  iw(1:n) = 0
!
!  The main loop.
!
  do ii = 1, n

    js = ju0
!
!  Generating row II of L and U.
!
    do j = ia(ii), ia(ii+1)-1
!
!  Copy row II of A, JA, IA into row II of ALU, JLU (L/U) matrix.
!
      jcol = ja(j)

      if ( jcol == ii ) then
        alu(ii) = a(j)
        iw(jcol) = ii
        ju(ii) = ju0
      else
        alu(ju0) = a(j)
        jlu(ju0) = ja(j)
        iw(jcol) = ju0
        ju0 = ju0 + 1
      end if

    end do

    jlu(ii+1) = ju0
    jf = ju0 - 1
    jm = ju(ii) - 1
!
!  Exit if the diagonal element is reached.
!
    do j = js, jm

      jrow = jlu(j)
      tl = alu(j) * alu(jrow)
      alu(j) = tl
!
!  Perform linear combination.
!
      do jj = ju(jrow), jlu(jrow+1)-1
        jw = iw(jlu(jj))
        if ( jw /= 0 ) then
          alu(jw) = alu(jw) - tl * alu(jj)
        end if
      end do

    end do
!
!  Invert and store the diagonal element.
!
    if ( abs(alu(ii)) == 0.0D+00 ) then
      ierr = ii
      return
    end if

    alu(ii) = 1.0D+00 / alu(ii)
!
!  Reset pointer IW to zero.
!
    iw(ii) = 0
    do i = js, jf
      iw(jlu(i)) = 0
    end do

    call updateprc(ii)
    !prc = 100*ii/n
    !if (prc/2 /= prca/2) then
    !    print*, inttostr(prc) // ' %'
    !    prca = prc
    !end if

  end do

  ierr = 0
  return
end subroutine


subroutine coocsrC ( nrow, nnz, a, ir, jc, ao, jao, iao )

!*****************************************************************************80
!
!! COOCSR converts COO to CSR.
!
!  Discussion:
!
!    This routine converts a matrix that is stored in COO coordinate format
!    a, ir, jc into a CSR row general sparse ao, jao, iao format.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer ( kind = il ) NROW, the row dimension of the matrix.
!
!    Input, integer ( kind = il ) NNZ, the number of nonzero elements.
!
! a,
! ir,
! jc    = matrix in coordinate format. a(k), ir(k), jc(k) store the nnz
!         nonzero elements of the matrix with a(k) = actual real value of
!         the elements, ir(k) = its row number and jc(k) = its column
!        number. The order of the elements is arbitrary.
!
! on return:
!
! ir       is destroyed
!
!    Output, real AO(*), JAO(*), IAO(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!coocsrC ( nrow ,     nnz,        a,       ir,     jc,         ao,   jao, iao )
!coocsrC ( num_e, n_esparcida, a_znear, i_znear, j_znear, Zesparcida, ja, ia )
  implicit none

  integer ( kind = il ) nrow

  complex ( kind = dp ) a(*)
  complex ( kind = dp_near ) ao(*)
  integer ( kind = il ) i
  integer ( kind = il ) iad
  integer ( kind = il ) iao(nrow+1)
  integer ( kind = il ) ir(*)
  integer ( kind = il ) j
  integer ( kind = il ) jao(*)
  integer ( kind = il ) jc(*)
  integer ( kind = il ) k
  integer ( kind = il ) k0
  integer ( kind = il ) nnz
  complex ( kind = dp ) x

  iao(1:nrow+1) = 0
!
!  Determine the row lengths.
!

  do k = 1, nnz
    iao(ir(k)) = iao(ir(k)) + 1
  end do

!
!  The starting position of each row.
!

  k = 1
  do j = 1, nrow+1
     k0 = iao(j)
     iao(j) = k
     k = k + k0
  end do

!
!  Go through the structure once more.  Fill in output matrix.
!

  do k = 1, nnz
     i = ir(k)
     j = jc(k)
     x = a(k)
     iad = iao(i)
     ao(iad) = x
     jao(iad) = j
     iao(i) = iad + 1
  end do

!
!  Shift back IAO.
!
  do j = nrow, 1, -1
    iao(j+1) = iao(j)
  end do
  iao(1) = 1


  return
end subroutine



subroutine csrcooC ( nrow, job, nzmax, a, ja, ia, nnz, ao, ir, jc, ierr )

!*****************************************************************************80
!
!! CSRCOO converts Compressed Sparse Row to Coordinate format.
!
!  Discussion:
!
!   This routine converts a matrix that is stored in row general sparse
!   A, JA, IA format into coordinate format AO, IR, JC.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer ( kind = il ) NROW, the row dimension of the matrix.
!
! job   = integer ( kind = il ) serving as a job indicator.
!         if job = 1 fill in only the array ir, ignore jc, and ao.
!         if job = 2 fill in ir, and jc but not ao
!         if job = 3 fill in everything.
!         The reason why these options are provided is that on return
!         ao and jc are the same as a, ja. So when job = 3, a and ja are
!         simply copied into ao, jc.  When job=2, only jc and ir are
!         returned. With job=1 only the array ir is returned. Moreover,
!         the algorithm is in place:
!           call csrcoo (nrow,1,nzmax,a,ja,ia,nnz,a,ia,ja,ierr)
!         will write the output matrix in coordinate format on a, ja,ia.
!         (Important: note the order in the output arrays a, ja, ia. )
!         i.e., ao can be the same as a, ir can be the same as ia
!         and jc can be the same as ja.
!
!    Input, real A(*), integer ( kind = il ) JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
! nzmax = length of space available in ao, ir, jc.
!         the code will stop immediatly if the number of
!         nonzero elements found in input matrix exceeds nzmax.
!
! on return:
!-
! ao, ir, jc = matrix in coordinate format.
!
! nnz        = number of nonzero elements in matrix.
!
! ierr       = integer ( kind = il ) error indicator.
!         ierr == 0 means normal retur
!         ierr == 1 means that the the code stopped
!         because there was no space in ao, ir, jc
!         (according to the value of  nzmax).
!
  implicit none

  integer ( kind = il ) nrow

  complex ( kind = dp_near ) a(*)
  complex ( kind = dp_near ) ao(*)
  integer ( kind = il ) i
  integer ( kind = il ) ia(nrow+1)
  integer ( kind = il ) ierr
  integer ( kind = il ) ir(*)
  integer ( kind = il ) ja(*)
  integer ( kind = il ) jc(*)
  integer ( kind = il ) job
  integer ( kind = il ) k
  integer ( kind = il ) k1
  integer ( kind = il ) k2
  integer ( kind = il ) nnz
  integer ( kind = il ) nzmax

  ierr = 0
  nnz = ia(nrow+1)-1

  if ( nzmax < nnz ) then
    ierr = 1
    return
  end if

  if ( 3 <= job ) then
    ao(1:nnz) = a(1:nnz)
  end if

  if ( 2 <= job ) then
    jc(1:nnz) = ja(1:nnz)
  end if
!
!  Copy backward.
!
  do i = nrow, 1, -1
    k1 = ia(i+1) - 1
    k2 = ia(i)
    do k = k1, k2, -1
      ir(k) = i
    end do
  end do

  return
end subroutine


subroutine msrcsrC ( n, a, ja, ao, jao, iao, wk )

!*****************************************************************************80
!
!! MSRCSR converts Modified Sparse Row to Compressed Sparse Row.
!
!  Discussion:
!
!    This routine converts a compressed matrix using a separated diagonal
!    (modified sparse row format) in the Compressed Sparse Row format.
!
!    does not check for zero elements in the diagonal.
!
!    This is an "in place" algorithm (see a, ja, ia).
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer ( kind = il ) N, the row dimension of the matrix.
!
! ao, jao  = sparse matrix in msr sparse storage format
!           see routine csrmsr for details
!
! on return :
!
! a, ja, ia = matrix in csr format. note that the
!           algorithm is in place: ao, jao can be the same
!            as a, ja, in which case it will be overwritten on it
!            upon return.
!
!             here nnz = number of nonzero elements+1
!
!    Workspace, real WK(N).
!
  implicit none

  integer ( kind = il ) n

  complex ( kind = dp_prec ) a(*)
  logical added
  complex ( kind = dp ) ao(*)
  integer ( kind = il ) iao(n+1)
  integer ( kind = il ) idiag
  integer ( kind = il ) ii
  integer ( kind = il ) iptr
  integer ( kind = il ) j
  integer ( kind = il ) ja(*)
  integer ( kind = il ) jao(*)
  integer ( kind = il ) k
  complex ( kind = dp ) wk(n)

  wk(1:n) = a(1:n)

  iao(1) = 1
  iptr = 1

  do ii = 1, n

    added = .false.
    idiag = iptr + ( ja(ii+1) - ja(ii) )

    do k = ja(ii), ja(ii+1)-1

      j = ja(k)

      if ( j < ii ) then
        ao(iptr) = a(k)
        jao(iptr) = j
        iptr = iptr + 1
      else if ( added ) then
        ao(iptr) = a(k)
        jao(iptr) = j
        iptr = iptr + 1
      else
!
!  Add diagonal element.  Only reserve a position for it.
!
        idiag = iptr
        iptr = iptr + 1
        added = .true.
!
!  Then other elements.
!
        ao(iptr) = a(k)
        jao(iptr) = j
        iptr = iptr + 1
      end if

    end do

    ao(idiag) = wk(ii)
    jao(idiag) = ii
    if ( .not. added ) then
      iptr = iptr + 1
    end if
    iao(ii+1) = iptr

  end do

  return
end subroutine



  subroutine csrmsrC ( n, a, ja, ia, ao, jao, wk, iwk )

  !*****************************************************************************80
  !
  !! CSRMSR converts Compressed Sparse Row to Modified Sparse Row.
  !
  !  Discussion:
  !
  !    This routine converts a general sparse matrix a, ja, ia into
  !    a compressed matrix using a separated diagonal (referred to as
  !    the bell-labs format as it is used by bell labs semi conductor
  !    group. We refer to it here as the modified sparse row format.
  !
  !    This has been coded in such a way that one can overwrite
  !    the output matrix onto the input matrix if desired by a call of
  !    the form
  !
  !     call csrmsr (n, a, ja, ia, a, ja, wk,iwk)
  !
  !    In case ao, jao, are different from a, ja, then one can
  !    use ao, jao as the work arrays in the calling sequence:
  !
  !     call csrmsr (n, a, ja, ia, ao, jao, ao,jao)
  !
  !    Algorithm is in place.  i.e. both:
  !
  !          call csrmsr (n, a, ja, ia, ao, jao, ao,jao)
  !          (in which  ao, jao, are different from a, ja)
  !           and
  !          call csrmsr (n, a, ja, ia, a, ja, wk,iwk)
  !          (in which  wk, jwk, are different from a, ja)
  !        are OK.
  !
  !  Modified:
  !
  !    07 January 2004
  !
  !  Author:
  !
  !    Youcef Saad
  !
  !  Parameters:
  !
  !    Input, integer ( kind = il ) N, the order of the matrix.
  !
  !    Input, real A(*), integer ( kind = il ) JA(*), IA(N+1), the matrix in CSR
  !    Compressed Sparse Row format.
  !
  ! on return :
  !
  ! ao, jao  = sparse matrix in modified sparse row storage format:
  !         +  ao(1:n) contains the diagonal of the matrix.
  !         +  ao(n+2:nnz) contains the nondiagonal elements of the
  !             matrix, stored rowwise.
  !         +  jao(n+2:nnz) : their column indices
  !         +  jao(1:n+1) contains the pointer array for the nondiagonal
  !             elements in ao(n+1:nnz) and jao(n+2:nnz).
  !             i.e., for i <= n+1 jao(i) points to beginning of row i
  !            in arrays ao, jao.
  !             here nnz = number of nonzero elements+1
  !
  !    Work array, real WK(N).
  !
  !    Work array, integer ( kind = il ) IWK(N+1).
  !
      implicit none

      integer ( kind = il ) n

      complex ( kind = dp ) a(*)
      complex ( kind = dp ) ao(*)
      integer ( kind = il ) i
      integer ( kind = il ) ia(n+1)
      integer ( kind = il ) icount
      integer ( kind = il ) ii
      integer ( kind = il ) iptr
      integer ( kind = il ) iwk(n+1)
      integer ( kind = il ) j
      integer ( kind = il ) ja(*)
      integer ( kind = il ) jao(*)
      integer ( kind = il ) k
      complex ( kind = dp ) wk(n)

      icount = 0
    !
    !  Store away diagonal elements and count nonzero diagonal elements.
    !
      do i = 1, n
        wk(i) = 0.0D+00
        iwk(i+1) = ia(i+1) - ia(i)
        do k = ia(i), ia(i+1)-1
          if ( ja(k) == i ) then
            wk(i) = a(k)
            icount = icount + 1
            iwk(i+1) = iwk(i+1) - 1
          end if
        end do
      end do
    !
    !  Compute total length.
    !
      iptr = n + ia(n+1) - icount
    !
    !  Copy backwards, to avoid collisions.
    !
      do ii = n, 1, -1
        do k = ia(ii+1)-1, ia(ii), -1
          j = ja(k)
          if ( j /= ii ) then
            ao(iptr) = a(k)
            jao(iptr) = j
            iptr = iptr - 1
          end if
        end do
      end do
    !
    !  Compute the pointer values and copy WK.
    !
      jao(1) = n + 2
      do i = 1, n
        ao(i) = wk(i)
        jao(i+1) = jao(i) + iwk(i+1)
      end do

      return
  end subroutine



  subroutine ilutC ( n, a, ja, ia, lfil, tol, alu, jlu, ju, iwk, wu, wl, jr, jwl, jwu, ierr )
!ilutC (num_e, a_near, ja_near, ia_near, fill_in, drop_tol, precond_alu, precond_jlu, precond_ju, tamano_esperado, wu, wl, jr, jwl, jwu, status )
  !*****************************************************************************80
  !
  !! ILUT is an ILUT preconditioner.
  !
  !  Discussion:
  !
  !    This routine carries ouot incomplete LU factorization with dual
  !    truncation mechanism.  Sorting is done for both L and U.
  !
  !    The dual drop-off strategy works as follows:
  !
  !    1) Theresholding in L and U as set by TOL.  Any element whose size
  !       is less than some tolerance (relative to the norm of current
  !       row in u) is dropped.
  !
  !    2) Keeping only the largest lenl0+lfil elements in L and the
  !       largest lenu0+lfil elements in U, where lenl0=initial number
  !       of nonzero elements in a given row of lower part of A
  !       and lenlu0 is similarly defined.
  !
  !    Flexibility: one can use tol=0 to get a strategy based on keeping the
  !    largest elements in each row of L and U. Taking tol /= 0 but lfil=n
  !    will give the usual threshold strategy (however, fill-in is then
  !    unpredictible).
  !
  !    A must have all nonzero diagonal elements.
  !
  !  Modified:
  !
  !    07 January 2004
  !
  !  Author:
  !
  !    Youcef Saad
  !
  !  Parameters:
  !
  !  Parameters:
  !
  !    Input, integer ( kind = il ) N, the order of the matrix.
  !
  !    Input, real A(*), integer ( kind = il ) JA(*), IA(N+1), the matrix in CSR
  !    Compressed Sparse Row format.
  !
  ! lfil    = integer ( kind = il ). The fill-in parameter. Each row of L and
  !           each row of U will have a maximum of lfil elements
  !           in addition to the original number of nonzero elements.
  !           Thus storage can be determined beforehand.
  !           lfil must be >= 0.
  !
  ! iwk     = integer ( kind = il ). The minimum length of arrays alu and jlu
  !
  ! On return:
  !
  ! alu,jlu = matrix stored in Modified Sparse Row (MSR) format containing
  !           the L and U factors together. The diagonal (stored in
  !           alu(1:n) ) is inverted. Each i-th row of the alu,jlu matrix
  !           contains the i-th row of L (excluding the diagonal entry=1)
  !           followed by the i-th row of U.
  !
  ! ju      = integer ( kind = il ) array of length n containing the pointers to
  !           the beginning of each row of U in the matrix alu,jlu.
  !
  ! ierr    = integer ( kind = il ). Error message with the following meaning.
  !           ierr  = 0    --> successful return.
  !           ierr > 0  --> zero pivot encountered at step number ierr.
  !           ierr  = -1   --> Error. input matrix may be wrong.
  !                            (The elimination process has generated a
  !                            row in L or U whose length is >  n.)
  !           ierr  = -2   --> The matrix L overflows the array al.
  !           ierr  = -3   --> The matrix U overflows the array alu.
  !           ierr  = -4   --> Illegal value for lfil.
  !           ierr  = -5   --> zero pivot encountered.
  !
  ! work arrays:
  !
  ! jr,jwu,jwl, integer ( kind = il ) work arrays of length n.
  ! wu, wl, real work arrays of length n+1, and n resp.
  !
  !ilutC ( n, a, ja, ia, lfil, tol, alu, jlu, ju, iwk, wu, wl, jr, jwl, jwu, ierr )
  !
      implicit none

      integer ( kind = il ) n

      complex ( kind = dp_near ) a(*)
      complex ( kind = dp_prec ) alu(*)
      complex ( kind = dp ) fact
      integer ( kind = ild ) ia(n+1)
      integer ( kind = ild ) idiag
      integer ( kind = il ) ierr
      integer ( kind = il ) ii
      integer ( kind = ild ) iwk
      integer ( kind = ild ) j
      integer ( kind = ild ) j1
      integer ( kind = ild ) j2
      integer ( kind = il ) ja(*)
      integer ( kind = ild ) jjj
      integer ( kind = il ) jlu(*)
      integer ( kind = ild ) jpos
      integer ( kind = ild ) jr(*)
      integer ( kind = ild ) jrow
      integer ( kind = ild ) ju(*)
      integer ( kind = ild ) ju0
      integer ( kind = ild ) jwl(n)
      integer ( kind = ild ) jwu(n)
      integer ( kind = ild ) k
      integer ( kind = ild ) len
      integer ( kind = ild ) lenl
      integer ( kind = ild ) lenl0
      integer ( kind = ild ) lenu
      integer ( kind = ild ) lenu0
      integer ( kind = il ) lfil
      integer ( kind = ild ) nl
      complex ( kind = dp ) s!
      complex ( kind = dp ) t!
      real ( kind = dp ) tnorm!
      real ( kind = dp ) tol!
      complex ( kind = dp ) wl(n)!
      complex ( kind = dp ) wu(n+1)!
      !integer (kind = il) :: prc, prca

      !prca = -1
      call setprc(n)
      if ( lfil < 0 ) then
        ierr = -4
        return
      end if
    !
    !  Initialize JU0 (points to next element to be added to ALU, JLU)
    !  and pointer.
    !
      ju0 = n + 2
      jlu(1) = ju0
    !
    !  integer ( kind = il ) double pointer array.
    !
      jr(1:n) = 0
    !
    !  The main loop.
    !
      do ii = 1, n
!print*, 'ii: ', ii
        j1 = ia(ii)
        j2 = ia(ii+1) - 1
        lenu = 0
        lenl = 0

        tnorm = 0.0D+00
        do k = j1, j2
          tnorm = tnorm + abs ( a(k) )
        end do
        tnorm = tnorm / real ( j2-j1+1, kind = dp )
    !
    !  Unpack L-part and U-part of row of A in arrays WL, WU.
    !
        do j = j1, j2
!print*, 'j: ', j
          k = ja(j)
          t = a(j)

          if ( tol * tnorm <= abs ( t ) ) then

            if ( k < ii ) then
              lenl = lenl + 1
              jwl(lenl) = k
              wl(lenl) = t
              jr(k) = lenl
            else
              lenu = lenu+1
              jwu(lenu) = k
              wu(lenu) = t
!print*, 'lenu: ', lenu
!print*, 'wu(lenu):', wu(lenu)
              jr(k) = lenu
            end if

          end if

        end do
!print*, '111todo wu:'
!print*, wu
        lenl0 = lenl
        lenu0 = lenu
        jjj = 0
        nl = 0
    !
    !  Eliminate previous rows.
    !

!print*, 'before continue'
    150 continue
!print*, 'after continue'
        jjj = jjj + 1
!print*, 'jjj: ', jjj
        if ( lenl < jjj ) then
          go to 160
        end if
    !
    !  In order to do the elimination in the correct order we need to
    !  exchange the current row number with the one that has
    !  smallest column number, among jjj, jjj+1, ..., LENL.
    !
        jrow = jwl(jjj)
        k = jjj
    !
    !  Determine the smallest column index.
    !
        do j = jjj+1, lenl
           if ( jwl(j) < jrow ) then
              jrow = jwl(j)
              k = j
           end if
        end do
    !
    !  Exchange in JWL.
    !
        j = jwl(jjj)
        jwl(jjj) = jrow
        jwl(k) = j
    !
    !  Exchange in JR.
    !
        jr(jrow) = jjj
        jr(j) = k
    !
    !  Exchange in WL.
    !
        s = wl(k)
        wl(k) = wl(jjj)
        wl(jjj) = s

        if ( ii <= jrow ) then
          go to 160
        end if
    !
    !  Get the multiplier for row to be eliminated: JROW.
    !
        fact = wl(jjj) * alu(jrow)
        jr(jrow) = 0
!print*, 'todo wu:'
!print*, wu
!print*, 'size: ', size(wu,1)
!print*, 'wu: ', wu(n+2-jrow)
!print*, 'n+2-jrow', n+2-jrow
!print*, 'n: ', n
!stop "fin"
        if ( abs ( fact ) * real(wu(n+2-jrow)) <= tol * tnorm ) then !ATENTION: if ( abs ( fact ) * wu(n+2-jrow) <= tol * tnorm )
            if ( abs ( fact ) * aimag(wu(n+2-jrow)) <= tol * tnorm ) then
                go to 150
            end if
        end if
    !
    !  Combine current row and row JROW.
    !
        do k = ju(jrow), jlu(jrow+1)-1
           s = fact * alu(k)
           j = jlu(k)
           jpos = jr(j)
    !
    !  If fill-in element and small disregard.
    !
           if ( abs ( s ) < tol * tnorm .and. jpos == 0 ) then
             cycle
           end if

           if ( ii <= j ) then
    !
    !  Dealing with upper part.
    !
              if ( jpos == 0 ) then
    !
    !  This is a fill-in element.
    !
                 lenu = lenu + 1

                 if ( n < lenu ) then
                   go to 995
                 end if

                 jwu(lenu) = j
                 jr(j) = lenu
!print*, '1lenu-n: ', lenu, n
                 wu(lenu) = - s
              else
    !
    !  No fill-in element.
    !
!print*, '1jpos-n: ', jpos, n
                 wu(jpos) = wu(jpos) - s
              end if
           else
    !
    !  Dealing with lower part.
    !
              if ( jpos == 0 ) then
    !
    !  This is a fill-in element.
    !
                 lenl = lenl + 1

                 if ( n < lenl ) then
                   go to 995
                 end if

                 jwl(lenl) = j
                 jr(j) = lenl
                 wl(lenl) = -s
              else
    !
    !  No fill-in element.
    !
                 wl(jpos) = wl(jpos) - s
              end if
           end if

      end do

        nl = nl + 1
        wl(nl) = fact
        jwl(nl) = jrow
      go to 150
    !
    !  Update the L matrix.
    !
     160 continue

        len = min ( nl, lenl0 + lfil )

        call bsort2C ( wl, jwl, nl, len )

        do k = 1, len

           if ( iwk < ju0 ) then
             ierr = -2
             return
           end if

           alu(ju0) =  wl(k)
           jlu(ju0) =  jwl(k)
           ju0 = ju0 + 1

        end do
    !
    !  Save pointer to beginning of row II of U.
    !
        ju(ii) = ju0
    !
    !  Reset double pointer JR to zero (L-part - except first
    !  jjj-1 elements which have already been reset).
    !
      do k = jjj, lenl
        jr(jwl(k)) = 0
      end do
    !
    !  Be sure that the diagonal element is first in W and JW.
    !
        idiag = jr(ii)

        if ( idiag == 0 ) then
          go to 900
        end if
!print*, 'pero si llego aqui'
        if ( idiag /= 1 ) then
!print*, '1j-n: ', j, n
!!print*, '1idiag-n: ', idiag, n
           s = wu(1)
!print*, j , ' ', j1, ' ', j2          
           wu(j) = wu(idiag)
           wu(idiag) = s

           j = jwu(1)
           jwu(1) = jwu(idiag)
           jwu(idiag) = j

        end if

        len = min ( lenu, lenu0 + lfil )

        call bsort2C ( wu(2), jwu(2), lenu-1, len )
    !
    ! Update the U-matrix.
    !
        t = 0.0D+00

        do k = 2, len

           if ( iwk < ju0 ) then
             ierr = -3
             return
           end if

           jlu(ju0) = jwu(k)
           alu(ju0) = wu(k)
           t = t + abs ( wu(k) )
           ju0 = ju0 + 1

        end do
    !
    !  Save norm in WU (backwards). Norm is in fact average absolute value.
    !

!print*, 'n+2-ii-n: ', n+2-ii, n    
        wu(n+2-ii) = t / real ( len + 1, kind = dp )
!print*, 'wu(n+2-ii)', wu(n+2-ii)
    !
    !  Store inverse of diagonal element of U.
    !
        if ( abs(wu(1)) == 0.0D+00 ) then
          ierr = -5
          return
        end if

        alu(ii) = 1.0D+00 / wu(1)
    !
    !  Update pointer to beginning of next row of U.
    !
      jlu(ii+1) = ju0
    !
    !  Reset double pointer JR to zero (U-part).
    !
      do k = 1, lenu
        jr(jwu(k)) = 0
      end do

        call updateprc(ii)
        !prc = 100*ii/n
        !if (prc/2 /= prca/2) then
        !    print*, inttostr(prc) // ' %'
        !    prca = prc
        !end if

        !print*, ii, ' de ', n
      end do

      ierr = 0

      return
    !
    !  Zero pivot :
    !
     900    ierr = ii
        return
    !
    !  Incomprehensible error. Matrix must be wrong.
    !
     995    ierr = -1


        return
  end subroutine


subroutine amuxC ( n, x, y, a, ja, ia )

!*****************************************************************************80
!
!! AMUX multiplies a CSR matrix A times a vector.
!
!  Discussion:
!
!    This routine multiplies a matrix by a vector using the dot product form.
!    Matrix A is stored in compressed sparse row storage.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the row dimension of the matrix.
!
!    Input, real X(*), and array of length equal to the column dimension 
!    of A.
!
!    Input, real A(*), integer ( kind = 4 ) JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!    Output, real Y(N), the product A * X.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) a(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia(*)
  integer ( kind = 4 ) ja(*)
  integer ( kind = 4 ) k
  complex ( kind = 8 ) t
  complex ( kind = 8 ) x(*)
  complex ( kind = 8 ) y(n)

  do i = 1, n
!
!  Compute the inner product of row I with vector X.
!
    t = 0.0D+00
    do k = ia(i), ia(i+1)-1
      t = t + a(k) * x(ja(k))
    end do

    y(i) = t

  end do

  return
end subroutine





end module modIterativo
