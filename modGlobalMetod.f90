module modGlobalMetod
    use modGlobalParam
    implicit none

public
private :: prev_prc, prc_max
    integer (kind = il) :: prev_prc
    integer (kind = il) :: prc_max
    character (len = 100), allocatable, dimension(:) :: logreg
    integer (kind = il) :: logcount = 0
    integer (kind = il) :: loginc = 1000
contains

    subroutine printconsole(texto, textoGUI, repeat)
        character (len = *), intent(in) :: texto
        character (len = *), optional, intent(in) :: textoGUI
        logical, optional, intent(in) :: repeat

        if (conf%optionGUI) then
            if (present(textoGUI)) then
                if (present(repeat)) then
                    if (repeat) then
                        print*, '#' // textoGUI // '#' // texto
                    else
                        print*, '#' // textoGUI
                    end if
                else
                    print*, '#' // textoGUI
                end if
                
            else
                print*, texto
            end if
        else
            print*, texto
        end if
        
    end subroutine

    subroutine savelog(fname)
        !dummy
        character (len = *), intent(in) :: fname
        !
        !local
        integer (kind = il) :: i, iunit
        !

    if (.false.) then
        call getunit ( iunit )
        open(UNIT=iunit, FILE = fname, ACTION="write", STATUS="replace")


        do i=1,logcount
            write(iunit,'(1X,A50)')  logreg(i)
            !101 FORMAT('1X,A50')
        end do
        close(UNIT = iunit)
        deallocate(logreg)
    end if
    end subroutine

    subroutine savetolog(texto, impr)
        !dummy
        character (len = *), intent(in) :: texto
        logical, optional, intent(in) :: impr
        !
        !local
        character (len = 100), allocatable, dimension(:) :: logregtemp
        !

        if (present(impr)) then
            if (impr) then
                print*, texto
            end if
        end if

    if (.false.) then
        if (allocated(logreg)) then
        else
            allocate(logreg(loginc))
        end if

        if (logcount >= size(logreg)) then
            allocate(logregtemp(size(logreg)))
            logregtemp = logreg
            deallocate(logreg)
            allocate(logreg(size(logregtemp)+loginc))
            logreg = logregtemp
            deallocate(logregtemp)
        end if

        logcount = logcount + 1
        !logreg(logcount) = '--------------------------------------------------'
        logreg(logcount)(1:len(trim(adjustl(texto)))) = trim(adjustl(texto))
        logreg(logcount)((len(trim(adjustl(texto))) + 1):) = ' '
    end if
    end subroutine

    subroutine printmsj(texto)
        character (len = *), intent(in) :: texto
        print*, texto
        call savetolog(texto)
    end subroutine


    subroutine setprc(maxim)
        !dummy
        integer (kind = il) :: maxim
        !
        prc_max = maxim
        call updateprc(0)
    end subroutine

    subroutine updateprc(valor)
        !dummy
        integer (kind = il) :: valor
        !

        
        !write(*,'(1h+,a,f5.1,a)') 'Progress: ', valor/prc_max*100, '%'


        OPEN (UNIT=6,FORM='FORMATTED')
        if (conf%optionGUI) then
            WRITE(6, 101) '#%/', REAL(valor)/REAL(prc_max)*100
        else
            WRITE(6, 100) REAL(valor)/REAL(prc_max)*100, ' %'
        end if
        100 FORMAT('+', F10.2, A3)
        101 FORMAT('+', A4, F10.2)

    end subroutine

    subroutine printInLine(texto)
        character(len=*), intent(in) :: texto

        OPEN (UNIT=6,FORM='FORMATTED')
        WRITE(6, 100) texto
        100 FORMAT('+', A)        

    end subroutine


    subroutine getcommand_console(comm, val)
        !dummy
        character (len = 256), intent(out) :: comm, val
        !
        !local
        character (len = 256) :: console
        character (len = 4) :: a
        !
        a = '>> '
        write(*,101,advance='no') a
        101 format(a)
        read(*,'(A256)'), console

        call getcommand(console,comm,val)

    end subroutine



    subroutine getcommand(str, comm, val)
        !dummy
        character (len = 256), intent(out) :: comm, val
        character (len = *), intent(in) :: str
        !
        !local
        integer (kind = ilo) :: ind
        character (len = 256) :: strtrim
        !
        strtrim = trim(adjustl(str))
        comm = trim(str)
        val = 'null'
        ind = index(strtrim, '=')
        if (ind /= 0) then
            comm = trim(adjustl(strtrim(1:(ind-1))))
            val = trim(adjustl(strtrim(ind+1:)))
        end if
        
    end subroutine


    function streq(str1,str2) result (res)
        !dummy
        character (len = *), intent(in) :: str1, str2
        logical :: res
        !
        if (to_upper(trim(str1)) == to_upper(trim(str2))) then
            res = .true.
        else
            res = .false.
        end if
    end function


        function to_upper(strIn) result(strOut)
        ! Adapted from http://www.star.le.ac.uk/~cgp/fortran.html (25 May 2012)

             implicit none

             character(len=*), intent(in) :: strIn
             character(len=len(strIn)) :: strOut
             integer :: i,j

             do i = 1, len(strIn)
                  j = iachar(strIn(i:i))
                  if (j>= iachar("a") .and. j<=iachar("z") ) then
                       strOut(i:i) = achar(iachar(strIn(i:i))-32)
                  else
                       strOut(i:i) = strIn(i:i)
                  end if
             end do

        end function to_upper



    !*  ************************************************************
    !*                          guardarVectorI
    !*  Guarda en un archivo (nombre_archivo) la matriz bidimensional
    !*  (matriz) de tipo real y con formato 'cuadrado' en el archivo
    !*  ************************************************************
    subroutine guardarVectorI(vector, nombre_archivo)
        !dummy
        integer (kind = il), intent(in), dimension(:) :: vector
        character (len = *),  intent(in) :: nombre_archivo
        !
        !local
        integer (kind = il) :: i
        integer (kind = il) :: iunit
        integer (kind = il) :: nfilas, ncolumnas
        !

        call getunit ( iunit )
        open(UNIT=iunit, FILE = nombre_archivo, ACTION="write", STATUS="replace")

        nfilas = size(vector,1)
        !ncolumnas = size(vector,2)

        do i=1,nfilas
            write(iunit,'(' // inttostr(1) // 'I10)') vector(i)
        end do
        close(UNIT = iunit)
    end subroutine



    !*  ************************************************************
    !*                          guardarVectorR
    !*  Guarda en un archivo (nombre_archivo) la matriz bidimensional
    !*  (matriz) de tipo real y con formato 'cuadrado' en el archivo
    !*  ************************************************************
    subroutine guardarVectorR(vector, nombre_archivo, fullprec)
        !dummy
        real (kind = dp), intent(in), dimension(:) :: vector
        character (len = *),  intent(in) :: nombre_archivo
        logical, intent(in), optional :: fullprec
        !
        !local
        integer (kind = il) :: i
        integer (kind = il) :: iunit
        integer (kind = il) :: nfilas, ncolumnas
        logical :: isFull
        !
        if (present(fullprec)) then
            if (fullprec) then
                isfull = .true.
            else
                isfull = .false.
            end if
        else
            isFull = .false.
        end if

        call getunit ( iunit )
        open(UNIT=iunit, FILE = nombre_archivo, ACTION="write", STATUS="replace")

        nfilas = size(vector,1)
        !ncolumnas = size(vector,2)

        do i=1,nfilas
            if (isfull) then
                write(iunit,'(*(G0.4,:))') vector(i)
            else
                write(iunit,'(' // inttostr(1) // 'f16.8)') vector(i)
            end if
        end do
        close(UNIT = iunit)
    end subroutine





    !*  ************************************************************
    !*							guardarMatrizR
    !*	Guarda en un archivo (nombre_archivo) la matriz bidimensional
    !*	(matriz) de tipo real y con formato 'cuadrado' en el archivo
    !*  ************************************************************
    subroutine guardarMatrizR(matriz, nombre_archivo, formatCSV)
        !dummy
        real (kind = dp), intent(in), dimension(:,:) :: matriz
        character (len = *),  intent(in) :: nombre_archivo
        logical, optional :: formatCSV
        !
        !local
        integer (kind = il) :: i
        integer (kind = il) :: iunit
        integer (kind = il) :: nfilas, ncolumnas
        logical :: isCSV
        !
      
        if (present(formatCSV)) then
            if (formatCSV) then
                isCSV = .true.
            else
                isCSV = .false.
            end if

        else
            isCSV = .false.
        end if
           
        call getunit ( iunit )
             
        open(UNIT=iunit, FILE = nombre_archivo, ACTION="write", STATUS="replace")


        nfilas = size(matriz,1)
    
        ncolumnas = size(matriz,2)

        do i=1,nfilas
            if (isCSV) then
                write(iunit,'(*(G0.4,:,","))') matriz(i,:)
            else
                write(iunit,'(' // inttostr(ncolumnas) // 'f16.8)') matriz(i,:)
            end if
        end do

        close(UNIT = iunit)
    end subroutine



    !*  ************************************************************
    !*                          guardarMatrizI
    !*  Guarda en un archivo (nombre_archivo) la matriz bidimensional
    !*  (matriz) de tipo entero y con formato 'cuadrado' en el archivo
    !*  ************************************************************
    subroutine guardarMatrizI(matriz, nombre_archivo, formatCSV)
        !dummy
        integer (kind = il), intent(in), dimension(:,:) :: matriz
        character (len = *),  intent(in) :: nombre_archivo
        logical, optional :: formatCSV
        !
        !local
        integer (kind = il) :: i
        integer (kind = il) :: iunit
        integer (kind = il) :: nfilas, ncolumnas
        logical :: isCSV
        !
        if (present(formatCSV)) then
            if (formatCSV) then
                isCSV = .true.
            else
                isCSV = .false.
            end if

        else
            isCSV = .false.
        end if

        call getunit ( iunit )
        open(UNIT=iunit, FILE = nombre_archivo, ACTION="write", STATUS="replace")

        nfilas = size(matriz,1)
        ncolumnas = size(matriz,2)

        do i=1,nfilas
            if (isCSV) then
                write(iunit,'(*(G0.4,:,","))') matriz(i,:)
            else
                write(iunit,'(' // inttostr(ncolumnas) // 'I6)') matriz(i,:)
            end if            
        end do

        close(UNIT = iunit)
    end subroutine



    !*  ************************************************************
    !*							guardarVectorC
    !*	Guarda en un archivo (nombre_archivo) el vector (vector) de
    !*  tipo complex, en formato dos columnas (real, imaginario)
    !*  ************************************************************
    subroutine guardarVectorC(vector, nombre_archivo)
        !dummy
        complex (kind = dp), intent(in), dimension(:) :: vector
        character (len = *),  intent(in) :: nombre_archivo
        !
        !local
        integer (kind = il) :: i
        integer (kind = il) :: iunit
        integer (kind = il) :: nfilas, ncolumnas
        !

        call getunit ( iunit )
        open(UNIT=iunit, FILE = nombre_archivo, ACTION="write", STATUS="replace")

        nfilas = size(vector,1)
        ncolumnas = 2

        do i=1,nfilas
            write(iunit,'(' // inttostr(ncolumnas) // 'f16.8)') vector(i)
        end do
        close(UNIT = iunit)
    end subroutine


    !*  ************************************************************
    !*							frecuenciaMaxMallado
    !*	Imprime l alongitud de onda minima de uso del dispersor segun
    !*	varios criterios
    !*  ************************************************************
    subroutine frecuenciaMaxMallado(t_area, e_long, num_t, num_e)
        !dummy
        real (kind = dp ),  dimension(:) :: e_long
        real ( kind = dp ), dimension(:) :: t_area
        integer ( kind = il ) :: num_e
        integer ( kind = il ) :: num_t
        !
        !local
        integer (kind = il) :: i
        real (kind = dp) :: lMax, lProm, lDesv, l
        real (kind = dp) :: aMax, aProm, aDesv, a
        real (kind = dp) :: aTotal
        real (kind = dp) :: factor
        !
        factor = 10.

        lMax = 0.
        lProm = 0.
        LDesv = 0.

        do i = 1, num_e
            l = e_long(i)
            if (l>lMax) then
                lMax = l
            end if
            lProm = lProm + l
        end do
        lProm = lProm/num_e
        do i = 1, num_e
            l = e_long(i)
            lDesv = lDesv + (l - lProm)**2
        end do
        lDesv = lDesv/num_e

print*, 'Lambda min q puede ser usado pa un mallado perfecto (criterio longitud RWG): ' // realtostr(factor*lMax)
print*, 'Lambda min q puede ser usado pa un mallado 95.44% correcto (criterio longitud RWG): ' // realtostr(factor*(lProm+2*lDesv))


        aMax = 0.
        aProm = 0.
        aDesv = 0.
        aTotal = 0.
        do i = 1, num_t
            a = t_area(i)
            if (a>aMax) then
                aMax = a
            end if
            aProm = aProm + a
        end do
        aTotal = aProm
        aProm = aProm/num_t
        do i = 1, num_t
            a = t_area(i)
            aDesv = aDesv + (a - aProm)**2
        end do
        aDesv = aDesv/num_t

print*, 'Lambda min pue c usado para un mallado perfecto (criterio area parches): ' // realtostr(factor*sqrt(aMax))
print*, 'Lambda min pue c usado para un mallado 95.44% correcto (criterio parches): ' // realtostr(factor*sqrt(aProm+2.*aDesv))



        print*, 'Lambda minimo que puede ser usado para un mallado correcto' // realtostr(factor*sqrt(aTotal/num_e))

    end subroutine


    !*  ************************************************************
    !*							limites_dispersor
    !*	Funcion que devuelve un vector con los 6 indices de los lados comunes
    !*	(edges) que tienen mayor y menor coordenada en X Y y Z
    !*  ************************************************************
    function limites_dispersor(e_centro, num_e) result (limites)
        ! devuelve un vector con los índices de los lados comunes responsables de los limites
        ! máximos en X Y y Z del dispersos: [maxX,minX,maxY,minY,maxZ,minZ]

        !dummy
        real ( kind = dp ), dimension(:,:), intent(in) :: e_centro
        integer (kind = il), dimension(6) :: limites
        integer ( kind = il ), intent(in) :: num_e
        !

        integer (kind = il) :: i

        real (kind = dp), dimension(3) :: centroi

        !inicializar
        limites(:) = 1

        do i = 2, num_e

            centroi = e_centro(:,i)

            !maxX
            if ( centroi(1) > e_centro(1,limites(1)) ) then
                limites(1) = i
            end if
            !minX
            if ( centroi(1) < e_centro(1,limites(2)) ) then
                limites(2) = i
            end if
            !maxY
            if ( centroi(2) > e_centro(2,limites(3)) ) then
                limites(3) = i
            end if
            !minY
            if ( centroi(2) < e_centro(2,limites(4)) ) then
                limites(4) = i
            end if
            !maxZ
            if ( centroi(3) > e_centro(3,limites(5)) ) then
                limites(5) = i
            end if
            !minZ
            if ( centroi(3) < e_centro(3,limites(6)) ) then
                limites(6) = i
            end if

        end do

    end function limites_dispersor

        !funciones matematicas


	function cross_productC(v1,v2) result(vout)

		complex (kind = dp), intent(in), dimension(3) :: v1,v2
		complex (kind = dp), dimension(3) :: vout

		vout(1) = v1(2)*v2(3) - v1(3)*v2(2)
		vout(2) = v1(3)*v2(1) - v1(1)*v2(3)
		vout(3) = v1(1)*v2(2) - v1(2)*v2(1)
	end function


    function cross_productR(v1,v2) result(vout)

        real (kind = dp), intent(in), dimension(3) :: v1,v2
        real (kind = dp), dimension(3) :: vout

        vout(1) = v1(2)*v2(3) - v1(3)*v2(2)
        vout(2) = v1(3)*v2(1) - v1(1)*v2(3)
        vout(3) = v1(1)*v2(2) - v1(2)*v2(1)
    end function


    function green (punto_m, punto_n) result (vg)
        !dummy
        real (kind = dp), dimension (3) :: punto_m, punto_n
        complex (kind = dp) :: vg
        !local
        real (kind = dp) :: r

        !r = modulov(punto_m-punto_n)
		r = sqrt(dot_product(punto_m-punto_n,punto_m-punto_n))

        !print*, punto_n


        vg = exp(-jj*w*(sqrt(mu*epsi))*r)/(r)

    end function green





    function hankel(n,x) result(resultado)
        !dumy
        real (kind = dp), intent(in) :: x
        integer (kind = il), intent(in) :: n
        complex (kind = dp) :: resultado

        !resultado = BESSEL_JN (n, x) - jj*BESSEL_YN (n, x)
        resultado = sphBesselJ (n, x) - jj*sphBesselY (n, x)

    end function hankel





    function sphBesselJ(n,x) result(valor)
        !dummy
        integer(kind = il), intent(in) :: n
        real (kind = dp), intent(in) :: x
        real (kind = dp) :: valor
        !
        !local
        real(kind = dp), dimension(n+1) :: SJ_,DJ_
        integer (kind = il) :: NM_
        !print*, n
        !print*, x
        !print*, 'pre llamo j'
        call sphj(n,x,NM_,SJ_,DJ_)
        !print*, 'llamo j'
        valor = SJ_(n+1)

    end function


    function sphBesselY(n,x) result(valor)
        !dummy
        integer(kind = il), intent(in) :: n
        real (kind = dp), intent(in) :: x
        real (kind = dp) :: valor
        !
        !local
        real(kind = dp), dimension(n+1) :: SY_,DY_
        integer (kind = il) :: NM_
        !print*, 'pre llamo y'
        call sphy(n,x,NM_,SY_,DY_)
        !print*, 'llamo y'
        valor = SY_(n+1)

    end function






    function envj ( n, x ) result(valueenvj)

        !*****************************************************************************80
        !
        !! ENVJ is a utility function used by MSTA1 and MSTA2.
        !
        !  Licensing:
        !
        !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However,
        !    they give permission to incorporate this routine into a user program
        !    provided that the copyright is acknowledged.
        !
        !  Modified:
        !
        !    14 March 2012
        !
        !  Author:
        !
        !    Shanjie Zhang, Jianming Jin
        !
        !  Reference:
        !
        !    Shanjie Zhang, Jianming Jin,
        !    Computation of Special Functions,
        !    Wiley, 1996,
        !    ISBN: 0-471-11963-6,
        !    LC: QA351.C45.
        !
        !  Parameters:
        !
        !    Input, integer ( kind = il ) N, ?
        !
        !    Input, real ( kind = dp ) X, ?
        !
        !    Output, real ( kind = dp ) ENVJ, ?
        !

        real ( kind = dp ) valueenvj
        integer ( kind = il ) n
        real ( kind = dp ) x

        valueenvj = 0.5D+00 * log10 ( 6.28D+00 * n ) - n * log10 ( 1.36D+00 * x / n )

      !return

    end function







    function msta1 ( x, mp ) result(valuemsta1)

        !*****************************************************************************80
        !
        !! MSTA1 determines a backward recurrence starting point for Jn(x).
        !
        !  Discussion:
        !
        !    This procedure determines the starting point for backward
        !    recurrence such that the magnitude of
        !    Jn(x) at that point is about 10^(-MP).
        !
        !  Licensing:
        !
        !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However,
        !    they give permission to incorporate this routine into a user program
        !    provided that the copyright is acknowledged.
        !
        !  Modified:
        !
        !    08 July 2012
        !
        !  Author:
        !
        !    Shanjie Zhang, Jianming Jin
        !
        !  Reference:
        !
        !    Shanjie Zhang, Jianming Jin,
        !    Computation of Special Functions,
        !    Wiley, 1996,
        !    ISBN: 0-471-11963-6,
        !    LC: QA351.C45.
        !
        !  Parameters:
        !
        !    Input, real ( kind = dp ) X, the argument.
        !
        !    Input, integer ( kind = il ) MP, the negative logarithm of the
        !    desired magnitude.
        !
        !    Output, integer ( kind = il ) MSTA1, the starting point.
        !

        real ( kind = dp ) a0
        !real ( kind = dp ) envj
        real ( kind = dp ) f
        real ( kind = dp ) f0
        real ( kind = dp ) f1
        integer ( kind = il ) it
        integer ( kind = il ) mp
        integer ( kind = il ) valuemsta1
        integer ( kind = il ) n0
        integer ( kind = il ) n1
        integer ( kind = il ) nn
        real ( kind = dp ) x

        a0 = abs ( x )
        n0 = int ( 1.1D+00 * a0 ) + 1
        f0 = envj ( n0, a0 ) - mp
        n1 = n0 + 5
        f1 = envj ( n1, a0 ) - mp
        do it = 1, 20
            nn = n1 - ( n1 - n0 ) / ( 1.0D+00 - f0 / f1 )
            f = envj ( nn, a0 ) - mp
            if ( abs ( nn - n1 ) < 1 ) then
                exit
            end if
            n0 = n1
            f0 = f1
            n1 = nn
            f1 = f
        end do

        valuemsta1 = nn

      !return
    end function
    function msta2 ( x, n, mp ) result(valuemsta2)

        !*****************************************************************************80
        !
        !! MSTA2 determines a backward recurrence starting point for Jn(x).
        !
        !  Discussion:
        !
        !    This procedure determines the starting point for a backward
        !    recurrence such that all Jn(x) has MP significant digits.
        !
        !  Licensing:
        !
        !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However,
        !    they give permission to incorporate this routine into a user program
        !    provided that the copyright is acknowledged.
        !
        !  Modified:
        !
        !    08 July 2012
        !
        !  Author:
        !
        !    Shanjie Zhang, Jianming Jin
        !
        !  Reference:
        !
        !    Shanjie Zhang, Jianming Jin,
        !    Computation of Special Functions,
        !    Wiley, 1996,
        !    ISBN: 0-471-11963-6,
        !    LC: QA351.C45.
        !
        !  Parameters:
        !
        !    Input, real ( kind = dp ) X, the argument of Jn(x).
        !
        !    Input, integer ( kind = il ) N, the order of Jn(x).
        !
        !    Input, integer ( kind = il ) MP, the number of significant digits.
        !
        !    Output, integer ( kind = il ) MSTA2, the starting point.
        !


        real ( kind = dp ) a0
        real ( kind = dp ) ejn
        !real ( kind = dp ) envj
        real ( kind = dp ) f
        real ( kind = dp ) f0
        real ( kind = dp ) f1
        real ( kind = dp ) hmp
        integer ( kind = il ) it
        integer ( kind = il ) mp
        integer ( kind = il ) valuemsta2
        integer ( kind = il ) n
        integer ( kind = il ) n0
        integer ( kind = il ) n1
        integer ( kind = il ) nn
        real ( kind = dp ) obj
        real ( kind = dp ) x

        a0 = abs ( x )
        hmp = 0.5D+00 * mp
        ejn = envj ( n, a0 )

        if ( ejn <= hmp ) then
            obj = mp
            n0 = int ( 1.1D+00 * a0 )
        else
            obj = hmp + ejn
            n0 = n
        end if

        f0 = envj ( n0, a0 ) - obj
        n1 = n0 + 5
        f1 = envj ( n1, a0 ) - obj

        do it = 1, 20
            nn = n1 - ( n1 - n0 ) / ( 1.0D+00 - f0 / f1 )
            f = envj ( nn, a0 ) - obj
            if ( abs ( nn - n1 ) < 1 ) then
                exit
            end if
            n0 = n1
            f0 = f1
            n1 = nn
            f1 = f
        end do

        valuemsta2 = nn + 10

      !return
    end function




    subroutine sphy ( n, x, nm, sy, dy )

        !*****************************************************************************80
        !
        !! SPHY computes spherical Bessel functions yn(x) and their derivatives.
        !
        !  Licensing:
        !
        !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However,
        !    they give permission to incorporate this routine into a user program
        !    provided that the copyright is acknowledged.
        !
        !  Modified:
        !
        !    15 July 2012
        !
        !  Author:
        !
        !    Shanjie Zhang, Jianming Jin
        !
        !  Reference:
        !
        !    Shanjie Zhang, Jianming Jin,
        !    Computation of Special Functions,
        !    Wiley, 1996,
        !    ISBN: 0-471-11963-6,
        !    LC: QA351.C45.
        !
        !  Parameters:
        !
        !    Input, integer ( kind = il ) N, the order.
        !
        !    Input, real ( kind = dp ) X, the argument.
        !
        !    Output, integer ( kind = il ) NM, the highest order computed.
        !
        !    Output, real ( kind = dp ) SY(0:N), DY(0:N), the values of yn(x) and yn'(x).
        !
        implicit none

        integer ( kind = il ) n

        real ( kind = dp ) dy(0:n)
        real ( kind = dp ) f
        real ( kind = dp ) f0
        real ( kind = dp ) f1
        integer ( kind = il ) k
        integer ( kind = il ) nm
        real ( kind = dp ) sy(0:n)
        real ( kind = dp ) x

        nm = n

        if ( x < 1.0D-60 ) then
            do k = 0, n
                sy(k) = -1.0D+300
                dy(k) = 1.0D+300
            end do
            return
        end if

        sy(0) = - cos ( x ) / x
        sy(1) = ( sy(0) - sin ( x ) ) / x
        f0 = sy(0)
        f1 = sy(1)
        do k = 2, n
            f = ( 2.0D+00 * k - 1.0D+00 ) * f1 / x - f0
            sy(k) = f
            if ( 1.0D+300 <= abs ( f ) ) then
                exit
            end if
            f0 = f1
            f1 = f
        end do

        nm = k - 1
        dy(0) = ( sin ( x ) + cos ( x ) / x ) / x
        do k = 1, nm
            dy(k) = sy(k-1) - ( k + 1.0D+00 ) * sy(k) / x
        end do

      !return
    end subroutine





    subroutine sphj ( n, x, nm, sj, dj )

        !*****************************************************************************80
        !
        !! SPHJ computes spherical Bessel functions jn(x) and their derivatives.
        !
        !  Licensing:
        !
        !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However,
        !    they give permission to incorporate this routine into a user program
        !    provided that the copyright is acknowledged.
        !
        !  Modified:
        !
        !    15 July 2012
        !
        !  Author:
        !
        !    Shanjie Zhang, Jianming Jin
        !
        !  Reference:
        !
        !    Shanjie Zhang, Jianming Jin,
        !    Computation of Special Functions,
        !    Wiley, 1996,
        !    ISBN: 0-471-11963-6,
        !    LC: QA351.C45.
        !
        !  Parameters:
        !
        !    Input, integer ( kind = il ) N, the order.
        !
        !    Input, real ( kind = dp ) X, the argument.
        !
        !    Output, integer ( kind = il ) NM, the highest order computed.
        !
        !    Output, real ( kind = dp ) SJ(0:N), the values of jn(x).
        !
        !    Output, real ( kind = dp ) DJ(0:N), the values of jn'(x).
        !
        implicit none

        integer ( kind = il ) n

        real ( kind = dp ) cs
        real ( kind = dp ) dj(0:n)
        real ( kind = dp ) f
        real ( kind = dp ) f0
        real ( kind = dp ) f1
        integer ( kind = il ) k
        integer ( kind = il ) m
        !integer ( kind = il ) msta1
        !integer ( kind = il ) msta2
        integer ( kind = il ) nm
        real ( kind = dp ) sa
        real ( kind = dp ) sb
        real ( kind = dp ) sj(0:n)
        real ( kind = dp ) x

        !allocate(sj(0:n))

        nm = n

        if ( abs ( x ) <= 1.0D-100 ) then
            do k = 0, n
                sj(k) = 0.0D+00
                dj(k) = 0.0D+00
            end do
            sj(0) = 1.0D+00
            dj(1) = 0.3333333333333333D+00
            return
        end if

        sj(0) = sin ( x ) / x
        sj(1) = ( sj(0) - cos ( x ) ) / x

        if ( 2 <= n ) then

            sa = sj(0)
            sb = sj(1)
            m = msta1 ( x, 200 )
            if ( m < n ) then
                nm = m
            else
                m = msta2 ( x, n, 15 )
            end if

            f0 = 0.0D+00
            f1 = 1.0D+00-100
            do k = m, 0, -1
                f = ( 2.0D+00 * k + 3.0D+00 ) * f1 / x - f0
                if ( k <= nm ) then
                    sj(k) = f
                end if
                f0 = f1
                f1 = f
            end do

            if ( abs ( sa ) <= abs ( sb ) ) then
                cs = sb / f0
            else
                cs = sa / f
            end if

            do k = 0, nm
                sj(k) = cs * sj(k)
            end do

        end if

        dj(0) = ( cos(x) - sin(x) / x ) / x
        do k = 1, nm
            dj(k) = sj(k-1) - ( k + 1.0D+00 ) * sj(k) / x
        end do

      !return
    end subroutine


    SUBROUTINE  gauleg(ngp, xabsc, weig)
              implicit none

       REAL (kind = 8) :: newv
       REAL(kind = 8)  :: EPS, M_PI
       PARAMETER (EPS=3.0d-15)          !EPS is the relative precision
       PARAMETER (M_PI=3.141592654d0)      ! Pi value
          INTEGER  i, j, m
          REAL(kind = 8)  p1, p2, p3, pp, z, z1
          INTEGER, INTENT(IN) :: ngp            ! # of Gauss Points
          REAL(kind = 8), INTENT(OUT) :: xabsc(ngp), weig(ngp)


           m = (ngp + 1) / 2
    !* Roots are symmetric in the interval - so only need to find half of them  */

           do i = 1, m              ! Loop over the desired roots */

                z = cos( M_PI * (i-0.25d0) / (ngp+0.5d0) )
    !*   Starting with the above approximation to the ith root,
    !*          we enter the main loop of refinement by NEWTON'S method   */
    100         p1 = 1.0d0
                p2 = 0.0d0
    !*  Loop up the recurrence relation to get the Legendre
    !*  polynomial evaluated at z                 */

                do j = 1, ngp
                p3 = p2
                p2 = p1
                p1 = ((2.0d0*j-1.0d0) * z * p2 - (j-1.0d0)*p3) / j
                enddo

    !* p1 is now the desired Legendre polynomial. We next compute pp,
    !* its derivative, by a standard relation involving also p2, the
    !* polynomial of one lower order.      */
                pp = ngp*(z*p1-p2)/(z*z-1.0d0)
                z1 = z
                z = z1 - p1/pp             ! Newton's Method  */

                if (dabs(z-z1) .gt. EPS) GOTO  100

            xabsc(i) =  - z                     ! Roots will be bewteen -1.0 & 1.0 */
            xabsc(ngp+1-i) =  + z                   ! and symmetric about the origin  */
            weig(i) = 2.0d0/((1.0d0-z*z)*pp*pp) ! Compute the weight and its       */
            weig(ngp+1-i) = weig(i)               ! symmetric counterpart         */

          end do     ! i loop
          xabsc = (-1.)*xabsc
    End subroutine gauleg





    function legendre(n,x) result(resultado)
        !dumy
        integer (kind = il), intent(in) :: n
        real (kind = dp), intent(in) :: x
        real (kind = dp) :: resultado
        !local
        integer (kind = il) :: i

        resultado = 0.
        do i = 0,n
            resultado = resultado + ((factorial(n)/(factorial(i)*factorial(n-i)))**2)*((x+1)**(n-i))*((x-1)**i)
        end do
        resultado = (0.5**n)*resultado

    end function legendre

    function factorial(x) result(resultado)
        !dumy
        integer (kind = il), intent(in) :: x
        real (kind = dp) :: resultado
        !local
        integer (kind = il) :: i

        resultado = 1.
        do i=1,x
            resultado = resultado*i
        end do

    end function factorial



    !funciones basicas


        function gettiempo() result(res)
            !dummy
            character (len = 15) :: res
            !
            !local
            integer ( kind= il ) values(8)
            integer ( kind= il ) h
            integer ( kind= il ) mm
            integer ( kind= il ) n
            integer ( kind= il ) s            
            !
            call date_and_time ( values = values )
            h = values(5)
            n = values(6)
            s = values(7)
            mm = values(8)

            res = inttostr(h) // ':' // inttostr(n) // ':' // inttostr(s) // ':' // inttostr(mm)
        end function

        subroutine msjtiempo ( )

        character ( len = 8 ) ampm
        integer ( kind= il ) d
        integer ( kind= il ) h
        integer ( kind= il ) m
        integer ( kind= il ) mm
        character ( len = 9 ), parameter, dimension(12) :: month = (/ &
            'January  ', 'February ', 'March    ', 'April    ', &
            'May      ', 'June     ', 'July     ', 'August   ', &
            'September', 'October  ', 'November ', 'December ' /)
        integer ( kind= il ) n
        integer ( kind= il ) s
        integer ( kind= il ) values(8)
        integer ( kind= il ) y

        call date_and_time ( values = values )

        y = values(1)
        m = values(2)
        d = values(3)
        h = values(5)
        n = values(6)
        s = values(7)
        mm = values(8)

        if ( h < 12 ) then
            ampm = 'AM'
        else if ( h == 12 ) then
            if ( n == 0 .and. s == 0 ) then
                ampm = 'Noon'
            else
                ampm = 'PM'
            end if
        else
            h = h - 12
            if ( h < 12 ) then
                ampm = 'PM'
            else if ( h == 12 ) then
                if ( n == 0 .and. s == 0 ) then
                    ampm = 'Midnight'
                else
                    ampm = 'AM'
                end if
            end if
        end if

        write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
            d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

        return
    end subroutine msjtiempo



    function realtostrE(num) result(text)
        !dummy
        character (len = 10) :: text
        real (kind = dp), intent(in) :: num
        !

        write (text,'(ES10.2)') num
    end function realtostrE

    function realtostrE2(num,n_dig) result(text)
        integer (kind = ilo), intent(in) :: n_dig
        character (len = 10 + n_dig) :: text
        
        real (kind = dp), intent(in) :: num
        !

        write (text,'(ES' // trim(inttostr(8+n_dig)) // '.' //  trim(inttostr(n_dig))  //  ')') num
    end function realtostrE2


    function dbllongtostr(num) result(text)

        !dummy
        integer(kind = ild), intent(in) :: num
        character (len = 16) text
        !

        Write(text, '(I14)') num

        text = trim(adjustl(text)) !// char(0)
    end function dbllongtostr

    function inttostr(num) result(text)

        !dummy
        integer(kind = il), intent(in) :: num
        character (len = Digitos(num)) text
        !
        !local
        character (len = 1) :: formatoDD
        character (len = Digitos(Digitos(num))):: formato
        Write(formatoDD, '(I1)') Digitos(Digitos(num))
        Write(formato, '(I' // formatoDD // ')') Digitos(num)
        Write(text, '(I' // formato // ')') num

        text = trim(adjustl(text)) !// char(0)
    end function inttostr
    function realtostr(num) result(text)
        !dummy
        real(kind = dp), intent(in) :: num
        character (len = (6+DigitosR(num))) text

        character (len = 1) :: formatoDD
        character (len = Digitos(6+DigitosR(num))):: formato


        Write(formatoDD, '(I1)') Digitos(6+DigitosR(num))
        Write(formato, '(I' // formatoDD // ')') (6+DigitosR(num))
        Write(text, '(f' // formato // '.4)') num

        !Write(text, '(f10.4)') num
        text = trim(adjustl(text))
    end function realtostr
    pure function Digitos(numero) result(d)
        !dummy
        integer(kind = il), intent(in) :: numero
        integer(kind = il) :: d
        !
        !local
        integer(kind = il) :: rest, cont
        rest = numero
        cont = 0
        do
            if (rest/10 == 0) then
                exit
            else
                cont = cont + 1
            end if
            rest = rest/10
        end do
        if (numero < 0 ) then
            cont = cont + 1
        end if
        d = cont + 1

    end function Digitos
    pure function DigitosR(numer) result(d)
        !dummy
        real(kind = dp), intent(in) :: numer
        integer(kind = il) :: d
        !
        !local
        integer(kind = il) :: rest, cont,numero
        numero = numer
        rest = numero
        cont = 0
        do
            if (rest/10 == 0) then
                exit
            else
                cont = cont + 1
            end if
            rest = rest/10
        end do
        if (numero < 0 ) then
            cont = cont + 1
        end if
        d = cont + 1 !+ 5
    end function DigitosR

    subroutine getunit ( iunit )
        integer ( kind = il ) i
        integer ( kind = il ) ios
        integer ( kind = il ) iunit
        logical lopen
        iunit = 0
        do i = 1, 99
            if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then
                inquire ( unit = i, opened = lopen, iostat = ios )
                if ( ios == 0 ) then
                    if ( .not. lopen ) then
                        iunit = i
                        return
                    end if
                end if
            end if
        end do
        return
    end subroutine getunit

    ! --------------------------------------------------------------------
    ! SUBROUTINE  Sort():
    !    This subroutine receives an array x() and sorts it into ascending
    ! order.
    ! --------------------------------------------------------------------

   SUBROUTINE  Sort(x, Size)
      IMPLICIT  NONE
      INTEGER, DIMENSION(1:), INTENT(INOUT) :: x
      INTEGER, INTENT(IN)                   :: Size
      INTEGER                               :: i
      INTEGER                               :: Location

      DO i = 1, Size-1          ! except for the last
         Location = FindMinimum(x, i, Size) ! find min from this to last
         CALL  Swap(x(i), x(Location))  ! swap this and the minimum
      END DO
   END SUBROUTINE  Sort

    ! --------------------------------------------------------------------
    ! INTEGER FUNCTION  FindMinimum():
    !    This function returns the location of the minimum in the section
    ! between Start and End.
    ! --------------------------------------------------------------------

   INTEGER FUNCTION  FindMinimum(x, Start, End)
      IMPLICIT  NONE
      INTEGER, DIMENSION(1:), INTENT(IN) :: x
      INTEGER, INTENT(IN)                :: Start, End
      INTEGER                            :: Minimum
      INTEGER                            :: Location
      INTEGER                            :: i

      Minimum  = x(Start)       ! assume the first is the min
      Location = Start          ! record its position
      DO i = Start+1, End       ! start with next elements
         IF (x(i) < Minimum) THEN   !   if x(i) less than the min?
            Minimum  = x(i)     !      Yes, a new minimum found
            Location = i                !      record its position
         END IF
      END DO
      FindMinimum = Location            ! return the position
   END FUNCTION  FindMinimum

    ! --------------------------------------------------------------------
    ! SUBROUTINE  Swap():
    !    This subroutine swaps the values of its two formal arguments.
    ! --------------------------------------------------------------------

   SUBROUTINE  Swap(a, b)
      IMPLICIT  NONE
      INTEGER, INTENT(INOUT) :: a, b
      INTEGER                :: Temp

      Temp = a
      a    = b
      b    = Temp
   END SUBROUTINE  Swap


end module modGlobalMetod
