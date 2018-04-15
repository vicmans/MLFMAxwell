program main
    use modGlobalParam
    use modGlobalMetod
    use modReadStl
    use modRWG
    use modMoM
    use modRCS
    use modIterativo
    use modFMM
    use modMLFMA
    implicit none


    complex (kind = dp_prec), allocatable, dimension(:) :: precond_alu
    integer (kind = il), allocatable, dimension(:) :: precond_jlu
    integer (kind = ild), allocatable, dimension(:) :: precond_ju

    complex (kind = dp), allocatable, dimension(:,:) :: Z!, Zinv
    complex ( kind = dp ),  pointer, dimension (:,:) :: Z_MatxVec

    complex (kind = dp), allocatable, dimension(:) :: I
    complex (kind = dp), allocatable, dimension(:) :: b
    complex (kind = dp), allocatable, dimension(:) :: guess
    real ( kind = dp ), allocatable, dimension(:,:) :: p_coord
    integer ( kind = il ), allocatable, dimension(:,:) :: t_p
    integer ( kind = il ), allocatable, dimension(:,:) :: e_p
    integer ( kind = il ), allocatable, dimension(:,:) :: e_t
    integer ( kind = il ), allocatable, dimension(:,:) :: e_po
    real (kind = dp ), allocatable, dimension(:) :: e_long
    real ( kind = dp ), allocatable, dimension(:) :: t_area
    real ( kind= dp ), allocatable, dimension(:,:) :: t_normal
    real ( kind = dp ), allocatable, dimension(:,:) :: e_centro
    real ( kind = dp ), allocatable, dimension(:,:) :: t_baric
    real ( kind = dp ), allocatable, dimension(:,:,:) :: t_baric_sub

    real (kind = dp) :: frequency

    integer ( kind = il ) :: num_e
    integer ( kind = il ) :: num_p
    integer ( kind = il ) :: num_t


    complex (kind = dp), dimension(3) :: pol_onda
    complex (kind = dp) :: ctte_onda

    real (kind = dp) :: error_solver
    real (kind = dp), dimension(3) :: dir_onda
    
    character (len = 256) :: comm, val
    character (len = 3), parameter :: strdefault= '/%\'
    character (len = 256) :: autocommand

    integer (kind = il) :: ent, numtest

    systemDelimiter = '/'
    !crea la carpeta para asegurarnos de que exista
    !call system('mkdir -p --parents' // '..' // trim(systemDelimiter) // 'Results')

    !Lectura de argumentos de la linea de comandos
    conf%argument1(1:len(conf%argument1)) = ' '
    conf%argument2(1:len(conf%argument2)) = ' '
    call getarg(1, conf%argument1)
    call getarg(2, conf%argument2)
    if (streq(conf%argument2, 'gui')) then
        conf%optionGUI = .true.
    else
        conf%optionGUI = .false.
    end if
    !

    !Inicializacion

    ctte_onda = 1.
    lambda = 1.
    call setlambda(lambda)
    !skipkD = .false.

    conf%efieOldMode = 0
    conf%test_name = 'lasttest'!strdefault
    conf%msj = .false.
    conf%numThetaRCS = 1
    conf%numPhiRCS = 360
    conf%errorSolverCGS = 0.001
    conf%Lforz = -1
    conf%dirTheta = 90.
    conf%dirPhi = 180.
    conf%rcsMonoFMin = 10
    conf%rcsMonoFMax = 275
    conf%rcsMonoSamples = 75
    conf%usar_precond = 'ILUT'
    conf%prec_activo = .true.
    conf%fill_in_prc = 5.
    conf%drop_tolerance = 0.001
    conf%analisis_rcs = 'MONO'
    conf%n_pasos_rcs_mono = 180
    conf%tipo_polarizac = 'VV'
    conf%variable_espec = 'FREC'
    conf%file_name = strdefault
    conf%file_script = strdefault
    conf%criterio_ladomax = 0.5_dp
    conf%CFIEfactor = 1._dp
    conf%force_ladomax = .false.
    conf%mult_maxprecision = 6
    conf%mult_minprecision = 0
    conf%isurf_phase = .false.
    conf%save_incrementals = .true.
    conf%rcsmonostepi = 1
    conf%rcsmonostepf = conf%n_pasos_rcs_mono

    alphaCFIE = conf%CFIEfactor      

    !conf%polVH = 1
    error_solver = conf%errorSolverCGS
    dir_onda(1) = sin(conf%dirTheta*pi/180.)*cos(conf%dirPhi*pi/180.)
    dir_onda(2) = sin(conf%dirTheta*pi/180.)*sin(conf%dirPhi*pi/180.)
    dir_onda(3) = cos(conf%dirTheta*pi/180.)
    call ajustar_polariz(dir_onda,pol_onda)


    if (streq(conf%argument1, '')) then
        call waitcommands()
    else
        autocommand = 'loadscript'
        call waitcommands(.true.,autocommand ,(adjustl(conf%argument1)))
    end if

contains


subroutine destruir_todo()

    call destruir_MLFMA()
    call destruir_FMM()

    if (allocated(p_coord)) then
        deallocate(p_coord)
    end if

    if (allocated(guess)) then
        deallocate(guess)
    end if

    if (allocated(b)) then
        deallocate(b)
    end if

    if (allocated(I)) then
        deallocate(I)
    end if

    if (allocated(Z)) then
        deallocate(Z)
    end if

    if (allocated(precond_ju)) then
        deallocate(precond_ju)
    end if

    if (allocated(precond_jlu)) then
        deallocate(precond_jlu)
    end if

    if (allocated(precond_alu)) then
        deallocate(precond_alu)
    end if

    if (allocated(t_p)) then
        deallocate(t_p)
    end if

    if (allocated(e_p)) then
        deallocate(e_p)
    end if

    if (allocated(e_t)) then
        deallocate(e_t)
    end if

    if (allocated(e_po)) then
        deallocate(e_po)
    end if

    if (allocated(e_long)) then
        deallocate(e_long)
    end if

    if (allocated(t_area)) then
        deallocate(t_area)
    end if

    if (allocated(t_normal)) then
        deallocate(t_normal)
    end if

    if (allocated(e_centro)) then
        deallocate(e_centro)
    end if

    if (allocated(t_baric)) then
        deallocate(t_baric)
    end if

    if (allocated(t_baric_sub)) then
        deallocate(t_baric_sub)
    end if


end subroutine

subroutine loadscript()

    !local
    integer (kind = il) :: iunit, i
    character (len = 256) :: console
    character (len = 256) :: argcomm2, argval2
    !
    call getunit ( iunit )
    open(UNIT=iunit, FILE = trim(conf%file_script))
    call printconsole ('Corriendo script...', 'running_script')
    do i = 1, 9999999
        read(iunit,'(A256)', end=10) console
        call printconsole(inttostr(i) // '  >> ' // trim(console) // ' - - -' // gettiempo(), 'command/' // trim(console))
        call getcommand(console, argcomm2, argval2)
        call waitcommands(.true.,argcomm2,argval2)
    end do
    
    


10 call printconsole('Script finalizado','end_script')
close(iunit)
call savelog(trim(conf%test_name)// '_LOG.dat') 
end subroutine



subroutine waitcommands(justone, argcomm, argval)
    !dummy
    logical, optional, intent(in) :: justone
    character (len = 256), optional, intent(in) :: argcomm, argval
    !
    !local
    logical :: flag_script

    if (present(justone)) then
        flag_script = justone
    else
        flag_script = .false.
    end if

    if (present(argcomm)) then
    else
        flag_script = .false.
    end if
    if (present(argval)) then
    else
        flag_script = .false.
    end if

    if (flag_script .eqv..false.) then
        call printconsole( '* * * * * * *MLFMAxwell* * * * * * *', 'welcome')
        call printconsole( ' ', 'welcome')
        call printconsole( 'ejecute help para ayuda', 'welcome')
    end if


    do
        if (flag_script) then
            comm = argcomm
            val = argval

        else
        call getcommand_console(comm, val)

        end if

        if (streq(comm,'run')) then
            if (streq(val,'mlfma')) then
                call testMLFMA()
            elseif (streq(val,'fmm')) then
                call testfmm()
            elseif (streq(val,'mom')) then
                call testMoM()
            elseif (streq(val,'wizard')) then
                call wizard
            else
                print*, 'Valor invalido para el comando ' // trim(comm)
            end if
        elseif  (streq(comm,'help')) then
            print*, 'Debe introducir comandos con la siguiente sintaxis:'
            print*, 'comando=valor'
            print*, 'Ejecute el comando [config] para mostrar las configuraciones y las posibles simulaciones a ejecutar'
        elseif  (streq(comm,'loadscript') .or. streq(comm,'load')) then
            conf%file_script = trim(adjustl(val))
            call loadscript()
        elseif  (streq(comm,'config')) then
            call showconfig()
        elseif  (streq(comm,'filename')) then
            conf%file_name = trim(val)
        elseif  (streq(comm,'testname')) then
            conf%test_name = trim(val)
        elseif  (streq(comm,'ladomax')) then
            read(val,*) conf%criterio_ladomax   
        elseif  (streq(comm,'errsolver')) then
            read(val,*) conf%errorSolverCGS
            error_solver = conf%errorSolverCGS
        elseif  (streq(comm,'rcsbi.ntheta')) then
            read(val,*) conf%numThetaRCS
        elseif  (streq(comm,'rcsbi.nphi')) then
            read(val,*) conf%numPhiRCS
        elseif  (streq(comm,'efieOldMode')) then
            read(val,*) conf%efieOldMode            
        elseif  (streq(comm,'force_multipoles')) then
            read(val,*) conf%Lforz
        elseif  (streq(comm,'rcsbi.dirtheta')) then
            read(val,*) conf%dirTheta
            dir_onda(1) = sin(conf%dirTheta*pi/180.)*cos(conf%dirPhi*pi/180.)
            dir_onda(2) = sin(conf%dirTheta*pi/180.)*sin(conf%dirPhi*pi/180.)
            dir_onda(3) = cos(conf%dirTheta*pi/180.)
            call ajustar_polariz(dir_onda,pol_onda)
        elseif  (streq(comm,'rcsbi.dirphi')) then
            read(val,*) conf%dirPhi
            dir_onda(1) = sin(conf%dirTheta*pi/180.)*cos(conf%dirPhi*pi/180.)
            dir_onda(2) = sin(conf%dirTheta*pi/180.)*sin(conf%dirPhi*pi/180.)
            dir_onda(3) = cos(conf%dirTheta*pi/180.)
            call ajustar_polariz(dir_onda,pol_onda)
        elseif  (streq(comm,'rcsmonof.fmin')) then
            read(val,*) conf%rcsMonoFMin
        elseif  (streq(comm,'rcsmonof.fmax')) then
            read(val,*) conf%rcsMonoFMax
        elseif  (streq(comm,'rcsmonof.samples')) then
            read(val,*) conf%rcsMonoSamples
        elseif  (streq(comm,'precond')) then
            if (streq(val,'no')) then
                conf%usar_precond = 'NO'
            elseif (streq(val,'ilu0')) then
                conf%usar_precond = 'ILU0'
            elseif (streq(val,'ilut')) then
                conf%usar_precond = 'ILUT'
            else
                print*, 'Valor invalido para el comando ' // trim(comm)
            end if
        elseif (streq(comm,'coefcfie')) then
            read(val,*) conf%CFIEfactor
            alphaCFIE = conf%CFIEfactor           
        elseif (streq(comm,'ilut.fillin')) then
            read(val,*) conf%fill_in_prc
        elseif (streq(comm,'ilut.droptol')) then
            read(val,*) conf%drop_tolerance
        elseif (streq(comm,'ilut.droptol')) then
            read(val,*) conf%drop_tolerance
        elseif (streq(comm,'rcs.type')) then
            if (streq(val,'mono')) then
                conf%analisis_rcs = 'MONO'
            elseif (streq(val,'bi')) then
                conf%analisis_rcs = 'BI'
            elseif (streq(val,'monof')) then
                conf%analisis_rcs = 'MONOF'            
            else
                print*, 'Valor invalido para el comando ' // trim(comm)
            end if
        elseif (streq(comm,'rcsmono.samples')) then
            read(val,*) conf%n_pasos_rcs_mono
        elseif (streq(comm,'rcsmono.stepi')) then
            read(val,*) conf%rcsmonostepi
        elseif (streq(comm,'rcsmono.stepf')) then
            read(val,*) conf%rcsmonostepf
        elseif (streq(comm,'polariz')) then
            if (streq(val,'VV')) then
                conf%tipo_polarizac = 'VV'
            elseif (streq(val,'HH')) then
                conf%tipo_polarizac = 'HH'
            elseif (streq(val,'VH') .or. streq(val,'HV')) then
                conf%tipo_polarizac = 'VH'
            else
                print*, 'Valor invalido para el comando ' // trim(comm)
            end if
        elseif (streq(comm,'lambda')) then
            conf%variable_espec = 'LAMBDA'
            read(val,*) lambda
            call setlambda(lambda)
        elseif (streq(comm,'frec')) then
            conf%variable_espec = 'FREC'
            read(val,*) frequency
            call setf(frequency*1000000000_dp)
        elseif (streq(comm,'force_ladomax')) then
            if (streq(val,'true') .or. streq(val,'t') .or. streq(val,'yes') .or. streq(val,'y')) then
                conf%force_ladomax = .true.
            elseif (streq(val,'false') .or. streq(val,'f') .or. streq(val,'no') .or. streq(val,'n')) then
                conf%force_ladomax = .false.
            else
                print*, 'Valor invalido para el comando ' // trim(comm)
            end if
        elseif (streq(comm,'isurf_phase')) then
            if (streq(val,'true') .or. streq(val,'t') .or. streq(val,'yes') .or. streq(val,'y')) then
                conf%isurf_phase = .true.
            elseif (streq(val,'false') .or. streq(val,'f') .or. streq(val,'no') .or. streq(val,'n')) then
                conf%isurf_phase = .false.
            else
                print*, 'Valor invalido para el comando ' // trim(comm)
            end if
        elseif (streq(comm,'save_incrementals')) then
            if (streq(val,'true') .or. streq(val,'t') .or. streq(val,'yes') .or. streq(val,'y')) then
                conf%save_incrementals = .true.
            elseif (streq(val,'false') .or. streq(val,'f') .or. streq(val,'no') .or. streq(val,'n')) then
                conf%save_incrementals = .false.
            else
                print*, 'Valor invalido para el comando ' // trim(comm)
            end if
        elseif (streq(comm,'mult_maxprecision')) then
            read(val,*) conf%mult_maxprecision
        elseif (streq(comm,'mult_minprecision')) then
            read(val,*) conf%mult_minprecision
        else
            print*, 'Comando no valido'
        end if

        if (flag_script) then
            exit
        end if

    end do


end subroutine


subroutine showconfig()

character (len = 30) :: charerrsolver, chardroptol, charlambda, charfrec, charcoefcfie

write (charerrsolver,'(ES10.2)') conf%errorSolverCGS
write (chardroptol,'(ES10.2)') conf%drop_tolerance
write (charlambda,'(ES10.2)') getlambda()
write (charfrec,'(ES10.2)') getf()
write (charcoefcfie,'(ES10.2)') conf%CFIEfactor


print*, 'Configuracion actual del simulador...'
print*, ' '
        print*, '> Nombre del test (guarda results automaticamente)  [testname=' , trim(conf%test_name) , ']'
        print*, '> Archivo a cargar                                  [filename=' , trim(conf%file_name) , ']'
        print*, '> # de muestras en THETA de RCS Biestatica          [rcsbi.ntheta=' // inttostr(conf%numThetaRCS) //']'
        print*, '> # de muestras en PHI de RCS Biestatica            [rcsbi.nphi=' // inttostr(conf%numPhiRCS) //']'
        !print*, '> Error del solver CGS               [errsolver=' , conf%errorSolverCGS , ']'
        print*, '> Error del solver CGS                              [errsolver=' , trim(charerrsolver) , ']'
        print*, '> Lado maximo de los cubos mas pequeños             [ladomax=' // realtostrE(conf%criterio_ladomax) // ']'
        print*, '> Numero forzado de multipolos (-1 para automatico) [force_multipoles=' // inttostr(conf%Lforz) // ']'
        print*, '> Ang(°) THETA de incidenc de onda (rcs biestatica) [rcsbi.dirtheta=' // realtostr(conf%dirTheta) //']'
        print*, '> Ang(°) PHI de incidenc de onda (rcs biestatica)   [rcsbi.dirphi=' // realtostr(conf%dirPhi) //']'
        print*, '> Frecuencia inicial barrido RCS mono frec (MHz)    [rcsmonof.fmin=' // realtostr(conf%rcsMonoFMin) //']'
        print*, '> Frecuencia final barrido RCS mono frec (MHz)      [rcsmonof.fmax=' // realtostr(conf%rcsMonoFMax) //']'
        print*, '> # de muestras para RCS mono frec                  [rcsmonof.samples=' // inttostr(conf%rcsMonoSamples) //']'
        print*, '> Tipo precondicionador (NO/ILU0/ILUT)              [precond=' // trim(conf%usar_precond) //']'
        print*, '> Fill in de ILUT (%)                               [ilut.fillin=' // realtostr(conf%fill_in_prc) //']'
        print*, '> Drop tolerance de ILUT                            [ilut.droptol=' , trim(chardroptol) , ']'
        print*, '> Tipo de RCS a calcular                            [rcs.type=' // trim(conf%analisis_rcs) // ']'
        print*, '> # de muestras RCS monostatica (0°<=phi<180°)      [rcsmono.samples=' // inttostr(conf%n_pasos_rcs_mono) //']'
        print*, '> Tipo de polarizacion (VV/HH/VH(solo rcsmono))     [polariz=' , trim(conf%tipo_polarizac) , ']'
        print*, '> Longitud de onda de trabajo (m)                   [lambda=' , trim(charlambda) , ']'
        print*, '> Frecuencia de trabajo (GHz)                       [frec=' , trim(charfrec) , ']'
        print*, '> Factor CFIE (1 = EFIE)                            [coefCFIE=' , trim(charcoefcfie) , ']'
print*, '- - - - - - - - - - - - - - - - - - - -'
print*, 'Comandos de ejecucion de simulaciones:'
print*, 'run=mlfma (Resolver con MLFMA)'
print*, 'run=fmm (Resolver con FMM)'
print*, 'run=mom (Resolver con MoM)'
print*, 'run=wizard (Ejecuta el asistente de Configuracion)'
!print*, '[inhabilitado] run=rcsmonomom (Barrido de frec para rcs monostatica con MoM)'
print*, 'loadscript=nombre_archivo (ejecuta un script de comandos del simulador'
print*, ' '
print*, ' '

end subroutine

subroutine wizard()


    numtest = 0

    print*, '#Test a correr:'
    print*, '[0] Default'
    print*, '[1] testComparaImatxvecMLFMA (Comparacion I MoM vs MLFMA puro)'
    print*, '[2] testComparaImatxvec (comparacion I MoM vs FMM puro)'
    print*, '[3] testMLFMA (Realizar solo MLFMA puro)'
    print*, '[4] testFMM (Realizar solo FMM puro)'
    print*, '[5] testMoM (Realizar solo MoM puro)'
    print*, '[6] testcomparaZMLFMA (Comparacion Zfar MLFMA vs Zfar MoM)'
    print*, '[7] testComparaZfar (Comparacion Zfar FMM vs Zfar MoM)'
    print*, '[8] testRCSMoM_monof (RCS monostatica barrido frecuencia)'
    read*, numtest

    do
        print*, ' '
        print*, ' '
        print*, '###Seleccione correr [1], o configure la simulacion'
        print*, '[1] OK, correr'
        print*, '   [2] # de muestras en THETA de RCS Biestatica [default = ' // inttostr(conf%numThetaRCS) //']'
        print*, '   [3] # de muestras en PHI de RCS Biestatica [default = ' // inttostr(conf%numPhiRCS) //']'
        print*, '   [4] Error del solver CGS [default = ' , conf%errorSolverCGS , ']'
        print*, '   [5] Numero forzado de multipolos (-1 para automatico) [default = ' // inttostr(conf%Lforz) // ']'
        print*, '   [6] Ang(°) THETA de incidenc de onda (rcs biestatica) [default = ' // realtostr(conf%dirTheta) //']'
        print*, '   [7] Ang(°) PHI de incidenc de onda (rcs biestatica) [default = ' // realtostr(conf%dirPhi) //']'
        print*, '   [8] Frecuencia inicial barrido RCS mono frec (MHz) [default = ' // realtostr(conf%rcsMonoFMin) //']'
        print*, '   [9] Frecuencia final barrido RCS mono frec (MHz) [default = ' // realtostr(conf%rcsMonoFMax) //']'
        print*, '   [10] # de muestras para RCS mono frec [default = ' // inttostr(conf%rcsMonoSamples) //']'
        print*, '   [11] Tipo precondicionador (NO/ILU0/ILUT) [default = ' // trim(conf%usar_precond) //']'
        print*, '   [12] Fill in de ILUT (%) [default = ' // realtostr(conf%fill_in_prc) //']'
        print*, '   [13] Drop tolerance de ILUT [default = ' , conf%drop_tolerance , ']'
        print*, '   [14] Tipo de analisis RCS (mono/bi/monof) [default = ' , trim(conf%analisis_rcs) , ']'
        print*, '   [15] # de muestras RCS monostatica (0°<=phi<180°) [default = ' // inttostr(conf%n_pasos_rcs_mono) //']'
        print*, '   [16] Tipo de polarizacion (VV/HH) [default = ' , trim(conf%tipo_polarizac) , ']'
        print*, '   [17] Variable a especificar (LAMBDA/FREC) [default = ' , trim(conf%variable_espec) , ']'

        read*, ent
        if (ent == 1) then
            exit
        else if (ent == 2) then

            print*, 'Num de muestras en THETA al calcular RCS:'
            read*, conf%numThetaRCS
        else if (ent == 3) then
            print*, 'Num de muestras en PHI al calcular RCS:'
            read*, conf%numPhiRCS
        elseif (ent == 4) then
            print*, 'Error del solver CGS:'
            read*, conf%errorSolverCGS
        elseif (ent == 5) then
            print*, 'Numero forzado de multipolos (-1 para automatico:'
            read*, conf%Lforz
        elseif (ent == 6) then
            print*, 'Angulo THETA de incidencia de la onda (Grados °)'
            read*, conf%dirTheta
        elseif (ent == 7) then
            print*, 'Angulo PHI de incidencia de la onda (Grados °)'
            read*, conf%dirPhi
        elseif (ent == 8) then
            print*, 'Frecuencia inicial de barrido RCS'
            read*, conf%rcsMonoFMin
        elseif (ent == 9) then
            print*, 'Frecuencia final de barrido RCS'
            read*, conf%rcsMonoFMax
        elseif (ent == 10) then
            print*, 'Numero de muestras RCS barrido frecuencia'
            read*, conf%rcsMonoSamples
        elseif (ent == 11) then
            print*, 'Usar precondicionador (NO/ILU0/ILUT)'
            read*, conf%usar_precond
            if (to_upper(trim(conf%usar_precond)) == 'NO') then
                conf%usar_precond = 'NO'
            else if (to_upper(trim(conf%usar_precond)) == 'ILUT') then
                conf%usar_precond = 'ILUT'
            else
                conf%usar_precond = 'ILU0'
            end if
        elseif (ent == 12) then
            print*, 'Fill in de ILUT (%)'
            read*, conf%fill_in_prc
        elseif (ent == 13) then
            print*, 'Drop tolerance de ILUT'
            read*, conf%drop_tolerance            
        elseif (ent == 14) then
            print*, 'Tipo de analisis RCS (mono/bi/monof)'
            read*, conf%analisis_rcs
            if (to_upper(trim(conf%analisis_rcs)) == 'MONO') then
                conf%analisis_rcs = 'MONO'
            elseif (to_upper(trim(conf%analisis_rcs)) == 'MONOF') then
                conf%analisis_rcs = 'MONOF'            
            else
                conf%analisis_rcs = 'BI'
            end if
        elseif (ent == 15) then
            print*, 'Numero de muestras para RCS monostatica'
            read*, conf%n_pasos_rcs_mono
        elseif (ent == 16) then
            print*, 'Tipo de polarizacion (VV/HH)'
            read*, conf%tipo_polarizac
            if (to_upper(trim(conf%tipo_polarizac)) == 'HH') then
                conf%tipo_polarizac = 'HH'
            else
                conf%tipo_polarizac = 'VV'
            end if
        elseif (ent == 17) then
            print*, 'Variable a especificar (LAMBDA/FREC)'
            read*, conf%variable_espec
            if (to_upper(trim(conf%variable_espec)) == 'LAMBDA') then
                conf%variable_espec = 'LAMBDA'
            else
                conf%variable_espec = 'FREC'
            end if            
        end if

    end do

    error_solver = conf%errorSolverCGS
    dir_onda(1) = sin(conf%dirTheta*pi/180.)*cos(conf%dirPhi*pi/180.)
    dir_onda(2) = sin(conf%dirTheta*pi/180.)*sin(conf%dirPhi*pi/180.)
    dir_onda(3) = cos(conf%dirTheta*pi/180.)

    call ajustar_polariz(dir_onda,pol_onda)

    print*, 'Ingrese ' // trim(conf%variable_espec) // ' (GHz o m):'

    if (to_upper(trim(conf%variable_espec)) == 'LAMBDA') then
        read*, lambda
        call setlambda(lambda)
    else
        read*, frequency
        call setf(frequency*1000000000_dp)
    end if


    if (numtest == 0) then
        !call testtest()
    elseif (numtest == 1) then
    !        call testComparaImatxvecMLFMA()
    elseif (numtest == 2) then
    !        call testComparaImatxvec()
    elseif (numtest == 3) then
        call testMLFMA()
    elseif (numtest == 4) then
        call testFMM()
    elseif (numtest == 5) then
        call testMoM()
    elseif (numtest == 6) then
    !        call testcomparaZMLFMA()
    elseif (numtest == 7) then
    !        call testComparaZfar()
    elseif (numtest == 8) then
    !        call testRCSMoM_monof(conf%rcsMonoFMin*1000000, conf%rcsMonoFMax*1000000, conf%rcsMonoSamples)
    end if


end subroutine

    subroutine ajustar_polariz(vdirecc, vpolariz)
        !dummy
        real (kind = dp), dimension(3), intent(inout) :: vdirecc
        complex (kind = dp), dimension(3), intent(inout) :: vpolariz
        !
        !local
        real (kind = dp), dimension(3) :: vector_vertical
        !

        if (streq(conf%tipo_polarizac,'VV')) then
            vpolariz(1) = 0.
            vpolariz(2) = 0.
            vpolariz(3) = 1.
        elseif (streq(conf%tipo_polarizac, 'HH')) then
            vector_vertical(1) = 0.
            vector_vertical(2) = 0.
            vector_vertical(3) = 1.
            vpolariz = cross_productR(vdirecc, vector_vertical)
            vpolariz = vpolariz/sqrt(sum(vpolariz*conjg(vpolariz)))
        end if

    end subroutine

    subroutine crear_precond(n_nozero, a_near, ja_near, ia_near, status)
        complex (kind = dp_near), dimension(:) :: a_near
        integer (kind = il), dimension(:) :: ja_near
        integer (kind = ild), dimension(:) :: ia_near
        integer (kind = ild), allocatable, dimension(:) :: iw
        integer ( kind = ild ), allocatable, dimension(:) :: jr
        integer ( kind = ild ), allocatable, dimension(:) :: jwl
        integer ( kind = ild ), allocatable, dimension(:) :: jwu
        complex ( kind = dp ), allocatable, dimension(:) :: wl
        complex ( kind = dp ), allocatable, dimension(:) :: wu

        integer (kind = ild) :: n_nozero, tamano_esperado
        integer (kind = ilo) :: status

        integer (kind = il) :: fill_in
        real (kind = dp) :: drop_tol

        if (allocated(precond_jlu)) then
            deallocate(precond_jlu)
        end if 
        if (allocated(precond_alu)) then
            deallocate(precond_alu)
        end if
        if (allocated(precond_ju)) then
            deallocate(precond_ju)
        end if
        if (allocated(iw)) then
            deallocate(iw)
        end if

        call printconsole('Creando precondicionador a partir de la Z cercana la cual tiene ' // &
            dbllongtostr(n_nozero) // ' elementos', 'creating_precond/' // dbllongtostr(n_nozero))
        print*, ' - - - ' // gettiempo()
        if (trim(conf%usar_precond) == 'ILU0') then    
            allocate(precond_jlu(n_nozero + 1))
            allocate(precond_alu(n_nozero + 1))
            allocate(precond_ju(num_e))
            allocate(iw(num_e))
            precond_jlu(:) = 0
            precond_alu(:) = 0.
            precond_ju(:) = 0
            iw(:) = 0

            call ilu0C(num_e, a_near, ja_near, ia_near, precond_alu, precond_jlu, precond_ju, iw, status)
            deallocate(iw)
            if (status == 0) then
               call printconsole('Hecho', 'done')
            else
               call printconsole('ERROR en precondicionador, se encontro pivote nulo en indice'// &
                inttostr(status), 'error_precond_nullpivot')
            end if

        else if (trim(conf%usar_precond) == 'ILUT') then

            fill_in = ceiling(0.01*conf%fill_in_prc*num_e)

            call printconsole('Numero de elementos de fill_in: ' // inttostr(fill_in), 'num_fillin/' // inttostr(fill_in) )
            drop_tol = conf%drop_tolerance
            tamano_esperado = 2*num_e*fill_in + n_nozero + 1

            call printconsole('Tamaño de la matriz asignado para el precondicionador :' // dbllongtostr(tamano_esperado), &
                'precond_size/' // dbllongtostr(tamano_esperado))
            allocate(precond_jlu(tamano_esperado))
            allocate(precond_alu(tamano_esperado))
            allocate(precond_ju(num_e))
            allocate(jwl(num_e))
            allocate(jwu(num_e))
            allocate(jr(num_e))
            allocate(wl(num_e))
            allocate(wu(num_e + 1))
            precond_jlu(:) = 0
            precond_alu(:) = 0.
            precond_ju(:) = 0
            jwl(:) = 0
            jwu(:) = 0
            jr(:) = 0
            wl(:) = 0.
            wu(:) = 0.


            call ilutC (num_e, a_near, ja_near, ia_near, fill_in, drop_tol, precond_alu, precond_jlu, precond_ju, &
                tamano_esperado, wu, wl, jr, jwl, jwu, status )
            deallocate(wu)
            deallocate(wl)
            deallocate(jr)
            deallocate(jwu)
            deallocate(jwl)

            if (status == 0) then
                call printconsole( 'Hecho', 'done')
            else if (status > 0) then
                call printconsole( 'ERROR en precondicionador, se encontro pivote nulo en indice' // inttostr(status), &
                    'error_precond_nullpivot')
            else if (status == -1) then
                call printconsole( 'ERROR en precondicionador: input matrix may be wrong. (The elimination process has &
                   & generated a row in L or U whose length is >  n.)', 'error_precond_LU')
            else if (status == -2) then
                call printconsole( 'ERROR en precondicionador: The matrix L overflows the array alu.', 'error_precond_L')
            else if (status == -3) then
                call printconsole( 'ERROR en precondicionador: The matrix U overflows the array alu.', 'error_precond_U')
            else if (status == -4) then
                call printconsole( 'ERROR en precondicionador: Illegal value for lfil.', 'error_precond_lfill')
            else if (status == -5) then
                call printconsole( 'ERROR en precondicionador: zero pivot encountered.', 'error_precond_zero')
            end if

        end if
        call printconsole( 'Memoria ocupada por el precondicionador: ', '-')
        call printconsole( inttostr( int((sizeof(precond_alu) + sizeof(precond_ju) + sizeof(precond_jlu))/1000000.)) // &
            ' MB', 'precond_mem_', .true.)
    end subroutine


    subroutine gen_rcs_monostatica(test_no, prec_activo)
        !dummy
        integer (kind = ilo), intent(in) :: test_no
        integer (kind = ilo), intent(in) :: prec_activo
        !
        !local
        integer (kind = il) :: ang
        real (kind = dp) :: angphi, dephi
        real (kind = dp), dimension(3) :: direccion, prev_direccion, deltaR
        complex (kind = dp), dimension(3) :: polariz
        complex (kind = dp), allocatable, dimension(:) :: b_mono
        real (kind = dp), allocatable, dimension(:,:) :: mat_rcs_mono, mat_currXtriang, mat_currXedgePART
        complex (kind = dp), allocatable, dimension(:,:) :: mat_currXedge
        character (len = 60) :: fname
        logical :: flag_both
        character (len = 60) :: local_pol
        complex (kind = dp), allocatable, dimension(:) :: phaseCorrection
        integer (kind = dp) :: ibase

        allocate(phaseCorrection(num_e))

        local_pol = conf%tipo_polarizac
        if (streq(conf%tipo_polarizac,'VH')) then
            flag_both = .true.
            conf%tipo_polarizac = 'VV'
        else
            flag_both=.false.
        end if

        allocate(mat_currXtriang(num_t + 1, conf%n_pasos_rcs_mono))
        allocate(mat_currXedge(num_e, conf%n_pasos_rcs_mono))
        allocate(mat_currXedgePART(num_e, conf%n_pasos_rcs_mono))
    
        do
            mat_currXtriang(:,:) = 0.
            mat_currXedge(:,:) = 0.
            allocate(b_mono(num_e))
            allocate(mat_rcs_mono(conf%n_pasos_rcs_mono,2))

            mat_rcs_mono(:,:) = 0.
            guess(:) = 0.

            dephi = 1.*pi/conf%n_pasos_rcs_mono

            prev_direccion(3) = 0.
            prev_direccion(1) = -1.
            prev_direccion(2) = 0.
            
            do ang = conf%rcsmonostepi, conf%rcsmonostepf
                call printconsole( ' ', '-')
                call printconsole( '*****Calculando RCS monostatica... Polarizacion: ' // conf%tipo_polarizac, 'rcsmono_pol/'&
                 // conf%tipo_polarizac)
                call printconsole( '*****Angulo de incidencia ' // inttostr(ang) // ' de ' // inttostr(conf%n_pasos_rcs_mono), &
                    'rcsangle/' // inttostr(ang))
                call printconsole( '*****Angulo de incidencia ' // inttostr(ang) // ' de ' // inttostr(conf%n_pasos_rcs_mono), &
                    'rcs_total_angles/' //  inttostr(conf%n_pasos_rcs_mono))
                angphi = (ang - 1)*dephi
                direccion(3) = 0.
                direccion(1) = -cos(angphi)
                direccion(2) = -sin(angphi)

                !Ajuste de fase
                deltaR = direccion - prev_direccion
                do ibase = 1, num_e
                    phaseCorrection(ibase) = exp( (-jj)*getkappa()*dot_product(deltaR, e_centro(:,ibase) )  )
                end do
                guess = guess*phaseCorrection
                !

                call ajustar_polariz(direccion, polariz)

                b_mono = generar_onda(direccion, ctte_onda, polariz)

                if (test_no == 3) then
                    if (prec_activo == 0) then
                        call CGSold(MatxVec_MLFMA,b_mono, guess, error_solver, num_e, I, 1, precond_alu, precond_jlu, precond_ju)
                    else
                        call CGSold(MatxVec_MLFMA,b_mono, guess, error_solver, num_e, I, 1)
                    end if
                elseif (test_no == 4) then
                    if (prec_activo == 0) then
                        call CGSold(MatxVec_FMM,b_mono, guess, error_solver, num_e, I, 1, precond_alu, precond_jlu, precond_ju)
                    else
                        call CGSold(MatxVec_FMM,b_mono, guess, error_solver, num_e, I, 1)
                    end if
                elseif (test_no == 5) then
                    if (prec_activo == 0) then
                        call CGSold(MatxVec_MoM,b_mono, guess, error_solver, num_e, I, 1, precond_alu, precond_jlu, precond_ju)
                    else
                        call CGSold(MatxVec_MoM,b_mono, guess, error_solver, num_e, I, 1)
                    end if
                end if

                call inicializar_RCS (e_centro, I, e_long, t_baric, e_t, num_e, num_t, polariz, ctte_onda, p_coord, t_p, e_po)
                mat_rcs_mono(ang,1) = angphi*180./pi
                mat_rcs_mono(ang,2) = 10*log10(sigmaRCS(pi/2.,angphi))
                guess = I
                prev_direccion = direccion
                call guardarMatrizR(mat_rcs_mono, trim(getResDir()) // trim(conf%tipo_polarizac) // 'RCSmo_azm_incremental.dat', &
                    .true.)
                
                mat_currXtriang(1,ang) = angphi*180./pi
                mat_currXtriang(2:num_t + 1, ang) = gen_currentspertriangle(I)

                mat_currXedge(:, ang) = I

                if (conf%save_incrementals) then
                    call guardarMatrizR(mat_currXtriang, trim(getResDir()) // trim(conf%tipo_polarizac) //&
                     'Isurf_incremental.dat', .true.)
                    mat_currXedgePART = aimag(mat_currXedge)
                    call guardarMatrizR(mat_currXedgePART, trim(getResDir()) // trim(conf%tipo_polarizac) // &
                        'IedgeIM_incremental.dat',.true.)
                    mat_currXedgePART = real(mat_currXedge)                
                    call guardarMatrizR(mat_currXedgePART, trim(getResDir()) // trim(conf%tipo_polarizac) // &
                        'IedgeRE_incremental.dat',.true.)
                end if
                print*, ' '

            end do

            if (streq(conf%test_name,strdefault)) then
                print*, 'Introduzca el nombre <base> de archivos con resultados: '
                read*, fname
            else
                fname = trim(conf%test_name)
            end if

            call guardarMatrizR(mat_rcs_mono, trim(getResDir()) //  'RCSmo_azm' // trim(conf%tipo_polarizac) // '.dat', .true.)
            call guardarMatrizR(mat_currXtriang, trim(getResDir()) //  'Isurf_azm' // trim(conf%tipo_polarizac) // '.dat', .true.)
            mat_currXedgePART = real(mat_currXedge)
            call guardarMatrizR(mat_currXedgePART, trim(getResDir()) //  'Iedge_azmRE' // trim(conf%tipo_polarizac) // &
                '.dat',.true.)
            mat_currXedgePART = aimag(mat_currXedge)            
            call guardarMatrizR(mat_currXedgePART, trim(getResDir()) //  'Iedge_azmIM' // trim(conf%tipo_polarizac) // &
                '.dat',.true.)

            deallocate(mat_rcs_mono)
            deallocate(b_mono)

            if (flag_both .eqv..false.) then
                exit
            else
                flag_both = .false.
                conf%tipo_polarizac = 'HH'
            end if

        end do

        if (prec_activo == 0) then
            deallocate(precond_alu)
            deallocate(precond_jlu)
            deallocate(precond_ju)                        
        end if        

        conf%tipo_polarizac = local_pol
        deallocate(phaseCorrection)
    end subroutine


    subroutine testMLFMA()
        !local
        character (len = 50) :: file_i_out, ent
        integer(kind = il) :: ierr
        logical :: isOriginalVH
        !

        ierr = -1
        call ajustar_polariz(dir_onda,pol_onda)
        call cap_geom(.true.)

        call printconsole( '********testMLFMA*********', 'testMLFMA')
        call timestamp

        call inicializar_MLFMA(e_centro, num_e, e_t, e_po, e_long, t_area, t_baric, t_baric_sub, p_coord, t_p, e_p, num_p,&
         num_t, t_normal)
        
        if (streq(conf%tipo_polarizac, 'VH')) then
            isOriginalVH = .true.
        else
            isOriginalVH = .false.
        end if

        if (trim(conf%usar_precond) /= 'NO') then
            !crear precondicionador a partir de la Z cercana (Zespacida) de MLFMA
            call crear_precond(n_esparcida_MLFMA, Zesparcida_MLFMA, ja_MLFMA, ia_MLFMA, ierr)
            !
        end if

        if (streq(conf%analisis_rcs, 'BI')) then
            call gen_rcs_bistatic(3, conf%numThetaRCS, conf%numPhiRCS, ierr)
        elseif (streq(conf%analisis_rcs, 'MONO')) then
            if (ierr == 0) then
                call gen_rcs_monostatica(3,0)
            else
                call gen_rcs_monostatica(3,-1)
            end if
        elseif (streq(conf%analisis_rcs, 'MONOF')) then
            print*, 'La rcs monostatica de barrido en frecuencia no esta disponible aun'
        end if

        call timestamp

        if (isOriginalVH) then
            conf%tipo_polarizac = 'VH'
        end if

        call printconsole( '***FINALIZADO***', 'end_test')
        call destruir_todo()
    end subroutine


    subroutine testFMM()
        !local
        character (len = 50) :: file_i_out, ent
        integer(kind = il) :: ierr
        logical :: isOriginalVH
        !
        call ajustar_polariz(dir_onda,pol_onda)
        ierr = -1
        call cap_geom(.true.)
        call printconsole( '********testFMM*********', 'testFMM')
        call timestamp

        !guess(:) = 0.
        call inicializar_FMM(e_centro, num_e, e_t, e_po, e_long, t_area, t_baric, t_baric_sub, p_coord, t_p, e_p, num_p,&
         num_t, t_normal)

        if (streq(conf%tipo_polarizac, 'VH')) then
            isOriginalVH = .true.
        else
            isOriginalVH = .false.
        end if        

        if (trim(conf%usar_precond) /= 'NO') then
            !crear precondicionador a partir de la Z cercana (Zespacida) de MLFMA
            call crear_precond(n_esparcida_FMM, Zesparcida_FMM, ja_FMM, ia_FMM, ierr)
            !
        end if

        if (streq(conf%analisis_rcs, 'BI')) then
            call gen_rcs_bistatic(4, conf%numThetaRCS, conf%numPhiRCS, ierr)
        elseif (streq(conf%analisis_rcs, 'MONO')) then
            if (ierr == 0) then
                call gen_rcs_monostatica(4,0)
            else
                call gen_rcs_monostatica(4,-1)
            end if
        elseif (streq(conf%analisis_rcs, 'MONOF')) then
            print*, 'La rcs monostatica de barrido en frecuencia no esta disponible aun'
        end if


        call timestamp

        if (isOriginalVH) then
            conf%tipo_polarizac = 'VH'
        end if
        
        call printconsole( '***FINALIZADO***', 'end_test')
        call destruir_todo()
    end subroutine

    subroutine testMoM()
        !local
        character (len = 50) :: file_i_out, ent
        complex (kind = dp), allocatable, dimension(:) :: Zmom_sparse
        complex (kind = dp_near), allocatable, dimension(:) :: Zmom_csr
        integer (kind = il), allocatable, dimension(:) :: i_Zmom_sparse, j_Zmom_sparse, i_Zmom_csr, j_Zmom_csr
        integer(kind = il), allocatable, dimension(:) :: iw
        integer(kind = il) :: ierr
        integer (kind = il) :: fil, col, cont
        logical :: isOriginalVH
        !
        call ajustar_polariz(dir_onda,pol_onda)
        call printconsole( '********testMoM*********', 'testMoM')
        ierr = -1

        call cap_geom()

        call timestamp

        call MoM

        guess(:) = 0.

        if (streq(conf%tipo_polarizac, 'VH')) then
            isOriginalVH = .true.
        else
            isOriginalVH = .false.
        end if        

        call setZ_MatxVec_MoM(Z)

        if (streq(conf%analisis_rcs, 'BI')) then
            call gen_rcs_bistatic(5, conf%numThetaRCS, conf%numPhiRCS, ierr)
        elseif (streq(conf%analisis_rcs, 'MONO')) then
            if (ierr == 0) then
                call gen_rcs_monostatica(5,0)
            else
                call gen_rcs_monostatica(5,-1)
            end if
        elseif (streq(conf%analisis_rcs, 'MONOF')) then
            print*, 'La rcs monostatica de barrido en frecuencia no esta disponible aun'
        end if

        call timestamp

        if (isOriginalVH) then
            conf%tipo_polarizac = 'VH'
        end if
        call printconsole( '***FINALIZADO***', 'end_test')
    end subroutine


    subroutine ComparaRCS(I1, I2, numTheta, numPhi)
        !dummy
        complex (kind = dp), intent(in), dimension(:) :: I1,I2
        integer (kind = il), intent(in) :: numTheta, numPhi
        !
        !local
        real (kind = dp), allocatable, dimension(:,:) :: rcsI1, rcsI2, rcsIDif, rcslegend
        real (kind = dp) :: normI1, normI2, normIRest, normIDif
        character (len = 30) :: nArch
        !

        allocate(rcsI1(numTheta, numPhi))
        allocate(rcsI2(numTheta, numPhi))
        allocate(rcsIDif(numTheta, numPhi))
        allocate(rcslegend(numTheta + 1, numPhi + 1))

        call printconsole(' ', 'calculating_rcs')
        print*, 'Calculando RCS de las corrientes dadas por MoM en todas las direcciones...'
        rcslegend = RCS_bi(I1,numTheta,numPhi)
        rcsI1 = rcslegend(2:(numTheta+1),2:(numPhi+1))
        normI1 = norm2(rcsI1)
        print*, 'Calculando RCS de las corrientes dadas por FMM/MLFMA en todas las direcciones...'
        rcslegend = RCS_bi(I2,numtheta,numphi)
        rcsI2 = rcslegend(2:(numTheta+1),2:(numPhi+1))
        normI2 = norm2(rcsI2)
        normIRest = norm2(rcsI1-rcsI2)
        print*, 'Calculando RCS de las corrientes de diferencia en todas las direcciones...'
        rcslegend = RCS_bi(I1-I2,numtheta,numphi)
        rcsIDif = rcslegend(2:(numTheta+1),2:(numPhi+1))
        normIDif = norm2(rcsIDif)


        print*, 'Error relativo porcentual de matriz RCS de FMM/MLFMA respecto a MoM:'
        print*, 100*normIRest/normI1
        print*, 'Error relativo porcentual de matriz RCS de corrientes diferencia respecto a MoM: '
        print*, 100*normIDif/normI1


        if (streq(conf%test_name,strdefault)) then
            print*, 'Desea guardar los resultados de la RCS en archivos? (S/N)'
            read*, nArch
        else
            nArch = 'S'
        end if

        if (nArch == 'N' .or. nArch == 'n') then
        else

            if (streq(conf%test_name,strdefault)) then
                print*, 'Introduzca el nombre <base> de archivos con resultados: '
                read*, nArch
            else
                nArch = trim(conf%test_name) // '_RCSBI'
            end if



            call printconsole( 'Guardando archivos', 'saving_rcs')
            call guardarMatrizR(rcsI1, trim(nArch) // '_MoM.dat')
            call guardarMatrizR(rcsI2, trim(nArch) // '_FMM-MLFMA.dat')
            call guardarMatrizR(rcsIDif, trim(nArch) // '_I-DIF.dat')

            call guardarVectorR(gen_currentspertriangle(I2), trim(nArch) // '_Ixtriang.dat', .true.)

        end if

    end subroutine


    subroutine gen_rcs_bistatic(test_no, numTheta, numPhi, isPreconditioned, I1)
        !dummy
        integer (kind = ilo), intent(in) :: test_no
        complex (kind = dp), optional, intent(in), dimension(:) :: I1
        integer (kind = il), intent(in) :: numTheta, numPhi
        integer (kind = ilo), intent(in) :: isPreconditioned
        !
        !local
        complex (kind = dp), allocatable,  dimension(:) :: Iwork, bwork
        real (kind = dp), allocatable, dimension(:,:) :: rcsI1, radPattern
        real (kind = dp), allocatable, dimension(:,:,:) :: polarCoord
        real (kind = dp) :: normI1
        character (len = 30) :: nArch
        logical :: isCurrentAvailable
        integer ( kind = ilo) :: itrcs, startrcs, endrcs, itheta, iphi
        real (kind = dp) :: vtheta, vphi
        !

        isCurrentAvailable = present(I1)
        allocate(Iwork(num_e))
        allocate(bwork(num_e))
        allocate(rcsI1(numTheta + 1, numPhi + 1))
        allocate(radPattern(numTheta, numPhi))
        allocate(polarCoord(numTheta, numPhi, 3))
        
        if (streq(conf%tipo_polarizac,'VV')) then
            startrcs = 1
            endrcs = 1
        elseif (streq(conf%tipo_polarizac,'HH')) then
            startrcs = 2
            endrcs = 2
        elseif (streq(conf%tipo_polarizac,'VH')) then
            startrcs = 1
            endrcs = 2
        end if

        do itrcs = startrcs, endrcs

            if (itrcs == 1) then
                conf%tipo_polarizac = 'VV'
            elseif (itrcs == 2) then
                conf%tipo_polarizac = 'HH'
            end if

            call printconsole( '*****Calculando RCS biestatica... Polarizacion: ' // trim(conf%tipo_polarizac), 'rcsbi_pol/'&
             // trim(conf%tipo_polarizac))

            guess(:) = 0.
            call ajustar_polariz(dir_onda,pol_onda)
            bwork = generar_onda(dir_onda, ctte_onda, pol_onda)

            if (isCurrentAvailable) then
                Iwork = I1
            else
                if (test_no == 3) then
                    if (isPreconditioned == 0) then
                        call CGSold(MatxVec_MLFMA,bwork, guess, error_solver, num_e, Iwork, 1,precond_alu,precond_jlu,precond_ju)
                    else
                        call CGSold(MatxVec_MLFMA,bwork, guess, error_solver, num_e, Iwork, 1)
                    end if
                elseif (test_no == 4) then
                    if (isPreconditioned == 0) then
                        call CGSold(MatxVec_FMM,bwork, guess, error_solver, num_e, Iwork, 1, precond_alu, precond_jlu, precond_ju)
                    else
                        call CGSold(MatxVec_FMM,bwork, guess, error_solver, num_e, Iwork, 1)
                    end if
                elseif (test_no == 5) then
                    if (isPreconditioned == 0) then
                        call CGSold(MatxVec_MoM,bwork, guess, error_solver, num_e, Iwork, 1, precond_alu, precond_jlu, precond_ju)
                    else
                        call CGSold(MatxVec_MoM,bwork, guess, error_solver, num_e, Iwork, 1)
                    end if
                end if

            end if

            call printconsole(' - - - ' // gettiempo(), 'calculating_rcs')
            print*, 'Calculando RCS de las corrientes dadas por FMM/MLFMA en todas las direcciones...'
            rcsI1 = RCS_bi(Iwork,numtheta,numphi)

            if (streq(conf%test_name,strdefault)) then
                print*, 'Introduzca el nombre <base> de los resultados de salida: '
                read*, nArch
            else
                nArch = trim(conf%test_name)
            end if

            call printconsole( 'Guardando archivos', 'saving_rcs')

            radPattern(:,:) = rcsI1(2:(numTheta + 1), 2:(numPhi + 1))

            radPattern =  radPattern - maxval(radPattern) + 50

            do itheta = 1,numTheta
                vtheta = rcsI1(itheta + 1, 1)
                do iphi = 1,numPhi
                    vphi = rcsI1(1,iphi + 1)
                    if (radPattern(itheta,iphi) < 0) then
                        radPattern(itheta, iphi) = 0.
                    end if
                    polarCoord(itheta, iphi, 1) = radPattern(itheta,iphi)*sin(vtheta*pi/180.)*cos(vphi*pi/180.)
                    polarCoord(itheta, iphi, 2) = radPattern(itheta,iphi)*sin(vtheta*pi/180.)*sin(vphi*pi/180.)
                    polarCoord(itheta, iphi, 3) = radPattern(itheta,iphi)*cos(vtheta*pi/180.)
                end do
            end do

            call guardarMatrizR(polarCoord(:,:,1), trim(getResDir()) //  'RadPatternX' // trim(conf%tipo_polarizac) // &
                '.dat' , .true.)
            call guardarMatrizR(polarCoord(:,:,2), trim(getResDir()) //  'RadPatternY' // trim(conf%tipo_polarizac) // &
                '.dat' , .true.)
            call guardarMatrizR(polarCoord(:,:,3), trim(getResDir()) //  'RadPatternZ' // trim(conf%tipo_polarizac) // &
                '.dat' , .true.)
            call guardarMatrizR(rcsI1, trim(getResDir()) //  'RCSbi' // trim(conf%tipo_polarizac) // '.dat' , .true.)
            call guardarVectorR(gen_currentspertriangle(Iwork), trim(getResDir()) //  'Isurfbi' // trim(conf%tipo_polarizac)&
             // '.dat' , .true.)
            call guardarVectorC(Iwork, trim(getResDir()) // 'Iedge' // trim(conf%tipo_polarizac) // '.dat' )


        end do

        if (isPreconditioned == 0) then
            deallocate(precond_alu)
            deallocate(precond_jlu)
            deallocate(precond_ju)
        end if

    end subroutine    

    function RCS_bi(Corriente, numtheta, numphi) result(MRCS)
        !dummy
        complex (kind = dp), intent(in), dimension(:) :: Corriente
        integer (kind = il), intent(in) :: numtheta, numphi
        real (kind = dp), allocatable, dimension(:,:) :: MRCS
        !
        !local
        integer (kind = il) :: p,t
        real (kind = dp) :: detheta, dephi, fi, tita
        integer (kind = il) :: porc, prevPorc
        !
        prevPorc = -1.
        detheta = pi/numtheta
        dephi = 2.*pi/numphi
        allocate(MRCS(numtheta+1,numphi+1))


        call inicializar_RCS (e_centro, Corriente, e_long, t_baric, e_t, num_e, num_t, pol_onda, ctte_onda, p_coord, t_p, e_po)
        call setprc(numphi)
        do p = 1, numphi
            fi = (2*p-1)*dephi/2.
            MRCS(1,p+1) = fi*180./pi
            do t = 1, numtheta
                tita = (2*t-1)*detheta/2.
                MRCS(t+1,p+1) = 10*log10(sigmaRCS(tita,fi))
                MRCS(t+1,1) = tita*180./pi
            end do
            call updateprc(p)
        end do
        MRCS(1,1) = 0.

    end function


    function gen_currentspertriangle(resultI) result(curr)
        
        complex (kind = dp), dimension(:), intent(in) :: resultI
        real (kind = dp), allocatable, dimension(:) :: curr
        !local
        complex (kind = dp) :: coefI
        complex (kind = dp), dimension(3, num_t) :: currVector
        integer (kind = il) :: bf, tri_plus, tri_minus, tri
        real (kind = dp), dimension(3) :: vopp_plus, vopp_minus, baric_plus, baric_minus, rho_plus, rho_minus
        real (kind = dp) :: edge_long, area_plus, area_minus
        !
        allocate(curr(num_t))
        currVector(:,:) = 0.

        do bf = 1, num_e
            coefI = resultI(bf)
            tri_plus = e_t(1, bf)
            tri_minus = e_t(2, bf)
            vopp_plus = p_coord(:, e_po(1,bf))
            vopp_minus = p_coord(:, e_po(2,bf))
            baric_plus = t_baric(:, tri_plus)
            baric_minus = t_baric(:, tri_minus)
            rho_plus = baric_plus - vopp_plus
            rho_minus =  vopp_minus - baric_minus
            area_plus = t_area(tri_plus)
            area_minus = t_area(tri_minus)
            edge_long = e_long(bf)


            !suma en triangulo +
            currVector(:, tri_plus) = currVector(:, tri_plus) + coefI*rho_plus*edge_long/(2*area_plus)
            !suma en triangulo -
            currVector(:, tri_minus) = currVector(:, tri_minus) + coefI*rho_minus*edge_long/(2*area_minus)
        end do

        if (conf%isurf_phase) then
            !parte real de Js
            do tri = 1, num_t
                curr(tri) = sqrt(sum(   real(currVector(:,tri))*real(currVector(:,tri))   ))
            end do
            !!!!!!!!!!!!!!!
        else
            !magnitud de Js
            do tri = 1, num_t
                curr(tri) = sqrt(sum(currVector(:,tri)*conjg(currVector(:,tri))))
            end do
            !!!!!!!!!!!!!!!
        end if


    end function


    subroutine cap_geom(skipZ)
        !dummy
        logical, optional :: skipZ
        character (len = 30) :: fname
        !


        ! Capturar Geometría en base a RWG
        if (streq(conf%file_name,strdefault)) then
            print*, 'Archivo a cargar:'
            read *, conf%file_name
        else
            call printconsole( 'Archivo a cargar: ' // conf%file_name, 'filename/' // conf%file_name)
        end if

        call printconsole( 'Capturando geometria...', 'capturing_geometry')
        call printconsole( 'Leyendo STL...', 'reading_STL')

        call obtener_p_t(conf%file_name, p_coord, t_p, num_p, num_t,t_normal)

        call printconsole( 'Generando parametros de la geometria...', 'generating_geom_parameters')

        call parametros_rwg (p_coord, t_p, num_t, e_p, e_long, e_t, e_po, num_e, e_centro)
        call hallar_geometria_rwg (p_coord,t_p,num_p,num_t,t_area,t_baric,t_baric_sub)
call pruebaNormal()
        if (streq(conf%test_name,strdefault)) then
            fname = 'lasttest'
        else
            fname = trim(conf%test_name)
        end if

        call guardarMatrizR(transpose(p_coord), trim(getResDir()) //  'points.dat', .true.)
        call guardarMatrizI(transpose(t_p), trim(getResDir()) //  'faces.dat', .true.)

        if (present(skipZ)) then
            if (skipZ) then
            else
                allocate (Z(num_e,num_e))
            end if
        else
            allocate (Z(num_e,num_e))
        end if

        allocate(I(num_e))
        allocate(guess(num_e))
        allocate(b(num_e)) 


        call printconsole( 'Geometria Capturada. ' // inttostr(num_e) // ' incognitas', 'end_capture_geom/' // inttostr(num_e))

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    end subroutine

subroutine pruebaNormal()

integer ( kind = il) :: ti
real (kind = dp) :: pe

integer (kind = il), dimension(6) :: limits
real (kind = dp) :: D, maxX, minX, maxY, minY, maxZ, minZ, rzl, lambda
real (kind = dp), dimension (3) :: roo
logical :: therewaserrnormal

therewaserrnormal = .false.

limits = limites_dispersor(e_centro, num_e)
maxX = e_centro(1,limits(1))
minX = e_centro(1,limits(2))
maxY = e_centro(2,limits(3))
minY = e_centro(2,limits(4))
maxZ = e_centro(3,limits(5))
minZ = e_centro(3,limits(6))
roo(1) = 0.5*(maxX+minX)
roo(2) = 0.5*(maxY+minY)
roo(3) = 0.5*(maxZ+minZ)
!print*, 'THIS IS ROO: ', roo
!print*, minX, maxX, minY, maxY, minZ, maxZ
!print*, ' '


do ti = 1, num_t
    pe = dot_product(t_normal(:,ti), t_baric(:,ti))
    if (pe <= 0.) then
        therewaserrnormal = .true.
t_normal(:, ti) = -t_normal(:,ti)
    end if
end do
if (therewaserrnormal) then
    print*, 'Hubo normales equivocadas en el STL y se han cambiado'
end if

end subroutine


    subroutine MoM()

        !Usar modulo MoM
        call inicializar_MoM(p_coord, t_baric, t_baric_sub, t_area, e_long, t_p, e_p, e_t, e_po, num_e, num_p, num_t, t_normal)
        call printconsole( 'Llenando matriz Z...', 'filling_Z')
        call llenarZ_MoM(Z,1)
        call printconsole( 'Hecho', 'done')
        call printconsole( 'Generando onda...', 'generating_wave')

        b = generar_onda(dir_onda, ctte_onda, pol_onda)
        call printconsole( 'Hecho', 'done')
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    end subroutine


    function generar_onda_old(direccion, ctte, polarizacion) result (res)
        !subroutine generar_onda(direccion, ctte, polarizacion, b)
        !dummy
        real (kind = dp), dimension(3), intent(in) :: direccion
        complex (kind = dp), dimension(3), intent(in) :: polarizacion
        complex (kind = dp), intent(in) :: ctte
        complex (kind = dp), dimension(:), allocatable :: res
        !complex (kind = dp), dimension(num_e), intent(out) :: b
        !
        !local
        integer (kind = il) :: i
        real (kind = dp), dimension(3) :: baricentro
        complex (kind = dp), dimension(3) :: ondaE, ondaH
        real (kind = dp), dimension(3) :: Udir, normalPlus, normalMinus
        complex (kind = dp) :: UNO
        !
        allocate(res(num_e))
        UNO = jj/jj
        !Udir = direccion/modulov(direccion)
        Udir = direccion/sqrt(dot_product(direccion,direccion))
        do i = 1, num_e
            baricentro = t_baric(:,e_t(1,i))
            !onda = polarizacion*ctte*exp((-jj)*(w*sqrt(mu*epsi))*producto_escalarR(Udir,baricentro))
            ondaE = polarizacion*ctte*exp((-jj)*(w*sqrt(mu*epsi))*dot_product(Udir,baricentro))
            ondaH = cross_productC(Udir*UNO, ondaE)/sqrt(mu/epsi)
            !b(i) = producto_escalar(onda,baricentro-p_coord(:,e_po(1,i)))
            normalPlus = t_normal(:,e_t(1,i))/norm2(t_normal(:,e_t(1,i)))
            res(i) = alphaCFIE*(1./(jj*w*mu))*sum(ondaE*( baricentro-p_coord(:,e_po(1,i)) ) ) + (jj/( getkappa() ) )*(1 - &
                alphaCFIE)*sum(cross_productC(UNO*normalPlus,ondaH)*( baricentro-p_coord(:,e_po(1,i)) ) )

            baricentro = t_baric(:,e_t(2,i))
            !onda = polarizacion*ctte*exp((-jj)*(w*sqrt(mu*epsi))*producto_escalarR(Udir,baricentro))
            ondaE = polarizacion*ctte*exp((-jj)*(w*sqrt(mu*epsi))*dot_product(Udir,baricentro))
            ondaH = cross_productC(Udir*UNO, ondaE)/sqrt(mu/epsi)
            !b(i) = b(i) + producto_escalar(onda,p_coord(:,e_po(2,i))-baricentro)
            normalMinus = t_normal(:,e_t(2,i))/norm2(t_normal(:,e_t(2,i)))
            res(i) = res(i) + alphaCFIE*(1./(jj*w*mu))*sum(ondaE*(p_coord(:,e_po(2,i))-baricentro) ) + (jj/( getkappa() ) &
                )*(1 - alphaCFIE)*sum(cross_productC(UNO*normalMinus,ondaH)*( p_coord(:,e_po(2,i)) - baricentro ) )
            res(i) = res(i)*0.5*e_long(i)
        end do

    !end subroutine
    end function generar_onda_old

        function generar_onda(direccion, ctte, polarizacion) result (res)
        !dummy
        real (kind = dp), dimension(3), intent(in) :: direccion
        complex (kind = dp), dimension(3), intent(in) :: polarizacion
        complex (kind = dp), intent(in) :: ctte
        complex (kind = dp), dimension(:), allocatable :: res!, revisar
        !complex (kind = dp), dimension(num_e), intent(out) :: b
        !
        !local
        integer (kind = il) :: m, itest, i_m, tri_m, GL_numberpoints
        real (kind = dp), dimension(3) :: baricentro, Udir, Vm1, Vm2, Vm3, Vom, nm_uni, rtest, rhom
        complex (kind = dp), dimension(3) :: Ei, Hi
        real (kind = dp), dimension(:,:), allocatable :: GL_pointsweights_test
        complex (kind = dp) :: UNO
        !
        UNO = jj/jj
        !Udir = direccion/modulov(direccion)
        Udir = direccion/norm2(direccion)
        allocate(res(num_e))
        res(:) = 0._dp

        GL_numberpoints = 7
        allocate(GL_pointsweights_test(GL_numberpoints, 4))

         do m = 1, num_e
            !all basis function
            do i_m = 1, 2
                !trangle + and -
                tri_m = e_t(i_m, m)
                Vm1 = p_coord(:,t_p(1, tri_m))
                Vm2 = p_coord(:,t_p(2, tri_m))
                Vm3 = p_coord(:,t_p(3, tri_m))
                Vom = p_coord(:,e_po(i_m, m))

                nm_uni = t_normal(:, tri_m)/norm2(t_normal(:, tri_m)) 

                GL_pointsweights_test = getGLpointsWeights(Vm1, Vm2, Vm3)
            
                do itest = 1, GL_numberpoints
                    rtest = GL_pointsweights_test(itest, 1:3)
                    rhom = rtest - Vom
                    if (i_m == 2) then
                        rhom = rhom*(-1._dp)
                    end if
                    Ei = polarizacion*ctte*exp((-jj)*(w*sqrt(mu*epsi))*dot_product(Udir,rtest))
                    Hi = cross_productC(Udir*UNO, Ei)/sqrt(mu/epsi)
                    res(m) = res(m) + (e_long(m)/2._dp)*( alphaCFIE*sum(rhom*Ei)/(jj*getw()*mu) + (jj/( getkappa()) )*(1 &
                        - alphaCFIE)*sum(rhom*cross_productC(nm_uni*UNO,Hi)) )*GL_pointsweights_test(itest, 4)
 
                end do

            end do            

        end do


    end function generar_onda


end program main