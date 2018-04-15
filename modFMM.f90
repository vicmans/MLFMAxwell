module modFMM
    use modGlobalParam
    use modGlobalMetod
    use modMoM
    implicit none

    private
    !Variables globales
    public ::  inicializar_FMM, MatxVec_FMM, Zesparcida_FMM, ia_FMM, ja_FMM, n_esparcida_FMM, destruir_FMM

    !Variables puntero de geometria del dispersor
    real ( kind = dp ),  pointer, dimension (:,:) :: p_coord
    integer ( kind = il ), pointer, dimension (:,:) :: t_p, e_p
    integer ( kind = il ), pointer :: num_p, num_t
    real ( kind = dp ), pointer, dimension(:,:) :: e_centro
    integer (kind = il), pointer :: num_e
    real ( kind = dp ), pointer, dimension (:,:) :: t_baric
    real ( kind = dp ), pointer, dimension (:,:,:) :: t_baric_sub
    real ( kind = dp ), pointer, dimension(:) :: t_area, e_long
    integer ( kind = il ), pointer, dimension (:,:) :: e_t, e_po
    real ( kind= dp ), pointer, dimension(:,:) :: t_normal
    !
    !Paramtros de entrada del metodo
    integer (kind = il), parameter :: dig_multipolos = 6 !Precision del metodo (usado para el calculo del numero de multipolos)
    real (kind = dp) :: criterio_ladoW !criterio de tamaño de lados de cubos (usualmente lambda/2)
    logical :: onlyEFIE
    !

    integer (kind = ild) :: n_esparcida_FMM
    integer (kind = ild), allocatable, dimension(:) :: ia_FMM
    integer (kind = il), allocatable, dimension(:) :: ja_FMM
    complex (kind = dp_near), allocatable, dimension(:) :: Zesparcida_FMM

    real (kind = dp) :: lado_cubo !longitud del lado de los cubos
    integer (kind = il) :: n_cubos_eje !cantidad de cubos a lo largo del eje mas largo del dispersor
    integer (kind = il) :: n_cubos !numero de cubos llenos totales
    integer (kind = il) :: n_multipolos !cantidad de multipolos que arroja la formula de truncamiento
    integer (kind = il) :: n_kappas !numero de muestras sobre la esfera unitaria
    integer (kind = il) :: n_muestras_theta, n_muestras_phi !numero de muestras en theta y en phi sobre la esfera unitaria


    complex (kind = dp), allocatable, dimension(:,:) :: functransf  !Matriz de conjunto de todas las funciones de transferencia muestreadas.
																                                                !Cada columna es una funcion de transferencia entre un par de cubos,
																                                                !muestreada en todos los kappas [:,:] = [kappas,indice]

    integer (kind = il) :: n_functransf

    real (kind = dp), allocatable, dimension(:,:) :: kappas, ang_kappas !Matriz de conjunto de muestras (kappas sobre la esfera unitaria.
																                                !Cada columna de kappas, son las coordenadas (x,y,z) de esa muestra
																	                            !Cada columna de ang_kappas, son los angulos theta y phi de esa muestra
																	                            !las muestras (columnas de kappas y ang_kappas) son paralelas (correspondientes)
    real (kind = dp), allocatable, dimension(:,:) :: atheta, aphi   !Matrices de vectores unitarios atheta y aphi (coord esfericas)
																	                            !muestreados sobre la esfera unitaria. Cada columna es un vector
																	                            !(x,y,z) correspondiente a una muestra paralela al vector kappas

    real (kind = dp), allocatable, dimension(:) :: pesos_integ  !Vector de pesos de integracion sobre la esfera unitaria
															                                                !segun la regla de cubatura, para todas las muestras sobre
															                                                !la esfera unitaria. (paralelo al vector de muestras kappas)

    complex (kind = dp), allocatable, dimension(:,:,:) :: funcRAD, funcREC  !Matrices de funciones de radiacion o recepcion funcRAD y funcREC respectivamente
																		                                                !Es una matriz tridimensional (una lista de matrices cuadradas), en donde cada
																		                                                !matriz cuadrada corresponde a un cubo (paralelo al vector cubos de tipo cubo)
																		                                                !y en dicha matriz tiene n_kappas filas y dos columnas (componentes en atheta
																		                                                !y aphi) [:,:,:] = [kappas,componentes_theta-phi,cubo)

    complex (kind = dp), allocatable, dimension(:,:,:) :: campoC, campoB    !Matrices de campoC (de agregacion) y campoB (de agregacion por transferencias)
																		                                                !Es una matriz tridimensional (una lista de matrices cuadradas), en donde cada
																		                                                !matriz cuadrada corresponde a un cubo (paralelo al vector cubos de tipo cubo)
																		                                                !y en dicha matriz tiene n_kappas filas y dos columnas (componentes en atheta
																		                                                !y aphi) [:,:,:] = [kappas,componentes_theta-phi,cubo)


    !pruebas especiales
    integer (kind = il) :: Lverd
    integer (kind = il), allocatable, dimension (:) :: Lverdadero
    integer (kind = il) :: Lforzado
    logical :: skipkD
    !


    !****** tyoe cubo ******
    !* Estructura de datos con las propiedades que debe tener
    !* cada cubo lleno en la geometria del dispersor
    type :: cubo
        real (kind = dp), dimension(3) :: centro !vector (x,y,z) del centro del cubo
        integer (kind = il) :: n_propias !cantidad de funciones bases dentro del cubo
        integer (kind = il) :: n_cercanas !cantidad de funciones bases dentro del cubo más las que estan en cubos cercanos
        integer (kind = il) :: n_cubos_lejanos !cantidad de cubos lejanos
        integer (kind = il), allocatable, dimension(:) :: bases_propias, bases_cercanas, cubos_lejanos, ind_functransf
													                                                                                                !vectores de lista: bases_propias->  Lista de indices de funciones bases dentro del cubo
													                                                                                                !					bases_cercanas-> Lista de indices de funciones bases dentro del cubo
													                                                                                                !									 y dentro de cubos cercanos
													                                                                                                !					cubos_lejanos->  Lista de indices de cubos lejanos
													                                                                                                !					ind_functransf-> Lista de indices de las funciones de transferencia. La compenente j-esima
													                                                                                                !									 de esta lista corresponde con la funcion de transferencia del cubo en
													                                                                                                !									 cuestion con el cubo lejano j-esimo del vector cubos_lejanos. Y dicho indice
													                                                                                                !									 de la componente j-esima de ind_functransf se refiere una columna de functransf

        complex (kind = dp), allocatable, dimension(:,:) :: Zcerc !Matriz formato esparcido donde filas corresponden a bases_propias y columnas corresponden a bases_cercanas

    end type
    type(cubo), allocatable, dimension(:) :: cubos !vector de cubos llenos de la geometria del dispersor


contains



    subroutine destruir_FMM()
        if (allocated(cubos)) then
            deallocate(cubos)
        end if

        if (allocated(functransf)) then
            deallocate(functransf)
        end if

        if (allocated(kappas)) then
            deallocate(kappas)
        end if

        if (allocated(ang_kappas)) then
            deallocate(ang_kappas)
        end if

        if (allocated(atheta)) then
            deallocate(atheta)
        end if

        if (allocated(aphi)) then
            deallocate(aphi)
        end if

        if (allocated(pesos_integ)) then
            deallocate(pesos_integ)
        end if

        if (allocated(funcRAD)) then
            deallocate(funcRAD)
        end if

        if (allocated(funcREC)) then
            deallocate(funcREC)
        end if

        if (allocated(campoC)) then
            deallocate(campoC)
        end if

        if (allocated(campoB)) then
            deallocate(campoB)
        end if

        if (allocated(Lverdadero)) then
            deallocate(Lverdadero)
        end if

    end subroutine




    !***************************** inicializar_FMM ************************************
    !* Subrutina de inicializacion de FMM. Es aqui donde se le asignan valores a
    !* los punteros de las variables de geometria del dispersor. Luego de obtenidos
    !" los datos de la geometria en el modulo, se comienzan a hacer todos los procesos
    !* para preparar todo el metodo FMM, y que solo haga falta comenzar con el proceso
    !* iterativo de solucion. Es decir, se llenan las funciones de radiacion, recepcion
    !* y transferencia, y todos los datos necesarios por la funcion del producto matriz
    !* por vector que requerira el metodo iterativo
    !***********************************************************************************

    subroutine inicializar_FMM (arg_e_centro, arg_num_e, arg_e_t, arg_e_po, arg_e_long, arg_t_area, arg_t_baric, &
        arg_t_baric_sub, arg_p_coord, arg_t_p, arg_e_p, arg_num_p, arg_num_t, arg_t_normal)
        !dummy
        real ( kind = dp ),  dimension(:,:), target :: arg_e_centro
        integer (kind = il), target :: arg_num_e
        real ( kind = dp ),  dimension (:,:), intent(in), target :: arg_t_baric
        real ( kind = dp ),  dimension (:,:,:), intent(in), target :: arg_t_baric_sub
        real ( kind = dp ),  dimension(:), intent(in), target :: arg_t_area, arg_e_long
        integer ( kind = il ),  dimension (:,:), intent(in), target :: arg_e_t, arg_e_po
        real ( kind = dp ),  target, dimension (:,:) :: arg_p_coord
        integer ( kind = il ), target, dimension (:,:) :: arg_t_p, arg_e_p
        integer ( kind = il ), target :: arg_num_p, arg_num_t
        real ( kind= dp ), target, dimension(:,:) :: arg_t_normal
        !
        !local
        integer (kind = il), dimension(6) :: limites
        real (kind = dp) :: ladomayor, deltaX, deltaY, deltaZ, maxX, minX, maxY, minY, maxZ, minZ
        real (kind = dp), dimension (3) :: roo
        real (kind = dp) :: tol

        !Asignar los punteros de variables de geometria a este modulo
        num_e => arg_num_e
        e_centro => arg_e_centro
        e_t => arg_e_t
        e_po => arg_e_po
        e_long => arg_e_long
        t_area => arg_t_area
        t_baric => arg_t_baric
        t_baric_sub => arg_t_baric_sub
        p_coord => arg_p_coord
        t_p => arg_t_p
        e_p => arg_e_p
        num_p => arg_num_p
        num_t => arg_num_t
        t_normal => arg_t_normal
        !
onlyEFIE = (abs(alphaCFIE - 1)<=0.001_dp )
        ! Inicializar MoM para el posterior calculo de la Zcercana
        call inicializar_MoM(p_coord, t_baric, t_baric_sub, t_area, e_long, t_p, e_p, e_t, e_po, num_e, num_p, num_t, t_normal)
        !

        criterio_ladoW = conf%criterio_ladomax!getlambda()/2.

        limites = limites_dispersor(e_centro, num_e)
        maxX = e_centro(1,limites(1))
        minX = e_centro(1,limites(2))
        maxY = e_centro(2,limites(3))
        minY = e_centro(2,limites(4))
        maxZ = e_centro(3,limites(5))
        minZ = e_centro(3,limites(6))
        deltaX = maxX-minX
        deltaY = maxY-minY
        deltaZ = maxZ-minZ
        roo(1) = 0.5*(maxX+minX)
        roo(2) = 0.5*(maxY+minY)
        roo(3) = 0.5*(maxZ+minZ)

        if (  (deltaX >= deltaY) .and. (deltaX >= deltaZ) ) then
            ladomayor = deltaX
        else if (  (deltaY >= deltaX) .and. (deltaY >= deltaZ) ) then
            ladomayor = deltaY
        else
            ladomayor = deltaZ
        end if
        !Tolerancia para las funciones bases en los bordes
        tol = getlambda()/10000.
        ladomayor = ladomayor + tol
        tol = tol/2.
        print*, 'El cubo que contiene al dispersor mide ' // realtostr(ladomayor) // ' metros'
        print*, 'Calculando geometria'

        call geomFMM(ladomayor,maxX+tol,maxY+tol,maxZ+tol)

        print*, 'La longitud del lado de cada cubo resulta: ' // realtostr(lado_cubo) // ' metros'
        print*, 'Hecho'

        print*, 'Llenando muestras sobre la esfera unitaria'
        call hallar_n_multipolos_muestras()

        call llenar_kappas()

        print*, 'Se tienen ' // inttostr(n_multipolos) // ' multipolos y ' // inttostr(n_kappas) // ' muestras sobre la &
        esfera unitaria'
        print*, 'Hecho'
        print*, 'Calculando funciones de radiacion'

        allocate(funcRAD(n_kappas,2,num_e))
        allocate(funcREC(n_kappas,2,num_e))
        allocate(campoC(n_kappas,2,n_cubos))
        allocate(campoB(n_kappas,2,n_cubos))
        call llenar_funcRAD(.false.)
        funcREC = conjg(funcRAD)
        print*, 'Hecho'
        print*, 'Llenando funciones de transferencia'
        call llenar_functransf()

        print*, 'Hecho'
        print*, 'Llenando Z cercana'

        call calc_Zcercana()
        print*, 'Hecho'

    end subroutine

    !********************************* geomFMM **************************************
    !* Subrutina que llena el vector cubos de tipo cubo, con los datos del centro del
    !* cubo, la cantidad de funciones base dentro de el mismo, y el vector bases_propias
    !* que contiene los indices de las funciones base dentro del cubo. El resto de los
    !* datos de ese vector de cubos (como las funciones base cercanas) seran llenados
    !* posteriormente en otra subrutina al calcular las funciones de transferencia
    !********************************************************************************

    subroutine geomFMM(lCuboMayor, esquinaX, esquinaY, esquinaZ)
        !dummy
        real (kind = dp), intent(in) :: lCuboMayor, esquinaX, esquinaY, esquinaZ
        !
        !local
        integer (kind = il) :: n, i, j ,k, n_cubos_eje, cont, l
        real (kind = dp), dimension(3) :: esquina, centro
        real (kind = dp) :: deltaX, deltaY, deltaZ
        integer (kind = il),allocatable, dimension(:) :: vBusqueda, vEncontrados, vRestantes
        integer (kind = il) :: n_busq, n_encon, n_rest
        logical :: encontro
        type(cubo), allocatable, dimension(:) :: auxcubos, aux

        integer (kind = il) :: tamano_aux	 		!tamaño actual del vector auxcubos

        integer (kind = il) :: criterio_incremento  !criterio de incremento del vector auxcubos
                                                !cuando se necesite de un tamaño mayor al
                                                !que tiene (todo esto para evitar que auxcubos
                                                !se inicialice con un tamaño igual a n_cubos_eje**3
                                                !el cual puede ser muy grande y ocupar memoria innecesariamente
                !

        criterio_incremento = 100
        !se busca una cantidad de cubos que quepan en el eje mayor y que sean
        !de lado menor o igual a criterio_ladoW
        n_cubos_eje = ceiling(lCuboMayor/(criterio_ladoW*getlambda()))
        lado_cubo = lCuboMayor/n_cubos_eje
        !
        !se creara entonces un cubo de lado lCuboMayor cuya esquina superior corresponda
        !con la esquina superior del dispersor


        allocate(auxcubos(criterio_incremento))
        tamano_aux = criterio_incremento


        !Inicializar vBusqueda para que la primera vez busque todas las funciones base
        allocate(vBusqueda(num_e))
        n_busq = num_e
        do n = 1,n_busq
            vBusqueda(n) = n
        end do
        !Inicializar los deltas
        deltaX = 0.
        deltaY = 0.
        deltaZ = 0.
        cont = 0
        !    pos = 1

        do i = 1,n_cubos_eje !Se mueve a lo largo del eje Z
            do j = 1,n_cubos_eje !Se mueve a lo largo del eje Y
                do k = 1,n_cubos_eje !Se mueve a lo largo del eje X
                    !Coordenada X
                    esquina(1) = esquinaX - deltaX
                    !Coordenada Y
                    esquina(2) = esquinaY - deltaY
                    !Coordenada Z
                    esquina(3) = esquinaZ - deltaZ


                    centro(1) = esquina(1) - lado_cubo/2.
                    centro(2) = esquina(2) - lado_cubo/2.
                    centro(3) = esquina(3) - lado_cubo/2.

                    call buscarFuncBases(centro, lado_cubo, vBusqueda,n_busq,vEncontrados,n_encon,vRestantes,n_rest, encontro)
                    if (encontro) then

                        cont = cont + 1

                        if (cont > tamano_aux) then
                            !hay que redimensionar
                            allocate(aux(tamano_aux))
                            aux = auxcubos
                            deallocate(auxcubos)
                            tamano_aux = tamano_aux
                            allocate(auxcubos(tamano_aux + criterio_incremento))
                            auxcubos(1:tamano_aux) = aux(:)
                            tamano_aux = tamano_aux + criterio_incremento
                            deallocate(aux)
                        end if

                        !Guardo el Centro
                        auxcubos(cont)%centro(:) = centro!auxCentros(:,cont) = centro

                        n_busq = n_rest
                        deallocate(vBusqueda)
                        allocate(vBusqueda(n_busq))
                        vBusqueda = vRestantes
                        allocate(auxcubos(cont)%bases_propias(n_encon))
                        auxcubos(cont)%bases_propias = vEncontrados
                        auxcubos(cont)%n_propias = n_encon !funcBases(pos:(pos+n_encon-1)) = vEncontrados

                        !auxPos_Num(1,cont) = pos
                        !auxPos_Num(2,cont) = n_encon
                        !pos = pos + n_encon
                    end if
                    deltaX = deltax + lado_cubo
                end do
                deltaX = 0.
                deltaY = deltaY + lado_cubo
            end do
            deltaY = 0.
            deltaZ = deltaZ + lado_cubo
        end do


        print*, 'Se tienen ' // inttostr(cont) // ' cubos llenos'

        n_cubos = cont
        allocate(cubos(n_cubos))
        cubos(:) = auxcubos(1:n_cubos)
        deallocate(auxcubos)

    end subroutine

    !************************************** buscarFuncBases *******************************************
    !* Subrutina usada por geomFMM, para buscar la lista de funciones base que estan dentro de un cubo
    !**************************************************************************************************
    
    subroutine buscarFuncBases(centro, lcubo, vBusqueda,n_busq,vEncontrados,n_encon,vRestantes,n_rest, encontro)
        !el vector vBusqueda contiene los indices de las funciones bases contenidas (e_centro)
        !dumy
        real (kind = dp), dimension(3), intent(in) :: centro
        real (kind = dp), intent(in) :: lcubo
        integer (kind = il), intent(in), dimension(:) :: vBusqueda
        integer (kind = il), intent(in) :: n_busq
        integer (kind = il), allocatable, intent(out), dimension(:) :: vEncontrados, vRestantes
        integer (kind = il), intent(out) :: n_encon, n_rest
        logical, intent(out) :: encontro
        !local
        real (kind = dp), dimension(6) :: limites
        real (kind = dp) :: lq
        integer (kind = il) :: i
        integer (kind = il), allocatable, dimension(:) :: aux1, aux2
        integer (kind = il) :: n_aux
        logical :: fuera
        !calcular limites (Xmax,Xmin,Ymax,Ymin,Zmax,Zmin) del cubo
        lq = lcubo/2.
        !Xmax
        limites(1) = centro(1)+lq
        !Xmin
        limites(2) = centro(1)-lq
        !Ymax
        limites(3) = centro(2)+lq
        !Ymin
        limites(4) = centro(2)-lq
        !Zmax
        limites(5) = centro(3)+lq
        !Zmin
        limites(6) = centro(3)-lq

        allocate(aux1(n_busq))
        allocate(aux2(n_busq))
        n_encon= 0
        n_rest = 0
        do i=1,n_busq
            fuera = .true.
            if (e_centro(1,vBusqueda(i)) <= limites(1) .and. &
                e_centro(1,vBusqueda(i)) >= limites(2) ) then
                if (e_centro(2,vBusqueda(i)) <= limites(3) .and. &
                    e_centro(2,vBusqueda(i)) >= limites(4) ) then
                    if (e_centro(3,vBusqueda(i)) <= limites(5) .and. &
                        e_centro(3,vBusqueda(i)) >= limites(6) ) then
                        !si esta dentro del cubo
                        !numero de funciones base enonctradas
                        n_encon = n_encon + 1
                        aux1(n_encon) = vBusqueda(i)
                        fuera = .false.
                    end if
                end if
            end if
            if (fuera) then
                !si esta fuera del cubo
                n_rest = n_rest + 1
                aux2(n_rest) = vBusqueda(i)
            end if
        end do
        if (n_encon > 0) then
            encontro = .true.
            allocate(vEncontrados(n_encon))
            allocate(vRestantes(n_rest))
            vEncontrados(:) = aux1(1:n_encon)
            vRestantes(:) = aux2(1:n_rest)
        else
            encontro = .false.
            allocate(vRestantes(n_rest))
            vRestantes(:) = aux2(1:n_rest)
        end if
        deallocate(aux1)
        deallocate(aux2)

    end subroutine


    !******************* hallar_n_multipolos_muestras ************************
    !* Establece el numero de multipolos estandar que debe ser usado para el
    !* truncamiento de la serie infinita de la expansion multipolar, y
    !* a partir de ese numero, se calcula el numero de muestras en theta y phi
    !* y totales sobre la esfera unitaria
    !*************************************************************************
    
    subroutine hallar_n_multipolos_muestras()
        !local
        real (kind = dp) :: d
        !
        d = sqrt(3.)*lado_cubo
        n_multipolos = floor(((2.*pi/getlambda())*d) + dig_multipolos*(((2.*pi/getlambda())*d)**(1./3.)))

        n_muestras_theta = n_multipolos
        n_muestras_phi = 2*n_multipolos
        n_kappas = n_muestras_phi*n_muestras_theta

    end subroutine

    !********************** llenar_kappas *************************
    !* Llena el vector kappas de muestras sobre la esfera unitaria
    !* y los vectores atheta y aphi
    !**************************************************************
    
    subroutine llenar_kappas_old()

        !local
        integer (kind = il) :: i,j,k, cont
        real (kind = dp) :: deltatheta, deltaphi, theta, phi,  pesoPhi
        !

        allocate(kappas(3,n_kappas))
        allocate(atheta(3,n_kappas))
        allocate(aphi(3,n_kappas))
        allocate(ang_kappas(2,n_kappas))
        allocate(pesos_integ(n_kappas))
        cont = 0
        deltaphi = 2.*pi/n_muestras_phi!deltaphi = 2.*pi/(2.*n_multipolos + 1)!deltaphi = 2*pi/(2*n_multipolos + 1) ###
        deltatheta = pi/(n_muestras_theta)!deltatheta = pi/(n_multipolos)
        theta = 0.
        phi = 0.


        do j = 1, n_muestras_phi

            !theta = (deltatheta/2.)*((2*j)-1)
            phi = (deltaphi/2.)*((2*j)-1)
            do k = 1, n_muestras_theta

                theta = (deltatheta/2.)*((2*k)-1)
                !phi = (deltaphi/2.)*((2*k)-1)
                cont = cont + 1

                kappas(1,cont) = sin(theta)*cos(phi)
                kappas(2,cont) = sin(theta)*sin(phi)
                kappas(3,cont) = cos(theta)

                atheta(1,cont) = cos(theta)*cos(phi)
                atheta(2,cont) = cos(theta)*sin(phi)
                atheta(3,cont) = -sin(theta)

                aphi(1,cont) = -sin(phi)
                aphi(2,cont) = cos(phi)
                aphi(3,cont) = 0.

                ang_kappas(1,cont) = theta
                ang_kappas(2,cont) = phi
                pesos_integ(cont) = (cos((k-1)*deltatheta)-cos(k*deltatheta))*deltaphi

            end do
        end do

    end subroutine


    subroutine llenar_kappas()
        !local
        integer (kind = il) :: i,j,k, cont
        real (kind = dp) :: deltatheta, deltaphi, theta, phi
        real (kind = dp), allocatable, dimension(:) :: thetaAbscissas, thetaWeights
        !

        allocate(kappas(3,n_kappas))
        allocate(atheta(3,n_kappas))
        allocate(aphi(3,n_kappas))
        allocate(ang_kappas(2,n_kappas))
        allocate(pesos_integ(n_kappas))

        allocate(thetaAbscissas(n_muestras_theta))
        allocate(thetaWeights(n_muestras_theta))
        call gauleg(n_muestras_theta,thetaAbscissas,thetaWeights)

        cont = 0
        deltaphi = (2.*pi)/(n_muestras_phi)!2.*pi/octree(i)%n_muestras_phi
        !deltatheta = pi/octree(i)%n_muestras_theta

        do j = 1, n_muestras_phi

            phi = (j - 1)*deltaphi!(deltaphi/2.)*((2*j)-1)

            do k = 1,n_muestras_theta

                theta = acos(thetaAbscissas(k))!(deltatheta/2.)*((2*k)-1)

                cont = cont + 1

                kappas(1,cont) = sin(theta)*cos(phi)
                kappas(2,cont) = sin(theta)*sin(phi)
                kappas(3,cont) = cos(theta)

                atheta(1,cont) = cos(theta)*cos(phi)
                atheta(2,cont) = cos(theta)*sin(phi)
                atheta(3,cont) = -sin(theta)

                aphi(1,cont) = -sin(phi)
                aphi(2,cont) = cos(phi)
                aphi(3,cont) = 0.

                ang_kappas(1,cont) = theta
                ang_kappas(2,cont) = phi
                pesos_integ(cont) = thetaWeights(k)*deltaphi !(cos((k-1)*deltatheta)-cos(k*deltatheta))*deltaphi
            end do
        end do

        deallocate(thetaAbscissas)
        deallocate(thetaWeights)

    end subroutine


    !*********************** calcFuncRad ***************************
    !* Calcula las funciones de radiacion y llena el vector funcRAD
    !* Esta la opcion de hacerlo con maxima precision o no. Cuando
    !* se usa maxima precision, se subdividen los parches en nueve
    !* sub baricentros
    !***************************************************************
    
    subroutine llenar_funcRAD_old(maxPrecision)
        !dummy
        logical, optional, intent(in) :: maxPrecision
        !local
        integer (kind = il) :: i, j, ii,b!, nbases, pos

        real (kind = dp), dimension(3) :: baricentroplus, vopuestoplus, rhoplus, baricentrominus, vopuestominus, rhominus, &
        vkappa, centro
        real (kind = dp) :: areaplus,areaminus,longitud, ka
        complex (kind = dp), dimension(3) :: vRadplus, vRadminus, vRadtotal
        logical :: mP
        !


        ka = (2.*pi)/getlambda()
        if (present(maxPrecision)) then
            mP = maxPrecision
        else
            mP = .false.
        end if


        do i = 1,n_cubos
            centro(:) = cubos(i)%centro(:)

            do b = 1,cubos(i)%n_propias
                if (mP) then
                    !triangulo mas
                    vopuestoplus = p_coord(:,e_po(1,cubos(i)%bases_propias(b)))
                    areaplus = t_area(e_t(1, cubos(i)%bases_propias(b)))
                    longitud = e_long(cubos(i)%bases_propias(b))
                    !triangulo menos
                    vopuestominus = p_coord(:,e_po(2,cubos(i)%bases_propias(b)))
                    areaminus = t_area(e_t(2, cubos(i)%bases_propias(b)))
                    !vector radiacion
                    do j=1,n_kappas
                        vkappa = kappas(:,j)
                        vRadplus=0.
                        vRadminus=0.
                        do ii = 1, 9
                            baricentroplus = t_baric_sub(:, ii, e_t(1, cubos(i)%bases_propias(b)))
                            rhoplus = baricentroplus - vopuestoplus
                            baricentrominus = t_baric_sub(:, ii, e_t(2, cubos(i)%bases_propias(b)))
                            rhominus = vopuestominus - baricentrominus

                            vRadplus = vRadplus + exp(jj*ka*dot_product(vkappa,baricentroplus-centro))*rhoplus
                            vRadminus = vRadminus + exp(jj*ka*dot_product(vkappa,baricentrominus-centro))*rhominus
                        end do

                        vRadtotal = (vRadplus + vRadminus)*longitud*(1./18.)

                        funcRAD(j,1,cubos(i)%bases_propias(b)) = sum(vRadtotal*atheta(:,j))
                        funcRAD(j,2,cubos(i)%bases_propias(b)) = sum(vRadtotal*aphi(:,j))
                    end do

                else

                    baricentroplus = t_baric(:, e_t(1, cubos(i)%bases_propias(b)))

                    vopuestoplus = p_coord(:,e_po(1,cubos(i)%bases_propias(b)))
                    areaplus = t_area(e_t(1, cubos(i)%bases_propias(b)))
                    longitud = e_long(cubos(i)%bases_propias(b))
                    rhoplus = baricentroplus - vopuestoplus

                    !triangulo menos
                    baricentrominus = t_baric(:, e_t(2, cubos(i)%bases_propias(b)))
                    vopuestominus = p_coord(:,e_po(2,cubos(i)%bases_propias(b)))
                    areaminus = t_area(e_t(2, cubos(i)%bases_propias(b)))
                    rhominus = vopuestominus - baricentrominus
                    !print*, 't-'
                    !vector radiacion
                    do j=1,n_kappas
                        !print*, 'entro kappa'
                        vkappa = kappas(:,j)
                        vRadplus = exp(jj*ka*dot_product(vkappa,baricentroplus-centro))*rhoplus!/(areaplus)
                        vRadminus = exp(jj*ka*dot_product(vkappa,baricentrominus-centro))*rhominus!/(areaminus)
                        vRadtotal = (vRadplus + vRadminus)*longitud*0.5

                        funcRAD(j,1,cubos(i)%bases_propias(b)) = sum(vRadtotal*atheta(:,j))
                        funcRAD(j,2,cubos(i)%bases_propias(b)) = sum(vRadtotal*aphi(:,j))

                    end do

                end if

            end do
        end do

    end subroutine



    subroutine llenar_funcRAD(maxPrecision)
        !dummy
        logical, optional :: maxPrecision
        !
        !local
        real (kind = dp), dimension(3) :: centro
        integer (kind = il) :: nBases , niveles

        integer (kind = il) :: n, k, fb, i_t, tri_i, GL_numberpoints, itest, i, b
        real (kind = dp), dimension(3) :: V1, V2, V3, Vo, rho, vkappa, vnormal, rtest
        real (kind = dp) :: ln, ka, area
        complex (kind = dp), dimension(3) :: vRad, vRadtotal, vRec
        complex (kind = dp) :: UNO
        real (kind = dp), dimension(:,:), allocatable :: GL_pointsweights_test
        !
        UNO = (jj)/(jj)
        ka = (2.*pi)/getlambda()
        GL_numberpoints = 7
       
        allocate(GL_pointsweights_test(GL_numberpoints, 4))


        do i = 1,n_cubos
            centro(:) = cubos(i)%centro(:)

            do b = 1,cubos(i)%n_propias
                fb = cubos(i)%bases_propias(b)
                ln = e_long(fb)

                !vector radiacion
                do k = 1, n_kappas
                    vkappa = kappas(:,k)                
                    vRec(:) = 0._dp
                    vRadtotal(:) = 0._dp

                    !triangle + and -
                    do i_t = 1, 2

                        vRad(:) = 0._dp                    
                        
                        !triangle ind
                        tri_i = e_t(i_t, fb)
                        !area
                        area = t_area(tri_i)
                        !triangle vertex
                        V1 = p_coord(:, t_p(1, tri_i))
                        V2 = p_coord(:, t_p(2, tri_i))
                        V3 = p_coord(:, t_p(3, tri_i))
                        !oposite vertex
                        Vo = p_coord(:, e_po(i_t,fb))
                        
                        
                        GL_pointsweights_test = getGLpointsWeights(V1, V2, V3)
                        
                        do itest = 1, GL_numberpoints

                            rtest = GL_pointsweights_test(itest, 1:3)
                            rho = rtest - Vo
                            if ( i_t == 2 ) then
                                !triangle -
                                rho = (-1._dp)*rho
                            end if
                            rho = rho*0.5_dp*ln
                            vRad = vRad + exp(jj*ka*dot_product(vkappa,rtest-centro))*rho*GL_pointsweights_test(itest, 4)

                        end do

                        vnormal = t_normal(:,tri_i)
                        vnormal = vnormal/norm2(vnormal)

                        vRec = vRec + cross_productC ( vkappa*UNO , cross_productC( conjg(vRad) , vnormal*UNO) )

                        vRadtotal = vRadtotal + vRad
                        
                    end do

                    funcRAD(k,1,fb) = sum(vRadtotal*atheta(:,k))
                    funcRAD(k,2,fb) = sum(vRadtotal*aphi(:,k))

    !funcRAD2(k,:,fb) = vRadtotal - vkappa*sum(vRadtotal*vkappa)
                    if (.not. onlyEFIE) then
                        funcREC(k,1,fb) = alphaCFIE*conjg(funcRAD(k,1,fb)) + (1._dp - alphaCFIE)*sum(vRec*atheta(:,k))
                        funcREC(k,2,fb) = alphaCFIE*conjg(funcRAD(k,2,fb)) + (1._dp - alphaCFIE)*sum(vRec*aphi(:,k))
                    end if
    !print*, '>>>', funcRAD(k, 1, fb)
                end do

            end do
        end do        

        deallocate(GL_pointsweights_test) 

    end subroutine llenar_funcRAD



    !************************* llenar_functransf ***************************
    !* Llena todo lo concerniente a las funciones de transferencia, esto es:
    !* functransf, n_functransf; de la estructura cubos: bases_cercanas,
    !* cubos_lejanos, ind_functransf, n_cercanas y n_cubos_lejanos
    !***********************************************************************
   
    subroutine llenar_functransf ()

        !local
        integer (kind = il) :: i, j, k, cont, clej, contcer, conttotaltrans
        real (kind = dp) :: tol
        real (kind = dp), dimension(3) :: rab, resta
        real (kind = dp), allocatable, dimension(:,:) :: rabTrans, auxR
        complex (kind = dp), allocatable, dimension(:,:) :: auxFuncTransf, auxF
        integer (kind = il), allocatable, dimension(:) ::  auxLejanos, auxIndFuncTransf
        integer (kind = il), allocatable, dimension(:) :: auxfuncCercanas
        integer (kind = il) :: prc, aprc


        integer (kind = il) :: tamano_aux	 		!tamaño actual de los vectores auxFuncTransf, rabTransf

        integer (kind = il) :: criterio_incremento  !criterio de incremento del vector auxcubos
                                                            !cuando se necesite de un tamaño mayor al
                                                            !que tiene (todo esto para evitar que auxcubos
                                                            !se inicialice con un tamaño igual a n_cubos_eje**3
                                                            !el cual puede ser muy grande y ocupar memoria innecesariamente

                !

        n_esparcida_FMM = 0
        criterio_incremento = 100

        aprc = -1
        tol = lambda/100.
        cont = 1
        clej = 0
        conttotaltrans = 0

        tamano_aux = criterio_incremento
        allocate(auxFuncTransf(n_kappas,tamano_aux))
        allocate(rabTrans(3,tamano_aux))

        allocate(auxIndFuncTransf(n_cubos))
        allocate(auxLejanos(n_cubos))
        allocate(auxfuncCercanas(num_e))

        rabTrans(:,:) = 0.
        call setprc(n_cubos)
        do i = 1,n_cubos

            contcer = 0
            clej = 0
            do j = 1,n_cubos
                rab(:) = cubos(i)%centro(:) - cubos(j)%centro(:)
                if ( norm2(rab) >= 2*lado_cubo-tol ) then
                    !Cubos lejanos
                    clej = clej + 1
                    auxLejanos(clej) = j
                    conttotaltrans = conttotaltrans + 1
                    do k = 1,cont
                        resta(:) = rabTrans(:,k) - rab(:)!resta(:) = rabTrans(:,cont) - rab(:)
                        if ((abs(resta(1))<=tol).and.(abs(resta(2))<=tol).and.(abs(resta(3))<=tol)) then
                            !Ya existe la funcion de transferencia
                            auxIndFuncTransf(clej) = k
                            exit
                        else if (k == cont) then
                            !Calcular funcion transferencia
                            rabTrans(:,cont) = rab(:)
                            auxFuncTransf(:,cont) = calc_functransf(rab)
                            auxIndFuncTransf(clej) = cont

                            cont = cont + 1

                            if (cont == tamano_aux) then
                                !redimensionar auxFuncTransf y auxrabTrans
                                allocate(auxF(n_kappas,tamano_aux-1))
                                allocate(auxR(3,tamano_aux-1))

                                auxF(:,1:(tamano_aux-1)) = auxFuncTransf(:,1:(tamano_aux-1))
                                auxR(:,1:(tamano_aux-1)) = rabTrans(:,1:(tamano_aux-1))
                                deallocate(auxFuncTransf)
                                deallocate(rabTrans)
                                allocate(auxFuncTransf(n_kappas,tamano_aux + criterio_incremento))
                                allocate(rabTrans(3,tamano_aux + criterio_incremento))
                                auxFuncTransf(:,1:tamano_aux-1) = auxF(:,1:tamano_aux-1)
                                rabTrans(:,1:tamano_aux-1) = auxR(:,1:tamano_aux-1)
								tamano_aux = tamano_aux + criterio_incremento
                                deallocate(auxF)
                                deallocate(auxR)
                            end if

                        end if
                    end do
                else
                    !cubos cercanos
                    auxfuncCercanas((contcer+1):(contcer+cubos(j)%n_propias)) = cubos(j)%bases_propias(:)
                    contcer = contcer + cubos(j)%n_propias

                end if

            end do
            call updateprc(i)
            !prc = i*100/n_cubos
            !if ((prc) /= (aprc)) then
            !    aprc = prc
            !    print*, inttostr(prc) // ' %'
            !end if

            !Guarda las propiedades del cubo i
            allocate(cubos(i)%bases_cercanas(contcer))
            cubos(i)%bases_cercanas(:) = auxfuncCercanas(1:contcer)
            cubos(i)%n_cercanas = contcer
            n_esparcida_FMM = n_esparcida_FMM + cubos(i)%n_propias*(cubos(i)%n_cercanas)
            allocate(cubos(i)%cubos_lejanos(clej))
            cubos(i)%cubos_lejanos(:) = auxLejanos(1:clej)
            cubos(i)%n_cubos_lejanos = clej
            allocate(cubos(i)%ind_functransf(clej))
            cubos(i)%ind_functransf(:) = auxIndFuncTransf(1:clej)

        end do

        allocate(functransf(n_kappas,(cont-1)))

        functransf(:,:) = auxFuncTransf(:,1:(cont-1))
        n_functransf = cont - 1
        print*, 'Se tienen ' // inttostr(n_functransf) // ' funciones de transferencia unicas'
        print*, 'Numero de funciones de transferencia totales: ' // inttostr(conttotaltrans)

        deallocate(auxIndFuncTransf)
        deallocate(auxFuncTransf)
        deallocate(rabTrans)
        deallocate(auxLejanos)
        deallocate(auxfuncCercanas)

    end subroutine


    !******************** calc_functransf ****************************
    !* Calculo de una funcion de transferencia (Serie de la expansion
    !* multipolar truncada en n_multipolos multipolos) y muestreada
    !* en todos los kappas de la esfera unitaria
    !*****************************************************************
    
    function calc_functransf(rab) result(v)
        !dumy
        real (kind = dp), dimension(3), intent(in) :: rab
        complex (kind =8), dimension(n_kappas) :: v
        !
        !local
        integer (kind = il) :: i, j, Ltrue, kDint
        real (kind = dp), dimension(3) :: vkappa
        real (kind = dp) :: kD, modul
        complex(kind = dp) :: valor
        !
        Ltrue = n_multipolos
        modul = norm2(rab)
        kD = (2.*pi/getlambda())*modul
        kDint = int(kD)
kDint = ceiling(kD)
        if (kDint < Ltrue) then
            Ltrue = kDint! - 1
        end if

        do i=1,n_kappas
            vkappa = kappas(:,i)
            valor = 0.
            do j=0,Ltrue
                valor = valor + (((-jj)**(j+1))*((2.*j)+1)*hankel(j,kD)*legendre(j,dot_product(vkappa,rab/modul)))
            end do
            valor = valor*(2.*pi/getlambda())/(4.*pi)!/**/
            v(i) = valor

        end do

    end function


    !****************************** calc_Zcercana **************************
    !" Llena la matriz Zcerc de cada cubo en formato esparcido para ahorro
    !* memoria y de computo
    !***********************************************************************
   
    

    subroutine calc_Zcercana()
        !local
        integer (kind = il) :: i,j,c,bp,bc,nc,m,n, prc, prca
        integer (kind = ild) :: fila, cont
        integer (kind = ild) :: total_matrix_elements
        real (kind = dp) :: near_percent
        !

        prca = -1
        fila = 1
        cont = 0

        allocate(Zesparcida_FMM(n_esparcida_FMM))
        allocate(ja_FMM(n_esparcida_FMM))
        allocate(ia_FMM(num_e + 1))

total_matrix_elements = num_e**2
near_percent = 100*(n_esparcida_FMM/(1._16*total_matrix_elements))
print*, 'Un ', near_percent, '% de las interacciones son cercanas'
!Zesparcida_FMM(:) = 0.
!ja_FMM(:) = 0
!ia_FMM(:) = 0
        call setprc(num_e)
        do m = 1, num_e
            Busq: do c = 1, n_cubos
                do bp = 1, Cubos(c)%n_propias
                    if (Cubos(c)%bases_propias(bp) == m) then
                        nc = Cubos(c)%n_cercanas
                        ia_FMM(m) = fila
                        do j =1, nc
                            cont = cont + 1
                            n = Cubos(c)%bases_cercanas(j)
                            Zesparcida_FMM(cont) = Zmn_MoM(m,n)
                            ja_FMM(cont) = n
                            fila = fila + 1
                        end do
                        exit Busq
                    end if
                end do
            end do Busq
            call updateprc(m)
            !prc = 100*m/num_e
            !if (prc/10 /= prca/10) then
            !    print*, inttostr(prc) // ' %'
            !    prca = prc
            !end if
        end do
        ia_FMM(num_e + 1) = fila      

!do c = 1, n_cubos
!    print*, '°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°'                    
!    print*, 'cubo: ', c, ' ' , 'n_propias: ', Cubos(c)%n_propias, ' n_cercanas: ', Cubos(c)%n_cercanas
!end do


print*, 'La Znear tiene ', n_esparcida_FMM, ' elementos'

!print*, 'ia:########################################'        
!print*, ia_FMM
!print*, 'ja:########################################'        
!print*, ja_FMM
!print*, 'Zesparcida:########################################'        
!print*, Zesparcida_FMM
    end subroutine calc_Zcercana


    subroutine calc_Zcercana_old()
        !local
        integer (kind = il) :: i,j,nprop,ncerc,c,ni,nj,posi,posj, prc, prca
		prca = -1
        call setprc(n_cubos)
        do c = 1,n_cubos

            allocate(cubos(c)%Zcerc(cubos(c)%n_propias,cubos(c)%n_cercanas))


            do i = 1,cubos(c)%n_propias

                do j = 1,cubos(c)%n_cercanas
                    nprop = cubos(c)%bases_propias(i)
                    ncerc = cubos(c)%bases_cercanas(j)
                    cubos(c)%Zcerc(i,j) = Zmn_MoM(nprop,ncerc)

                end do

            end do
            call updateprc(c)
            !prc = 100*c/n_cubos
            !if (prc/10 /= prca/10) then
            !	prca = prc
            !	print*, inttostr(prc) // ' %'
            !end if
        end do

    end subroutine

    function MatxVec_FMM(vecX) result(res)
    	!dummy
        complex (kind = dp), intent(in), dimension(:) :: vecX
        complex (kind = dp), dimension(size(vecX,1)) :: res
		!

        !local
        integer (kind = il) :: c,b,pos,posi,posj, nbases, indi,indj,m,n, cl, posL, cubL, cont, nf
        complex(kind = dp) :: intlej
		!

        !Agregacion
        campoC(:,:,:) = 0.
        res(:) = 0.
        campoB(:,:,:) = 0.

        do c = 1, n_cubos
            do b = 1,cubos(c)%n_propias
                campoC(:,:,c) = campoC(:,:,c) + funcRad(:,:,cubos(c)%bases_propias(b))*vecX(cubos(c)%bases_propias(b))
            end do
        end do
        !Fin Agregacion

        !Calculo campo B
        do c = 1, n_cubos
            do cl = 1, cubos(c)%n_cubos_lejanos
                campoB(:,1,c) = campoB(:,1,c) + functransf(:,cubos(c)%ind_functransf(cl))*campoC(:,1,cubos(c)%cubos_lejanos(cl))
                campoB(:,2,c) = campoB(:,2,c) + functransf(:,cubos(c)%ind_functransf(cl))*campoC(:,2,cubos(c)%cubos_lejanos(cl))
            end do
        end do
        !Fin campo B

        !MatxVec
        do c = 1, n_cubos
            do m = 1,cubos(c)%n_propias
                indi = cubos(c)%bases_propias(m)
                if (.not. onlyEFIE) then
                    res(indi) = res(indi) + sum(pesos_integ(:)*sum(campoB(:,:,c)*funcREC(:,:,indi),DIM = 2))
                else
                    res(indi) = res(indi) + sum(pesos_integ(:)*sum(campoB(:,:,c)*conjg(funcRAD(:,:,indi)),DIM = 2))
                end if
                !res(indi) = res(indi)*(jj*w*mu/(4.*pi))
                res(indi) = res(indi)*(1/(4.*pi))

                !do n = 1,cubos(c)%n_cercanas
                !    indj = cubos(c)%bases_cercanas(n)
                !    res(indi) = res(indi) + (  vecX(indj)*cubos(c)%Zcerc(m,n)  )
                !end do

            end do

        end do

        !MxV con Zcercana
        cont = 0
        do m = 1, num_e
            do nf = 1, ia_FMM(m+1)-ia_FMM(m)
                cont = cont + 1
                res(m) = res(m) +  vecX(ja_FMM(cont))*Zesparcida_FMM(cont)
            end do
        end do        



    end function

end module modFMM
