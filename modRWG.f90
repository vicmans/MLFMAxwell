module modRWG

    use modGlobalParam
    implicit none

contains

    subroutine llenar_lados (t, num_t, l_todos)

        integer ( kind = il ), dimension(:,:), intent(in) :: t
        integer ( kind = il ), intent(in) :: num_t
        integer ( kind = il ), dimension(:,:), intent(out) :: l_todos

        integer (kind = il) it



        do it = 1, num_t

            l_todos(1,3*it-2)=t(1,it)
            l_todos(2,3*it-2)=t(2,it)
            l_todos(1,3*it-1)=t(1,it)
            l_todos(2,3*it-1)=t(3,it)
            l_todos(1,3*it)=t(2,it)
            l_todos(2,3*it)=t(3,it)

        end do

    end subroutine llenar_lados


    function Buscar_L_Comun(indice, ntodos, l_todos) result (salida)

        integer ( kind = il ) :: salida
        integer ( kind = il ), intent(in) :: ntodos, indice
        integer ( kind = il ), dimension (:,:), intent(in) :: l_todos
        integer i


        salida = -1

        do i = indice+1, ntodos

            if ((  (l_todos(1,i) == l_todos(1,indice)) .and. &
                (l_todos(2,i) == l_todos(2,indice)) ) .or. &
                (     (l_todos(1,i) == l_todos(2,indice)) .and.  &
                (l_todos(2,i) == l_todos(1,indice)) )) then

                salida = i
                exit

            end if

        end do

    end function Buscar_L_Comun


    subroutine llenar_todo(l_todos, ntodos, e_p, e_t, t, e_po, num_e)

        integer (kind = il ), dimension (:,:), intent(inout) :: l_todos
        integer (kind = il ), dimension (:,:), intent(in) :: t
        integer (kind = il ), intent (in) :: ntodos
        integer (kind = il ), allocatable, dimension (:,:), intent (out) :: e_p
        integer (kind = il ), allocatable, dimension (:,:), intent (out) :: e_t
        integer (kind = il ), allocatable, dimension (:,:), intent (out) :: e_po
        integer (kind = il ), intent (out) :: num_e

        integer (kind = il) :: busqueda, x, y, i
        allocate (e_p(2,ntodos))
        allocate (e_t(2,ntodos))
        allocate (e_po(2,ntodos))
        num_e=0
        Do i = 1, ntodos
            if (l_todos(1,i) /= -1) then
                busqueda = Buscar_L_Comun(i,ntodos,l_todos)
                if (busqueda/=-1) then
                    num_e = num_e + 1
                    !llenado de lados comunes
                    e_p(:,num_e) = l_todos(:,i)
                    !borrado del lado repetido
                    l_todos(1,busqueda)=-1
                    l_todos(2,busqueda)=-1
                    !l_todos(3,busqueda)=-1
                    if (mod(i,3) == 0) then
                        x = (i/3)
                    else
                        x = (i/3)+1
                    end if
                    if (mod(busqueda,3)==0) then
                        y = busqueda/3
                    else
                        y = (busqueda/3)+1
                    end if
                    !llenado de triangulos adyacentes
                    e_t(1,num_e) = x
                    e_t(2,num_e) = y
                    !llenado de vertices opuestos
                    !busqueda de opuesto del triangulo x
                    if ( (t(1,x) /= e_p(1,num_e))  .and. (t(1,x) /= e_p(2,num_e)) ) then
                        !entonces el opuesto es el punto 1 del triangulo x
                        e_po(1,num_e) = t(1,x)
                    else if ( (t(2,x) /= e_p(1,num_e)) .and. (t(2,x) /= e_p(2,num_e)) ) then
                        !entonces el opuesto es el punto 2 del triangulo x
                        e_po(1,num_e) = t(2,x)
                    else
                        !entonces el opuesto es el punto 3 del triangulo x
                        e_po(1,num_e) = t(3,x)
                    end if
                    !busqueda de opuesto del triangulo y
                    if ( (t(1,y) /= e_p(1,num_e))  .and. (t(1,y) /= e_p(2,num_e)) ) then
                        !entonces el opuesto es el punto 1 del triangulo y
                        e_po(2,num_e) = t(1,y)
                    else if ( (t(2,y) /= e_p(1,num_e)) .and. (t(2,y) /= e_p(2,num_e)) ) then
                        !entonces el opuesto es el punto 2 del triangulo y
                        e_po(2,num_e) = t(2,y)
                    else
                        !entonces el opuesto es el punto 3 del triangulo y
                        e_po(2,num_e) = t(3,y)
                    end if
                else
                    l_todos(1,i)=-1
                    l_todos(2,i)=-1
                    !l_todos(3,i)=-1
                end if

            end if

        end do

    end subroutine llenar_todo

    subroutine parametros_rwg (p_coord, t_p, num_t, e_p, e_long, e_t, e_po, num_e, e_centro)

        !dummy
        real ( kind = dp ), allocatable, dimension(:,:), intent(in) :: p_coord
        integer (kind = il ), dimension (:,:), intent(in) :: t_p
        integer (kind = il ), intent (in) :: num_t
        integer (kind = il ), intent (out) :: num_e
        integer (kind = il ), allocatable, dimension (:,:), intent(out) :: e_p
        integer (kind = il ), allocatable, dimension (:,:), intent(out) :: e_t
        integer (kind = il ), allocatable, dimension (:,:), intent(out) :: e_po
        real (kind = dp ), allocatable, dimension(:), intent (out) :: e_long
        real ( kind = dp ), allocatable, dimension(:,:),intent (out) :: e_centro
        !local
        integer (kind = il ), dimension (2,3*num_t) :: l_todos
        integer (kind = il) :: i

        !print*, 't_p:'
        !print*, t_p
        call llenar_lados (t_p, num_t, l_todos)
        !print*, 'l_todos:'
        !print*, l_todos
        call llenar_todo (l_todos,3*num_t, e_p, e_t, t_p, e_po, num_e)
        call ResizeArray(e_p,num_e,3*num_t)
        call ResizeArray(e_t,num_e,3*num_t)
        call ResizeArray(e_po,num_e, 3*num_t)
        allocate (e_long(num_e))
        allocate (e_centro(3,num_e))

        do i = 1, num_e
            e_long(i) = modulo_v( p_coord(:,e_p(1,i)) - p_coord(:,e_p(2,i)) )
            e_centro(:,i) = 0.5*( p_coord(:,e_p(1,i))+p_coord(:,e_p(2,i)) )
        end do
    !call stla_lcomun_print ( num_e, e_p)
    !call stla_tcomun_print ( num_e, e_t )
    !call stla_vo_print ( num_e, e_po)

    !deallocate ( e_p )
    !deallocate ( e_t )
    !deallocate ( e_po )


    end subroutine parametros_rwg

    function modulo_v(v1) result(md)

        real (kind = dp), dimension(3), intent(in) :: v1
        real (kind = dp) :: md
        md = sqrt(producto_escalar(v1,v1))

    end function modulo_v

    SUBROUTINE ResizeArray(A, newSize, sizeA)

        IMPLICIT NONE

        INTEGER, DIMENSION(:,:),allocatable, INTENT(INOUT) :: A
        INTEGER, INTENT(IN) :: newSize
        INTEGER, INTENT(IN) :: sizeA

        INTEGER, DIMENSION(:,:), ALLOCATABLE :: B

        !ALLOCATE(B(2,LBOUND(A):UBOUND(A))
        ALLOCATE(B(2,sizeA))

        B = A

        DEALLOCATE(A)

        ALLOCATE(A(2,newSize))

        A = B

        DEALLOCATE(B)
    END SUBROUTINE ResizeArray

    subroutine stla_lcomun_print ( node_num, node_xyz )
        integer ( kind = il ) node_num
        integer ( kind = il ) node
        integer (kind = il ) node_xyz(2,node_num)

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '    Lados         Puntos'
        write ( *, '(a)' ) ' '

        do node = 1, node_num
            write ( *, '(2x,i8,3(2x,i8))' ) node, node_xyz(1:2,node)
        end do
        return
    end subroutine stla_lcomun_print

    subroutine stla_tcomun_print ( node_num, node_xyz )
        integer ( kind = il ) node_num
        integer ( kind = il ) node
        integer (kind = il ) node_xyz(2,node_num)

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '    Lados        Triangulos Comunes'
        write ( *, '(a)' ) ' '

        do node = 1, node_num
            write ( *, '(2x,i8,3(2x,i8))' ) node, node_xyz(1:2,node)
        end do
        return
    end subroutine stla_tcomun_print

    subroutine stla_vo_print ( node_num, node_xyz )
        integer ( kind = il ) node_num
        integer ( kind = il ) node
        integer (kind = il ) node_xyz(2,node_num)

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '    Lados         Vertices'
        write ( *, '(a)' ) ' '

        do node = 1, node_num
            write ( *, '(2x,i8,3(2x,i8))' ) node, node_xyz(1:2,node)
        end do
        return
    end subroutine stla_vo_print

    function producto_cruz ( v1, v2 ) result (v3 )

        implicit none

        integer ( kind = il ), parameter :: dim_num = 3

        real ( kind = dp ) v1(dim_num)
        real ( kind = dp ) v2(dim_num)
        real ( kind = dp ) v3(dim_num)

        v3(1) = v1(2) * v2(3) - v1(3) * v2(2)
        v3(2) = v1(3) * v2(1) - v1(1) * v2(3)
        v3(3) = v1(1) * v2(2) - v1(2) * v2(1)

        return
    end function producto_cruz

    subroutine hallar_geometria_rwg (p_coord,t_p,num_p,num_t,t_area,t_baric,t_baric_sub)
        !dummy
        real ( kind = dp ), dimension (:,:), intent(in) :: p_coord
        integer ( kind = il ), dimension (:,:), intent(in) ::t_p
        integer ( kind = il ) , intent(in) :: num_p, num_t
        real ( kind = dp ), allocatable, dimension(:), intent(out) :: t_area
        real ( kind = dp ), allocatable, dimension(:,:), intent(out) :: t_baric
        real ( kind = dp ), allocatable, dimension(:,:,:), intent(out) :: t_baric_sub
        !local
        integer ( kind =  4 ) i
        !real ( kind = dp ), dimension(3,10,num_t) :: vert
        real ( kind = dp ), allocatable, dimension(:,:,:) :: vert
        real ( kind = dp ), dimension(3) :: ba, ca

        allocate(vert(3,10,num_t))
        allocate(t_baric(3,num_t))
        allocate(t_baric_sub(3,9,num_t))
        allocate(t_area(num_t))

        do i = 1, num_t
            t_baric(:,i)=(1./3.)*(p_coord(:,t_p(1,i))+p_coord(:,t_p(2,i))+p_coord(:,t_p(3,i)))
            !punto a (punto 1)
            vert(:,1,i) = p_coord(:,t_p(1,i))
            !punto b (punto 2)
            vert(:,2,i) = p_coord(:,t_p(2,i))
            !punto c (punto 3)
            vert(:,3,i) = p_coord(:,t_p(3,i))
            !punto 4
            vert(:,4,i) = p_coord(:,t_p(1,i))+(1./3.)*(p_coord(:,t_p(2,i))-p_coord(:,t_p(1,i)))
            !punto 5
            vert(:,5,i) = p_coord(:,t_p(1,i))+(2./3.)*(p_coord(:,t_p(2,i))-p_coord(:,t_p(1,i)))
            !punto 6
            vert(:,6,i) = p_coord(:,t_p(2,i))+(1./3.)*(p_coord(:,t_p(3,i))-p_coord(:,t_p(2,i)))
            !punto 7
            vert(:,7,i) = p_coord(:,t_p(2,i))+(2./3.)*(p_coord(:,t_p(3,i))-p_coord(:,t_p(2,i)))
            !punto 8
            vert(:,8,i) = p_coord(:,t_p(3,i))+(1./3.)*(p_coord(:,t_p(1,i))-p_coord(:,t_p(3,i)))
            !punto 9
            vert(:,9,i) = p_coord(:,t_p(3,i))+(2./3.)*(p_coord(:,t_p(1,i))-p_coord(:,t_p(3,i)))
            !punto 10
            vert(:,10,i) = t_baric(:,i)

            !calculo baricentro de los subtrianguos
            t_baric_sub(:,1,i) = (1./3.)*(vert(:,1,i)+vert(:,4,i)+vert(:,9,i))
            t_baric_sub(:,2,i) = (1./3.)*(vert(:,4,i)+vert(:,5,i)+vert(:,10,i))
            t_baric_sub(:,3,i) = (1./3.)*(vert(:,5,i)+vert(:,2,i)+vert(:,6,i))
            t_baric_sub(:,4,i) = (1./3.)*(vert(:,6,i)+vert(:,10,i)+vert(:,7,i))
            t_baric_sub(:,5,i) = (1./3.)*(vert(:,7,i)+vert(:,8,i)+vert(:,3,i))
            t_baric_sub(:,6,i) = (1./3.)*(vert(:,7,i)+vert(:,8,i)+vert(:,10,i))
            t_baric_sub(:,7,i) = (1./3.)*(vert(:,8,i)+vert(:,9,i)+vert(:,10,i))
            t_baric_sub(:,8,i) = (1./3.)*(vert(:,9,i)+vert(:,4,i)+vert(:,10,i))
            t_baric_sub(:,9,i) = (1./3.)*(vert(:,5,i)+vert(:,6,i)+vert(:,10,i))

            !area de los triangulos
            ba = p_coord(:,t_p(1,i)) - p_coord(:,t_p(2,i))
            ca = p_coord(:,t_p(1,i)) - p_coord(:,t_p(3,i))
            t_area (i) = 0.5*sqrt(producto_escalar(producto_cruz(ba,ca),producto_cruz(ba,ca)))

        end do

    end subroutine hallar_geometria_rwg

    function producto_escalar(v1,v2) result(valor)

        real ( kind = dp ), dimension(3), intent(in) :: v1, v2
        real ( kind = dp ) valor

        valor = (v1(1)*v2(1)) + (v1(2)*v2(2)) + (v1(3)*v2(3))


    end function producto_escalar









end module modRWG
