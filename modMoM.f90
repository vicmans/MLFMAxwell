module modMoM
    use modGlobalParam
    use modGlobalMetod
    implicit none
    private
    public :: inicializar_MoM, llenarZ_MoM, Zmn_MoM, MatxVec_MoM, setZ_MatxVec_MoM, getGLpointsWeights, &
    getNGLpointsWeights, Zmn_MoM_old
    complex (kind = dp), pointer, dimension (:,:) :: ZMoM_MxV !puntero copia de Z para la funcion de MatxVec_MoM
    real ( kind = dp ),  pointer, dimension (:,:) :: p_coord
    real ( kind = dp ), pointer, dimension (:,:) :: t_baric
    real ( kind = dp ), pointer, dimension (:,:,:) :: t_baric_sub
    real ( kind = dp ), pointer, dimension(:) :: t_area, e_long
    integer ( kind = il ), pointer, dimension (:,:) :: t_p, e_p, e_t, e_po
    integer ( kind = il ), pointer :: num_e, num_p, num_t
    real ( kind= dp ), pointer, dimension(:,:) :: t_normal

logical bandmn

contains



    subroutine inicializar_MoM (arg_p_xyz, arg_baric_t, arg_baric_sub_t, &
        arg_area_t, arg_long_l, arg_triang_p, arg_l_comun, arg_t_comun, arg_vo, arg_num_lcom, arg_num_p, arg_num_t, arg_t_normal)
        !dummy

        real ( kind = dp ),  dimension (:,:), intent(in), target :: arg_p_xyz
        real ( kind = dp ),  dimension (:,:), intent(in), target :: arg_baric_t
        real ( kind = dp ),  dimension (:,:,:), intent(in), target :: arg_baric_sub_t
        real ( kind = dp ),  dimension(:), intent(in), target :: arg_area_t, arg_long_l
        integer ( kind = il ),  dimension (:,:), intent(in), target :: arg_triang_p, arg_l_comun, arg_t_comun, arg_vo
        integer ( kind = il ), intent(in), target :: arg_num_lcom, arg_num_p, arg_num_t
        real ( kind= dp ), dimension(:,:), intent(in), target :: arg_t_normal
        !real ( kind = dp ), intent(in) :: arg_w

        !allocate(p_coord(3,arg_num_p))
        p_coord => arg_p_xyz
        !allocate(t_baric(3,arg_num_t))
        t_baric => arg_baric_t
        !allocate(t_baric_sub(3,9,arg_num_t))
        t_baric_sub => arg_baric_sub_t
        !allocate(t_area(arg_num_t))
        t_area => arg_area_t
        !allocate(e_long(arg_num_lcom))
        e_long => arg_long_l
        !allocate(t_p(3,arg_num_t))
        t_p => arg_triang_p
        !allocate(e_p(2,arg_num_lcom))
        e_p => arg_l_comun
        !allocate(e_t(2,arg_num_lcom))
        e_t => arg_t_comun
        !allocate(e_po(2,arg_num_lcom))
        e_po => arg_vo

        t_normal => arg_t_normal

        num_e => arg_num_lcom
        num_p => arg_num_p
        num_t => arg_num_t
        !w = arg_w

    end subroutine


    subroutine llenarZ_MoM(Z, mensaje)
        !dummy
        complex (kind = dp), dimension(num_e,num_e), intent(out) :: Z
        integer (kind = il), optional, intent (in) :: mensaje
        !local
        integer (kind = il) :: m, n
        !integer (kind = il) :: prevPrc, Prc
        !prevPrc=-1
        call setprc(num_e)
        do m = 1, num_e  !punto pbservación
            !    print*, m
            do n = 1, num_e  !punto fuente
                Z(m,n) = Zmn_MoM(m,n)
            end do

            if (present(mensaje)) then
                call updateprc(m)
                !Prc = (100*m/num_e)
                !if (prevPrc < Prc) then
                !    !Write(strprint, '(I5)') Prc
                !    print*, inttostr(Prc) // ' %'
                !    !call imprPRC(Prc)

                !    prevPrc = Prc
                !end if
            end if

        end do
		call setZ_MatxVec_MoM(Z)
    end subroutine llenarZ_MoM



    function Zmn_MoM_old(m,n) result (valorz)
        !dummy
        integer (kind = il), intent (in) :: m, n
        complex (kind = dp) :: valorz, valueZMFIE
        !local
        integer (kind = il) :: i
        complex (kind = dp) :: PhiP, PhiM, PintA, PintPhi
        complex (kind = dp), dimension(3) :: AmnP, AmnM, Spp, Spm, Smp, Smm
        complex (kind = dp) :: g11, g12, g21, g22 !mn: n=fuente y m=observacion del argumento de la funcion de green
        real (kind = dp), dimension(3) :: rho1, rho2, rhomPlus, rhomMinus, normalPlus, normalMinus, robsPlus, rsourcePlus, &
        robsMinus, rsourceMinus, VectorA, VectorB
        real (kind = dp) :: AreaPlus, AreaMinus, AreaPlusN, AreaMinusN, SolidAngle, signValue, NormalAngle
        complex (kind = dp) :: UNO
        !
        UNO = (jj)/(jj)

        AmnP = 0
        AmnM = 0
        PhiP = 0
        PhiM = 0

        do i=1,9
            g11 = green(t_baric(:,e_t(1,m)),t_baric_sub(:,i,e_t(1,n)))
            g12 = green(t_baric(:,e_t(1,m)),t_baric_sub(:,i,e_t(2,n)))
            g21 = green(t_baric(:,e_t(2,m)),t_baric_sub(:,i,e_t(1,n)))
            g22 = green(t_baric(:,e_t(2,m)),t_baric_sub(:,i,e_t(2,n)))
            !print*, g11
            !print*, g12
            !print*, g21
            !print*, g22
            rho1 = t_baric_sub(:,i,e_t(1,n))-p_coord(:,e_po(1,n)) !T+
            rho2 = p_coord(:,e_po(2,n))-t_baric_sub(:,i,e_t(2,n)) !T-
            !
            AmnP = AmnP + (g11*rho1) + (g12*rho2)
            AmnM = AmnM + (g21*rho1) + (g22*rho2)
            !
            PhiP = PhiP + (g11-g12)
            PhiM = PhiM + (g21-g22)
        end do
        !Constantes multiplicativas faltantes
        AmnP = AmnP*(mu/pi)*(1./72.)*e_long(n)
        AmnM = AmnM*(mu/pi)*(1./72.)*e_long(n)
        !
        PhiP = (-1.)*PhiP*e_long(n)/(jj*36*pi*w*epsi)
        PhiM = (-1.)*PhiM*e_long(n)/(jj*36*pi*w*epsi)

        !Productos internos
        !PintA = (jj*0.5*w*e_long(m))*(producto_escalar(AmnP, t_baric(:,e_t(1,m))-p_coord(:,e_po(1,m)))+producto_escalar(AmnM, &
        !-t_baric(:,e_t(2,m))+p_coord(:,e_po(2,m))))
        PintA = (jj*0.5*w*e_long(m))*(sum(AmnP*( t_baric(:,e_t(1,m))- p_coord(:,e_po(1,m)))) + sum(AmnM*(-t_baric(:,e_t(2,m))+ &
            p_coord(:,e_po(2,m)) )  ))
        PintPhi = e_long(m)*(PhiM-PhiP)

        valorz = PintPhi+PintA
        valorz = valorz*(1./(jj*w*mu))

        valueZMFIE = 0.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
        if (abs(alphaCFIE - 1) > 0.00001_dp ) then
            !Ahora vamos con CFIE
            
            AreaPlus = t_area(e_t(1,m))
            AreaMinus = t_area(e_t(2,m))
            AreaPlusN = t_area(e_t(1,n))
            AreaMinusN = t_area(e_t(2,n))  
            
            rhomPlus = t_baric(:,e_t(1,m))- p_coord(:,e_po(1,m))
            rhomMinus = p_coord(:,e_po(2,m)) - t_baric(:,e_t(2,m))        

            if (m == n) then
                !RWG overlaps
                SolidAngle = 2.*pi !Por defecto
                VectorA = ( ( p_coord(:,e_po(1,m)) - t_baric(:,e_t(1,m)) )/norm2(p_coord(:,e_po(1,m)) - t_baric(:,e_t(1,m))) )+ &
                ( (  p_coord(:,e_po(2,m)) - t_baric(:,e_t(2,m))  )/norm2(p_coord(:,e_po(2,m)) - t_baric(:,e_t(2,m)))   )
                VectorB = ( t_normal(:,e_t(1,m))/norm2(t_normal(:,e_t(1,m)))  )  +   ( t_normal(:,e_t(2,m))/ &
                    norm2(t_normal(:,e_t(2,m)))  )
                if (abs(dot_product(VectorA, VectorB)) == 0.) then
                    signValue = -1.
                else
                    signValue = dot_product(VectorA, VectorB)/abs(dot_product(VectorA, VectorB))
                end if
                signValue = (-1.)*signValue
                NormalAngle = acos(0.999999999999_dp*dot_product( t_normal(:,e_t(1,m))/norm2(t_normal(:,e_t(1,m))) , &
                    t_normal(:,e_t(2,m))/norm2(t_normal(:,e_t(2,m))) ))
                SolidAngle = 2.*(pi + (signValue*NormalAngle))
                !
                if (isnan(SolidAngle)) then
                    print*, 'SolidAngle: '            , SolidAngle, ' signo: ', signValue, ' dot: ', &
                    dot_product( t_normal(:,e_t(1,m))/norm2(t_normal(:,e_t(1,m))) , t_normal(:,e_t(2,m))/ &
                        norm2(t_normal(:,e_t(2,m))) )
                end if

                valueZMFIE = 2.*(1. - ( SolidAngle/(4.*pi))  )*( (e_long(m)**2)/8. )*(   ( dot_product(rhomPlus,rhomPlus) &
                    /AreaPlus ) &
                  +  ( dot_product(rhomMinus,rhomMinus)/AreaMinus )       )
                !print*, 'valor EFIE: ', valorz, '  valor MFIE: ', valueZMFIE            
            else
                !RWG don't overlaps


                robsPlus = t_baric(:,e_t(1,m))
                robsMinus = t_baric(:,e_t(2,m))

                normalPlus = t_normal(:,e_t(1,m))/norm2(t_normal(:,e_t(1,m)))
                normalMinus = t_normal(:,e_t(2,m))/norm2(t_normal(:,e_t(2,m)))

                Spp = 0.
                Spm = 0.
                Smp = 0.
                Smm = 0.

                do i=1,9

                    rho1 = t_baric_sub(:,i,e_t(1,n))-p_coord(:,e_po(1,n)) !T+
                    rho2 = p_coord(:,e_po(2,n))-t_baric_sub(:,i,e_t(2,n)) !T-
                    rsourcePlus = t_baric_sub(:,i,e_t(1,n))
                    rsourceMinus = t_baric_sub(:,i,e_t(2,n))


                    Spp = Spp + cross_productC(gradGreen(robsPlus, rsourcePlus), UNO*rho1)
                    Spm = Spm + cross_productC(gradGreen(robsPlus, rsourceMinus), UNO*rho2)

                    Smp = Smp + cross_productC(gradGreen(robsMinus, rsourcePlus), UNO*rho1)
                    Smm = Smm + cross_productC(gradGreen(robsMinus, rsourceMinus), UNO*rho2)


                end do
                Spp = cross_productC(normalPlus*UNO, Spp)
                Spm = cross_productC(normalPlus*UNO, Spm)

                Smp = cross_productC(normalMinus*UNO, Smp)
                Smm = cross_productC(normalMinus*UNO, Smm)

                valueZMFIE = (  e_long(m)*e_long(n)/(144.*pi)  )*(   sum(rhomPlus*(Spp + Spm))  +  sum(rhomMinus*(Smp + Smm)))
                
            end if
        end if
        !
        !valueZMFIE = valueZMFIE*(-1.)
        valorz = alphaCFIE*valorz + (jj/( getkappa() ) )*(1-alphaCFIE)*valueZMFIE
    end function Zmn_MoM_old



    function Zmn_MoM(m,n) result (valorz)
        !dummy
        integer (kind = il), intent (in) :: m, n
        complex (kind = dp) :: valorz
        !local
        integer (kind = il) :: i_m, i_n, tri_m, tri_n
        real (kind = dp) :: threshold, value_overlapping, value_near
        complex (kind = dp) :: value_soft, value_total, resultZefie, resultZMfie, value_soft2
        logical justEFIE, justMFIE
        !
        justEFIE = (abs(alphaCFIE - 1)<=0.001_dp )
        justMFIE = (abs(alphaCFIE - 0)<=0.001_dp )

        threshold = 0.3_dp*getlambda()
        resultZefie = 0.
        resultZMfie = 0.

if (conf%efieOldMode == 1) then
    valorz = Zmn_MoM_old(m,n)

else

        do i_m = 1, 2
            do i_n = 1, 2
                tri_m = e_t(i_m, m)
                tri_n = e_t(i_n, n)

                if (( norm2(t_baric(:,tri_m) - t_baric(:, tri_n)) <= threshold ) ) then
                    !If triangles are near
                    if (( tri_m == tri_n ) )  then
                        !If triangles overlaps


                        if (.not. justEFIE) then
                            !MFIE for overlapping
                            value_total = integrationSoftMFIE( e_po(i_m, m), e_po(i_n, n) , tri_m, tri_n, .true., .true. )
                            if (i_m /= i_n) then
                                !multiply by minus
                                value_total = -1.*value_total
                            end if  
                            resultZMfie = resultZMfie + value_total*(e_long(m)*e_long(n))/(t_area(tri_m)*t_area(tri_n))
                        end if

                    else
                        !near but not overlapping

                        if (.not. justEFIE) then
                            !MFIE for near
                            value_soft = integrationSoftMFIE( e_po(i_m, m), e_po(i_n, n) , tri_m, tri_n, .false., .false. )
                            value_near = integrationNearMFIE( e_po(i_m, m), e_po(i_n, n) , tri_m, tri_n)

                            if (i_m /= i_n) then
                                !multiply by minus
                                value_soft = -1.*value_soft
                                value_near = -1.*value_near
                            end if

                            value_total = value_soft + value_near
                            resultZMfie = resultZMfie + value_total*(e_long(m)*e_long(n))/(t_area(tri_m)*t_area(tri_n))
                        end if

                    end if

                    if (.not. justMFIE) then
                        !EFIE for near
                        value_near = integrationNear(e_po(i_m, m), e_po(i_n, n) , tri_m, tri_n)
                        value_soft = integrationSoft(e_po(i_m, m), e_po(i_n, n) , tri_m, tri_n, .false. )
                        if (i_m /= i_n) then
                            !multiply by minus
                            value_near = -1.*value_near
                            value_soft = -1.*value_soft
                        end if                     
                        value_total = value_near + value_soft
                        value_total = ((e_long(m)*e_long(n))/(t_area(tri_m)*t_area(tri_n)))*value_total
                        resultZefie = resultZefie + value_total
                    end if

                else
                    !If triangles are far
                    if (.not. justMFIE) then
                        !EFIE for far
                        value_total = integrationSoft(e_po(i_m, m), e_po(i_n, n) , tri_m, tri_n, .true. )
                        if (i_m /= i_n) then
                            value_total = -1.*value_total
                        end if
                        value_total = ((e_long(m)*e_long(n))/(t_area(tri_m)*t_area(tri_n)))*value_total
                        resultZefie = resultZefie + value_total
                    end if
                    
                    if (.not. justEFIE) then
                        !MFIE for far
                        value_total = integrationSoftMFIE(e_po(i_m, m), e_po(i_n, n) , tri_m, tri_n, .true., .false. )

                        if (i_m /= i_n) then
                            value_total = -1.*value_total
                        end if
                        resultZMfie = resultZMfie + value_total*(e_long(m)*e_long(n))/(t_area(tri_m)*t_area(tri_n))
                    end if

                end if
            end do 
        end do 

        valorz = alphaCFIE*resultZefie*(1._dp/(4._dp*pi))  +  (jj/( getkappa() ) )*(1-alphaCFIE)*resultZMfie

end if

    end function Zmn_MoM

    subroutine analytical_sourceIntegrals(tri_n, rtest_in, I1_out, I2_out, I3_out, I4_out, onlyI1I2)
        !dummy
        integer (kind = il), intent(in) ::  tri_n
        real (kind = dp), dimension(3), intent(in) :: rtest_in
        real (kind = dp), intent(out), dimension(3) :: I1_out, I3_out
        real (kind = dp), intent(out) :: I2_out, I4_out
        logical, intent(in) :: onlyI1I2
        !
        !local
        integer (kind = il) :: i
        real (kind = dp), dimension(3) :: Vn1, Vn2, Vn3
        real (kind = dp), dimension(3) :: rtest, n_uni, r_minus, r_plus, rho, rho_plus, rho_minus, l_uni, u_uni, Po_uni, &
        rhom, nm_uni
        real (kind = dp), dimension(3) :: I1, I3, IntegA, IntegB
        real (kind = dp) :: I2, I4
        real (kind = dp), dimension(3,3) :: vertices_n
        real (kind = dp) :: cc, rest1, rest2, ee
        real (kind = dp) :: l_plus, l_minus, P_o, P_plus, P_minus, d, R_o, RR_minus, RR_Plus
        !
        rtest = rtest_in

        Vn1 = p_coord(:, t_p(1, tri_n))
        Vn2 = p_coord(:, t_p(2, tri_n))
        Vn3 = p_coord(:, t_p(3, tri_n))      

        vertices_n (:,1) = Vn1
        vertices_n (:,2) = Vn2
        vertices_n (:,3) = Vn3
        
        n_uni = cross_productR(Vn2 - Vn1, Vn3 - Vn1)
        n_uni = n_uni/norm2(n_uni)

        ee = 0.0001_dp*getlambda()


        do i = 1,3
            do
                if ( i == 3 ) then
                    r_minus(:) = vertices_n(:,i)
                    r_plus(:) = vertices_n(:,1)
                else
                    r_minus(:) = vertices_n(:,i)
                    r_plus(:) = vertices_n(:, i+1)
                end if

                rho = rtest - n_uni*dot_product(n_uni,rtest)

                rho_plus = r_plus - n_uni*dot_product(n_uni,r_plus)
                rho_minus = r_minus - n_uni*dot_product(n_uni,r_minus)
                l_uni = (rho_plus - rho_minus)/norm2(rho_plus - rho_minus)
                u_uni = cross_productR(l_uni,n_uni)

                P_o = abs(dot_product((rho_plus-rho),u_uni))
                d = dot_product(n_uni, (rtest - r_plus))
                R_o = sqrt( (P_o**2._dp) + (d**2._dp) )

                if ( R_o <= 0.00001_dp*getlambda() ) then
                    rtest = rtest + ee*u_uni
                    cycle
                else
                    exit
                end if
            end do
        end do


        I1(:) = 0._dp
        I2 = 0._dp
        I3(:) = 0._dp
        I4 = 0._dp

        rho = rtest - n_uni*dot_product(n_uni,rtest)

        do i = 1, 3

            !end point of line segment C
            if ( i == 3 ) then
                r_minus(:) = vertices_n(:,i)
                r_plus(:) = vertices_n(:,1)
            else
                r_minus(:) = vertices_n(:,i)
                r_plus(:) = vertices_n(:, i+1)
            end if

            rho_plus = r_plus - n_uni*dot_product(n_uni,r_plus)
            rho_minus = r_minus - n_uni*dot_product(n_uni,r_minus)
            l_uni = (rho_plus - rho_minus)/norm2(rho_plus - rho_minus)
            u_uni = cross_productR(l_uni,n_uni)
            l_plus = dot_product((rho_plus-rho),l_uni)
            l_minus = dot_product((rho_minus-rho),l_uni)
            P_o = abs(dot_product((rho_plus-rho),u_uni))
            P_plus = norm2(rho_plus - rho)
            P_minus = norm2(rho_minus - rho)
            if (P_o < 0.0000000001_dp*getlambda()) then
                Po_uni(:) = 0.

            else
                Po_uni = ((rho_plus - rho) - (l_plus*l_uni))/P_o
            end if

            d = dot_product(n_uni, (rtest - r_plus))
            R_o = sqrt( (P_o**2._dp) + (d**2._dp) )

            RR_Plus = sqrt( (P_plus**2._dp) + (d**2._dp) )
            RR_minus = sqrt( (P_minus**2._dp) + (d**2._dp) )

            !I1
            cc = 0.5_dp*( (R_o**2)*log( (RR_Plus+l_plus)/(RR_minus+l_minus) ) + l_plus*RR_Plus - l_minus*RR_minus  )
            I1 = I1 + u_uni*cc
            !I2
            I2 = I2 + dot_product(Po_uni,u_uni)*( P_o*log( (RR_Plus+l_plus)/(RR_minus+l_minus) ) - &
                abs(d)*( atan(P_o*l_plus/((R_o**2) &
                + abs(d)*RR_Plus)) - atan(P_o*l_minus/((R_o**2) + abs(d)*RR_minus)) ) )
            
            if (.not. onlyI1I2) then
                I3 = I3 + (-1._dp)*u_uni*log((RR_Plus + l_plus)/(RR_minus + l_minus))
                if (abs(d) <= 0.0001_dp*getlambda()) then


                    if (P_o < 0.00000000001_dp*getlambda()) then
                        I4 = I4 + 0._dp

                    else
                        cc = (-1._dp)*dot_product(Po_uni, u_uni)*( (l_plus/(P_o*RR_Plus)) - (l_minus/(P_o*RR_minus)) )
                        I4 = I4 + cc                   
                    end if
                else
                    rest1 = (R_o**2._dp)+ (abs(d)*RR_Plus)
                    rest2 = (R_o**2._dp)+ (abs(d)*RR_minus)

                    I4 = I4 + dot_product(Po_uni, u_uni)*(1._dp/abs(d))*( atan( (P_o*l_plus)/( rest1 ) )  - &
                    atan( (P_o*l_minus)/( rest2 ) ) )

                end if
            end if
        end do

        I1_out = I1
        I2_out = I2
        I3_out = I3
        I4_out = I4



    end subroutine analytical_sourceIntegrals


    function integrationNearMFIE(Vom_ind, Von_ind, tri_m, tri_n) result(valor)
        !dummy
        integer (kind = il), intent(in) :: Vom_ind, Von_ind, tri_m, tri_n
        real (kind = dp) :: valor
        !local
        real (kind = dp), dimension(3) :: Vm1, Vm2, Vm3, Vom, Von
        real (kind = dp), dimension(3,3) :: vertices_n
        real (kind = dp), dimension(3) :: n_uni, r_minus, r_plus, rho, rho_plus, rho_minus, l_uni, u_uni, Po_uni, I1, rtest, &
         rhom, Von_proy, IntegA, IntegB, I3, nm_uni, Rn, Iinternal
        real (kind = dp) :: l_plus, l_minus, P_o, P_plus, P_minus, d, R_o, v, RR_minus, RR_Plus, I2, cc, ee, IntegralI, I4
        integer( kind = il) :: GL_numberpoints, itest
        real (kind = dp), dimension(:,:), allocatable :: GL_pointsweights_test
        !

        Vm1 = p_coord(:, t_p(1, tri_m))
        Vm2 = p_coord(:, t_p(2, tri_m))
        Vm3 = p_coord(:, t_p(3, tri_m))
        Vom = p_coord(:, Vom_ind)

        Von = p_coord(:, Von_ind)


        nm_uni = t_normal(:, tri_m)/norm2(t_normal(:, tri_m))
        n_uni = t_normal(:, tri_n)/norm2(t_normal(:, tri_n))

        GL_numberpoints = 7
        allocate(GL_pointsweights_test(GL_numberpoints, 4))
        GL_pointsweights_test = getGLpointsWeights(Vm1, Vm2, Vm3)

        IntegralI = 0._dp
        
        do itest = 1, GL_numberpoints
            rtest = GL_pointsweights_test(itest, 1:3)
            rho = rtest - n_uni*dot_product(n_uni,rtest)
            rhom = rtest - Vom
            Rn = rtest - Von
            Von_proy = Von - n_uni*dot_product(n_uni,Von)

            call analytical_sourceIntegrals(tri_n, rtest, I1, I2, I3, I4, .false.)

            IntegB = I1 + (rho - Von_proy)*I2
            IntegA = I3 + (rho - Von_proy)*I4

            Iinternal = IntegA + 0.5_dp*(getkappa()**2)*IntegB

            IntegralI = IntegralI + sum(rhom*cross_productR( nm_uni, cross_productR(Rn, Iinternal) ))* &
             GL_pointsweights_test(itest, 4)*t_area(tri_m)
        end do

        valor = IntegralI/(16._dp*pi)
    end function integrationNearMFIE


    function integrationSoftMFIE( Vom_ind, Von_ind, tri_m, tri_n, isCompleteGreen, isOverlapped ) result(valint)
        !dummy
        integer (kind = il), intent(in) :: Vom_ind, Von_ind, tri_m, tri_n
        complex (kind = dp) :: valint
        logical :: isCompleteGreen, isOverlapped
        !
        !local
        integer( kind = il) :: GL_numberpoints, itest, isour
        real (kind = dp), dimension(:,:), allocatable :: GL_pointsweights_test, GL_pointsweights_sour
        real (kind = dp), dimension(3) :: Vm1, Vm2, Vm3, Vn1, Vn2, Vn3, Vom, Von, rtest, rsour, rhom, rhon, Rn, n_uni
        complex (kind = dp), dimension(3) :: Ii
        complex (kind = dp) :: UNO
        !
        UNO = jj/jj
        GL_numberpoints = 7
        allocate(GL_pointsweights_test(GL_numberpoints, 4))
        allocate(GL_pointsweights_sour(GL_numberpoints, 4))

        Vm1 = p_coord(:, t_p(1, tri_m))
        Vm2 = p_coord(:, t_p(2, tri_m))
        Vm3 = p_coord(:, t_p(3, tri_m))
        Vom = p_coord(:, Vom_ind)

        Vn1 = p_coord(:, t_p(1, tri_n))
        Vn2 = p_coord(:, t_p(2, tri_n))
        Vn3 = p_coord(:, t_p(3, tri_n))
        Von = p_coord(:, Von_ind)

        n_uni = t_normal(:, tri_m)/norm2(t_normal(:, tri_m))

        GL_pointsweights_test = getGLpointsWeights(Vm1, Vm2, Vm3)
        GL_pointsweights_sour = getGLpointsWeights(Vn1, Vn2, Vn3)
        valint = 0._dp
        Ii(:) = 0._dp

        if (isOverlapped) then

            do itest = 1, GL_numberpoints
                rtest = GL_pointsweights_test(itest, 1:3)
                rhom = rtest - Vom

                    rsour = rtest!GL_pointsweights_sour(isour, 1:3)
                    rhon = rsour - Von
                    valint = valint + dot_product(rhom,rhon)*GL_pointsweights_test(itest, 4)

            end do
            valint = valint*t_area(tri_m)/(8._dp)

        else

            do itest = 1, GL_numberpoints
                rtest = GL_pointsweights_test(itest, 1:3)
                rhom = rtest - Vom
                Rn = rtest - Von
                Ii(:) = 0._dp
                do isour = 1, GL_numberpoints
                    rsour = GL_pointsweights_sour(isour, 1:3)
                    rhon = rsour - Von
                    if (isCompleteGreen) then
                        Ii = Ii + rhon*greenNormalMFIE(rtest, rsour)*GL_pointsweights_test(itest, 4)*GL_pointsweights_sour(isour,4)
                    else
                        Ii = Ii + rhon*greenSoftMFIE(rtest, rsour)*GL_pointsweights_test(itest, 4)*GL_pointsweights_sour(isour, 4)
                    end if
                end do
                Ii = cross_productC(Rn*UNO, Ii)
                valint = valint + sum(rhom*cross_productC(n_uni*UNO, Ii))
            end do
            valint = valint*t_area(tri_m)*t_area(tri_n)/(16._dp*pi)

        end if

    end function integrationSoftMFIE


    function integrationOverlapping(V1, V2, V3, Vm, Vn, tri_m) result(valor)
        !dummy
        real (kind = dp), dimension(3), intent(in) :: V1, V2, V3, Vm, Vn
        integer (kind = il), intent(in) :: tri_m
        real (kind = dp) :: valor
        !
        !local
        real (kind = dp), dimension(3) :: valABC
        real (kind = dp) :: I1, I2
        !
        I1 = 0.
        I2 = 0.
        valor = 0.

        !I1
        I1 = I1 + dot_product((V1-Vm),(V1-Vn))
        valABC = getABC(V1,V2,V3)
        I1 = I1 + dot_product((V2-V1),(2.*V1-Vm-Vn))*eqnAnalytic1(valABC(1), valABC(2), valABC(3))
        !change V2->V3
        valABC = getABC(V1,V3,V2)
        I1 = I1 + dot_product((V3-V1),(2.*V1-Vm-Vn))*eqnAnalytic1(valABC(1), valABC(2), valABC(3))
        valABC = getABC(V1,V2,V3)
        I1 = I1 + dot_product((V2-V1),(V2-V1))*eqnAnalytic11(valABC(1), valABC(2), valABC(3))
        !change V2->V3
        valABC = getABC(V1,V3,V2)
        I1 = I1 + dot_product((V3-V1),(V3-V1))*eqnAnalytic11(valABC(1), valABC(2), valABC(3))
        valABC = getABC(V1,V2,V3)
        I1 = I1 + dot_product((V2-V1),(V3-V1))*eqnAnalytic12(valABC(1), valABC(2), valABC(3))
        !change V2->V3
        valABC = getABC(V1,V3,V2)
        I1 = I1 + dot_product((V3-V1),(V2-V1))*eqnAnalytic12(valABC(1), valABC(2), valABC(3))

        !I2
        valABC = getABC(V1,V2,V3)
        I2 = eqnAnalytic0(valABC(1), valABC(2), valABC(3))

        !value
        valor = 0.25_dp*I1 - ( (1._dp/(getkappa()**2._dp))*I2 )
        !fix value
        valor = valor*4._dp*(t_area(tri_m)**2)

    end function integrationOverlapping



    function integrationNear(Vom_ind, Von_ind, tri_m, tri_n) result(valor)
        !dummy
        integer (kind = il), intent(in) :: Vom_ind, Von_ind, tri_m, tri_n
        real (kind = dp) :: valor
        !local
        real (kind = dp), dimension(3) :: Vm1, Vm2, Vm3, Vom, Von
        real (kind = dp), dimension(3,3) :: vertices_n
        real (kind = dp), dimension(3) :: n_uni, r_minus, r_plus, rho, rho_plus, rho_minus, l_uni, u_uni, Po_uni, I1, rtest, &
         rhom, Von_proy, IntegA, IntegB, I3, nm_uni, Rn, Iinternal
        real (kind = dp) :: l_plus, l_minus, P_o, P_plus, P_minus, d, R_o, v, RR_minus, RR_Plus, I2, cc, ee, IntegralI, I4
        integer( kind = il) :: GL_numberpoints, itest
        real (kind = dp), dimension(:,:), allocatable :: GL_pointsweights_test
        !

        Vm1 = p_coord(:, t_p(1, tri_m))
        Vm2 = p_coord(:, t_p(2, tri_m))
        Vm3 = p_coord(:, t_p(3, tri_m))
        Vom = p_coord(:, Vom_ind)

        Von = p_coord(:, Von_ind)


        nm_uni = t_normal(:, tri_m)/norm2(t_normal(:, tri_m))
        n_uni = t_normal(:, tri_n)/norm2(t_normal(:, tri_n))

        GL_numberpoints = 7
        allocate(GL_pointsweights_test(GL_numberpoints, 4))
        GL_pointsweights_test = getGLpointsWeights(Vm1, Vm2, Vm3)

        IntegralI = 0._dp
        
        do itest = 1, GL_numberpoints
            rtest = GL_pointsweights_test(itest, 1:3)
            rho = rtest - n_uni*dot_product(n_uni,rtest)
            rhom = rtest - Vom
            Rn = rtest - Von
            Von_proy = Von - n_uni*dot_product(n_uni,Von)

            call analytical_sourceIntegrals(tri_n, rtest, I1, I2, I3, I4, .true.)

        IntegralI = IntegralI + ( dot_product(rhom, 0.25_dp*(I1 + (rho - Von_proy)*I2)) - (1._dp/(getkappa()**2))*I2 )* &
        GL_pointsweights_test(itest, 4)*t_area(tri_m)
        end do

        valor = IntegralI
    end function integrationNear


    function integrationSoft( Vom_ind, Von_ind, tri_m, tri_n, isCompleteGreen ) result(valint)
        !dummy
        integer (kind = il), intent(in) :: Vom_ind, Von_ind, tri_m, tri_n
        complex (kind = dp) :: valint
        logical :: isCompleteGreen
        !
        !local
        integer( kind = il) :: GL_numberpoints, itest, isour
        real (kind = dp), dimension(:,:), allocatable :: GL_pointsweights_test, GL_pointsweights_sour
        real (kind = dp), dimension(3) :: Vm1, Vm2, Vm3, Vn1, Vn2, Vn3, Vom, Von, rtest, rsour, rhom, rhon
        !
        GL_numberpoints = 7
        allocate(GL_pointsweights_test(GL_numberpoints, 4))
        allocate(GL_pointsweights_sour(GL_numberpoints, 4))

        Vm1 = p_coord(:, t_p(1, tri_m))
        Vm2 = p_coord(:, t_p(2, tri_m))
        Vm3 = p_coord(:, t_p(3, tri_m))
        Vom = p_coord(:, Vom_ind)

        Vn1 = p_coord(:, t_p(1, tri_n))
        Vn2 = p_coord(:, t_p(2, tri_n))
        Vn3 = p_coord(:, t_p(3, tri_n))
        Von = p_coord(:, Von_ind)

        GL_pointsweights_test = getGLpointsWeights(Vm1, Vm2, Vm3)
        GL_pointsweights_sour = getGLpointsWeights(Vn1, Vn2, Vn3)
        valint = 0._dp
        do itest = 1, GL_numberpoints
            rtest = GL_pointsweights_test(itest, 1:3)
            rhom = rtest - Vom
            do isour = 1, GL_numberpoints
                rsour = GL_pointsweights_sour(isour, 1:3)
                rhon = rsour - Von
                if (isCompleteGreen) then
                    valint = valint + ( (0.25_dp)*dot_product(rhom, rhon) - (1._dp/(getkappa()**2._dp)) )*greenNormal(rtest, &
                     rsour)*GL_pointsweights_test(itest, 4)*GL_pointsweights_sour(isour, 4)
                else
                    valint = valint + ( (0.25_dp)*dot_product(rhom, rhon) - (1._dp/(getkappa()**2._dp)) )*greenSoft(rtest, &
                        rsour)*GL_pointsweights_test(itest, 4)*GL_pointsweights_sour(isour, 4)
                end if
            end do
        end do
        valint = valint*t_area(tri_m)*t_area(tri_n)
    end function integrationSoft


    function greenNormal(rtest, rsour) result(valgreen)
        !dummy
        real (kind = dp), dimension(3), intent(in) :: rtest, rsour
        complex(kind = dp) :: valgreen
        !
        !local
        real (kind = dp) :: rdif
        !
        rdif = norm2(rtest - rsour)

        valgreen = exp(-jj*getkappa()*rdif)/(rdif) 

    end function

    function greenSoft(rtest, rsour) result(valgreen)
        !dummy
        real (kind = dp), dimension(3), intent(in) :: rtest, rsour
        complex(kind = dp) :: valgreen
        !
        !local
        real (kind = dp) :: rdif
        !
        rdif = norm2(rtest - rsour)

        valgreen = ( exp(-jj*getkappa()*rdif)/(rdif) ) - (1./rdif)
        if (rdif == 0._dp) then
            !Evaluar el limite
            valgreen = -jj*getkappa()            
        end if

    end function

    function greenNormalMFIE(rtest, rsour) result(valgreen)
        !dummy
        real (kind = dp), dimension(3), intent(in) :: rtest, rsour
        complex(kind = dp) :: valgreen
        !
        !local
        real (kind = dp) :: rdif
        !
        rdif = norm2(rtest - rsour)

        valgreen = (1._dp + jj*getkappa()*rdif)*exp(-jj*getkappa()*rdif)/(rdif**3._dp)

    end function

    function greenSoftMFIE(rtest, rsour) result(valgreen)
        !dummy
        real (kind = dp), dimension(3), intent(in) :: rtest, rsour
        complex(kind = dp) :: valgreen
        !
        !local
        real (kind = dp) :: rdif
        !
        rdif = norm2(rtest - rsour)

        valgreen = (   (1._dp + jj*getkappa()*rdif)*exp(-jj*getkappa()*rdif) - (1._dp + 0.5_dp*(getkappa()**2._dp)* &
            (rdif**2._dp) ))/(rdif**3._dp)
        if (rdif == 0._dp) then
            !Evaluar el limite
            valgreen = -jj*(getkappa()**3._dp)/3._dp
        end if

    end function    

    function getGLpointsWeights(V1, V2, V3) result(points_weights)
        !dummy
        real (kind = dp), dimension(3), intent(in) :: V1, V2, V3
        real (kind = dp), dimension(7,4) :: points_weights ! points_weights(:,1:3)->points; points_weights(:,4)->weights;
        !
        !local
        real (kind = dp), dimension(7,4) :: pw_table
        !

        pw_table(1,:) = (/0.333333333333333_dp, 0.333333333333333_dp, 0.333333333333333_dp, 0.225000000000000_dp /)
        pw_table(2,:) = (/0.059715871789770_dp, 0.470142064105115_dp, 0.470142064105115_dp, 0.132394152788506_dp /)
        pw_table(3,:) = (/0.470142064105115_dp, 0.059715871789770_dp, 0.470142064105115_dp, 0.132394152788506_dp /)
        pw_table(4,:) = (/0.470142064105115_dp, 0.470142064105115_dp, 0.059715871789770_dp, 0.132394152788506_dp /)
        pw_table(5,:) = (/0.797426985353087_dp, 0.101286507323456_dp, 0.101286507323456_dp, 0.125939180544827_dp /)
        pw_table(6,:) = (/0.101286507323456_dp, 0.797426985353087_dp, 0.101286507323456_dp, 0.125939180544827_dp /)
        pw_table(7,:) = (/0.101286507323456_dp, 0.101286507323456_dp, 0.797426985353087_dp, 0.125939180544827_dp /)


        points_weights(:, 1) = pw_table(:, 1)*V2(1) + pw_table(:, 2)*V3(1) + pw_table(:, 3)*V1(1)
        points_weights(:, 2) = pw_table(:, 1)*V2(2) + pw_table(:, 2)*V3(2) + pw_table(:, 3)*V1(2)
        points_weights(:, 3) = pw_table(:, 1)*V2(3) + pw_table(:, 2)*V3(3) + pw_table(:, 3)*V1(3)

        points_weights(:, 4) = pw_table(:, 4)

    end function getGLpointsWeights

    function getNGLpointsWeights(V1, V2, V3, n_points) result(points_weights)
        !dummy
        real (kind = dp), dimension(3), intent(in) :: V1, V2, V3
        real (kind = dp), dimension(n_points,4) :: points_weights ! points_weights(:,1:3)->points; points_weights(:,4)->weights;
        integer (kind = il) :: n_points
        !
        !local
        real (kind = dp), dimension(n_points,4) :: pw_table
        !


        if (n_points == 7) then
            pw_table(1,:) = (/0.333333333333333_dp, 0.333333333333333_dp, 0.333333333333333_dp, 0.225000000000000_dp /)
            pw_table(2,:) = (/0.059715871789770_dp, 0.470142064105115_dp, 0.470142064105115_dp, 0.132394152788506_dp /)
            pw_table(3,:) = (/0.470142064105115_dp, 0.059715871789770_dp, 0.470142064105115_dp, 0.132394152788506_dp /)
            pw_table(4,:) = (/0.470142064105115_dp, 0.470142064105115_dp, 0.059715871789770_dp, 0.132394152788506_dp /)
            pw_table(5,:) = (/0.797426985353087_dp, 0.101286507323456_dp, 0.101286507323456_dp, 0.125939180544827_dp /)
            pw_table(6,:) = (/0.101286507323456_dp, 0.797426985353087_dp, 0.101286507323456_dp, 0.125939180544827_dp /)
            pw_table(7,:) = (/0.101286507323456_dp, 0.101286507323456_dp, 0.797426985353087_dp, 0.125939180544827_dp /)
        elseif (n_points == 16) then
            pw_table(1,:) = (/0.333333333333333_dp, 0.333333333333333_dp, 0.333333333333333_dp, 0.144315607677787_dp /)
            pw_table(2,:) = (/0.081414823414554_dp, 0.459292588292723_dp, 0.459292588292723_dp, 0.095091634267285_dp /)
            pw_table(3,:) = (/0.459292588292723_dp, 0.459292588292723_dp, 0.081414823414554_dp, 0.095091634267285_dp /)
            pw_table(4,:) = (/0.459292588292723_dp, 0.081414823414554_dp, 0.459292588292723_dp, 0.095091634267285_dp /)
            pw_table(5,:) = (/0.658861384496480_dp, 0.170569307751760_dp, 0.170569307751760_dp, 0.103217370534718_dp /)
            pw_table(6,:) = (/0.170569307751760_dp, 0.170569307751760_dp, 0.658861384496480_dp, 0.103217370534718_dp /)
            pw_table(7,:) = (/0.170569307751760_dp, 0.658861384496480_dp, 0.170569307751760_dp, 0.103217370534718_dp /)
            pw_table(8,:) = (/0.898905543365938_dp, 0.050547228317031_dp, 0.050547228317031_dp, 0.032458497623198_dp /)
            pw_table(9,:) = (/0.050547228317031_dp, 0.898905543365938_dp, 0.050547228317031_dp, 0.032458497623198_dp /)
            pw_table(10,:) = (/0.050547228317031_dp, 0.050547228317031_dp, 0.898905543365938_dp, 0.032458497623198_dp /)

            pw_table(11,:) = (/0.008394777409958_dp, 0.263112829634638_dp, 0.728492392955404_dp, 0.027230314174435_dp /)
            pw_table(12,:) = (/0.263112829634638_dp, 0.728492392955404_dp, 0.008394777409958_dp, 0.027230314174435_dp /)
            pw_table(13,:) = (/0.728492392955404_dp, 0.008394777409958_dp, 0.263112829634638_dp, 0.027230314174435_dp /)
            pw_table(14,:) = (/0.728492392955404_dp, 0.263112829634638_dp, 0.008394777409958_dp, 0.027230314174435_dp /)
            pw_table(15,:) = (/0.263112829634638_dp, 0.008394777409958_dp, 0.728492392955404_dp, 0.027230314174435_dp /)
            pw_table(16,:) = (/0.008394777409958_dp, 0.728492392955404_dp, 0.263112829634638_dp, 0.027230314174435_dp /)
        end if

            points_weights(:, 1) = pw_table(:, 1)*V2(1) + pw_table(:, 2)*V3(1) + pw_table(:, 3)*V1(1)
            points_weights(:, 2) = pw_table(:, 1)*V2(2) + pw_table(:, 2)*V3(2) + pw_table(:, 3)*V1(2)
            points_weights(:, 3) = pw_table(:, 1)*V2(3) + pw_table(:, 2)*V3(3) + pw_table(:, 3)*V1(3)

            points_weights(:, 4) = pw_table(:, 4)            

    end function getNGLpointsWeights


    function getABC(V1, V2, V3) result(valABC)
        !dummy
        real (kind = dp), dimension(3), intent(in) :: V1, V2, V3
        real (kind = dp), dimension(3) :: valABC
        !local
        real (kind = dp) :: valA, valB, valC

        valA = dot_product((V1-V3),(V1-V3))    
        valB = dot_product((V1-V3),(V1-V2)) 
        valC = dot_product((V1-V2),(V1-V2))

        valABC(1) = valA
        valABC(2) = valB
        valABC(3) = valC

    end function getABC


    function eqnAnalytic11(valA, valB, valC) result(valor)
        !dummy
        real (kind = dp), intent(in) :: valA, valB, valC
        real (kind = dp) :: valor
        !

        !local
        real (kind = dp) :: sqabc, sqa, sqb, sqc
        !
        sqabc = sqrt(valA - (2.*valB) + valC)
        sqa = sqrt(valA)
        sqb = sqrt(valB)
        sqc = sqrt(valC)
        valor = 0.

        valor = valor + &
                log( (valB + (sqa*sqc))/(valB - valC + (sqc*sqabc)) ) / (40.*sqc)

        valor = valor + &
                log( (-valB + valC + (sqc*sqabc))/(-valB + (sqa*sqc)) ) / (40.*sqc)           

        valor = valor + ((sqa*sqabc) - (sqc*sqabc)) &
                        / (120.*(sqabc**3))

        valor = valor + ( (2.*valA) - (5.*valB) + (3.*valC) ) & 
                        *log( ( (valA - valB + (sqa*sqabc))*(valc - valB + (sqc*sqabc)) )/( (valB - valA + (sqa*sqabc))* &
                            (valB - valC + (sqc*sqabc)) ) ) &
                        /(120.*(sqabc**3))

        valor = valor + ((-sqa*sqc) + (sqa*sqabc)) &
                        /(120.*(sqa**3))

        valor = valor + ((2.*valA) + valB)*log( ((valB + sqa*sqc)*(valA - valB + sqa*sqabc))/((-valB + sqa*sqc)*(-valA + &
            valB + sqa*sqabc)) )/(120.*(sqa**3))
                

    end function eqnAnalytic11


    function eqnAnalytic12(valA, valB, valC) result(valor)
        !dummy
        real (kind = dp), intent(in) :: valA, valB, valC
        real (kind = dp) :: valor
        !

        !local
        real (kind = dp) :: sqabc, sqa, sqb, sqc
        !
        sqabc = sqrt(valA - (2.*valB) + valC)
        sqa = sqrt(valA)
        sqb = sqrt(valB)
        sqc = sqrt(valC)
        valor = 0.

        valor = valor + &
                log( (valB + (sqa*sqc))/(valB - valC + (sqc*sqabc)) ) / (120.*sqc)

        valor = valor + &
                log( (valA - valB + (sqa*sqabc))/(-valB + (sqa*sqc)) ) / (120.*sqa)           

        valor = valor + ((-sqa*sqabc) + (sqc*sqabc)) &
                        / (120.*(sqabc**3))

        valor = valor + ( (2.*valA) - (3.*valB) + (valC) ) & 
                        *log( (valA - valB + (sqa*sqabc))/( (valB - valC + (sqc*sqabc)) ) ) &
                        /(120.*(sqabc**3))

        valor = valor + ((sqa*sqabc) - (sqc*sqabc)) &
                        /(120.*(sqabc**3))

        valor = valor + ( (valA) - (3.*valB) + (2.*valC) ) & 
                        *log( (-valB + valC + (sqc*sqabc))/( (-valA + valB + (sqa*sqabc)) ) ) &
                        /(120.*(sqabc**3))

        valor = valor + ((-3.*sqa*sqc) + (3.*sqc*sqabc)) &
                        /(120.*(sqc**3))

        valor = valor + ( (3.*valB) + (2.*valC) ) & 
                        *log( (-valB + valC + (sqc*sqabc))/( (-valB + (sqa*sqc)) ) ) &
                        /(120.*(sqc**3))

        valor = valor + ((-3.*sqa*sqc) + (3.*sqa*sqabc)) &
                        /(120.*(sqa**3))

        valor = valor + ((2.*valA) + (3.*valB))*log( (valB + sqa*sqc)/(-valA + valB + sqa*sqabc) )/(120.*(sqa**3))
                

    end function eqnAnalytic12


    function eqnAnalytic1(valA, valB, valC) result(valor)
        !dummy
        real (kind = dp), intent(in) :: valA, valB, valC
        real (kind = dp) :: valor
        !

        !local
        real (kind = dp) :: sqabc, sqa, sqb, sqc
        !
        sqabc = sqrt(valA - (2.*valB) + valC)
        sqa = sqrt(valA)
        sqb = sqrt(valB)
        sqc = sqrt(valC)
        valor = 0.

        valor = valor + &
                -log( (-valB + (sqa*sqc))/(valA - valB + (sqa*sqabc)) ) / (24.*sqa)

        valor = valor + &
                log( (valB + (sqa*sqc))/(valB - valC + (sqc*sqabc)) ) / (24.*sqc)         

        valor = valor + ((-sqa*sqc) + (sqa*sqabc)) &
                        / (24.*(sqa**3))

        valor = valor + ( (valA) + (valB) ) & 
                        *log( ( valB + (sqa*sqc) )/( -valA + valB + (sqa*sqabc) ) ) &
                        /(24.*(sqa**3))

        valor = valor + log( ( valA - valB + (sqa*sqabc) )/( valB - valC + (sqc*sqabc) ) ) &
                        /(24.*(sqabc))
                        
        valor = valor - log( ( valB + (sqa*sqc) )/( -valB + valC + (sqc*sqabc) ) ) &
                        /(12.*(sqc))


        valor = valor + ((sqa*sqabc) - (sqc*sqabc)) &
                        /(24.*(sqabc**3))

        valor = valor + ( (valA) - (3.*valB) + (2.*valC) ) & 
                        *log( ( -valB + valC + (sqc*sqabc) )/( -valA + valB + (sqa*sqabc) ) ) &
                        /(24.*(sqabc**3))
                

    end function eqnAnalytic1


    function eqnAnalytic0(valA, valB, valC) result(valor)
        !dummy
        real (kind = dp), intent(in) :: valA, valB, valC
        real (kind = dp) :: valor
        !

        !local
        real (kind = dp) :: sqabc, sqa, sqb, sqc
        !
        sqabc = sqrt(valA - (2.*valB) + valC)
        sqa = sqrt(valA)
        sqb = sqrt(valB)
        sqc = sqrt(valC)
        valor = 0.


        valor = valor + log( ( (valA - valB + (sqa*sqabc))*(valB + (sqa*sqc)) )/( (-valB + (sqa*sqc))*(-valA + valB +  &
            (sqa*sqabc)) ) ) &
                        /(6.*(sqa))

        valor = valor + log( ( (valB + (sqa*sqc))*(-valB + valC + (sqc*sqabc)) )/( (valB - valC + (sqc*sqabc))*(-valB +  &
            (sqa*sqc)) ) ) &
                        /(6.*(sqc))

        valor = valor + log( ( (valA - valB + (sqa*sqabc))*(-valB + valC + (sqc*sqabc)) )/( (valB - valC + (sqc*sqabc))* &
            (-valA + valB + (sqa*sqabc)) ) ) &
                        /(6.*(sqabc))
                

    end function eqnAnalytic0




    function gradGreen(robs, rsource) result (gg)
        !dummy
        complex (kind = dp), dimension(3) :: gg
        real (kind = dp), dimension(3), intent(in) :: robs, rsource
        real (kind = dp) :: dist
        !
        !local
        complex (kind = dp) :: UNO
        !
        UNO = (jj)/(jj)
        dist = sqrt(dot_product(robs - rsource, robs - rsource))

        gg = (robs - rsource)*(1 + jj*getkappa()*dist)*(exp(-jj*getkappa()*dist))/(dist**3)

    end function gradGreen

	subroutine setZ_MatxVec_MoM(matZ)
		complex ( kind = dp ), intent(in), target, dimension(:,:) :: matZ

        ZMoM_MxV => matZ
	end subroutine

    function MatxVec_MoM(vecX) result(res)
    	complex (kind = dp), intent(in), dimension(:) :: vecX
		complex (kind = dp), dimension(size(vecX,1)) :: res
		res = matmul(ZMoM_MxV,vecX)
    end function MatxVec_MoM

end module modMoM
