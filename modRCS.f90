module modRCS
    use modGlobalParam
    use modGlobalMetod
    use modMoM
    implicit none

    private
    !Variables globales
    public :: inicializar_RCS, sigmaRCS
    !real ( kind = dp ), dimension (3) :: r_obs
    complex ( kind = dp ), pointer,  dimension (:) :: Vi
    real ( kind = dp ), pointer, dimension(:) :: e_long
    real ( kind = dp ), pointer, dimension(:,:) :: e_centro
    real ( kind = dp ), pointer, dimension (:,:) :: t_baric
    integer ( kind = il ), pointer, dimension(:,:) :: e_t
    integer ( kind = il ), pointer :: num_e
    integer ( kind = il ), pointer :: num_t

    complex (kind = dp), dimension(3) :: polariz
    complex (kind = dp) :: ctteOnda

    real ( kind = dp ),  pointer, dimension (:,:) :: p_coord
    integer ( kind = il ), pointer, dimension (:,:) :: t_p, e_po
contains



    subroutine inicializar_RCS (arg_centros_comun, arg_Vi, arg_long_l, arg_baric_t, arg_t_comun, arg_num_lcom, arg_num_t,&
     arg_polariz, arg_ctteOnda, arg_p_coord, arg_t_p, arg_e_po)

        !dummy
        real ( kind = dp ),  dimension(:,:), target :: arg_centros_comun
        !real ( kind = dp ), dimension (3), intent(in) :: arg_r_obs
        complex ( kind = dp ),  dimension (:), intent(in), target :: arg_Vi
        real ( kind = dp ),  dimension(:), intent(in), target :: arg_long_l
        real ( kind = dp ),  dimension (:,:),intent(in), target :: arg_baric_t
        integer ( kind = il ),  dimension(:,:),intent(in), target :: arg_t_comun
        integer ( kind = il ), intent(in), target :: arg_num_lcom
        integer ( kind = il ), intent(in), target :: arg_num_t

        complex (kind = dp), dimension(3) :: arg_polariz
        complex (kind = dp) :: arg_ctteOnda

        real ( kind = dp ),  target, dimension (:,:) :: arg_p_coord
        integer ( kind = il ), target, dimension (:,:) :: arg_t_p, arg_e_po
        !!!!!!!

        !r_obs = arg_r_obs
        num_e => arg_num_lcom
        num_t => arg_num_t
        !allocate (Vi(num_e))
        Vi => arg_Vi
        !allocate (e_centro(3, num_e))
        e_centro => arg_centros_comun
        !allocate (e_long(num_e))
        e_long => arg_long_l
        !allocate (t_baric(3,num_t))
        t_baric => arg_baric_t
        !allocate (e_t(2,num_e))
        e_t => arg_t_comun

        p_coord => arg_p_coord
        t_p => arg_t_p
        e_po => arg_e_po

        polariz = arg_polariz
        ctteOnda = arg_ctteOnda

    end subroutine


    function calc_E_old (r_obs, roo) result (E)

        !dummy
        real ( kind = dp ), dimension (3), intent(in) :: r_obs
        complex (kind = dp), dimension(3) :: E
        real (kind = dp), dimension (3), intent(in) :: roo
        !!!!!!!!
        complex (kind = dp), dimension(3) :: acumE
        complex (kind = dp), dimension(3) :: m
        complex (kind = dp), dimension(3) :: MM
        real (kind = dp), dimension(3) :: rm
        real (kind = dp), dimension(3) :: rrm
        integer (kind = il) :: i

        acumE = 0


        do i = 1,num_e

            rm = 0.5*(  t_baric(:, e_t(1,i)) + t_baric(:, e_t(2,i))   )
            m  = (  t_baric(:, e_t(2,i)) - t_baric(:, e_t(1,i))   )*Vi(i)*e_long(i)
            rrm = r_obs - rm
            !MM = producto_escalar(m,rrm)*rrm/(modulov(rrm)**2)
            MM = sum(m*rrm)*rrm/(dot_product(rrm,rrm))
            acumE=acumE + ( green(r_obs,rm)*(MM(:)-m(:)) )
        !print*, acumE
        end do

        E = acumE*jj*w*mu/(4*pi)
        !print*, E
    end function calc_E_old



    function calc_E_noAprox (r_obs, roo) result (E)

        !dummy
        real ( kind = dp ), dimension (3), intent(in) :: r_obs
        complex (kind = dp), dimension(3) :: E
        real (kind = dp), dimension (3), intent(in) :: roo
        !!!!!!!!
        complex (kind = dp), dimension(3) :: acumE, potA
        integer (kind = il) :: fb, i_t, tri_i
        real (kind = dp), dimension(3) :: V1, V2, V3, Vo, rho, rtest, ln, AR
        real (kind = dp), dimension(:,:), allocatable :: GL_pointsweights_test

        integer (kind = il) :: GL_numberpoints, itest

        GL_numberpoints = 7        
        
        allocate(GL_pointsweights_test(GL_numberpoints, 4))

        acumE(:) = 0.


        do fb = 1,num_e
            ln = e_long(fb)
            do i_t = 1, 2
                tri_i = e_t(i_t, fb)
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
                    AR = r_obs - rtest
                    AR = AR/norm2(AR)
                    if ( i_t == 2 ) then
                        !triangle -
                        rho = (-1._dp)*rho
                    end if
                    rho = rho*0.5_dp*ln

                    potA = rho*Vi(fb)*green(r_obs, rtest)*GL_pointsweights_test(itest,4)
                    acumE = acumE + (potA - AR*sum(AR*potA))

                end do                
            end do
        end do
        
        E = acumE*(-jj)*w*mu/(4*pi)

        deallocate(GL_pointsweights_test)
        !print*, E
    end function calc_E_noAprox


!     function calc_E_analytical (r_obs, roo) result (E)

!         !dummy
!         real ( kind = dp ), dimension (3), intent(in) :: r_obs, roo
!         complex (kind = dp), dimension(3) :: E
!         !local
!         complex (kind = dp), dimension(3) :: Integral, acumE, potA
!         real (kind = dp), dimension(3) :: r_vector, r_uni, vector_s
!         integer (kind = il) :: fb, i_t, tri_i, iVo
!         real (kind = dp) :: area, ln
!         !
!         acumE = 0.
!         potA = 0.
!         !r uni
!         r_vector = r_obs - roo
!         r_uni = r_obs/norm2(r_obs)

!         do fb = 1,num_e
!             ln = e_long(fb)
!             do i_t = 1, 2                
!                 !triangle ind
!                 tri_i = e_t(i_t, fb)
!                 !area
!                 area = t_area(tri_i)
!                 !opposite vertex
!                 iVo = e_po(i_t,fb)

!                 Integral = radIntegral(r_uni, tri_i, iVo)
!                 if ( i_t == 2 ) then
!                     Integral = (-1._dp)*Integral
!                 endif

!                 potA =  potA + Integral*Vi(fb)*ln/(2._dp*area)

!             end do      
!         end do        

!         potA = potA*exp(-jj*getkappa()*norm2(r_vector)) /norm2(r_vector)

!         acumE = potA - r_uni*sum(r_uni*potA)

!         E = acumE*(-jj)*w*mu/(4._dp*pi)
!         !print*, E
!     end function calc_E_analytical


    function calc_E (r_obs, roo) result (E)

        !dummy
        real ( kind = dp ), dimension (3), intent(in) :: r_obs
        complex (kind = dp), dimension(3) :: E
        real (kind = dp), dimension (3), intent(in) :: roo
        !!!!!!!!
        complex (kind = dp), dimension(3) :: acumE, potA
        integer (kind = il) :: fb, i_t, tri_i
        real (kind = dp), dimension(3) :: V1, V2, V3, Vo, rho, rtest, ln, AR, r_source
        real (kind = dp), dimension(:,:), allocatable :: GL_pointsweights_test
        real (kind = dp) :: Rmag

        integer (kind = il) :: GL_numberpoints, itest

        GL_numberpoints = 7        
        
        allocate(GL_pointsweights_test(GL_numberpoints, 4))

        acumE(:) = 0.
        potA(:) = 0.
        AR = r_obs - roo
        Rmag = norm2(AR)
        AR = AR/Rmag

        do fb = 1,num_e
            ln = e_long(fb)
            do i_t = 1, 2
                tri_i = e_t(i_t, fb)
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
                    r_source = rtest - roo
                    potA = potA + rho*Vi(fb)*exp(jj*getkappa()*dot_product(r_source,AR))*GL_pointsweights_test(itest,4)
                    

                end do                
            end do
        end do
        
        potA = potA*exp(-jj*getkappa()*Rmag)/Rmag

        acumE = (potA - AR*sum(AR*potA))

        E = acumE*(-jj)*w*mu/(4*pi)

        deallocate(GL_pointsweights_test)
        !print*, E
    end function calc_E


    function sigmaRCS(theta, phi) result (sigma)
        !dummy
        real (kind = dp), intent(in) :: theta, phi
        real (kind = dp) :: sigma
        !
        !local
        real (kind = dp) :: r_zl
        complex (kind = dp), dimension(3) :: E
        E = calc_RCS(theta, phi, r_zl)

        sigma = dot_product(E,E)

        sigma = sigma*4*pi*(r_zl**2)/(dot_product(polariz*ctteOnda,polariz*ctteOnda))

    end function sigmaRCS



    function calc_RCS(theta, phi, r_zl) result(ERCS)
        !dummy
        real (kind = dp), intent(in) :: theta, phi
        complex (kind = dp), dimension(3) :: ERCS
        real (kind = dp), intent(out) :: r_zl
        !!!!!
        !local
        real (kind = dp), dimension (3) :: roo
        integer (kind = il), dimension(6) :: limits
        real (kind = dp) :: D, maxX, minX, maxY, minY, maxZ, minZ, rzl, lambda
        real (kind = dp), dimension(3) :: r_obs
        !

        limits = limites_dispersor(e_centro, num_e)

        maxX = e_centro(1,limits(1))
        minX = e_centro(1,limits(2))
        maxY = e_centro(2,limits(3))
        minY = e_centro(2,limits(4))
        maxZ = e_centro(3,limits(5))
        minZ = e_centro(3,limits(6))
        lambda = 2*pi/(w*sqrt(mu*epsi))
        D = sqrt(    ((maxX-minX)**2) + ((maxY-minY)**2) + ((maxZ-minZ)**2)   )

        rzl = 2*(D**2)/lambda
        if (50*D>rzl) then
            rzl = 50*D
        end if
        if (20*lambda>rzl) then
            rzl = 20*lambda
        end if
        rzl = rzl*10 !multiplicador de zona lejana
        r_zl = rzl
        roo(1) = 0.5*(maxX+minX)
        roo(2) = 0.5*(maxY+minY)
        roo(3) = 0.5*(maxZ+minZ)

        r_obs(1) = roo(1)+ ( rzl*sin(theta)*cos(phi) )
        r_obs(2) = roo(2)+ ( rzl*sin(theta)*sin(phi) )
        r_obs(3) = roo(3)+ ( rzl*cos(theta) )

        !return
        ERCS = calc_E(r_obs, roo)

    end function calc_RCS







end module modRCS
