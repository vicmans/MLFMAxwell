module modReadStl
use modGlobalParam
    implicit none

contains

    subroutine ch_cap ( c )
        character c
        integer ( kind= il ) itemp
        itemp = ichar ( c )
        if ( 97 <= itemp .and. itemp <= 122 ) then
            c = char ( itemp - 32 )
        end if
        return
    end subroutine ch_cap

    function ch_eqi ( c1, c2 ) result (ch_eqi_r)
        character c1
        character c1_cap
        character c2
        character c2_cap
        logical ch_eqi_r
        c1_cap = c1
        c2_cap = c2
        call ch_cap ( c1_cap )
        call ch_cap ( c2_cap )
        if ( c1_cap == c2_cap ) then
            ch_eqi_r = .true.
        else
            ch_eqi_r = .false.
        end if
        return
    end function ch_eqi

    subroutine ch_to_digit ( c, digit )
        character c
        integer ( kind= il ) digit
        if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
            digit = ichar ( c ) - 48
        else if ( c == ' ' ) then
            digit = 0
        else
            digit = -1
        end if
        return
    end subroutine ch_to_digit

    subroutine get_unit ( iunit )
        integer ( kind= il ) i
        integer ( kind= il ) ios
        integer ( kind= il ) iunit
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
    end subroutine get_unit

    subroutine r8vec_cross_3d ( v1, v2, v3 )
        integer ( kind= il ), parameter :: dim_num = 3
        real ( kind= dp ) v1(dim_num)
        real ( kind= dp ) v2(dim_num)
        real ( kind= dp ) v3(dim_num)
        v3(1) = v1(2) * v2(3) - v1(3) * v2(2)
        v3(2) = v1(3) * v2(1) - v1(1) * v2(3)
        v3(3) = v1(1) * v2(2) - v1(2) * v2(1)
        return
    end subroutine r8vec_cross_3d

    subroutine s_cat ( s1, s2, s3 )
        character ( len = * ) s1
        character ( len = * ) s2
        character ( len = * ) s3
        if ( s1 == ' ' .and. s2 == ' ' ) then
            s3 = ' '
        else if ( s1 == ' ' ) then
            s3 = s2
        else if ( s2 == ' ' ) then
            s3 = s1
        else
            s3 = trim ( s1 ) // trim ( s2 )
        end if
        return
    end subroutine s_cat

    function s_eqi ( s1, s2 ) result (s_eqi_r)
        logical  s_eqi_r
        character c1
        character c2
        integer ( kind= il ) i
        integer ( kind= il ) len1
        integer ( kind= il ) len2
        integer ( kind= il ) lenc

        character ( len = * ) s1
        character ( len = * ) s2
        len1 = len ( s1 )
        len2 = len ( s2 )
        lenc = min ( len1, len2 )
        s_eqi_r = .false.
        do i = 1, lenc
            c1 = s1(i:i)
            c2 = s2(i:i)
            call ch_cap ( c1 )
            call ch_cap ( c2 )
            if ( c1 /= c2 ) then
                return
            end if
        end do
        do i = lenc + 1, len1
            if ( s1(i:i) /= ' ' ) then
                return
            end if
        end do
        do i = lenc + 1, len2
            if ( s2(i:i) /= ' ' ) then
                return
            end if
        end do
        s_eqi_r = .true.
        return
    end function s_eqi

    subroutine s_to_r8 ( s, dval, ierror, length )
        logical ch_eqi_r
        character c
        real ( kind= dp ) dval
        integer ( kind= il ) ierror
        integer ( kind= il ) ihave
        integer ( kind= il ) isgn
        integer ( kind= il ) iterm
        integer ( kind= il ) jbot
        integer ( kind= il ) jsgn
        integer ( kind= il ) jtop
        integer ( kind= il ) length
        integer ( kind= il ) nchar
        integer ( kind= il ) ndig
        real ( kind= dp ) rbot
        real ( kind= dp ) rexp
        real ( kind= dp ) rtop
        character ( len = * ) s
        nchar = len_trim ( s )
        ierror = 0
        dval = 0.0D+00
        length = -1
        isgn = 1
        rtop = 0
        rbot = 1
        jsgn = 1
        jtop = 0
        jbot = 1
        ihave = 1
        iterm = 0
        do
            length = length + 1
            if ( nchar < length+1 ) then
                exit
            end if
            c = s(length+1:length+1)
            if ( c == ' ' ) then
                if ( ihave == 2 ) then
                else if ( ihave == 6 .or. ihave == 7 ) then
                    iterm = 1
                else if ( 1 < ihave ) then
                    ihave = 11
                end if
            else if ( c == ',' .or. c == ';' ) then
                if ( ihave /= 1 ) then
                    iterm = 1
                    ihave = 12
                    length = length + 1
                end if
            else if ( c == '-' ) then
                if ( ihave == 1 ) then
                    ihave = 2
                    isgn = -1
                else if ( ihave == 6 ) then
                    ihave = 7
                    jsgn = -1
                else
                    iterm = 1
                end if
            else if ( c == '+' ) then
                if ( ihave == 1 ) then
                    ihave = 2
                else if ( ihave == 6 ) then
                    ihave = 7
                else
                    iterm = 1
                end if
            else if ( c == '.' ) then
                if ( ihave < 4 ) then
                    ihave = 4
                else if ( 6 <= ihave .and. ihave <= 8 ) then
                    ihave = 9
                else
                    iterm = 1
                end if
            else if ( ch_eqi ( c, 'E' ) .or. ch_eqi ( c, 'D' ) ) then
                if ( ihave < 6 ) then
                    ihave = 6
                else
                    iterm = 1
                end if
            else if (  ihave < 11 .and. lle ( '0', c ) .and. lle ( c, '9' ) ) then
                if ( ihave <= 2 ) then
                    ihave = 3
                else if ( ihave == 4 ) then
                    ihave = 5
                else if ( ihave == 6 .or. ihave == 7 ) then
                    ihave = 8
                else if ( ihave == 9 ) then
                    ihave = 10
                end if
                call ch_to_digit ( c, ndig )
                if ( ihave == 3 ) then
                    rtop = 10.0D+00 * rtop + real ( ndig, kind= dp )
                else if ( ihave == 5 ) then
                    rtop = 10.0D+00 * rtop + real ( ndig, kind= dp )
                    rbot = 10.0D+00 * rbot
                else if ( ihave == 8 ) then
                    jtop = 10 * jtop + ndig
                else if ( ihave == 10 ) then
                    jtop = 10 * jtop + ndig
                    jbot = 10 * jbot
                end if
            else
                iterm = 1
            end if
            if ( iterm == 1 ) then
                exit
            end if
        end do
        if ( iterm /= 1 .and. length+1 == nchar ) then
            length = nchar
        end if
        if ( ihave == 1 .or. ihave == 2 .or. ihave == 6 .or. ihave == 7 ) then
            ierror = ihave
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'S_TO_R8 - Serious error!'
            write ( *, '(a)' ) '  Illegal or nonnumeric input:'
            write ( *, '(a)' ) '    ' // trim ( s )
            return
        end if
        if ( jtop == 0 ) then
            rexp = 1.0D+00
        else
            if ( jbot == 1 ) then
                rexp = 10.0D+00 ** ( jsgn * jtop )
            else
                rexp = 10.0D+00 ** ( real ( jsgn * jtop, kind= dp ) &
                    / real ( jbot, kind= dp ) )
            end if
        end if
        dval = real ( isgn, kind= dp ) * rexp * rtop / rbot
        return
    end subroutine s_to_r8

    function stla_check ( input_file_name ) result (stla_check_r)
        logical check
        logical done
        real ( kind= dp ) dval
        integer ( kind= il ) i
        integer ( kind= il ) ierror
        character ( len = * ) input_file_name
        integer ( kind= il ) ios
        integer ( kind= il ) iunit
        integer ( kind= il ) lchar
        logical s_eqi_r
        integer ( kind= il ) state
        logical stla_check_r
        character ( len = 255 ) text
        integer ( kind= il ) text_num
        integer ( kind= il ) vertex
        character ( len = 255 ) word1
        character ( len = 255 ) word2
        state = 0
        text_num = 0
        call get_unit ( iunit )
        open ( unit = iunit, file = input_file_name, status = 'old', iostat = ios )
        if ( ios /= 0 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'STLA_CHECK - Fatal error!'
            write ( *, '(a)' ) '  Could not open the file "' &
                // trim ( input_file_name ) // '".'
            stla_check_r = .false.
            return
        end if
        do
            read ( iunit, '(a)', iostat = ios ) text
            if ( ios /= 0 ) then
                if ( state /= 0 .and. state /= 1 ) then
                    write ( *, '(a)' ) ' '
                    write ( *, '(a)' ) 'STLA_CHECK - Fatal error!'
                    write ( *, '(a,i8)' ) '  File line number = ', text_num
                    write ( *, '(a)' ) '  End-of-file, but model not finished.'
                    stla_check_r = .false.
                    return
                end if
                exit
            end if
            text_num = text_num + 1
            done = .true.
            call word_next_read ( text, word1, done )
            if ( done ) then
                write ( *, '(a)' ) ' '
                write ( *, '(a)' ) 'STLA_CHECK - Fatal error!'
                write ( *, '(a,i8)' ) '  File line number = ', text_num
                write ( *, '(a)' ) '  No information on line.'
                stla_check_r = .false.
                return
            end if
            if ( s_eqi ( word1, 'END' ) ) then
                call word_next_read ( text, word2, done )
                if ( .not. s_eqi ( word2, 'FACET' ) .and. &
                    .not. s_eqi ( word2, 'LOOP' ) .and. &
                    .not. s_eqi ( word2, 'SOLID' ) ) then
                    write ( *, '(a)' ) ' '
                    write ( *, '(a)' ) 'STLA_CHECK - Fatal error!'
                    write ( *, '(a,i8)' ) '  File line number = ', text_num
                    write ( *, '(a)' ) '  The tag END was followed by an illegal '
                    write ( *, '(a)' ) '  word: "' // trim ( word2 ) // '", when expecting'
                    write ( *, '(a)' ) '  "FACET", "LOOP", or "SOLID".'
                    stla_check_r = .false.
                    return
                end if
                call s_cat ( word1, word2, word1 )
            else if ( s_eqi ( word1, 'FACET' ) ) then
                call word_next_read ( text, word2, done )
                if ( .not. s_eqi ( word2, 'NORMAL' ) ) then
                    write ( *, '(a)' ) ' '
                    write ( *, '(a)' ) 'STLA_CHECK - Fatal error!'
                    write ( *, '(a,i8)' ) '  File line number = ', text_num
                    write ( *, '(a)' ) '  The tag FACET was followed by an illegal '
                    write ( *, '(a)' ) '  word: "' // trim ( word2 ) // '", when expecting'
                    write ( *, '(a)' ) '  "NORMAL".'
                    stla_check_r = .false.
                    return
                end if
                call s_cat ( word1, word2, word1 )
            else if ( s_eqi ( word1, 'OUTER' ) ) then
                call word_next_read ( text, word2, done )
                if ( .not. s_eqi ( word2, 'LOOP' ) ) then
                    write ( *, '(a)' ) ' '
                    write ( *, '(a)' ) 'STLA_CHECK - Fatal error!'
                    write ( *, '(a,i8)' ) '  File line number = ', text_num
                    write ( *, '(a)' ) '  The tag OUTER was followed by an illegal '
                    write ( *, '(a)' ) '  word: "' // trim ( word2 ) // '", when expecting'
                    write ( *, '(a)' ) '  "LOOP".'
                    stla_check_r = .false.
                    return
                end if
                call s_cat ( word1, word2, word1 )
            end if
            if ( s_eqi ( word1, 'SOLID' ) ) then
                if ( state /= 0 ) then
                    write ( *, '(a)' ) ' '
                    write ( *, '(a)' ) 'STLA_CHECK - Fatal error!'
                    write ( *, '(a,i8)' ) '  File line number = ', text_num
                    write ( *, '(a)' ) '  A new SOLID statement was encountered, but we'
                    write ( *, '(a)' ) '  have not finished processing the current solid.'
                    stla_check_r = .false.
                    return
                end if
                state = 1
            else if ( s_eqi ( word1, 'ENDSOLID' ) ) then
                if ( state /= 1 ) then
                    write ( *, '(a)' ) ' '
                    write ( *, '(a)' ) 'STLA_CHECK - Fatal error!'
                    write ( *, '(a,i8)' ) '  File line number = ', text_num
                    write ( *, '(a)' ) '  An END SOLID statement was encountered, but'
                    write ( *, '(a)' ) '  either we have not begun a solid at all, or we'
                    write ( *, '(a)' ) '  are not at an appropriate point to finish the'
                    write ( *, '(a)' ) '  current solid.'
                    stla_check_r = .false.
                    return
                end if
                state = 0
            else if ( s_eqi ( word1, 'FACETNORMAL' ) ) then
                if ( state /= 0 .and. state /= 1 ) then
                    write ( *, '(a)' ) ' '
                    write ( *, '(a)' ) 'STLA_CHECK - Fatal error!'
                    write ( *, '(a,i8)' ) '  File line number = ', text_num
                    write ( *, '(a)' ) '  Model not in right state for FACET.'
                    stla_check_r = .false.
                    return
                end if
                state = 2
                do i = 1, 3
                    call word_next_read ( text, word2, done )
                    if ( done ) then
                        write ( *, '(a)' ) ' '
                        write ( *, '(a)' ) 'STLA_CHECK - Fatal error!'
                        write ( *, '(a,i8)' ) '  File line number = ', text_num
                        write ( *, '(a)' ) '  End of information while reading a component'
                        write ( *, '(a)' ) '  of the normal vector.'
                        stla_check_r = .false.
                        return
                    end if
                    call s_to_r8 ( word2, dval, ierror, lchar )
                    if ( ierror /= 0 ) then
                        write ( *, '(a)' ) ' '
                        write ( *, '(a)' ) 'STLA_CHECK - Fatal error!'
                        write ( *, '(a,i8)' ) '  File line number = ', text_num
                        write ( *, '(a)' ) &
                            '  Error while reading a component of the normal vector.'
                        stla_check_r = .false.
                        return
                    end if
                end do
            else if ( s_eqi ( word1, 'ENDFACET' ) ) then
                if ( state /= 2 ) then
                    write ( *, '(a)' ) ' '
                    write ( *, '(a)' ) 'STLA_CHECK - Fatal error!'
                    write ( *, '(a,i8)' ) '  File line number = ', text_num
                    write ( *, '(a)' ) '  Model not in right state for ENDFACET.'
                    stla_check_r = .false.
                    return
                end if
                state = 1
            else if ( s_eqi ( word1, 'OUTERLOOP' ) ) then
                if ( state /= 2 ) then
                    write ( *, '(a)' ) ' '
                    write ( *, '(a)' ) 'STLA_CHECK - Fatal error!'
                    write ( *, '(a,i8)' ) '  File line number = ', text_num
                    write ( *, '(a)' ) '  Model not in right state for OUTERLOOP.'
                    stla_check_r = .false.
                    return
                end if
                state = 3
                vertex = 0
            else if ( s_eqi ( word1, 'ENDLOOP' ) ) then
                if ( state /= 3 ) then
                    write ( *, '(a)' ) ' '
                    write ( *, '(a)' ) 'STLA_CHECK - Fatal error!'
                    write ( *, '(a,i8)' ) '  File line number = ', text_num
                    write ( *, '(a)' ) '  Model not in right state for ENDLOOP.'
                    stla_check_r = .false.
                    return
                end if
                state = 2
            else if ( s_eqi ( word1, 'VERTEX' ) ) then
                if ( state /= 3 ) then
                    write ( *, '(a)' ) ' '
                    write ( *, '(a)' ) 'STLA_CHECK - Fatal error!'
                    write ( *, '(a,i8)' ) '  File line number = ', text_num
                    write ( *, '(a)' ) '  Model not in right state for VERTEX.'
                    stla_check_r = .false.
                    return
                end if
                if ( 3 <= vertex ) then
                    write ( *, '(a)' ) ' '
                    write ( *, '(a)' ) 'STLA_CHECK - Fatal error!'
                    write ( *, '(a,i8)' ) '  File line number = ', text_num
                    write ( *, '(a)' ) '  More than 3 vertices specified for a face.'
                    stla_check_r = .false.
                    return
                end if
                do i = 1, 3
                    call word_next_read ( text, word2, done )
                    if ( done ) then
                        write ( *, '(a)' ) ' '
                        write ( *, '(a)' ) 'STLA_CHECK - Fatal error!'
                        write ( *, '(a,i8)' ) '  File line number = ', text_num
                        write ( *, '(a)' ) '  The value of a vertex coordinate is missing.'
                        stla_check_r = .false.
                        return
                    end if

                    call s_to_r8 ( word2, dval, ierror, lchar )

                    if ( ierror /= 0 ) then
                        write ( *, '(a)' ) ' '
                        write ( *, '(a)' ) 'STLA_CHECK - Fatal error!'
                        write ( *, '(a,i8)' ) '  File line number = ', text_num
                        write ( *, '(a)' ) '  The value of a vertex coordinate makes'
                        write ( *, '(a)' ) '  no sense.'
                        stla_check_r = .false.
                        return
                    end if

                end do

                vertex = vertex + 1

            else

                write ( *, '(a)' ) ' '
                write ( *, '(a)' ) 'STLA_CHECK - Fatal error!'
                write ( *, '(a,i8)' ) '  File line number = ', text_num
                write ( *, '(a)' ) '  Unrecognized line in file.'
                stla_check_r = .false.
                return

            end if

        end do
        !
        !  Close the file.
        !
        close ( unit = iunit )

        stla_check_r = .true.

        return
    end function stla_check

    subroutine stla_face_node_print ( face_num, face_node )
        integer ( kind= il ) face_num
        integer ( kind= il ) face
        integer ( kind= il ) face_node(3,face_num)

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '    Face         Nodes'
        write ( *, '(a)' ) ' '

        do face = 1, face_num

            write ( *, '(2x,i8,3(2x,i8))' ) face, face_node(1:3,face)

        end do

        return
    end subroutine stla_face_node_print

    subroutine stla_face_normal_compute ( node_num, face_num, node_xyz, &
        face_node, face_normal )

        integer ( kind= il ) face_num
        integer ( kind= il ) node_num

        integer ( kind= il ) face
        integer ( kind= il ) face_node(3,face_num)
        real ( kind= dp ) face_normal(3,face_num)
        integer ( kind= il ) n1
        integer ( kind= il ) n2
        integer ( kind= il ) n3
        real ( kind= dp ) node_xyz(3,node_num)
        real ( kind= dp ) norm
        real ( kind= dp ) v1(3)
        real ( kind= dp ) v2(3)

        do face = 1, face_num

            n1 = face_node(1,face)
            n2 = face_node(2,face)
            n3 = face_node(3,face)

            v1(1:3) = node_xyz(1:3,n2) - node_xyz(1:3,n1)
            v2(1:3) = node_xyz(1:3,n3) - node_xyz(1:3,n1)

            call r8vec_cross_3d ( v1, v2, face_normal(1:3,face ) )

            norm = sqrt ( sum ( ( face_normal(1:3,face) )**2 ) )

            if ( norm /= 0.0D+00 ) then
                face_normal(1:3,face) = face_normal(1:3,face) / norm
            end if

        end do

        return
    end subroutine stla_face_normal_compute

    subroutine stla_face_normal_print ( face_num, face_normal )

        integer ( kind= il ) face_num

        integer ( kind= il ) face
        real ( kind= dp ) face_normal(3,face_num)

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '    Face         Normal Vectors'
        write ( *, '(a)' ) ' '

        do face = 1, face_num

            write ( *, '(2x,i8,3(2x,g14.6))' ) face, face_normal(1:3,face)

        end do

        return
    end subroutine stla_face_normal_print

    subroutine stla_node_xyz_print ( node_num, node_xyz )

        integer ( kind= il ) node_num

        integer ( kind= il ) node
        real ( kind= dp ) node_xyz(3,node_num)

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '    Node         Coordinates'
        write ( *, '(a)' ) ' '

        do node = 1, node_num

            write ( *, '(2x,i8,3(2x,g14.6))' ) node, node_xyz(1:3,node)

        end do

        return
    end subroutine stla_node_xyz_print

    subroutine stla_read ( input_file_name, node_num, face_num, node_xyz, &
        face_node, face_normal, ierror )

        integer ( kind= il ) face_num
        integer ( kind= il ) node_num

        logical done
        real ( kind= dp ) dval
        integer ( kind= il ) face
        integer ( kind= il ) face_node(3,face_num)
        real ( kind= dp ) face_normal(3,face_num)
        integer ( kind= il ) i
        integer ( kind= il ) ierror
        character ( len = * ) input_file_name
        integer ( kind= il ) ios
        integer ( kind= il ) iunit
        integer ( kind= il ) lchar
        integer ( kind= il ) node
        real ( kind= dp ) node_xyz(3,node_num)
        logical s_eqi_r
        integer ( kind= il ) state
        real ( kind= dp ) temp(3)
        character ( len = 255 ) text
        integer ( kind= il ) text_num
        integer ( kind= il ) vertex
        character ( len = 255 ) word1
        character ( len = 255 ) word2

        ierror = 0
        state = 0
        text_num = 0
        face = 0
        node = 0
        !
        !  Open the file.
        !
        call get_unit ( iunit )

        open ( unit = iunit, file = input_file_name, status = 'old', iostat = ios )

        if ( ios /= 0 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'STLA_READ - Fatal error!'
            write ( *, '(a)' ) '  Could not open the file "' // &
                trim ( input_file_name ) // '".'
            ierror = 1
            return
        end if
        !
        !  Read the next line of text.
        !
        do

            read ( iunit, '(a)', iostat = ios ) text

            if ( ios /= 0 ) then
                if ( state /= 0 .and. state /= 1 ) then
                    write ( *, '(a)' ) ' '
                    write ( *, '(a)' ) 'STLA_READ - Warning.'
                    write ( *, '(a)' ) '  End-of-file, but model not finished.'
                    write ( *, '(a,i8)' ) '  File line number = ', text_num
                    ierror = 1
                    return
                end if
                exit
            end if

            text_num = text_num + 1
            done = .true.
            !
            !  Read the first word in the line.
            !
            call word_next_read ( text, word1, done )

            if ( s_eqi ( word1, 'END' ) .or. &
                s_eqi ( word1, 'FACET' ) .or. &
                s_eqi ( word1, 'OUTER' ) ) then

                call word_next_read ( text, word2, done )
                call s_cat ( word1, word2, word1 )

            end if

            if ( s_eqi ( word1, 'SOLID' ) ) then

                if ( state /= 0 ) then
                    write ( *, '(a)' ) ' '
                    write ( *, '(a)' ) 'STLA_READ - Warning!'
                    write ( *, '(a)' ) '  Model not in right state for SOLID.'
                    write ( *, '(a,i8)' ) '  File line number = ', text_num
                    ierror = 1
                    return
                end if

                state = 1

            else if ( s_eqi ( word1, 'ENDSOLID' ) ) then

                if ( state /= 1 ) then
                    write ( *, '(a)' ) ' '
                    write ( *, '(a)' ) 'STLA_READ - Warning!'
                    write ( *, '(a)' ) '  Model not in right state for ENDSOLID.'
                    write ( *, '(a,i8)' ) '  File line number = ', text_num
                    ierror = 1
                    return
                end if

                state = 0

            else if ( s_eqi ( word1, 'FACETNORMAL' ) ) then

                if ( state /= 0 .and. state /= 1 ) then
                    write ( *, '(a)' ) ' '
                    write ( *, '(a)' ) 'STLA_READ - Warning!'
                    write ( *, '(a)' ) '  Model not in right state for FACET.'
                    write ( *, '(a,i8)' ) '  File line number = ', text_num
                    ierror = 1
                    return
                end if

                state = 2
                face = face + 1

                if ( face_num < face ) then
                    write ( *, '(a)' ) ' '
                    write ( *, '(a)' ) 'STLA_READ - Warning!'
                    write ( *, '(a)' ) '  More faces being read than expected.'
                    write ( *, '(a,i8)' ) '  File line number = ', text_num
                    ierror = 1
                    return
                end if

                do i = 1, 3
                    face_normal(i,face) = 0.0D+00
                    call word_next_read ( text, word2, done )
                    if ( .not. done ) then
                        call s_to_r8 ( word2, dval, ierror, lchar )
                        if ( ierror == 0 ) then
                            face_normal(i,face) = dval
                        end if
                    end if
                end do

            else if ( s_eqi ( word1, 'ENDFACET' ) ) then

                if ( state /= 2 ) then
                    write ( *, '(a)' ) ' '
                    write ( *, '(a)' ) 'STLA_READ - Warning!'
                    write ( *, '(a)' ) '  Model not in right state for ENDFACET.'
                    write ( *, '(a,i8)' ) '  File line number = ', text_num
                    ierror = 1
                    return
                end if

                state = 1

            else if ( s_eqi ( word1, 'OUTERLOOP' ) ) then

                if ( state /= 2 ) then
                    write ( *, '(a)' ) ' '
                    write ( *, '(a)' ) 'STLA_READ - Warning!'
                    write ( *, '(a)' ) '  Model not in right state for OUTERLOOP.'
                    write ( *, '(a,i8)' ) '  File line number = ', text_num
                    ierror = 1
                    return
                end if

                state = 3
                vertex = 0

            else if ( s_eqi ( word1, 'ENDLOOP' ) ) then

                if ( state /= 3 ) then
                    write ( *, '(a)' ) ' '
                    write ( *, '(a)' ) 'STLA_READ - Warning!'
                    write ( *, '(a)' ) '  Model not in right state for ENDLOOP.'
                    write ( *, '(a,i8)' ) '  File line number = ', text_num
                    ierror = 1
                    return
                end if

                state = 2

            else if ( s_eqi ( word1, 'VERTEX' ) ) then

                if ( state /= 3 ) then
                    write ( *, '(a)' ) ' '
                    write ( *, '(a)' ) 'STLA_READ - Warning!'
                    write ( *, '(a)' ) '  Model not in right state for VERTEX.'
                    write ( *, '(a,i8)' ) '  File line number = ', text_num
                    ierror = 1
                    return
                end if

                if ( 3 <= vertex ) then
                    write ( *, '(a)' ) ' '
                    write ( *, '(a)' ) 'STLA_READ - Warning!'
                    write ( *, '(a)' ) '  Too many vertices for face.'
                    write ( *, '(a,i8)' ) '  File line number = ', text_num
                    ierror = 1
                    return
                end if

                do i = 1, 3
                    call word_next_read ( text, word2, done )
                    call s_to_r8 ( word2, dval, ierror, lchar )
                    temp(i) = dval
                end do

                if ( node_num <= node ) then
                    write ( *, '(a)' ) ' '
                    write ( *, '(a)' ) 'STLA_READ - Warning!'
                    write ( *, '(a)' ) '  More nodes being read than expected.'
                    write ( *, '(a,i8)' ) '  File line number = ', text_num
                    ierror = 1
                    return
                end if

                node = node + 1
                node_xyz(1:3,node) = temp(1:3)

                vertex = vertex + 1
                face_node(vertex,face) = node

            else

                write ( *, '(a)' ) ' '
                write ( *, '(a)' ) 'STLA_READ - Warning!'
                write ( *, '(a)' ) '  Unrecognized line in file.'
                write ( *, '(a,i8)' ) '  File line number = ', text_num
                ierror = 1
                return

            end if

        end do
        !
        !  Close the file.
        !
        close ( unit = iunit )

        return
    end subroutine stla_read

    subroutine stla_size ( input_file_name, solid_num, node_num, face_num, &
        text_num )
        logical check
        logical done
        real ( kind= dp ) dval
        integer ( kind= il ) face_num
        integer ( kind= il ) i
        integer ( kind= il ) ierror
        character ( len = * ) input_file_name
        integer ( kind= il ) ios
        integer ( kind= il ) iunit
        integer ( kind= il ) lchar
        integer ( kind= il ) node_num
        logical s_eqi_r
        integer ( kind= il ) solid_num
        integer ( kind= il ) state
        character ( len = 255 ) text
        integer ( kind= il ) text_num
        integer ( kind= il ) vertex
        character ( len = 255 ) word1
        character ( len = 255 ) word2

        ierror = 0

        state = 0
        text_num = 0

        solid_num = 0
        node_num = 0
        face_num = 0
        !
        !  Open the file.
        !
        call get_unit ( iunit )

        open ( unit = iunit, file = input_file_name, status = 'old', iostat = ios )

        if ( ios /= 0 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'STLA_SIZE - Fatal error!'
            write ( *, '(a)' ) '  Could not open the file "' &
                // trim ( input_file_name ) // '".'
            ierror = 1
            return
        end if
        !
        !  Read the next line of text.
        !
        do

            read ( iunit, '(a)', iostat = ios ) text

            if ( ios /= 0 ) then
                if ( state /= 0 .and. state /= 1 ) then
                    return
                end if
                exit
            end if

            text_num = text_num + 1

            done = .true.
            !
            !  Read the first word in the line.
            !
            call word_next_read ( text, word1, done )

            if ( done ) then
                return
            end if

            if ( s_eqi ( word1, 'END' ) ) then

                call word_next_read ( text, word2, done )

                if ( .not. s_eqi ( word2, 'FACET' ) .and. &
                    .not. s_eqi ( word2, 'LOOP' ) .and. &
                    .not. s_eqi ( word2, 'SOLID' ) ) then
                    return
                end if

                call s_cat ( word1, word2, word1 )

            else if ( s_eqi ( word1, 'FACET' ) ) then

                call word_next_read ( text, word2, done )

                if ( .not. s_eqi ( word2, 'NORMAL' ) ) then
                    return
                end if

                call s_cat ( word1, word2, word1 )

            else if ( s_eqi ( word1, 'OUTER' ) ) then

                call word_next_read ( text, word2, done )

                if ( .not. s_eqi ( word2, 'LOOP' ) ) then
                    return
                end if

                call s_cat ( word1, word2, word1 )

            end if

            if ( s_eqi ( word1, 'SOLID' ) ) then

                if ( state /= 0 ) then
                    return
                end if

                state = 1

            else if ( s_eqi ( word1, 'ENDSOLID' ) ) then

                if ( state /= 1 ) then
                    return
                end if

                state = 0

                solid_num = solid_num + 1

            else if ( s_eqi ( word1, 'FACETNORMAL' ) ) then

                if ( state /= 0 .and. state /= 1 ) then
                    return
                end if

                state = 2

                do i = 1, 3

                    call word_next_read ( text, word2, done )

                    if ( done ) then
                        return
                    end if

                    call s_to_r8 ( word2, dval, ierror, lchar )

                    if ( ierror /= 0 ) then
                        return
                    end if

                end do

            else if ( s_eqi ( word1, 'ENDFACET' ) ) then

                if ( state /= 2 ) then
                    return
                end if

                state = 1

                face_num = face_num + 1

            else if ( s_eqi ( word1, 'OUTERLOOP' ) ) then

                if ( state /= 2 ) then
                    return
                end if

                state = 3
                vertex = 0

            else if ( s_eqi ( word1, 'ENDLOOP' ) ) then

                if ( state /= 3 ) then
                    return
                end if

                state = 2

            else if ( s_eqi ( word1, 'VERTEX' ) ) then

                if ( state /= 3 ) then
                    return
                end if

                if ( 3 <= vertex ) then
                    write ( *, '(a)' ) ' '
                    write ( *, '(a)' ) 'STLA_SIZE - Fatal error!'
                    write ( *, '(a)' ) '  Too many vertices for a face.'
                    ierror = 1
                    return
                end if

                do i = 1, 3

                    call word_next_read ( text, word2, done )

                    if ( done ) then
                        return
                    end if

                    call s_to_r8 ( word2, dval, ierror, lchar )

                    if ( ierror /= 0 ) then
                        return
                    end if

                end do

                vertex = vertex + 1
                node_num = node_num + 1

            else

                return

            end if

        end do
        !
        !  Close the file.
        !
        close ( unit = iunit )

        return
    end subroutine stla_size

    subroutine stla_size_print ( input_file_name, solid_num, node_num, face_num, &
        text_num )

        integer ( kind= il ) face_num
        character ( len = * ) input_file_name
        integer ( kind= il ) node_num
        integer ( kind= il ) solid_num
        integer ( kind= il ) text_num

        write ( *, '(a)'    ) ' '
        write ( *, '(a)'    ) '  Object sizes for STLA file "' // &
            trim ( input_file_name ) // '":'
        write ( *, '(a)'    ) ' '
        write ( *, '(a,i8)' ) '  Solids =                   ', solid_num
        write ( *, '(a,i8)' ) '  Nodes (may be repeated) =  ', node_num
        write ( *, '(a,i8)' ) '  Faces (triangular only) =  ', face_num
        write ( *, '(a)'    ) ' '
        write ( *, '(a,i8)' ) '  Number of lines of text =  ', text_num

        return
    end subroutine stla_size_print

    subroutine stla_write ( output_file_name, node_num, face_num, node_xyz, &
        face_node, face_normal )

        integer ( kind= il ) face_num
        integer ( kind= il ) node_num

        integer ( kind= il ) face
        integer ( kind= il ) face_node(3,face_num)
        real ( kind= dp ) face_normal(3,face_num)
        integer ( kind= il ) ios
        integer ( kind= il ) iunit
        integer ( kind= il ) node
        real ( kind= dp ) node_xyz(3,node_num)
        character ( len = * ) output_file_name
        integer ( kind= il ) text_num
        integer ( kind= il ) vertex

        text_num = 0
        !
        !  Open the file.
        !
        call get_unit ( iunit )

        open ( unit = iunit, file = output_file_name, status = 'replace', &
            iostat = ios )

        if ( ios /= 0 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'STLA_WRITE - Fatal error!'
            write ( *, '(a)' ) '  Could not open the file "' &
                // trim ( output_file_name ) // '".'
            stop
        end if

        write ( iunit, '(a)' ) 'solid MYSOLID'
        text_num = text_num + 1

        do face = 1, face_num

            write ( iunit, '(a,3(2x,g14.6))' ) '  facet normal', face_normal(1:3,face)
            text_num = text_num + 1

            write ( iunit, '(a)' ) '    outer loop'
            text_num = text_num + 1

            do vertex = 1, 3

                node = face_node(vertex,face)

                write ( iunit, '(a,2x,3(2x,g14.6))' ) '      vertex', node_xyz(1:3,node)
                text_num = text_num + 1

            end do

            write ( iunit, '(a)' ) '    end loop'
            text_num = text_num + 1
            write ( iunit, '(a)' ) '  end facet'
            text_num = text_num + 1

        end do

        write ( iunit, '(a)' ) 'end solid MYSOLID'
        text_num = text_num + 1

        close ( unit = iunit )

        return
    end subroutine stla_write

    subroutine timestamp ( )

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
    end subroutine timestamp

    subroutine word_next_read ( s, word, done )

        logical stla_check

        logical done
        integer ( kind= il ) ilo
        integer ( kind= il ), save :: lenc = 0
        integer ( kind= il ), save :: next = 1
        character ( len = * ) s
        character, parameter :: TAB = char ( 9 )
        character ( len = * ) word
        !
        !  We "remember" LENC and NEXT from the previous call.
        !
        !  An input value of DONE = TRUE signals a new line of text to examine.
        !
        if ( done ) then

            next = 1
            done = .false.
            lenc = len_trim ( s )

            if ( lenc <= 0 ) then
                done = .true.
                word = ' '
                return
            end if

        end if
        !
        !  Beginning at index NEXT, search the string for the next nonblank,
        !  which signals the beginning of a word.
        !
        ilo = next
        !
        !  ...S(NEXT:) is blank.  Return with WORD = ' ' and DONE = TRUE.
        !
        do

            if ( lenc < ilo ) then
                word = ' '
                done = .true.
                next = lenc + 1
                return
            end if
            !
            !  If the current character is blank, skip to the next one.
            !
            if ( s(ilo:ilo) /= ' ' .and. s(ilo:ilo) /= TAB ) then
                exit
            end if

            ilo = ilo + 1

        end do

        if ( s(ilo:ilo) == '"' .or. &
            s(ilo:ilo) == '(' .or. &
            s(ilo:ilo) == ')' .or. &
            s(ilo:ilo) == '{' .or. &
            s(ilo:ilo) == '}' .or. &
            s(ilo:ilo) == '[' .or. &
            s(ilo:ilo) == ']' ) then

            word = s(ilo:ilo)
            next = ilo + 1
            return

        end if
        !
        !  Now search for the last contiguous character that is not a
        !  blank, TAB, or special character.
        !
        next = ilo + 1

        do while ( next <= lenc )

            if ( s(next:next) == ' ' ) then
                exit
            else if ( s(next:next) == TAB ) then
                exit
            else if ( s(next:next) == '"' ) then
                exit
            else if ( s(next:next) == '(' ) then
                exit
            else if ( s(next:next) == ')' ) then
                exit
            else if ( s(next:next) == '{' ) then
                exit
            else if ( s(next:next) == '}' ) then
                exit
            else if ( s(next:next) == '[' ) then
                exit
            else if ( s(next:next) == ']' ) then
                exit
            end if

            next = next + 1

        end do

        if ( s(next-1:next-1) == ',' ) then
            word = s(ilo:next-2)
        else
            word = s(ilo:next-1)
        end if

        return
    end subroutine word_next_read


    subroutine test01 ( input_file_name )

        character ( len = * ) input_file_name
        logical stla_check_r

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TEST01'
        write ( *, '(a)' ) '  STLA_CHECK makes some simple checks on a file.'

        write ( *, '(a)' ) ' '
        if ( stla_check ( input_file_name ) ) then
            write ( *, '(a)' ) '  The file "' // trim ( input_file_name ) // &
                '" seems to be a legal ASCII STL file.'
        else
            write ( *, '(a)' ) '  The file "' // trim ( input_file_name ) // &
                '" does NOT seem to be a legal ASCII STL file.'
        end if

        return
    end subroutine test01

    subroutine test02 ( input_file_name )

        integer ( kind= il ) face_num
        character ( len = * ) input_file_name
        integer ( kind= il ) node_num
        integer ( kind= il ) solid_num
        integer ( kind= il ) text_num

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TEST02'
        write ( *, '(a)' ) '  STLA_SIZE determines the size of various objects'
        write ( *, '(a)' ) '  in an ASCII STL file.'

        call stla_size ( input_file_name, solid_num, node_num, face_num, text_num )

        call stla_size_print ( input_file_name, solid_num, node_num, face_num, &
            text_num )

        return
    end subroutine test02

    !subroutine test03 ( input_file_name )
    subroutine obtener_p_t(input_file_name, p_xyz, triang_p, num_p, num_t, t_normal)

        !dummy
        real ( kind= dp ), allocatable, dimension(:,:), intent(out) :: p_xyz
        integer ( kind= il ), allocatable, dimension(:,:), intent(out) :: triang_p
        real ( kind= dp ), allocatable, dimension(:,:), intent(out) :: t_normal
        integer ( kind= il ), intent(out) :: num_p
        integer ( kind= il ), intent(out) :: num_t

        integer ( kind= il ), allocatable, dimension(:,:) :: face_node
        real ( kind= dp ), allocatable, dimension(:,:) :: face_normal
        integer ( kind= il ) face_num
        integer ( kind= il ) ierror
        character ( len = * ) input_file_name
        integer ( kind= il ) node_num
        real ( kind= dp ), allocatable, dimension(:,:) :: node_xyz
        integer ( kind= il ) solid_num
        integer ( kind= il ) text_num
        num_p=0
        num_t=0
        !write ( *, '(a)' ) ' '
        !write ( *, '(a)' ) 'TEST03'
        !write ( *, '(a)' ) '  STLA_READ reads an object in an ASCII STL file.'

        call stla_size ( input_file_name, solid_num, node_num, face_num, text_num )

        allocate ( face_node(3,face_num) )
        allocate ( face_normal(3,face_num) )
        allocate ( node_xyz(3,node_num) )

        call stla_read ( input_file_name, node_num, face_num, node_xyz, &
            face_node, face_normal, ierror )

        if ( ierror /= 0 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a,i8)' ) '  STLA_READ returned IERROR = ', ierror
            return
        end if

        !call stla_size_print ( input_file_name, solid_num, node_num, face_num, &
        !  text_num )

        !call stla_face_node_print ( face_num, face_node)!a eliminar
        !call stla_node_xyz_print ( node_num, node_xyz)!a eliminar

        !print*, 'despues de limpiar nodos'

        call limpiarnodos (node_xyz, face_node, node_num, face_num)

        allocate ( triang_p(3,face_num) )
        allocate ( p_xyz(3,node_num) )
		allocate (t_normal(3,face_num))
        p_xyz = node_xyz
        triang_p = face_node
        num_p = node_num
        num_t = face_num
		t_normal = face_normal
        !call stla_face_node_print ( num_t, triang_p)
        !call stla_node_xyz_print ( num_p, p_xyz )

        deallocate ( face_node )
        deallocate ( face_normal )
        deallocate ( node_xyz )

    end subroutine obtener_p_t

    subroutine limpiarnodos ( p_xyz, triang_p, num_p, num_t)
        ! Limpia los nodos repetidos de p_xyz y actualiza (refresh) la matriz
        ! vertices de triangulos p_xyz en base al p_xyz limpio
        ! GABRIEL
        !Dummy
        real ( kind= dp ), dimension (3,num_p), intent (inout) :: p_xyz
        integer ( kind= il ), dimension (3,num_t), intent (inout) :: triang_p
        integer (kind= il), intent (inout) :: num_p
        integer (kind= il), intent (inout) :: num_t
        !Local
        integer (kind= il) :: i
        integer res

        do i = 1, num_p

            !        if (  (p_xyz(1,i)==p_xyz(1,num_p)) .and. &
            !        (p_xyz(2,i)==p_xyz(2,num_p)) .and. &
            !        (p_xyz(3,i)==p_xyz(3,num_p)) ) then
            !            p_xyz(1,num_p)=0.0
            !            p_xyz(2,num_p)=0.0
            !            p_xyz(3,num_p)=0.0
            !            call actualizar (num_p,i,triang_p,num_t)
            !            num_p=num_p-1
            !            !i=i-1
            !        else
            busqueda: do

                call buscar(p_xyz,i,num_p, res)
                if (res == -1) then
                    exit busqueda
                else

                    p_xyz(:,res) = p_xyz(:,num_p)
                    p_xyz(1,num_p)=0.0
                    p_xyz(2,num_p)=0.0
                    p_xyz(3,num_p)=0.0
                    call actualizar(res,i,triang_p,num_t)
                    call actualizar(num_p,res,triang_p,num_t)
                    num_p=num_p-1
                end if

            end do busqueda

        !        end if

        end do

    end subroutine limpiarnodos

    subroutine buscar(P,indice,np,resu)
        !Dummy
        real ( kind= dp ), dimension (3,np), intent (in) :: P
        integer (kind= il), intent (in) :: indice
        integer (kind= il), intent (in) :: np
        integer (kind= il), intent (out) :: resu
        !Local
        integer (kind= il) :: i

        resu = -1
        do i = indice+1,np
            if ( (P(1,i)==P(1,indice)) .and. &
                (P(2,i)==P(2,indice)) .and. &
                (P(3,i)==P(3,indice)) ) then
                resu = i
                exit
            end if
        end do

    end subroutine buscar

    subroutine actualizar(a,b,M,nt)
        !Dummy
        integer ( kind= il ), dimension (3,nt), intent (inout) :: M
        integer ( kind= il ), intent (in) :: a
        integer (kind= il), intent (in) :: b
        integer (kind= il), intent (in) :: nt
        !Local
        integer (kind= il) :: i

        do i = 1, nt
            if (M(1,i)==a) then
                M(1,i)=b
            elseif (M(2,i)==a) then
                M(2,i)=b
            elseif (M(3,i)==a) then
                M(3,i)=b
            end if
        end do

    end subroutine actualizar


end module modReadStl
