module m_string
    implicit none

    contains

        function to_string(i) result(s)
            implicit none

            character(:),allocatable :: s

            integer :: i
            
            character(32) :: ss
            integer :: l
            
            write(ss, *) i
            ss = adjustl(ss)

            l = len_trim(ss)

            allocate(character(l)::s)

            s = trim(ss)

        end function

end module

