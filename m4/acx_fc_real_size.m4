dnl acx_fc_real_size.m4 --- determine size of fortran real types
dnl
dnl Copyright  (C)  2010  Thomas Jahns <jahns@dkrz.de>
dnl
dnl Version: 1.0
dnl Keywords: configure configure.ac autotools fortran real size
dnl Author: Thomas Jahns <jahns@dkrz.de>
dnl Maintainer: Thomas Jahns <jahns@dkrz.de>
dnl URL: https://www.dkrz.de/redmine/projects/show/scales-ppm
dnl
dnl Redistribution and use in source and binary forms, with or without
dnl modification, are  permitted provided that the following conditions are
dnl met:
dnl
dnl Redistributions of source code must retain the above copyright notice,
dnl this list of conditions and the following disclaimer.
dnl
dnl Redistributions in binary form must reproduce the above copyright
dnl notice, this list of conditions and the following disclaimer in the
dnl documentation and/or other materials provided with the distribution.
dnl
dnl Neither the name of the DKRZ GmbH nor the names of its contributors
dnl may be used to endorse or promote products derived from this software
dnl without specific prior written permission.
dnl
dnl THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
dnl IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
dnl TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
dnl PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
dnl OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
dnl EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
dnl PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
dnl PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
dnl LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
dnl NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
dnl SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
dnl
dnl Commentary:
dnl
dnl
dnl
dnl Code:
dnl
dnl
dnl ACX_FORTRAN_RUN_CHECK_SIZEOF(var_for_size, )
dnl
dnl sets var_for_size to the size.  Ignores if the size cannot be determined
dnl
AC_DEFUN([ACX_FORTRAN_RUN_CHECK_SIZEOF],
  [AC_REQUIRE([ACX_C_CHAR_BITS])
   AS_VAR_PUSHDEF([acx_fortran_Type], [acx_cv_fortran_type_$1])dnl
   ACX_FORTRAN_CHECK_TYPE([$1])
   AS_IF([test x]AS_VAR_GET([acx_fortran_Type])[ = xyes],[
     AS_VAR_PUSHDEF([acx_fortran_Sizeof],[acx_cv_fortran_sizeof_$1])dnl
     AC_CACHE_CHECK([size of Fortran type $1], [acx_fortran_Sizeof],dnl
       [AC_LANG_PUSH([Fortran])
        AC_RUN_IFELSE(AC_LANG_PROGRAM([], [
    integer :: itest
    real    :: rtest
    integer :: integer_io_length
    integer :: real_io_length
    integer :: integer_bits
    integer :: real_bits
    inquire(iolength=integer_io_length) itest
    inquire(iolength=real_io_length) rtest
    integer_bits=bit_size(itest)
    real_bits = real_io_length * integer_bits / integer_io_length
    open(10,file="conftestval")
    write(10,*)real_bits
    close(10)
]),dnl
          [AS_VAR_SET([acx_fortran_Sizeof], `cat conftestval`)
           AS_VAR_SET([acx_fortran_Sizeof],[`expr ]AS_VAR_GET([acx_fortran_Sizeof])[ / $acx_cv_c_char_bits`])],
          [AS_VAR_SET([acx_fortran_Sizeof], [unavailable])])
        /bin/rm -f conftestval
        AC_LANG_POP([Fortran])
       ])
     AS_IF([test x"]AS_VAR_GET([acx_fortran_Sizeof])[" = xunavailable], [$2], [$3])
     m4_ifval([$2], $2=AS_VAR_GET([acx_fortran_Sizeof]))
     AS_VAR_POPDEF([acx_fortran_Sizeof])dnl
  ])dnl
  AS_VAR_POPDEF([acx_fortran_Type])dnl
])
dnl Local Variables:
dnl mode: autoconf
dnl license-project-url: "https://www.dkrz.de/redmine/projects/show/scales-ppm"
dnl license-default: "bsd"
dnl End:
