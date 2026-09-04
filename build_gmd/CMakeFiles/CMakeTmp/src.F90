      program main
      implicit none
#ifdef __NEC__
      integer a
#else
      choke me
#endif
#ifndef __NEC__
      choke me
#else
      integer b
#endif
      a = 4
      b = 2
      end
