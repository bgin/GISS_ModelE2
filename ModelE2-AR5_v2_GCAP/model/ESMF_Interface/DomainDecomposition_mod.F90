#include "rundeck_opts.h"

module DOMAIN_DECOMP_1D
!@sum  This module conveniently wraps (lumps together) various domain 
!@+ decomposition modules.
!@auth Tom Clune GSFC/SSSO/610.3
  use dist_grid_mod
  use Halo_mod
  use SpecialIO_mod
  use GatherScatter_mod
  use GlobalSum_mod
  implicit none
  public

end module DOMAIN_DECOMP_1D

#ifndef CUBED_SPHERE
#define DOMAIN_DECOMP_ATM_IS_1D
#endif

#ifdef DOMAIN_DECOMP_ATM_IS_1D
! If the atmosphere has a 1D domain decomposition, pass along the contents
! of DOMAIN_DECOMP_1D
      MODULE DOMAIN_DECOMP_ATM
        use dist_grid_mod
        use Halo_mod
        use SpecialIO_mod
        use GatherScatter_mod
        use GlobalSum_mod
      implicit none
      public
      END MODULE DOMAIN_DECOMP_ATM
#endif
