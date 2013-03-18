//Non-inline functions for Linear elasticity elements
#include "Thomo_lin_elasticity_elements.h"


namespace oomph
{



/////////////////////////////////////////////////////////////////////////
// TLinearElasticityElement
/////////////////////////////////////////////////////////////////////////


//====================================================================
// Force build of templates
//====================================================================
 template class THomogenisedLinearElasticityElement<2,2>;
 template class THomogenisedLinearElasticityElement<2,3>;
 template class THomogenisedLinearElasticityElement<2,4>;

 template class THomogenisedLinearElasticityElement<3,2>;
 template class THomogenisedLinearElasticityElement<3,3>;
}
