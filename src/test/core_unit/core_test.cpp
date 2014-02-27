#include "domain_test.h"
#include "shapeset_test.h"
#include "element_test.h"
#include "space_test.h"


// #include "../library_unit/laplace_singular_test.hpp"

int main(void)
{
    space_test();
    domain_test();
    shapeset_test();
    element_test();

//    laplace_singular_integrals_test();

    return 0;
}
