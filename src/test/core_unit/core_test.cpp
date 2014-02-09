#include "domain_test.h"
#include "shapeset_test.h"
#include "element_test.h"
#include "space_test.h"

int main(void)
{
    space_test();
    domain_test();
    shapeset_test();
    //element_test();
    return 0;
}
