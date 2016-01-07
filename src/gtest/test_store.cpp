#include <gtest/gtest.h>
#include "util/store_pattern.hpp"

// Dummy tester class that counts the number of instantiations
class Dummy {
private:
	static unsigned n_constr;

public:
	Dummy() { n_constr++; }
	static unsigned get_n_constr(void) { return n_constr; }
};

unsigned Dummy::n_constr = 0;


TEST(Store, Store1){
	auto d = store<Dummy>::get_data();
	EXPECT_EQ(d.get_n_constr(), 1);
	auto d2 = store<Dummy>::get_data();
	EXPECT_EQ(d2.get_n_constr(), 1);
}


