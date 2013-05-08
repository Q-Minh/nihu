#include "kernel.hpp"

const helmholtz_G_kernel::kernel_precision helmholtz_G_kernel::limits[] = {
	{5.0, 2},
	{2.0, 5},
	{0.0, 7}
/*
	{9.2, 1},
	{1.6, 3},
	{1.0, 5},
	{0.8, 7},
	{0.7, 9},
	{0.6, 11}
*/
	};

const helmholtz_H_kernel::kernel_precision helmholtz_H_kernel::limits[]  = {
	{5.0, 2},
	{2.0, 5},
	{0.0, 7}
	};


const helmholtz_HG_kernel::kernel_precision helmholtz_HG_kernel::limits[]  = {
	{5.0, 2},
	{2.0, 5},
	{0.0, 7}
	};

