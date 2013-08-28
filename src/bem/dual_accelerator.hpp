#ifndef DUAL_ACCELERATOR_HPP_INCLUDED
#define DUAL_ACCELERATOR_HPP_INCLUDED

#include "../util/crtp_base.hpp"

template <class Derived>
class dual_accelerator_iterator_base
{
public:
	NIHU_CRTP_HELPERS

	bool xi_new(void) const
    {
		return derived().xi_new();
    }

	bool eta_new(void) const
    {
		return derived().eta_new();
    }

	decltype(const_crtp_ptr<Derived>()->get_xi())
    get_xi(void) const
    {
		return derived().get_xi();
    }

	decltype(const_crtp_ptr<Derived>()->get_eta())
    get_eta(void) const
    {
		return derived().get_eta();
    }

	decltype(const_crtp_ptr<Derived>()->get_Nxi())
    get_Nxi(void) const
    {
		return derived().get_Nxi();
    }

	decltype(const_crtp_ptr<Derived>()->get_Neta())
    get_Neta(void) const
    {
		return derived().get_Neta();
    }

	decltype(const_crtp_ptr<Derived>()->get_w())
    get_w(void) const
    {
		return derived().get_w();
    }

    dual_accelerator_iterator_base &operator++(void)
    {
		return derived().operator++();
    }
};


#endif
