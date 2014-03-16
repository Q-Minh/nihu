// This file is a part of NiHu, a C++ BEM template library.
// 
// Copyright (C) 2012-2014  Peter Fiala <fiala@hit.bme.hu>
// Copyright (C) 2012-2014  Peter Rucz <rucz@hit.bme.hu>
// 
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#ifndef LIB_DOMAIN_HPP_INCLUDED
#define LIB_DOMAIN_HPP_INCLUDED

#include "../core/domain.hpp"

class line_domain;
class tria_domain;
class quad_domain;
class brick_domain;

namespace domain_traits
{
    template <> struct space_type<line_domain> : space_1d {};
    template <> struct space_type<tria_domain> : space_2d {};
    template <> struct space_type<quad_domain> : space_2d {};
    template <> struct space_type<brick_domain> : space_3d {};

    template <> struct volume<line_domain> { static constexpr double value = 2.0; };
    template <> struct volume<tria_domain> { static constexpr double value = 0.5; };
    template <> struct volume<quad_domain> { static constexpr double value = 4.0; };
    template <> struct volume<brick_domain> { static constexpr double value = 8.0; };

    template <> struct num_corners<line_domain> { enum { value = 2 }; };
    template <> struct num_corners<tria_domain> { enum { value = 3 }; };
    template <> struct num_corners<quad_domain> { enum { value = 4 }; };
    template <> struct num_corners<brick_domain> { enum { value = 8 }; };
}

class line_domain : public domain_base<line_domain>
{
public:
    static xi_t const *get_corners_impl(void) { return m_corners; }
    static xi_t const &get_center_impl(void) { return m_center; }

private:
    static corners_t const m_corners;
    static xi_t const m_center;
};

class tria_domain : public domain_base<tria_domain>
{
public:
    static xi_t const *get_corners_impl(void) { return m_corners; }
    static xi_t const &get_center_impl(void) { return m_center; }

private:
    static corners_t const m_corners;
    static xi_t const m_center;
};

class quad_domain : public domain_base<quad_domain>
{
public:
    static xi_t const *get_corners_impl(void) { return m_corners; }
    static xi_t const &get_center_impl(void) { return m_center; }

private:
    static corners_t const m_corners;
    static xi_t const m_center;
};

class brick_domain : public domain_base<brick_domain>
{
public:
    static xi_t const *get_corners_impl(void) { return m_corners; }
    static xi_t const &get_center_impl(void) { return m_center; }

private:
    static corners_t const m_corners;
    static xi_t const m_center;
};

#endif // LIB_DOMAIN_HPP_INCLUDED
