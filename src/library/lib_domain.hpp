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
    static xi_t get_center_impl(void) { return m_center; }

private:
    static corners_t const m_corners;
    static xi_t const m_center;
};

class tria_domain : public domain_base<tria_domain>
{
public:
    static xi_t const *get_corners_impl(void) { return m_corners; }
    static xi_t get_center_impl(void) { return m_center; }

private:
    static corners_t const m_corners;
    static xi_t const m_center;
};

class quad_domain : public domain_base<quad_domain>
{
public:
    static xi_t const *get_corners_impl(void) { return m_corners; }
    static xi_t get_center_impl(void) { return m_center; }

private:
    static corners_t const m_corners;
    static xi_t const m_center;
};

class brick_domain : public domain_base<brick_domain>
{
public:
    static xi_t const *get_corners_impl(void) { return m_corners; }
    static xi_t get_center_impl(void) { return m_center; }

private:
    static corners_t const m_corners;
    static xi_t const m_center;
};

#endif // LIB_DOMAIN_HPP_INCLUDED
