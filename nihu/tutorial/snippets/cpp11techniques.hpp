/**
 * \file cpp11techniques.hpp
 * \brief Gives a brief overview of the C++11 techniques applied in NiHu.
 * \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
 */

//! [VariadicTemplateClass]
template <class...Args>
class VariadicTemplateClass 
{
public:
	static int a_function(Args...a) {return 0;}
	static int b_function(Args...b) {return a_function(b...);}
};
//! [VariadicTemplateClass]

//! [VariadicExample]
typename traits<Args>::type...
//! [VariadicExample]

//! [VariadicInterpretation]
typename traits<bool>::type, typename traits<char>::type
//! [VariadicInterpretation]

//! [VariadicSpec1]
template <bool...Args>
struct vari_or {
	static const bool value = false;
};
//! [VariadicSpec1]

//! [VariadicSpec2]
template <bool...Args>
struct vari_or<true, Args...> {
	static const bool value = true;
};

template <bool...Args>
struct vari_or<false, Args...> {
	static const bool value = vari_or<Args...>::value;
};
//! [VariadicSpec2]

//! [VariadicSpec3]
vari_or<false, true>::value;                               // 2 args
vari_or<false, true, false, false, true, false>::value;    // 6 args
// and so on with any number of args 
//! [VariadicSpec3]

//! [Auto first]
auto c = 'a';
//! [Auto first]

//! [Dirac wrapper class]
template <class FuncSpace>
class dirac_wrapper
{
...
public:
	// constructor from function space
	dirac_wrapper(FuncSpace const &fs)
	{
		...
	}
};
//! [Dirac wrapper class]

//! [Dirac wrapper factory]
// factory function
template <class FuncSpace>
dirac_wrapper<FuncSpace> dirac(FuncSpace const &fs)
{
	return dirac_wrapper<FuncSpace>(fs);
}
//! [Dirac wrapper factory]

//! [Dirac wrapper old declaration]
dirac_wrapper<func_space_t> dfs(func_space);
//! [Dirac wrapper old declaration]

//! [Dirac wrapper usage]
auto dfs = dirac(func_space);
//! [Dirac wrapper usage]

