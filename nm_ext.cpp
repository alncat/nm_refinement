#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <mmtbx/nm/nm.h>
#include <scitbx/array_family/boost_python/shared_wrapper.h>
#include <scitbx/boost_python/is_polymorphic_workaround.h>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>
#include <boost/python.hpp>

SCITBX_BOOST_IS_POLYMORPHIC_WORKAROUND(mmtbx::nm::common)

namespace mmtbx { namespace nm {
    namespace bp = boost::python;

namespace {
    void init_module()
    {
        using namespace boost::python;
        using boost::python::arg;
        
        class_<nm_init>("nm_init",
                    init<const char*,
                        std::size_t,
                        af::shared<iotbx::pdb::hierarchy::atom> const&,
                        bool,
                        bool>(
                            (arg("filename"),
                             arg("n_modes"),
                             arg("atoms"),
                             arg("zero_mode_input_flag"),
                             arg("zero_mode_flag"))))
        .def("return_modes", &nm_init::return_modes)
        .def("print_eigenvector", &nm_init::print_eigenvector)
        .def("gen_zero_modes", &nm_init::gen_zero_modes)
        .def("return_zero_modes", &nm_init::return_zero_modes)
        ;
    }
}
}}

BOOST_PYTHON_MODULE(mmtbx_nm_ext)
{
    mmtbx::nm::init_module();
}
