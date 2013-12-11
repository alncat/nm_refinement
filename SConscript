Import("env_base", "env_etc")
env = env_base.Clone(LIBS=["cctbx"]+env_etc.libm)
env_etc.include_registry.append(
    env=env,
    paths=env_etc.mmtbx_common_includes)
if (not env_etc.no_boost_python):
    Import("env_cctbx_boost_python_ext")
    env_bpl = env_cctbx_boost_python_ext.Clone()
    env_etc.include_registry.append(
        env=env_bpl,
        paths=env_etc.mmtbx_common_includes)
    env_bpl.Prepend(LIBS=["cctbx", "iotbx_pdb"])
    env_bpl.SharedLibrary(
        target="#lib/mmtbx_nm_ext",
        source=["nm_ext.cpp"])

