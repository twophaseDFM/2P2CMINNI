add_input_file_links(FILES params.input)

add_gstat_file_links()

# a gstat test (becaue it's a random permeability field we can't test against a reference solution)
dumux_add_test(NAME 2prandomK_box
              CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-alugrid_FOUND )"
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=InjectionBox ENABLECACHING=0
              CMAKE_GUARD HAVE_GSTAT
              COMMAND ./2prandomK
              CMD_ARGS params_gstat.input)
              
dumux_add_test(NAME 2prandomK_tpfa
              SOURCES main.cc
              CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-alugrid_FOUND )"
              COMPILE_DEFINITIONS TYPETAG=InjectionCCTpfa ENABLECACHING=0
              CMAKE_GUARD HAVE_GSTAT
              COMMAND ./2prandomK
              CMD_ARGS params_gstat.input)

