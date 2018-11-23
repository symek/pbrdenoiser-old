#!/usr/bin/env bash
g++ -v -c *.cpp -I$REZ_BOOST_ROOT/include -I../ext/eigen -I$REZ_OPENEXR_ROOT/include/OpenEXR -I$REZ_ILMBASE_ROOT/include/OpenEXR
g++ -v *.o -o test -pthread -L$REZ_OPENEXR_ROOT/lib -lIlmImf
