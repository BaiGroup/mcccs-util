#!/bin/sh
rm -f config.cache acconfig.h
touch INSTALL COPYING NEWS README AUTHORS ChangeLog
aclocal #-I m4
autoconf
#autoheader
automake -a
exit

#Create statically linked libraries
g++ -c molsim.cpp
ar rcsfvq libmolsim.a molsim.o

#Create dynamically linked libraries
#On Mac OS X
g++ -dynamiclib -flat_namespace -o libmolsim.so.1.0 molsim.o
g++ -o driver.exe driver.cpp
#On Linux
g++ -shared -Wl,-soname,libshared.so.1 -o libmolsim.so.1.0 molsim.o

ln -sf libmolsim.so.1.0 libmolsim.so
ln -sf libmolsim.so.1.0 libmolsim.so.1
