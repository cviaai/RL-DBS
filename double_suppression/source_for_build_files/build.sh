 swig -python gfg.i
 g++ -c -fpic gfg_wrap.c gfg.c -I/usr/include/python3.7
 g++ -shared gfg.o gfg_wrap.o -o _oscillator_cpp.so
