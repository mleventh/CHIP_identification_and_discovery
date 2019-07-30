import sys
import os
import shutil

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = [ Extension("blat_realign_cython",
        sources=["blat_realign_cython.pyx","blat_realign.1.cpp"],  # additional source file(s)
        language="c++",
		extra_compile_args=["-I../SeqLib","-I../SeqLib/htslib","-I../SeqLib/blat/inc","-I../SeqLib/blat/jkOwnLib","-W","-O2","-O3","-Wall","-pedantic","-std=c++11","-fPIC","-rdynamic", "-shared"],
		extra_link_args=["../SeqLib/src/libseqlib.a","../SeqLib/htslib/libhts.a","../SeqLib/htslib/libhts.a","../SeqLib/blat/lib/jkOwnLib.a","../SeqLib/blat/lib/jkweb.a","-lpthread","-lz"]
		)]
	)
