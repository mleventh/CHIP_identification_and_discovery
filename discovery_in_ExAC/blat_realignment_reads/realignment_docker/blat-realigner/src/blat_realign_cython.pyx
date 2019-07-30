from libcpp.string cimport string

cdef extern from "blat_realign.1.cpp":
    cdef cppclass RealignmentFilter:
        RealignmentFilter() except +
        RealignmentFilter(const string  genome, const string  ooc) except +
        int blat_seq(const string  seq, const string  contig,const int mut_position)

cdef class pyRealignmentFilter:
    cdef RealignmentFilter rf      # hold a C++ instance which we're wrapping
    def __cinit__(self,const string genome,const string ooc):
        self.rf = RealignmentFilter( genome,ooc)
    def blat_seq(self,seq,contig,mut_position):
        return self.rf.blat_seq(seq,contig,mut_position)
