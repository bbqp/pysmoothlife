import cffi
import os


clas s

if __name__ == '__main__':
    print("Building CFFI Module")
    ffi = cffi.FFI()
    
    this_dir = os.path.abspath('.')
    h_file_name = this_dir + '/' + 'transitionfunction.h'
    with open(h_file_name) as h_file:
        ffi.cdef(h_file.read())
        
    # tasks.py
    ffi.set_source(
        "cffi_example",
        # Since you're calling a fully-built library directly, no custom source
        # is necessary. You need to include the .h files, though, because behind
        # the scenes cffi generates a .c file that contains a Python-friendly
        # wrapper around each of the functions.
        '#include "transitionfunction.h"',
        # The important thing is to include the pre-built lib in the list of
        # libraries you're linking against:
        libraries=["transitionfunction"],
        library_dirs=[this_dir.as_posix()],
        extra_link_args=[''],
    )

    # tasks.py
    ffi.compile()
