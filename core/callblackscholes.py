# This file was automatically generated by SWIG (http://www.swig.org).
# Version 2.0.11
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.





from sys import version_info
if version_info >= (2,6,0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_callblackscholes', [dirname(__file__)])
        except ImportError:
            import _callblackscholes
            return _callblackscholes
        if fp is not None:
            try:
                _mod = imp.load_module('_callblackscholes', fp, pathname, description)
            finally:
                fp.close()
            return _mod
    _callblackscholes = swig_import_helper()
    del swig_import_helper
else:
    import _callblackscholes
del version_info
try:
    _swig_property = property
except NameError:
    pass # Python < 2.2 doesn't have 'property'.
def _swig_setattr_nondynamic(self,class_type,name,value,static=1):
    if (name == "thisown"): return self.this.own(value)
    if (name == "this"):
        if type(value).__name__ == 'SwigPyObject':
            self.__dict__[name] = value
            return
    method = class_type.__swig_setmethods__.get(name,None)
    if method: return method(self,value)
    if (not static):
        self.__dict__[name] = value
    else:
        raise AttributeError("You cannot add attributes to %s" % self)

def _swig_setattr(self,class_type,name,value):
    return _swig_setattr_nondynamic(self,class_type,name,value,0)

def _swig_getattr(self,class_type,name):
    if (name == "thisown"): return self.this.own()
    method = class_type.__swig_getmethods__.get(name,None)
    if method: return method(self)
    raise AttributeError(name)

def _swig_repr(self):
    try: strthis = "proxy of " + self.this.__repr__()
    except: strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)

try:
    _object = object
    _newclass = 1
except AttributeError:
    class _object : pass
    _newclass = 0


class SwigPyIterator(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, SwigPyIterator, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, SwigPyIterator, name)
    def __init__(self, *args, **kwargs): raise AttributeError("No constructor defined - class is abstract")
    __repr__ = _swig_repr
    __swig_destroy__ = _callblackscholes.delete_SwigPyIterator
    __del__ = lambda self : None;
    def value(self): return _callblackscholes.SwigPyIterator_value(self)
    def incr(self, n=1): return _callblackscholes.SwigPyIterator_incr(self, n)
    def decr(self, n=1): return _callblackscholes.SwigPyIterator_decr(self, n)
    def distance(self, *args): return _callblackscholes.SwigPyIterator_distance(self, *args)
    def equal(self, *args): return _callblackscholes.SwigPyIterator_equal(self, *args)
    def copy(self): return _callblackscholes.SwigPyIterator_copy(self)
    def next(self): return _callblackscholes.SwigPyIterator_next(self)
    def __next__(self): return _callblackscholes.SwigPyIterator___next__(self)
    def previous(self): return _callblackscholes.SwigPyIterator_previous(self)
    def advance(self, *args): return _callblackscholes.SwigPyIterator_advance(self, *args)
    def __eq__(self, *args): return _callblackscholes.SwigPyIterator___eq__(self, *args)
    def __ne__(self, *args): return _callblackscholes.SwigPyIterator___ne__(self, *args)
    def __iadd__(self, *args): return _callblackscholes.SwigPyIterator___iadd__(self, *args)
    def __isub__(self, *args): return _callblackscholes.SwigPyIterator___isub__(self, *args)
    def __add__(self, *args): return _callblackscholes.SwigPyIterator___add__(self, *args)
    def __sub__(self, *args): return _callblackscholes.SwigPyIterator___sub__(self, *args)
    def __iter__(self): return self
SwigPyIterator_swigregister = _callblackscholes.SwigPyIterator_swigregister
SwigPyIterator_swigregister(SwigPyIterator)

class IntVector(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, IntVector, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, IntVector, name)
    __repr__ = _swig_repr
    def iterator(self): return _callblackscholes.IntVector_iterator(self)
    def __iter__(self): return self.iterator()
    def __nonzero__(self): return _callblackscholes.IntVector___nonzero__(self)
    def __bool__(self): return _callblackscholes.IntVector___bool__(self)
    def __len__(self): return _callblackscholes.IntVector___len__(self)
    def pop(self): return _callblackscholes.IntVector_pop(self)
    def __getslice__(self, *args): return _callblackscholes.IntVector___getslice__(self, *args)
    def __setslice__(self, *args): return _callblackscholes.IntVector___setslice__(self, *args)
    def __delslice__(self, *args): return _callblackscholes.IntVector___delslice__(self, *args)
    def __delitem__(self, *args): return _callblackscholes.IntVector___delitem__(self, *args)
    def __getitem__(self, *args): return _callblackscholes.IntVector___getitem__(self, *args)
    def __setitem__(self, *args): return _callblackscholes.IntVector___setitem__(self, *args)
    def append(self, *args): return _callblackscholes.IntVector_append(self, *args)
    def empty(self): return _callblackscholes.IntVector_empty(self)
    def size(self): return _callblackscholes.IntVector_size(self)
    def clear(self): return _callblackscholes.IntVector_clear(self)
    def swap(self, *args): return _callblackscholes.IntVector_swap(self, *args)
    def get_allocator(self): return _callblackscholes.IntVector_get_allocator(self)
    def begin(self): return _callblackscholes.IntVector_begin(self)
    def end(self): return _callblackscholes.IntVector_end(self)
    def rbegin(self): return _callblackscholes.IntVector_rbegin(self)
    def rend(self): return _callblackscholes.IntVector_rend(self)
    def pop_back(self): return _callblackscholes.IntVector_pop_back(self)
    def erase(self, *args): return _callblackscholes.IntVector_erase(self, *args)
    def __init__(self, *args): 
        this = _callblackscholes.new_IntVector(*args)
        try: self.this.append(this)
        except: self.this = this
    def push_back(self, *args): return _callblackscholes.IntVector_push_back(self, *args)
    def front(self): return _callblackscholes.IntVector_front(self)
    def back(self): return _callblackscholes.IntVector_back(self)
    def assign(self, *args): return _callblackscholes.IntVector_assign(self, *args)
    def resize(self, *args): return _callblackscholes.IntVector_resize(self, *args)
    def insert(self, *args): return _callblackscholes.IntVector_insert(self, *args)
    def reserve(self, *args): return _callblackscholes.IntVector_reserve(self, *args)
    def capacity(self): return _callblackscholes.IntVector_capacity(self)
    __swig_destroy__ = _callblackscholes.delete_IntVector
    __del__ = lambda self : None;
IntVector_swigregister = _callblackscholes.IntVector_swigregister
IntVector_swigregister(IntVector)

class DoubleVector(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, DoubleVector, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, DoubleVector, name)
    __repr__ = _swig_repr
    def iterator(self): return _callblackscholes.DoubleVector_iterator(self)
    def __iter__(self): return self.iterator()
    def __nonzero__(self): return _callblackscholes.DoubleVector___nonzero__(self)
    def __bool__(self): return _callblackscholes.DoubleVector___bool__(self)
    def __len__(self): return _callblackscholes.DoubleVector___len__(self)
    def pop(self): return _callblackscholes.DoubleVector_pop(self)
    def __getslice__(self, *args): return _callblackscholes.DoubleVector___getslice__(self, *args)
    def __setslice__(self, *args): return _callblackscholes.DoubleVector___setslice__(self, *args)
    def __delslice__(self, *args): return _callblackscholes.DoubleVector___delslice__(self, *args)
    def __delitem__(self, *args): return _callblackscholes.DoubleVector___delitem__(self, *args)
    def __getitem__(self, *args): return _callblackscholes.DoubleVector___getitem__(self, *args)
    def __setitem__(self, *args): return _callblackscholes.DoubleVector___setitem__(self, *args)
    def append(self, *args): return _callblackscholes.DoubleVector_append(self, *args)
    def empty(self): return _callblackscholes.DoubleVector_empty(self)
    def size(self): return _callblackscholes.DoubleVector_size(self)
    def clear(self): return _callblackscholes.DoubleVector_clear(self)
    def swap(self, *args): return _callblackscholes.DoubleVector_swap(self, *args)
    def get_allocator(self): return _callblackscholes.DoubleVector_get_allocator(self)
    def begin(self): return _callblackscholes.DoubleVector_begin(self)
    def end(self): return _callblackscholes.DoubleVector_end(self)
    def rbegin(self): return _callblackscholes.DoubleVector_rbegin(self)
    def rend(self): return _callblackscholes.DoubleVector_rend(self)
    def pop_back(self): return _callblackscholes.DoubleVector_pop_back(self)
    def erase(self, *args): return _callblackscholes.DoubleVector_erase(self, *args)
    def __init__(self, *args): 
        this = _callblackscholes.new_DoubleVector(*args)
        try: self.this.append(this)
        except: self.this = this
    def push_back(self, *args): return _callblackscholes.DoubleVector_push_back(self, *args)
    def front(self): return _callblackscholes.DoubleVector_front(self)
    def back(self): return _callblackscholes.DoubleVector_back(self)
    def assign(self, *args): return _callblackscholes.DoubleVector_assign(self, *args)
    def resize(self, *args): return _callblackscholes.DoubleVector_resize(self, *args)
    def insert(self, *args): return _callblackscholes.DoubleVector_insert(self, *args)
    def reserve(self, *args): return _callblackscholes.DoubleVector_reserve(self, *args)
    def capacity(self): return _callblackscholes.DoubleVector_capacity(self)
    __swig_destroy__ = _callblackscholes.delete_DoubleVector
    __del__ = lambda self : None;
DoubleVector_swigregister = _callblackscholes.DoubleVector_swigregister
DoubleVector_swigregister(DoubleVector)


def callBlackScholes(x0=100, r=0.05, sigma=0.2, K=100, T=5, seed=0):
  return _callblackscholes.callBlackScholes(x0, r, sigma, K, T, seed)
callBlackScholes = _callblackscholes.callBlackScholes
# This file is compatible with both classic and new-style classes.


