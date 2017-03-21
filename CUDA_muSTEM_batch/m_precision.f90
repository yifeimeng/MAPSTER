module m_precision

integer, parameter, public :: Single = kind(0.0)
integer, parameter, public :: Double = kind(0.0d0)

#ifdef double_precision
integer, parameter, public :: fp_kind = Double
#elif single_precision
integer, parameter, public :: fp_kind = Single
#else
fail
#endif

end module m_precision
