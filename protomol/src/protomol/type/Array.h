/*  -*- c++ -*-  */
#ifndef ARRAY_H
#define ARRAY_H

#include <protomol/config.h>

// Wrapper header file for N-dim. Array class
//
// If your compiler does not support partial template
// specialization, define NO_PARTIAL_TEMPLATE_SPECIALIZATION

#ifdef NO_PARTIAL_TEMPLATE_SPECIALIZATION
#include <protomol/type/ArrayNoPartialSpecialization.h>
#else
#include <protomol/type/ArrayFastest.h>
#endif

#endif

