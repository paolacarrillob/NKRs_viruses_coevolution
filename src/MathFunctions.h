/*
 * MathFunctions.h
 *
 *  Created on: Apr 19, 2011
 *      Author: paola
 */

#ifndef MATHFUNCTIONS_H_
#define MATHFUNCTIONS_H_

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_smallint.hpp>
#include <boost/random/normal_distribution.hpp>

#include <vector>
#include <algorithm>
#include <numeric>

int RandomNumber(int min, int max);
int RandomNumberSmallInt(int min, int max);
double RandomNumberDouble();
double RandomNumberDouble(int min, int max);
int RandomNumber(double* probabilities, int count);
double ImmunityTime(int mean, double sigma);

#endif /* MATHFUNCTIONS_H_ */

