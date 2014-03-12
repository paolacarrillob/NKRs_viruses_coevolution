/*
 * MathFunctions.cpp
 *
 *  Created on: Apr 19, 2011
 *      Author: paola
 */

#include "MathFunctions.h"

//Prepare the randomizer

namespace MathFunctions{
   boost::mt19937 boostRandomSeed(static_cast<unsigned int>(time(NULL)));
   boost::uniform_real<> u01( 0., 1. );
   boost::variate_generator<boost::mt19937&, boost::uniform_real<> > s01(boostRandomSeed, u01);
   boost::variate_generator<boost::mt19937&, boost::uniform_real<> > RandomNumberDouble(boostRandomSeed, u01);
}

int RandomNumberSmallInt(int min, int max)
{
	boost::uniform_smallint<> dist(min,max);
	boost::variate_generator<boost::mt19937&, boost::uniform_smallint<> > die(MathFunctions::boostRandomSeed, dist);
	return die();
}

int RandomNumber(int min, int max)
{
   return (int)(MathFunctions::s01()*(max-min+1))+min;
	//boost::uniform_int<> dist(min,max);
	//boost::variate_generator<boost::mt19937&, boost::uniform_int<> > die(boostRandomSeed, dist);
	//return die();
}
double RandomNumberDouble(int min, int max)
{
   double r = MathFunctions::s01();
   return r*(max-min)+min;
	//boost::uniform_real<> a(min, max);
	//boost::variate_generator<boost::mt19937&, boost::uniform_real<> > bin(boostRandomSeed, a);
	//return bin();
}

double RandomNumberDouble()
{
   return MathFunctions::s01();
	//boost::uniform_real<> a(min, max);
	//boost::variate_generator<boost::mt19937&, boost::uniform_real<> > bin(boostRandomSeed, a);
	//return bin();
}

int RandomNumber(double* probabilities, int count)
{
	std::vector<double> cumulative;
	partial_sum(&probabilities[0],&probabilities[0] + count-1, back_inserter(cumulative));
	boost::uniform_real<> dist(0, cumulative.back());
	boost::variate_generator<boost::mt19937&, boost::uniform_real<> > die(MathFunctions::boostRandomSeed, dist);
	int i = (std::lower_bound(cumulative.begin(), cumulative.end(), die()) - cumulative.begin());
	return i;
}
double ImmunityTime(int mean, double sigma)
{
	boost::normal_distribution<> a(mean, sigma);
	boost::variate_generator<boost::mt19937&, boost::normal_distribution<> > bin(MathFunctions::boostRandomSeed, a);
	return bin();
}
