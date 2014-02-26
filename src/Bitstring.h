/*
 * Bitstring.h
 *
 *  Created on: Apr 19, 2011
 *      Author: paola
 */

#ifndef BITSTRING_H_
#define BITSTRING_H_

#include <iostream>
#include <bitset>
#include <string>
#include "MathFunctions.h"
#define MOLECULE_LENGTH 16
using namespace std;

class Bitstring {
public:
	Bitstring(); //works
	Bitstring(int molecule_number){Int2Bit(molecule_number);}
//	Bitstring(const string& mhc); // constructs  mhc molecules
//	Bitstring(const string& kir, Bitstring& mhc, int specificity); // constructs kir molecules
	virtual ~Bitstring(){};
	vector<int>& GetBits(){return bits;}
	unsigned int GetBitSize(){return bits.size();}
	bool operator==(Bitstring& rhs);
	int Bit2Int(); //works -
	void Int2Bit(int molecule_number);//
	int AdjacentMatch(Bitstring& kirBits); //works
	//int AdjacentMatch(Bitstring& otherBits); //works
	void PointMutation(); //works
	void PrintBits();
protected:
	vector<int> bits;


};

class BitstringMHC : public Bitstring
{
public:
	BitstringMHC();
	virtual ~BitstringMHC(){};
};

class BitstringKIR: public Bitstring
{
public:
	BitstringKIR(Bitstring& mhc, int specificity);
	virtual ~BitstringKIR(){};
\
};

#endif /* BITSTRING_H_ */
