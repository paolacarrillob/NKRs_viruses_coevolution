/*
 * Bitstring.cpp
 *
 *  Created on: Apr 19, 2011
 *      Author: paola
 */

#include "Bitstring.h"
/* Constructor: creates MHC molecules as a random bit string of length 16*/

Bitstring::Bitstring()
{
}


BitstringMHC::BitstringMHC()
{
	for(int i = 0; i< MOLECULE_LENGTH; i++)
	{
		char randmhc = int(RandomNumberDouble()<0.5);
		Bitstring::bits.push_back(randmhc);
	}
}

/* Constructor: creates KIR molecules as a random bit string that recognizes MHC molecules*/
BitstringKIR::BitstringKIR(Bitstring& mhc, int specificity)
{
	int max_start = MOLECULE_LENGTH - specificity - 1;
	int randomBit = RandomNumber(0,max_start); //pick a random position to start the complementary matching
	for (int i = 0; i < MOLECULE_LENGTH+1; i ++)
	{
		if (i == randomBit)
		{
			for(int j = randomBit; j< (randomBit + specificity); j++)
			{
				Bitstring::bits.push_back(1 - mhc.GetBits().at(j));
				i++;
			}
		}

		else
		{
			char randkir = int(RandomNumberDouble()<0.5);
			Bitstring::bits.push_back(randkir);
		}
	}
}



/* This functions checks whether two bit strings are equal*/
bool Bitstring:: operator == (Bitstring& rhs)
{
	if(rhs.GetBitSize() != bits.size())
		return false;
	for(unsigned int i = 0; i < bits.size(); i++)
	{
		if(bits.at(i) != rhs.bits.at(i))
			return false;
	}
	return true;
}

/* This function gives each bit string an "ID" corresponding to its binary value*/
int Bitstring :: Bit2Int()
{
	int basis = 1;
	int moleculeID = 0;
	vector <int>::iterator bit_it;
	for (bit_it = bits.begin(); bit_it != bits.end(); bit_it ++)
	{
		moleculeID += (*bit_it)*basis;
		basis*=2;
	}
	if(moleculeID <0)
	{
		cout << "ERROR! Int2Bit: gene ID cannot be negative!!! Exiting ... \n";
		exit(0);
	}
	return moleculeID;
}

/*This function builds a bit string from an "ID" */
void Bitstring :: Int2Bit(int moleculeID)
{
	if(moleculeID <0)
	{
		cout << "ERROR! Bit2Int: gene ID cannot be negative!!! Exiting ... \n";
		exit(0);
	}

	for(unsigned int i=0; i< MOLECULE_LENGTH; i ++)
	{
		bits.push_back(moleculeID%2);
		moleculeID = moleculeID/2;
	}
}

/* This function checks the complementary adjacent match between two bit strings */
int Bitstring :: AdjacentMatch(Bitstring& anotherBits)
{
	int adjLength = 0;
	int maxStrength = 0;

	// search for the maximal binding strength
	for(unsigned int i=0; i<anotherBits.GetBitSize(); i++)
	{
		if(anotherBits.bits.at(i) == bits.at(i))
		{
			if(adjLength > maxStrength)
				maxStrength = adjLength;
			adjLength = 0;
		}
		else
			adjLength ++;
	}
	if(adjLength > maxStrength)
		maxStrength = adjLength;
	return maxStrength;
}

/* This function checks the complementary adjacent match between two bit strings */
/*int Bitstring :: AdjacentMatch(Bitstring& otherBits)
{
	int counter = 0;
	int max = 0;
	for( int l = 1; l < pow(2.0,16); l = l*2)
	{
		if((bit_ID & l) ^ (otherBits.bit_ID & l)) //check with XOR operator whether the bits are  complementary
		{
			counter ++;
			if(counter > max)
				max = counter;
		}
		else
			counter = 0;
	}
	if(counter > max)
		max = counter;
	return max;
}
*/

void Bitstring :: PointMutation()
{
	unsigned int position = RandomNumber(0,bits.size()-1);
	for (unsigned int i=0; i<bits.size(); i++)
	{
		if(i == position)
		{
			bits.at(i)= 1-bits.at(i);
			cout <<"point mutation just happened\n";
		}
	}
}

void Bitstring:: PrintBits()
{
	for (unsigned int i=0; i<bits.size(); i++)
	{
		cout << bits.at(i) << " ";
	}
}
