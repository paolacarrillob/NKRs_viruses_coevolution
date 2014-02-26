/*
 * Genes.h
 *
 *  Created on: Apr 20, 2011
 *      Author: paola
 */

#ifndef GENES_H_
#define GENES_H_

#include "Bitstring.h"
#include "math.h"
#include <fstream>
#include <sstream>

#include <list>
#include <map>
#include <utility>
#include <bitset>

class Gene {
public:
	Gene();
	virtual ~Gene(){};
	bool operator ==(Gene& rhs);
	void SetGeneFunctionality(bool functionality);
	void SetGeneExpression(bool expression);
	void SetGeneSpecificity(); //works DO I NEED IT?
	void SetGeneSpecificity(int s){specificity = s;};//DO I NEED IT?
	void SetPseudogene(bool expression);
	void SetPseudogene(int M_id);
	int GetGeneMid(){return mID;};
	int GetGeneSpecificity(){return specificity;};
	bool IsFunctional(){return functional;};
	bool IsExpressed(){return isExpressed;};
	bool IsPseudogene(){return pseudoGene;};
	Gene& Copy(Gene& rhsGene);//works;
	void SetGeneID(int _ID){geneID = _ID;};
	int GetGeneID(){return geneID;};

	void Mutate(double mutationRate); //works IS OLD
	void PointMutation();//works
	void PrintBits();//works
	void MutateSpecificity();//works
	bool IsGeneUnique(vector<Gene>& genePool, int counter);//works
	int BindMolecule(Gene& anotherMolecule); //works NEW
	int CountHammingDistance(Gene& anotherMolecule);

	void SaveBackupGenes(fstream& backupFile);//works
	void SaveGenes (fstream& outfile);
	void PrintParameters();
	string RestoreGenes(stringstream& sgline);
	//pair <int,int> genePhenotype;

protected:
	bool functional;
	bool isExpressed;
	int geneID;
	int specificity;
	bool pseudoGene;
	unsigned long int mID;


};


class MHCGene: public Gene
{
public:
	MHCGene();
	virtual ~ MHCGene(){};
};

class KIRGene: public Gene
{
public:
	KIRGene(){};
	KIRGene(int L);
	virtual ~ KIRGene(){};
	enum state{inhibitory, activating};

	bool operator ==(KIRGene& rhs);
	KIRGene& Copy(KIRGene& rhsGene);//works;
	void SaveBackupGenes(fstream& backupFile); //works
	void SaveGenes (fstream& outfile);//works
	string RestoreGenes(stringstream& sgline); //works
	int GetGeneType(){return geneType;}; //works
	void SetGeneType(int type);
	void MutateReceptorType();//works
	bool IsInhibitory();
	bool IsActivating();
	void PrintGenes();


protected:
	state geneType;
};


class GenePool
{
public:
	GenePool(){};
	virtual ~GenePool(){};
	float GetAlleleFreq(int alleleIndex); //works
	bool GeneAlreadyInPool(int geneID); //works
	unsigned int GetPoolSize(){return genes.size();}
	int RandomlyPickGene(bool distribution); //works
	vector<int>& GetGenes(){return genes;};
	void WriteOutGenes();//works

protected:
	vector<int> genes;//a vector of integers
	fstream genePoolFile;
};

class MHCGenePool: public GenePool
{
public:
	MHCGenePool(){};
	void FillMHCGenePool(int size);//works
	void FillMHCGenePoolWithSimilarMHCs(int size);
	virtual ~MHCGenePool(){};
protected:
	int poolSize;
};

class KIRGenePool: public GenePool
{
public:
	KIRGenePool(GenePool& mhcPool, bool dist, int specificity);
	virtual ~KIRGenePool(){};

};


class Map : public GenePool
{
public:
	Map(){};
	virtual ~Map(){};
	void FillMap(GenePool& mhcpool, KIRGene& gene); //works
	bool IsGeneInMap(KIRGene& gene); //works
	unsigned int GetMapSize(){return mapGenes.size();} //works

	string RestoreMap(stringstream& smline); //works
	void SaveBackupMap(fstream& backupFile); //works

	multimap< pair< int,int>,  pair<int, int > >& GetMap() {return mapGenes;};

//protected:
	multimap <pair< int,int>, pair <int,int> > mapGenes; //the map has the design: (geneID,gene_type), (L, mID)
};


#endif /* GENES_H_ */

/*
 * 22.08.2012: new member variable int specificity.
 * Mutate Gene enhanced with an increase of the specificity value by one
 *
 * 30.08.2012:
 *
 * what i need to do is:
 * make a big table at the beginning, rather: make a map with:
 * KIR_id | L | M_id
 *
 *where M_id_i is the integer coding for which MHC molecules in the pool the KIR_id_i recognizes!
 *
 *05.09.2012:
 *
 *added new subclasses: KIRGene; MHCGene, to construct them differently (KIR with L = random, and MHC with L = 0)
 *added functions that directly stores the created bit string into an int.
 *BindMolecule() has been modified, and it calculates directly the longest adjacent match from the integers by using XOR
 *
 *06.09.2012
 *added new Subclass: Map. Creates a table with KIR_id | L | M_id
 * + the new functions of this new subclass
 *
 *05.10.2012
 *changes the constructor of MHCPool to an empty one, allocating only memory.
 *New function "FillMHCPool(size)" generates the pool according to the specified size
 *RandomlyPickGene was modified: the probabilites[] array was set to a max length of 50 alleles, but it will be filled, depending on the genePool.size()
 *
 *12.10.2012
 *PointMutation() implemented: it flips one random bit of the bitstring to simulate mutation of the aminoacid sequence
 *MutateSpecificity(): it changes the L value gradually to +- 1 (to simulate changes in phenotype)
 *
 *24.12.2012: Implemented KIRMap::SaveBackup and Restore functions!
 *
 */
