/*
 * Host.cpp
 *
 *  Created on: Apr 21, 2011
 *      Author: paola
 */

#include "Host.h"
/* FUNTIONS OF CLASS VIRUS
 * 	Constructs the virus*/
Virus::Virus()
{
	SetViralParameters(0.0, 0.0, 0.0, 0.0, 0, false, false,0,false);
}

void Virus::SetViralParameters(double _downregulation, double _decoy, double _viralLoad, double _lifeTimeVirus, int virus_type, bool only_acute, bool specific_mhc, int mhc_id, bool decoy_type)
{
	mutationRateDownregulation = _downregulation;
	mutationRateDecoy = _decoy;
	viralLoad = _viralLoad;
	originalViralLoad = _viralLoad;
	lifeTimeVirus = _lifeTimeVirus;
	switch(virus_type)
	{
	case 0:
		virusType = wildType;
	break;
	case 1:
		virusType = downregulating;
	break;
	}

	//every virus should express a random viral molecule
	int bla = RandomNumber(0,65000);
	viralMolecule.SetGeneID(bla);

	onlyAcute = only_acute;
	specificMHC = specific_mhc;
	downregulatedMHC.SetGeneID(mhc_id);
	isDecoyAnMHC = decoy_type;
}

/*This function sets the virus to downregulate the MHC expression*/
void Virus ::DownregulateMHC(int mhcID)
{
	if(RandomNumberDouble()<mutationRateDownregulation)
	{
		if(specificMHC)
		{
			downregulatedMHC.SetGeneID(mhcID);
		}
		virusType = downregulating;
		return;
	}
}

/*This function allows the virus to evolve an viral molecule*/
void Virus::BuildViralMolecule(int mhcID, int mutationType)
{
	if(isDecoyAnMHC) //the viral molecule could be an MHC molecule (a decoy) [we check for that here, regardless whether we passed an mhc as an argument or not]
	{
		if(RandomNumberDouble() <mutationRateDecoy)
		{
			viralMolecule.SetGeneID(mhcID);
		}
	}
	else
	{
		if(mutationType == 1)
		{
			if(RandomNumberDouble() <mutationRateDecoy) //or a random viral molecule
			{
				MHCGene random_decoy;
				viralMolecule.SetGeneID(random_decoy.GetGeneID());
			}
		}
		if(mutationType == 2)//point mutation
		{
			if(RandomNumberDouble() <mutationRateDecoy) //or a viral molecule one bit away!
			{
				viralMolecule.PointMutation();
			}
		}
	}
}

/*This function performs a deep copy*/
Virus & Virus::Copy(Virus& rhsVirus)
{
	// checking if it is a self-assignment
	if(this == &rhsVirus)
		return *this;
	//copying member variables
	this->mutationRateDownregulation = rhsVirus.mutationRateDownregulation;
	this->mutationRateDecoy = rhsVirus.mutationRateDecoy;
	this->lifeTimeVirus = rhsVirus.lifeTimeVirus;
	this->viralLoad = rhsVirus.viralLoad;
	this->originalViralLoad = rhsVirus.originalViralLoad;
	this->virusType = rhsVirus.virusType;
	this->viralMolecule.Copy(rhsVirus.viralMolecule);
	this->onlyAcute = rhsVirus.onlyAcute;
	this->downregulatedMHC.Copy(rhsVirus.downregulatedMHC);
	this->specificMHC = rhsVirus.specificMHC;
	this->isDecoyAnMHC = rhsVirus.isDecoyAnMHC;
	return *this;
}

bool Virus :: operator == (Virus& rhs)
{
	if(this->virusType == rhs.virusType)
		return true;
	else
		return false;
}

bool Virus :: IsDownregulatingMHC()
{
	if(virusType == downregulating)
		return true;
	else
		return false;
}

bool Virus :: IsWildType()
{
	if(virusType == wildType)
		return true;
	else
		return false;
}

bool Virus :: IsExpressingMolecules()
{
	if(viralMolecule.GetGeneID()) //if the virus has a viral molecule other than zero
		return true;
	else
		return false;
}

void Virus::SaveBackupVirus(fstream& file)
{
	file << mutationRateDownregulation << "\t" << mutationRateDecoy << "\t" << lifeTimeVirus << "\t" << viralLoad <<"\t" <<originalViralLoad<<"\t"<< virusType <<"\t" << viralMolecule.GetGeneID()<< "\t"<<onlyAcute <<"\t";
	file << downregulatedMHC.GetGeneID() << "\t" <<specificMHC <<"\t"<<isDecoyAnMHC << "\t";
}

void Virus::PrintParametersVirus()
{
	cout << mutationRateDownregulation << "\t" << mutationRateDecoy << "\t" << lifeTimeVirus << "\t" << viralLoad <<"\t" <<originalViralLoad<<"\t"<< virusType <<"\t" << viralMolecule.GetGeneID()<< "\t"<<onlyAcute <<"\t";
	cout << downregulatedMHC.GetGeneID() << "\t" <<specificMHC << "\t"<<isDecoyAnMHC << "\t";
}

string Virus::RestoreVirus(stringstream& svline)
{
	//cout <<"in the viral class: "<<svline.str()<<endl;
	int type;
	int id;
	int other_id;
	svline >> mutationRateDownregulation;
	svline >> mutationRateDecoy;
	svline >> lifeTimeVirus;
	svline >> viralLoad;
	svline >> originalViralLoad;
	svline >> type;
	switch(type)
	{
	case 0:
		virusType = wildType;
	break;
	case 1:
		virusType = downregulating;
	break;
	}
	svline >> id;
	viralMolecule.SetGeneID(id);
	svline >> onlyAcute;

	//out commented JUST for the version of 15max loci!!!! (it was an old version, and it didn't have these last parameters yet!
	svline >> other_id;
	downregulatedMHC.SetGeneID(other_id);
	svline >> specificMHC;
	svline >> isDecoyAnMHC;

	string vstring;
	svline.clear();
	getline(svline,vstring);
	//cout <<">>>>>>>>>>>>>>>VIRUS string:  "<<vstring<<endl;
	return vstring;
}

void Virus :: SaveParametersVirus(fstream& outfile)
{
	outfile << virusType << "\t" << viralLoad <<"\t" <<originalViralLoad << "\t"<< viralMolecule.GetGeneID() <<"\t"<<onlyAcute<<"\t";
	outfile << downregulatedMHC.GetGeneID() << "\t" <<specificMHC <<"\t" << isDecoyAnMHC << "\t";
}


/*FUNCTIONS OF CLASS HOST
 *Constructs a host for the initialization of the population: it fills the MHC genes with a randomly picked allele from the population
 * and creates KIRs that match their own MHC according to the specificity*/
Host::Host(int loci_kir, int loci_mhc, double _mutationRate, bool _tuning, int numberOfExtraKirs,Map& kirMap, MHCGenePool& mhcPool, bool hla) {

	InitializeHostParameters(_mutationRate,_tuning, loci_kir, loci_mhc);
	//fill the mhc Genes
	for(int i = 0; i <LOCI_MHC*TWO; i ++)
	{
		Gene firstGene;
		int mhc1 = mhcPool.RandomlyPickGene(hla);
		firstGene.SetGeneID(mhc1);
		mhcGenes.push_back(firstGene);
	}

/*
	Gene secondGene;
	secondGene.SetGeneID(mhc2);
	mhcGenes.push_back(secondGene);*/

	//create random KIRs with random specificity for ONE haplotype
	// (at the beginning of the simulation, the size of the map should equal the LOCI_NUMBER)
	map< pair<int,int>, pair <int, int> > ::iterator it;
	for(it = kirMap.GetMap().begin(); it != kirMap.GetMap().end(); it ++)
	{
		int id = it->first.first;
		int geneType = it->first.second;
		int L = it->second.first;
		int pseudo = it->second.second;
		KIRGene kir;
		kir.SetGeneSpecificity(L);
		kir.SetGeneID(id);
		kir.SetPseudogene(pseudo);
		kir.SetGeneType(geneType);
		kir.SetGeneFunctionality(false);
		kir.SetGeneExpression(false);
		kirGenes.push_back(kir);
		//cout << "printing genes in host first constructor" <<endl;
		//kir.PrintGenes();
	}

	int size = kirGenes.size();
	//make sure that the first haplotype has number LOCI_KIR (regardless of the Map size!)
	while(size < LOCI_KIR)
	{
		KIRGene bla = kirGenes.at(size-1);
		kirGenes.push_back(bla);
		size ++;
	}

	//copy the first kirs into the other haplotype (all individuals are homozygous!)
	for(int i=0; i<LOCI_KIR; i++)
	{
		KIRGene kir = kirGenes.at(i);
		kirGenes.push_back(kir);

	}

	if(tuning == true)
		EducateKIRs();
	ExpressKIRs(numberOfExtraKirs);
	//CountFunctionalKIRs();
	age = RandomNumber(1,70); //population initialized with a random age between 1 and 70
	//cout <<inhibitoryKIRs << "|" <<activatingKIRs <<"|"<<CountExpressedKIRs() <<endl;
}

/*Constructs a baby host out of two parents*/
Host::Host(int loci_kir, int loci_mhc, vector<Gene>& mhcGenesParent, GenePool& mhcPool, bool dist, vector<KIRGene>& kirGenesMother, vector<KIRGene>& kirGenesFather,double _mutationRate, bool _tuning, int numberOfExtraKirs, Map& kirMap, int mutationType, int gene_type, double simulation_time, double time_invation, bool invasion_analysis)
//Host::Host(int loci_kir, vector<Gene>& mhcGenesParent, GenePool& mhcPool, bool dist, GenePool& kirPool, vector<Gene>& kirGenesMother, vector<Gene>& kirGenesFather, int specificity, double _mutationRate, bool _tuning, int numberOfExtraKirs)
{	//to create a NEW host: the haplotypes of KIR of BOTH parents are needed. Besides one MHC haplotype of one parent plus one of the pool

	InitializeHostParameters(_mutationRate, _tuning, loci_kir, loci_mhc);
	 //cout<< "do I get stuck here?? | Host constructor!"<<endl;

	int KIR_init_mum = 0;
	int KIR_end_mum = 0;
	int MHC_init_mum = 0;
	int MHC_end_mum = 0;

	int KIR_init_dad = 0;
	int KIR_end_dad = 0;

	//pick haplotype 1/0 of each parent for the KIRs!
	bool hap_mum=(RandomNumberDouble()<0.5);
	bool hap_dad= (RandomNumberDouble()<0.5);
	//int k= 0;
	if(hap_mum)
	{
		KIR_init_mum = LOCI_KIR;
		KIR_end_mum = LOCI_KIR*TWO;
	}
	else
	{
		KIR_init_mum = 0;
 		KIR_end_mum = LOCI_KIR;
	}

	if(hap_dad)
	{
		KIR_init_dad = LOCI_KIR;
		KIR_end_dad = LOCI_KIR*TWO;
	}

	else
	{
		KIR_init_dad = 0;
		KIR_end_dad = LOCI_KIR;
	}

	//pick haplotype 1/0 of the parent for the MHCs!
	bool hap_mhc = (RandomNumberDouble()<0.5);
	if(hap_mhc)
	{
		MHC_init_mum = LOCI_MHC;
		MHC_end_mum = LOCI_MHC*TWO;
	}

	else
	{
		MHC_init_mum = 0;
		MHC_end_mum = LOCI_MHC;
	}

	//copy the KIR haplotype into the new host (mutation occurs!)
	bool hap_child = (RandomNumberDouble()<0.5);
	if(hap_child)
	{
		for(int i=KIR_init_dad; i<KIR_end_dad; i++)
		{
			KIRGene kir_hap2;
			kir_hap2.Copy(kirGenesFather.at(i));
			//MUTATION!
			if(RandomNumberDouble() < mutationRateHost)
			{
				if(!invasion_analysis)//is there a mutual invasion analysis?
					MutateGenes(mutationType, kir_hap2, kirMap, mhcPool, gene_type); //if no, mutate regularly
				else
					MutateGenesForMutualInvasion(mutationType, kir_hap2, kirMap, mhcPool,simulation_time, time_invation, gene_type);
			}
			kirGenes.push_back(kir_hap2);
		}
		for(int j=KIR_init_mum; j<KIR_end_mum; j++)
		{
			KIRGene kir_hap1;
			kir_hap1.Copy(kirGenesMother.at(j));
			if(RandomNumberDouble() < mutationRateHost)
			{
				if(!invasion_analysis)//is there a mutual invasion analysis?
					MutateGenes(mutationType, kir_hap1, kirMap, mhcPool, gene_type); //if no, mutate regularly
				else
					MutateGenesForMutualInvasion(mutationType, kir_hap1, kirMap, mhcPool,simulation_time, time_invation, gene_type);
			}
			kirGenes.push_back(kir_hap1);
		}

		//copy the MHC haplotype into the new host
		for(int m = MHC_init_mum; m < MHC_end_mum; m++)
		{
			Gene mhc1;
			mhc1.Copy(mhcGenesParent.at(m));
			mhcGenes.push_back(mhc1);
			//cout << mhc1.GetGeneID() << "_" << mhc1.GetGeneSpecificity() << " ";
		}

		//generate new mhc genes (from the pool) for the new born
		for(int mm = 0; mm<LOCI_MHC; mm++)
		{
			Gene mhc2;
			int mhcFromTheGenePool = mhcPool.RandomlyPickGene(dist);
			mhc2.SetGeneID(mhcFromTheGenePool);
			mhcGenes.push_back(mhc2);
		}
	}
	else
	{
		for(int i=KIR_init_mum; i<KIR_end_mum; i++)
		{
			KIRGene kir_hap1;
			kir_hap1.Copy(kirGenesMother.at(i));
			if(RandomNumberDouble() < mutationRateHost)
			{
				if(!invasion_analysis)//is there a mutual invasion analysis?
					MutateGenes(mutationType, kir_hap1, kirMap, mhcPool, gene_type); //if no, mutate regularly
				else
					MutateGenesForMutualInvasion(mutationType, kir_hap1, kirMap, mhcPool,simulation_time, time_invation, gene_type);
			}
			kirGenes.push_back(kir_hap1);
		}
		for(int j=KIR_init_dad; j<KIR_end_dad; j++)
		{
			KIRGene kir_hap2;
			kir_hap2.Copy(kirGenesFather.at(j));
			if(RandomNumberDouble() < mutationRateHost)
			{
				if(!invasion_analysis)//is there a mutual invasion analysis?
					MutateGenes(mutationType, kir_hap2, kirMap, mhcPool, gene_type); //if no, mutate regularly
				else
					MutateGenesForMutualInvasion(mutationType, kir_hap2, kirMap, mhcPool,simulation_time, time_invation, gene_type);
			}
			kirGenes.push_back(kir_hap2);
		}
		//copy the MHC haplotype into the new host
		//generate new mhc genes for the new born

		//generate new mhc genes (from the pool) for the new born
		for(int mm = 0; mm<LOCI_MHC; mm++)
		{
			Gene mhc2;
			int mhcFromTheGenePool = mhcPool.RandomlyPickGene(dist);
			mhc2.SetGeneID(mhcFromTheGenePool);
			mhcGenes.push_back(mhc2);
		}
		//copy the MHC haplotype into the new host
		for(int m = MHC_init_mum; m < MHC_end_mum; m++)
		{
			Gene mhc1;
			mhc1.Copy(mhcGenesParent.at(m));
			mhcGenes.push_back(mhc1);
			//cout << mhc1.GetGeneID() << "_" << mhc1.GetGeneSpecificity() << " ";
		}

	}
	ResetGenes(); //reset the expression and functionality back to zero and recalculate these parameters in the education function
	if(tuning == true)
	{
		EducateKIRs();
	}
	ExpressKIRs(numberOfExtraKirs);
	//cout <<inhibitoryKIRs << "|" <<activatingKIRs <<"|"<<CountExpressedKIRs() <<endl;
	age = 1.0; //newborns are given the age of 15
}

void Host::InitializeHostParameters(double mutationRate, bool _tuning, int loci_kir, int loci_mhc)
{
	dead = false;
	mutationRateHost = mutationRate;
	tuning = _tuning;
	mainInfectionType = susceptible;
	LOCI_KIR = loci_kir;
	LOCI_MHC = loci_mhc;
	viralDeathRate = 0.0;
	age = 0.0;
	inhibitoryKIRs = 0;
	activatingKIRs = 0;
	numberOfInfections = 0;
	//MHCsearchinPool = 0;

}

void Host::ResetGenes()
{
	vector<KIRGene>::iterator kirIt;
	for(kirIt = kirGenes.begin(); kirIt!=kirGenes.end(); kirIt++)
	{
		kirIt->SetGeneFunctionality(false);
		kirIt->SetGeneExpression(false);
	}
}
void Host :: MutateGenes(int mutationType, KIRGene& kir_hap2, Map& kirMap, GenePool& mhcPool, int gene_type)
{
	if(mutationType == 1)//pick another molecules as mutation
	{
		KIRGene newGene(RandomNumber(2,16));
		/*if(!kirMap.IsGeneInMap(newGene))
		{
			kirMap.FillMap(mhcPool, newGene);
			if(gene_type !=2)////force to have only one type of receptors, if the user wants it!
				newGene.SetGeneType(gene_type);
			kir_hap2.Copy(newGene);
		}//*/
		int M_id = 0;
		int mhcPoolSize = mhcPool.GetPoolSize();
		
		//calculate the value of M_id to determine whether the gene is pseudogene or not
		for(unsigned int i = 0; i < mhcPoolSize; i++)
		{
			Gene mhcGene;
			mhcGene.SetGeneID(mhcPool.GetGenes().at(i));
			int L = kir_hap2.BindMolecule(mhcGene);
			if(L >= kir_hap2.GetGeneSpecificity())
				M_id += (1<<i);
		}
		
		if(gene_type !=2)////force to have only one type of receptors, if the user wants it!
			newGene.SetGeneType(gene_type);
		kir_hap2.Copy(newGene);
	}

	if(mutationType == 2)//point mutation + L
	{

		if(RandomNumberDouble()<0.8) //pointmutation
		{
			kir_hap2.PointMutation();
			//kirMap.FillMap(mhcPool, kir_hap2);
			int M_id = 0;
			int mhcPoolSize = mhcPool.GetPoolSize();
			
			//calculate the value of M_id to determine whether the gene is pseudogene or not
			for(unsigned int i = 0; i < mhcPoolSize; i++)
			{
				Gene mhcGene;
				mhcGene.SetGeneID(mhcPool.GetGenes().at(i));
				int L = kir_hap2.BindMolecule(mhcGene);
				if(L >= kir_hap2.GetGeneSpecificity())
					M_id += (1<<i);
			}
		}

		if(RandomNumberDouble()<0.8) //mutate L
		{
			kir_hap2.MutateSpecificity();
			//kirMap.FillMap(mhcPool, kir_hap2);
			int M_id = 0;
			int mhcPoolSize = mhcPool.GetPoolSize();
			
			//calculate the value of M_id to determine whether the gene is pseudogene or not
			for(unsigned int i = 0; i < mhcPoolSize; i++)
			{
				Gene mhcGene;
				mhcGene.SetGeneID(mhcPool.GetGenes().at(i));
				int L = kir_hap2.BindMolecule(mhcGene);
				if(L >= kir_hap2.GetGeneSpecificity())
					M_id += (1<<i);
			}
		}
		if(RandomNumberDouble()<0.8)//mutate type
		{
			kir_hap2.MutateReceptorType();
		}
	}
}

void Host::MutateGenesForMutualInvasion(int mutationType, KIRGene& kir_hap2, Map& kirMap, GenePool& mhcPool, double simulationTime, double time_invasion, int gene_type)
{

	if(mutationType == 1)//pick another molecule as mutation
	{
		KIRGene newGene(RandomNumber(2,16));
		/*if(!kirMap.IsGeneInMap(newGene))
		{
			kirMap.FillMap(mhcPool, newGene);
			//cout <<simulationTime << "|" <<time_invasion <<endl;
			if(simulationTime < time_invasion)//only after the invasion time, the receptor type should be allowed to mutate
				newGene.SetGeneType(gene_type);
			kir_hap2.Copy(newGene);
		}//*/
		int M_id = 0;
		int mhcPoolSize = mhcPool.GetPoolSize();
		
		//calculate the value of M_id to determine whether the gene is pseudogene or not
		for(unsigned int i = 0; i < mhcPoolSize; i++)
		{
			Gene mhcGene;
			mhcGene.SetGeneID(mhcPool.GetGenes().at(i));
			int L = kir_hap2.BindMolecule(mhcGene);
			if(L >= kir_hap2.GetGeneSpecificity())
				M_id += (1<<i);
		}
		//set the Gene pseudo
		newGene.SetPseudogene(M_id);
		if(gene_type !=2)////force to have only one type of receptors, if the user wants it!
			newGene.SetGeneType(gene_type);
		kir_hap2.Copy(newGene);
		return;
	}

	if(mutationType == 2)//point mutation + L
	{
		if(RandomNumberDouble()<0.8)
		{
			kir_hap2.PointMutation();
			//kirMap.FillMap(mhcPool, kir_hap2);
			int M_id = 0;
			int mhcPoolSize = mhcPool.GetPoolSize();
			
			//calculate the value of M_id to determine whether the gene is pseudogene or not
			for(unsigned int i = 0; i < mhcPoolSize; i++)
			{
				Gene mhcGene;
				mhcGene.SetGeneID(mhcPool.GetGenes().at(i));
				int L = kir_hap2.BindMolecule(mhcGene);
				if(L >= kir_hap2.GetGeneSpecificity())
					M_id+= (1<<i);
			}
			//set the Gene pseudo
			kir_hap2.SetPseudogene(M_id);
		}

		if(RandomNumberDouble()<0.2)
		{
			kir_hap2.MutateSpecificity();
			//kirMap.FillMap(mhcPool, kir_hap2);
			int M_id = 0;
			int mhcPoolSize = mhcPool.GetPoolSize();
			
			//calculate the value of M_id to determine whether the gene is pseudogene or not
			for(unsigned int i = 0; i < mhcPoolSize; i++)
			{
				Gene mhcGene;
				mhcGene.SetGeneID(mhcPool.GetGenes().at(i));
				int L = kir_hap2.BindMolecule(mhcGene);
				if(L >= kir_hap2.GetGeneSpecificity())
					M_id+= (1<<i);
			}
			//set the Gene pseudo
			kir_hap2.SetPseudogene(M_id);
		}
		if(RandomNumberDouble()<0.2)
		{
			if(simulationTime >= time_invasion) //only after the invasion time, the receptor type should be allowed to mutate
				kir_hap2.MutateReceptorType();
		}
		return;
	}
}

/*This functions tunes the KIR repertoire according to the self MHC repertoire and whether they are inhibiting or activating*/
void Host :: EducateKIRs()
{
	vector<KIRGene>::iterator kirIt;
	vector<Gene>::iterator mhcIt;
	//cout <<"hello...educating the kirs..."<<endl;
	for(kirIt = kirGenes.begin(); kirIt !=kirGenes.end(); kirIt ++)
	{

		//kirIt->PrintGenes();
		if(kirIt->IsPseudogene()||kirIt->IsActivating()) //ignore pseudogenes and consider only iNKRs
			continue;
		//kirIt->PrintGenes();cout<<endl;
		for(mhcIt = mhcGenes.begin(); mhcIt !=mhcGenes.end(); mhcIt ++)
		{
			int bindingStrength = kirIt->BindMolecule(*mhcIt);
			//mhcIt->PrintBits();
			//kirIt->PrintBits();
			if(bindingStrength>=kirIt->GetGeneSpecificity()) //AND it binds to the MHC
			{
				//cout << "it binds!"<<endl;
				//kirIt->PrintGenes();
				kirIt->SetGeneFunctionality(true); //it is licensed!
				kirIt->SetGeneExpression(true);
				//cout<< "inhibitory |" <<kirIt->IsFunctional()<< kirIt->GetGeneType()<<kirIt->IsExpressed()<<endl;
			}
			else // but if it doesn't bind
			{
				//cout << "it doesn't bind!"<<endl;
				kirIt->SetGeneFunctionality(false); // it should not be licensed
				kirIt->SetGeneExpression(false);
				//cout<< "inhibitory |" <<kirIt->IsFunctional()<< kirIt->GetGeneType() <<kirIt->IsExpressed()<<endl;
			}
		}
	}
	CountFunctionalInhibitoryKIRs();
	//only if there is at least one licensed iNKRs, the NK cells are functional..
	// otherwise, we consider cell having only aNKRs to be anergic/hyporesponsive
	if(inhibitoryKIRs)
	{
		//cout <<"yeeeeeaaaaa, at least one licensed iNKR, now the NK cells are functional and we can start looking at the aNKR"<<endl;
		for(kirIt=kirGenes.begin(); kirIt!=kirGenes.end(); kirIt++)
		{
			if(kirIt->IsInhibitory())
				continue;
			if (kirIt->IsPseudogene())
				continue;
			int education_signal_activating = 0; // to keep track of the activating signal per receptor
			for(mhcIt=mhcGenes.begin(); mhcIt!=mhcGenes.end(); mhcIt++)
			{
				int bindingStrength = kirIt->BindMolecule(*mhcIt);
				//mhcIt->PrintBits();
				//kirIt->PrintBits();
				if(bindingStrength>=kirIt->GetGeneSpecificity()) // if an aNKR binds to the MHC
				{
					//cout << "it binds!"<<endl;
					kirIt->SetGeneFunctionality(false); //this gene is unlicensed
					kirIt->SetGeneExpression(false);
					//cout<< "activating |" <<kirIt->IsFunctional()<< kirIt->GetGeneType()<<endl;
				}
				else //BUT if it deosn't bind
				{
					//cout << "it doesn't bind!"<<endl;
					education_signal_activating ++; //keep track of how many MHC it doesn't bind
					//cout << education_signal_activating <<endl;
					//cout<< "activating |" <<kirIt->IsFunctional()<< kirIt->GetGeneType()<<endl;
				}
			}
			if(education_signal_activating == mhcGenes.size()) //if the activating KIR doesn't recognize ANY of the MHC
			{
				//cout << education_signal_activating << "|" << mhcGenes.size() << ": ";
				kirIt->SetGeneFunctionality(true); //the gene is licensed
				kirIt->SetGeneExpression(true);
				//cout<< "activating |"<<kirIt->IsFunctional()<< kirIt->GetGeneType()<<endl;
			}
		}
		CountFunctionalActivatingKIRs();
	}
	/*cout <<"after educating all of them"<<endl;
	for(kirIt=kirGenes.begin(); kirIt!=kirGenes.end();kirIt++)
	{
		kirIt->PrintGenes();
	}*/
	//exit(-1);
}

void Host :: ExpressKIRs(int numberOfExtraKirs)
{
	int counter = 0;
	vector<KIRGene>::iterator kirIt;
	if(tuning == false)
	{
		for(kirIt = kirGenes.begin(); kirIt !=kirGenes.end(); kirIt ++)
			kirIt->SetGeneExpression(true);
		// if there is no education, all KIRs should be expressed!
	}
	else //otherwise express more KIRs besides those which are already functional
	{
		for(kirIt = kirGenes.begin(); kirIt !=kirGenes.end(); kirIt ++)
		{
			if(kirIt->IsFunctional()) //KIRs that are functional are already expressed (see EducateKIRs!)
				continue;
			else
			{
				if(counter<numberOfExtraKirs)
				{
					kirIt->SetGeneExpression(true);
					counter++;
				}
				else
					return;
			}
		}
	}
}

/* This function counts the UNIQUE KIRs within one host that are functional, i.e. that are able to recognize their own MHC*/
void Host::CountFunctionalInhibitoryKIRs()
{
	//cout <<"before counting: "<<inhibitoryKIRs <<endl;
	for(unsigned int i= 0; i<kirGenes.size(); i++)
	{
		if(kirGenes.at(i).IsActivating()) //ignore non-functional KIRs, i.e. aNKRs and pseudo ones
			continue;
		if(kirGenes.at(i).IsPseudogene())
			continue;
		if(kirGenes.at(i).IsFunctional())
		{
			int flagNotUnique = 1;
			for(unsigned int j= 0; j<i; j++)
			{
				if(kirGenes.at(i).GetGeneID()==kirGenes.at(j).GetGeneID())
					flagNotUnique++;
			}
			if (flagNotUnique == 1)
				inhibitoryKIRs ++;
		}
	}
	//cout <<"after counting: "<<inhibitoryKIRs <<endl;
}

void Host::CountFunctionalActivatingKIRs()
{
	//cout <<"before counting aNKRs: "<<activatingKIRs <<endl;
	for(unsigned int i= 0; i<kirGenes.size(); i++)
	{
		if(kirGenes.at(i).IsInhibitory())//ignore non-functional KIRs
			continue;
		if(kirGenes.at(i).IsPseudogene())
			continue;
		if(kirGenes.at(i).IsFunctional())
		{
			int flagNotUnique = 1;
			for(unsigned int j= 0; j<i; j++)
			{
				if(kirGenes.at(i).GetGeneID()==kirGenes.at(j).GetGeneID())
					flagNotUnique++;
			}
			if (flagNotUnique == 1)
				activatingKIRs ++;
		}
	}
	//cout <<"after counting aNKRs: "<<activatingKIRs <<endl;
}

/* This function counts ALL KIRs within one host that are expressed*/
int Host::CountExpressedKIRs()
{
	int kirsThatAreExpressed = 0;
	vector<KIRGene>::iterator kirIt;
	for(kirIt = kirGenes.begin(); kirIt !=kirGenes.end(); kirIt ++)
	{
		//if(kirIt->IsPseudogene())
			//continue;
		if(kirIt->IsExpressed())
		{
			kirsThatAreExpressed ++;
		}
	}

	return kirsThatAreExpressed;
}
/*This function performs a deep copy*/
Host& Host:: Copy(Host& rightHandSideHost)
{
	// checking if it is a self-assignment
	if(this == &rightHandSideHost)
		return *this;

	//copying member variables
	this->dead = rightHandSideHost.dead;
	this->tuning = rightHandSideHost.tuning;
	this->age = rightHandSideHost.age;
	this->LOCI_KIR = rightHandSideHost.LOCI_KIR;
	this->LOCI_MHC = rightHandSideHost.LOCI_MHC;
	this->mainInfectionType = rightHandSideHost.mainInfectionType;
	this->mutationRateHost = rightHandSideHost.mutationRateHost;
	this->viralDeathRate = rightHandSideHost.viralDeathRate;
	this->inhibitoryKIRs = rightHandSideHost.inhibitoryKIRs;
	this->activatingKIRs = rightHandSideHost.activatingKIRs;
	this->numberOfInfections = rightHandSideHost.numberOfInfections;

	//copy genes
	for(unsigned int i=0; i<rightHandSideHost.mhcGenes.size(); i++)
	{
		this->mhcGenes.push_back(rightHandSideHost.mhcGenes.at(i));
		//cout << "mhc genes: \t"<<mhcGenes.at(i).GetGeneID() << "\t" << rightHandSideHost.mhcGenes.at(i).GetGeneID() <<endl;
	}
	for(unsigned int i=0; i<rightHandSideHost.kirGenes.size(); i++)
	{
		this->kirGenes.push_back(rightHandSideHost.kirGenes.at(i));
	}
	
	list<Infection>::iterator it;
	for (it = rightHandSideHost.infections.begin() ; it !=rightHandSideHost.infections.end(); it++)
	{
		Infection dummy;
		dummy.Copy((*it));
		this->infections.push_back(dummy);
	}
	//check whether it works! do i need a copy function for the infection class???copy infections TO DO!!!!!!!!!!!!!!!!!!1

	return *this; // returns self-reference so cascaded assignment works
}

/*this function sets the vector of infections*/
void Host::InfectWith(Virus& newVirus, double simulationTime, int maxNumberInfections) //function of the host receiving a new virus
{
	//cout << "do i get stuck here??? InfectWith()"<<endl;
	Infection newInfection;
	list<Infection>::iterator it;
	list<Infection>::iterator end = infections.end();

	if(!IsInfected()) // if the host is NOT infected yet
	{
		newInfection.TransmitInfection(newVirus,simulationTime);//set the parameters to the new infection!
		infections.push_back(newInfection);
	}
	else //but if the host IS infected with some viruses
	{
		if(infections.size()<maxNumberInfections) //check first whether it is already at its maximum
		{
			//*
			bool isInTheInfectionsList = false;
			for (it = infections.begin(); it != end; it++)
			{
				if(!it->IsPathogenNew(newVirus)) // if the virus is not new  (i.e. as incubating, acute, chronic or immune) ignore it
				{
					isInTheInfectionsList=true;
				}
			}
			if(!isInTheInfectionsList)
			{
				newInfection.TransmitInfection(newVirus,simulationTime);//set the parameters to the new infection!
				infections.push_back(newInfection);
			}
			/*
			 int howManyInfections = 0;
			 for (it = infections.begin(); it != end; it++)
			 {
			 if(!it->IsPathogenNew(newVirus)) // if the virus is already there (i.e. as incubating, acute, chronic or immune) ignore it
			 continue;
			 else //but if it's not the same, keep track of how many infections are different from the new one
			 howManyInfections++;
			 }
			 
			 if(howManyInfections == infections.size()) // if the new virus is different from ALL infections present in that host
			 {
			 newInfection.TransmitInfection(newVirus,simulationTime);//set the parameters to the new infection!
			 infections.push_back(newInfection);
			 }
			 //*/
		}
	}
}

void Host:: InfectWithMoreDownregulating(Virus& newVirus, double simulationTime, int maxNumberInfections) //function of the host receiving a new virus
{
	//cout << "do i get stuck here??? InfectWith()"<<endl;
	Infection newInfection;
	list<Infection>::iterator it;
	list<Infection>::iterator end = infections.end();

	if(!IsInfected()) // if the host is NOT infected yet
	{
		newInfection.TransmitInfection(newVirus,simulationTime);//set the parameters to the new infection!
		infections.push_back(newInfection);
	}
	else //but if the host IS infected with some viruses
	{
		if(infections.size()<maxNumberInfections) //check first whether it is already at its maximum
		{
			if(newVirus.IsDownregulatingMHC())
			{
				newInfection.TransmitInfection(newVirus,simulationTime);//set the parameters to the new infection!
				infections.push_back(newInfection);
			}
		}
	}
}

void Host:: InfectWithMoreWildType(Virus& newVirus, double simulationTime, int maxNumberInfections) //function of the host receiving a new virus
{
	//cout << "do i get stuck here??? InfectWith()"<<endl;
	Infection newInfection;
	list<Infection>::iterator it;
	list<Infection>::iterator end = infections.end();

	if(!IsInfected()) // if the host is NOT infected yet
	{
		newInfection.TransmitInfection(newVirus,simulationTime);//set the parameters to the new infection!
		infections.push_back(newInfection);
	}
	else //but if the host IS infected with some viruses
	{
		if(infections.size()<maxNumberInfections) //check first whether it is already at its maximum
		{
			if(newVirus.IsWildType())
			{
				newInfection.TransmitInfection(newVirus,simulationTime);//set the parameters to the new infection!
				infections.push_back(newInfection);
			}
		}
	}
}

int Host :: IsInfected() //this functions basically counts whether there is more than zero infections!
{
	return infections.size();
}

/*This functions sets each host with an infection type according to the "main "state of all his infections*/
int Host::GetMainInfectionType()
{
	list<Infection>::iterator it;
	int chronicInfections = 0;
	int recoveredInfections = 0;

	for(it = infections.begin(); it!= infections.end(); it++)
	{
		int type = it->GetInfectionType();
		switch(type)
		{
		case 0:
			mainInfectionType = susceptible;
		break; //is incubating... so is not really infectious yet!

		case 1:
			mainInfectionType = acute;
			return mainInfectionType;//if only one of them is acute return it as acute and leave the loop!
		break;

		case 2:
			chronicInfections++;
		break; // is chronic

		case 3:recoveredInfections++;
		break; //is memory

		case 4:mainInfectionType = susceptible;
		break;//is cleared
		}
	}
	if(chronicInfections > 0)
		mainInfectionType = chronic;
	else
	{
		if(recoveredInfections > 0)
			mainInfectionType = immune;
	}
	return mainInfectionType;
}


void Host::ClearInfection(double simulationTime, Infection& _infection)
{
	int virusType = _infection.pathogen.GetVirusType();
	switch (virusType)
	{
		case 0: //if it's a wild-type virus, check first whether the viral molecule engages an aNKR
		{
			vector<KIRGene> :: iterator firstIt;
			int aNKRs_bindingViralMol = 0;
			for(firstIt = kirGenes.begin(); firstIt !=kirGenes.end(); firstIt++)
			{
				if(!firstIt->IsExpressed())
					continue;

				int score = firstIt->BindMolecule(_infection.pathogen.viralMolecule);
				if(score>=firstIt->GetGeneSpecificity()) //check if the NKRs bind the decoy
				{
					if(firstIt->IsActivating()) //if yes, counte them only if they are activating
							aNKRs_bindingViralMol++;
				}
			}
			if(aNKRs_bindingViralMol)
			{
				if(RandomNumberDouble()<0.95)//if there is at least one aNKR that binds the viral protein
				{
					_infection.ResetInfection(simulationTime);//there should be maximal protection
					_infection.SetProtectionLevel("best");
				}
				//else
				//	_infection.SetProtectionLevel("zero");
			}
			else //if there not aNKR binding the viral protein
			{
				if(RandomNumberDouble()<0.8) //the protection is still good (like the wild-type virus)
				{
					_infection.ResetInfection(simulationTime);//because of responses of the adaptive immune system
					_infection.SetProtectionLevel("high");
				}
				//else
					//_infection.SetProtectionLevel("zero");
			}

		}break;
		case 1: //if it's an mHC downregulating virus, check first whether there is only one MHC that has been downregulated
		{
			if(_infection.pathogen.IsDownregulatingMHCSpecifically())
			{
				cout<<"specific MHC downregulation is still on construction......"<<endl; exit(-1);
			}
			else //if ALL MHCs are being downregulated
			{
				vector<KIRGene>::iterator nkrIt; //check whether any iKIR is recognizing it!
				int inhibiting_kirs_recognizing_decoy = 0;
				int activating_kirs_recognizing_decoy = 0;
				//go through each receptor and see whether it binds to the viral molecule!
				for(nkrIt = kirGenes.begin(); nkrIt !=kirGenes.end(); nkrIt ++)
				{
					if(!nkrIt->IsExpressed()) //ignore KIRs that are not functional / expressed!
						continue;
					else
					{
						int score = nkrIt->BindMolecule(_infection.pathogen.viralMolecule); //check if they bind to the decoy!
						if(score>=nkrIt->GetGeneSpecificity()) //receptors binds to decoy
						{
							//now check what kind of receptors we have
							if(nkrIt->IsInhibitory()) //if it is inhibiting
								inhibiting_kirs_recognizing_decoy++;//keep track of how many KIRs are recognizing decoys

							if(nkrIt->IsActivating()) //if it is activating
								activating_kirs_recognizing_decoy ++;//keep track of how many KIRs are recognizing decoys
						}
					}
				}
				if(inhibitoryKIRs && !activatingKIRs) //if that host has ONLY inhibitory receptors
					ClearDecoyWithInhibitoryOnly(inhibiting_kirs_recognizing_decoy, simulationTime, _infection);

				 //if(!inhibitoryKIRs && activatingKIRs) //if that host has ONLY activating receptors
					// ClearDecoyWithActivatingOnly(activating_kirs_recognizing_decoy, simulationTime, _infection);


				if(inhibitoryKIRs && activatingKIRs) //if that host has both types of receptors
					ClearDecoyWithActivatingAndInhibitory(inhibiting_kirs_recognizing_decoy, activating_kirs_recognizing_decoy, simulationTime, _infection);
			}
		}break;
	}
}

void Host :: ClearDecoyWithInhibitoryOnly(int inhibiting_signal, double simulationTime, Infection& _infection)
{
	if(!inhibiting_signal) //if there is no inhibiting signal, receptor didn't bind to the decoy: protection like MHC down
	{
		if(RandomNumberDouble()<0.45)
		{
			_infection.ResetInfection(simulationTime);
			_infection.SetProtectionLevel("medium");
			//cout <<"clearing infection with iNKRs.."<<endl;
		}
			
	}
}

void Host :: ClearDecoyWithActivatingOnly(int activating_signal, double simulationTime, Infection& _infection)
{
	//cout <<"clearing decoy: "<<inhibitoryKIRs << "|"<< activatingKIRs << endl;
	int counter = 0;
	if(activating_signal) //if activating receptor recognizes decoy, but there are no inhibitory receptors, the virus still escapes response of the T-cells
	{                       // protection as with an MHC-downregulating one
		if(RandomNumberDouble()<0.25)
		{
			_infection.ResetInfection(simulationTime);
			_infection.SetProtectionLevel("low");

		}
	}
}

void Host ::ClearDecoyWithActivatingAndInhibitory(int inhibiting_signal, int activating_signal, double simulationTime, Infection& _infection)
{
	//if functional KIRs recognize MHC down
	//cout <<"clearing decoy: "<<inhibitoryKIRs << "|"<< activatingKIRs << endl;
	if(inhibitoryKIRs > 0 && activating_signal) //and there is enough activating signal
	{
		if(!inhibiting_signal) //best protection
		{
			if(RandomNumberDouble()<0.7)
			{
				_infection.ResetInfection(simulationTime);
				//cout <<"having the best protection.."<<endl;
				_infection.SetProtectionLevel("normal");
				return;
			}
		}
		if(inhibiting_signal) //protection provided only by the aKIRs
		{
			if(RandomNumberDouble()<0.25)
			{
				//cout <<"good protection!" <<endl;
				_infection.ResetInfection(simulationTime);
				_infection.SetProtectionLevel("low");
				return;
			}
		}
	}
	//if functional KIRs recognize MHC down
	if(inhibitoryKIRs > 0 && !activating_signal)	//but there is not enough activation OR only inhibitory receptors!
	{
		//check whether the host is fooled or not
		if(!inhibiting_signal) // if no inhibitory receptor recognizes the decoy:
		{                                             //same as MHC downregulation p= 0.6
			if(RandomNumberDouble()<0.45)
			{
				_infection.ResetInfection(simulationTime);
				_infection.SetProtectionLevel("medium");
				return;
			}
		}
		else //if one iNKRs recognizes the decoy, then there is zero protection!
			_infection.SetProtectionLevel("zero");
	}
}

/*This functions returns the Virus of the acute infection*/
Virus& Host :: GetAcuteInfection(Virus& dummy)
{
	list<Infection>:: iterator it;
	vector<Virus> randomAcuteInfections;
	int inf_type = 0;
	for(it = infections.begin(); it!= infections.end(); it++)
	{
		inf_type = it->GetInfectionType();
		if(inf_type != 1) //if it's not an acute infection
			continue;
		else //if it is, return the acute virus
			randomAcuteInfections.push_back(it->pathogen);
	}
	
	int random = RandomNumber(0,randomAcuteInfections.size()-1);
	//randomAcuteInfections.at(random).PrintParametersVirus(); cout<<"getAcute()"<<endl;
	dummy.Copy(randomAcuteInfections.at(random));
	return dummy; // in case the list is empty, return an empty virus... shoudln't happen anyway! if the list is empty, this function should not be called!
}

/*This functions returns a randomly chose chronic virus infection  !!!!!!*/
Virus& Host :: GetChronicInfection(Virus &dummy)
{
	vector<Virus> randomChronicInfections;
	list<Infection>:: iterator it;
	int inf_type = 0;
	for(it = infections.begin(); it!= infections.end(); it++)
	{
		inf_type = it->GetInfectionType();
		//cout << inf_type << "infection type! "<<endl;
		if(inf_type!=2) //if it's not chronic
			continue;
		else
			randomChronicInfections.push_back(it->pathogen);
	}
	int random = RandomNumber(0,randomChronicInfections.size()-1);
	dummy.Copy(randomChronicInfections.at(random));
	return dummy;
}
/*This functions counts the number of chronic and acute infections (i.e ignores incubating individuals and obviously immune ones)*/
int Host :: CountInfections()
{
	int total_infections = 0;
	list <Infection>::iterator it;
	for(it = infections.begin(); it != infections.end(); it++)
	{

		if(it->IsIncubating() || it->IsImmune() || it->IsCured())
			continue;
		else
			total_infections++;
	}
	return total_infections;
}

void Host :: UpdateParameters(double timeStep, double simulationTime)
{
	age += (timeStep/YEAR);
	
	//update (get) all the information for each infection
	double vl = 0.0;
	double extra_vl = 0.0;
	list<Infection>::iterator it = infections.begin();
	while(it != infections.end())
	{
		it->SetInfectionType(simulationTime);
		if(it->IsCured())
		{//if the infection has been cleared again,remove it from the infections list!
			it = infections.erase(it);
			continue;
		}
		if(it->GetViralLoad() > vl)
			vl = it->GetViralLoad();
		if(it->GetDeadFlag()==true)
			SetDead();
		
		it++;
		//extra_vl += 0.0005;
	}
	viralDeathRate = vl+extra_vl;
	numberOfInfections = CountInfections();
}

//functions for SAVING PARAMETERS, BACKUP, etc
/*this function saves each locus -> to keep track of KIR and MHC diversity*/
void Host::SaveGenes(fstream& outfile)
{
	vector<KIRGene>::iterator kirIt;
	vector<Gene>::iterator mhcIt;
	for(mhcIt = mhcGenes.begin(); mhcIt!=mhcGenes.end(); mhcIt++)
	{
		mhcIt->SaveGenes(outfile);
		//outfile << mhcGenes.at(i).GetGeneID() << "\t";
	}
	for(kirIt = kirGenes.begin(); kirIt!=kirGenes.end(); kirIt++)
	{
		kirIt->SaveGenes(outfile);
		//outfile << kirGenes.at(i).GetGeneID() << "\t";
	}
	int functional_kirs = inhibitoryKIRs + activatingKIRs;
	outfile << functional_kirs << "\t"<< CountExpressedKIRs() << "\t" << Get_Host_ID()<<"\n";
}

void Host ::SaveParameters(fstream& outfile)
{
	outfile <<age << "\t"<<numberOfInfections<< "\t";
	list<Infection>::iterator infIt;
	for(infIt = infections.begin(); infIt!=infections.end(); infIt++)
	{
		 infIt->SaveParametersInfection(outfile);
	}
	outfile << Get_Host_ID() <<  "\n";
}

void Host :: PrintParametersHost()
{
	//cout <<age << "\t"<<infections.size()<< "\t"<< acute << chronic << immune << "\t";
	cout <<LOCI_KIR<<"\t"<<LOCI_MHC<<"\t"<<age << "\t"<< tuning << "\t"<< dead << "\t"<< mutationRateHost << "\t"<<inhibitoryKIRs<<"\t"<<activatingKIRs<<"\t"<< infections.size()<<"\t";
	list<Infection>::iterator infIt;
	vector<KIRGene>::iterator kirIt;
	vector<Gene>::iterator mhcIt;
	for(infIt = infections.begin(); infIt!=infections.end(); infIt++)
	{
		infIt->PrintParameters();
	}
	for(mhcIt = mhcGenes.begin(); mhcIt!=mhcGenes.end(); mhcIt++)
	{
		mhcIt->PrintParameters();
	}
	for(kirIt = kirGenes.begin(); kirIt!=kirGenes.end(); kirIt++)
	{
		kirIt->PrintParameters();
	}

	cout <<  "finished with this host ;-)"<<endl;
}

void Host :: SaveAgeDyingHost(fstream& outfile)
{
	outfile << age << "\t";
	vector<KIRGene>::iterator kirIt;
	for(kirIt = kirGenes.begin(); kirIt!=kirGenes.end(); kirIt++)
	{
		kirIt->SaveGenes(outfile);
	}
	//pathogen.SaveParametersVirus(outfile);
	//outfile << infectionType << "\t"<< ageInfection << "\t"<< ageClearance<<"\t";
}

void Host::SaveBackupHost(fstream& backupFile)
{
	//cout << "saving hosts..."<<endl;
	vector<KIRGene>::iterator kirIt;
	vector<Gene>::iterator mhcIt;
	list<Infection>::iterator infIt;

	backupFile << LOCI_KIR << "\t"<<LOCI_MHC<<"\t"<< age << "\t"<< tuning << "\t"<< dead << "\t"<< mutationRateHost << "\t"<<inhibitoryKIRs<<"\t"<<activatingKIRs<<"\t"<< infections.size()<<"\t";

	for(infIt = infections.begin(); infIt != infections.end(); infIt++)
	{
		infIt->SaveBackupInfection(backupFile);
	}

	for(mhcIt = mhcGenes.begin(); mhcIt!=mhcGenes.end(); mhcIt++)
	{
		mhcIt->SaveBackupGenes(backupFile);
	}
	for(kirIt = kirGenes.begin(); kirIt!=kirGenes.end(); kirIt++)
	{
		kirIt->SaveBackupGenes(backupFile);
	}
	backupFile << "\n";
}

void Host::RestoreHost(const string& sline)
{
	stringstream ssline(sline);

	ssline >> LOCI_KIR;
	ssline >> LOCI_MHC;
	ssline >>age;
	ssline >> tuning;
	ssline >> dead;
	ssline >> mutationRateHost;
	ssline >> inhibitoryKIRs;
	ssline >> activatingKIRs;
	int totalInfections;
	ssline >> totalInfections;

	/*use the streamstring from the position after "totalInfections"
	 *for that I copy the rest of the line into a string, and then reassign the streamstring with the content of
	 *for "restLine. For it to work, I have to clear the string first!
	 */
	string restLine;
	getline(ssline,restLine);
	ssline.clear();
	ssline.str(restLine);
	string infectionString;
	int type; double inf_time; double immune_time; double clearance_time;
	for(int i=0; i<totalInfections; i++)
	{
		Infection dummyInfection;
		infectionString = dummyInfection.RestoreInfection(ssline);
		infections.push_back(dummyInfection);
		ssline.clear();
		ssline.str(infectionString);
	}

	string geneString;
	//cout <<"mhc string>>>>>"<<ssline.str() <<endl;
	for(int i=0; i<LOCI_MHC*TWO; i++)
	{
		Gene mhc;
		geneString = mhc.RestoreGenes(ssline);
		mhcGenes.push_back(mhc);
		ssline.clear();
		ssline.str(geneString);
	}

	string kirString;
	for(int i=0; i<LOCI_KIR*TWO; i++)
	{
		//cout <<"kir string>>>>>"<<ssline.str() <<endl;
		KIRGene kir;
		kirString = kir.RestoreGenes(ssline);
		kirGenes.push_back(kir);
		ssline.clear();
		ssline.str(kirString);

	}
}


//FUNCTIONS OF CLASS INFECTION
Infection::Infection()
{
	SetInfectionParameters(cleared, dummy, 0.0, 0.0, 0.0);
}

Infection::Infection(Virus& _virus)
{
	pathogen.Copy(_virus);
	SetInfectionParameters(cleared, dummy, 0.0, 0.0, 0.0);

}

void Infection::SetProtectionLevel(const string protection)
{
	if(protection.compare("best") == 0)
		protectionLevel = best;
	if(protection.compare("high") == 0)
		protectionLevel = high;
	if(protection.compare("normal") == 0)
		protectionLevel = normal;
	if(protection.compare("medium") == 0)
		protectionLevel = medium;
	if(protection.compare("low") == 0)
		protectionLevel = low;
	if(protection.compare("zero") == 0)
		protectionLevel = zero;
	if(protection.compare("dummy") == 0)
		protectionLevel = dummy;
}

void Infection :: SetInfectionParameters(state _type, protection protection_type, double inf_time, double immune_time, double clearance_time)
{
	infectionType = _type;
	protectionLevel = protection_type;
	infectionTime = inf_time;
	immunityTime = immune_time;
	clearanceTime = clearance_time;
	deadFlagForHost = false;
	//cout << incubating  <<acute<< chronic<< memory<< cleared <<endl;
	//ageInfection = 0.0;
	//ageClearance = 0.0;
}

/*This function sets the pathogen parameters as infectious*/
void Infection::TransmitInfection(Virus& nastyInfection, double simulationTime)
{
	pathogen.Copy(nastyInfection); //infect it with the pathogen
	pathogen.SetViralLoad(pathogen.GetOriginalViralLoad());//get the original one for a new infection being transmission!
	infectionTime = simulationTime;
	if(pathogen.IsOnlyAcute())
		immunityTime = 80;
	else
		immunityTime = 10;
}

bool Infection :: IsPathogenNew(Virus& _newVirus)
{
	if(pathogen == _newVirus) //if the host is already infected with the same virus
		return false; //don't infect it
	else //but if it is a different one
		return true;
}

bool Infection::IsCured()
{
	if(infectionType == cleared)
		return true;
	else
		return false;
}

bool Infection :: IsAcute()
{
	if(infectionType == acute)
		return true;
	else
		return false;
}

bool Infection :: IsChronic()
{
	if(infectionType == chronic)
		return true;
	else
		return false;
}

bool Infection :: IsImmune()
{
	if(infectionType == memory)
		return true;
	else
		return false;
}

bool Infection :: IsIncubating()
{
	if(infectionType == incubating)
		return true;
	else
		return false;
}

bool Infection :: IsProtectionBest()
{
	if(protectionLevel == best)
		return true;
	else
		return false;
}

bool Infection :: IsProtectionHigh()
{
	if(protectionLevel == high)
		return true;
	else
		return false;
}
bool Infection :: IsProtectionNormal()
{
	if(protectionLevel == normal)
		return true;
	else
		return false;
}
bool Infection :: IsProtectionMedium()
{
	if(protectionLevel == medium)
		return true;
	else
		return false;
}
bool Infection :: IsProtectionLow()
{
	if(protectionLevel == low)
		return true;
	else
		return false;
}
bool Infection :: IsProtectionZero()
{
	if(protectionLevel == zero)
		return true;
	else
		return false;
}
void Infection::ResetInfection(double simulationTime)
{
	//set the viral load to zero, but keep the viral type to keep track of which infections is still present in the host!
	pathogen.SetViralLoad(0.0);
	infectionTime = 0.0;
	clearanceTime = simulationTime;
	deadFlagForHost = false;
}


/*This function changes the infection type of the host according to the time of infection*/
void Infection::SetInfectionType(double simulationTime)
{
	//infection has started but has not been cleared yet:
	if(infectionTime > 0.0 && clearanceTime == 0.0)
	{
		//check first if the host is in the infection period (i.e. between 1-4 weeks)
		if((simulationTime-infectionTime)<=1.0*WEEK)
			infectionType = incubating;

		if((simulationTime-infectionTime)>1.0*WEEK && (simulationTime-infectionTime)<1.0*WEEK+4.0*WEEK*pathogen.GetLifeTimeVirus())
		{
			//viralDeathRate = pathogen.GetViralLoad();STILL NEED TO IMPLEMENT THIS IN THE HOST!!!!!!!!
			infectionType = acute;
		}
		if((simulationTime -infectionTime) > 1.0*WEEK+ 4.0*WEEK*pathogen.GetLifeTimeVirus())
		{
			if(pathogen.IsOnlyAcute())
				deadFlagForHost = true;
			else
			{
				pathogen.SetViralLoad(0.6*pathogen.GetOriginalViralLoad());
				infectionType = chronic;
			}
		}
	}

	if(infectionTime == 0.0 && clearanceTime > 0.0)
	{

		if(simulationTime - clearanceTime < ImmunityTime(immunityTime,0.5)*YEAR)
		{
			infectionType = memory;
		}
		else
		{
			infectionType = cleared;
			clearanceTime = 0.0;
		}
	}

	if(infectionTime == 0.0 && clearanceTime == 0.0)
	{
		Virus deadVirus;
		pathogen.Copy(deadVirus);
		infectionType = cleared;
	}
}
void Infection::PrintParameters()
{
	//cout <<"infectionType infectionTime clearanceTime immunityTime and other viral parameters"<<endl;
	cout << infectionType<< "\t"<< infectionTime << "\t"<<clearanceTime << "\t"<<immunityTime << "\t";
	cout <<"protection:" <<protectionLevel <<"\t";
	pathogen.PrintParametersVirus();
}

void Infection ::SaveParametersInfection(fstream& outfile)
{
	outfile << infectionTime/YEAR <<"\t" << infectionType << "\t";
	pathogen.SaveParametersVirus(outfile);

}

void Infection::SaveBackupInfection(fstream& backupFile)
{
	backupFile<< infectionType<< "\t"<< infectionTime << "\t"<<clearanceTime << "\t"<<immunityTime << "\t";
	pathogen.SaveBackupVirus(backupFile);
}

// This function only works if the number of infections is one!
string Infection::RestoreInfection(stringstream& siline)
{
	int inf;
	siline >> inf;
	switch(inf)
	{
	case 4:infectionType = cleared;
	break;

	case 0:infectionType = incubating;
	break;

	case 1:infectionType = acute;
	break;

	case 2:infectionType = chronic;
	break;

	case 3:infectionType = memory;
	break;

	}
	siline >> infectionTime;
	siline >> clearanceTime;
	siline >> immunityTime;

	string restLine;
	getline(siline,restLine);
	siline.clear();
	siline.str(restLine);
	string virusString;
	//cout <<"passing to virus class:"<<siline.str()<<endl;
	virusString= pathogen.RestoreVirus(siline);
	//cout << "from the virus class:" <<virusString <<endl;
	//do I need this???
	//siline.clear();
	//siline.str(virusString);
	//string istring = siline.str();
	//return istring;

	return virusString;
}


/*This function performs a deep copy*/
Infection & Infection::Copy(Infection& rhsInfection)
{
	// checking if it is a self-assignment
	if(this == &rhsInfection)
		return *this;
	//copying member variables
	this->infectionType = rhsInfection.infectionType;
	this->infectionTime = rhsInfection.infectionTime;
	this->clearanceTime = rhsInfection.clearanceTime;
	this->immunityTime = rhsInfection.immunityTime;
	this->pathogen.Copy(rhsInfection.pathogen);
	return *this;
}


/*
 * to do:
 * write down in the file how many infections the host  has
 */
