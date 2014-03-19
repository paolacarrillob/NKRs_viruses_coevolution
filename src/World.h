/*
 * World.h
 *
 *  Created on: Apr 26, 2011
 *      Author: paola
 */

#ifndef WORLD_H_
#define WORLD_H_
#include "Host.h"//Host.h includes Virus.h and Gene.h which includes Bitstrings.h and Mathfunctions.h

class World {
public:
	World(); //Default constructor
	virtual ~World(){}; //destructor that so far does nothing at all!

	void WriteInfo();
	void LoadParameterFile(const string& fileName);//works
	bool Initialize(); //works
	void CreateBirthAndDeathRates();

	bool Birth(int index, unsigned long int next_id);//, Host& host); 	// Birth() function takes a host at index (index) and picks a random host at a random index and create a "child" that is copied to (host).
										// It returns true or false according if the child should be added to the population or not
	bool Death (int index);
	void Infect(int index);
	void EscapeOnlyDecoy(int index);


	void ShuffleHosts();
	void IntroduceVirus(Virus& secondVirus, double timeToIntroduceTheVirus, int _mutationTypeVirus);
	void DetermineSecondVirus(); //works
	void AddMoreViruses();
	void RemoveDeadHosts_HappyNewYear();
	void TrackInfectedIndividuals();
	
	double GetAgeDependentBirthRate(const double age)const;
	double GetIntrinsicDeathRate(const double age, const double viralrate)const;
	
	void CountProtectionType(Infection& _infection);

	void Simulate();

	void SaveGenes();//works
	void SavePopulationSize();//works
	void SaveParameters();//works
	void SaveMap();
	bool SaveAgeDyingHosts(double lastOutfileTime, double lastStopOutfileTime);

	void SaveBackupFile();//works
	void LoadBackupFile(const string& backupName);//works


	vector<Host> hosts; //vector of hosts
	MHCGenePool MHCPool; //vector of integers representing MHC pools
	//GenePool KIRPool;
	Map KIRGenesMap; //map of gene_type | KIR_id|L|M_id
	Virus nastyVirus;
	Virus downregulatingVirus;
	//vector<int> viralMolecules;
	//Virus decoyVirus;
	

	double mutationRate;
	bool education;
	int expressionExtraKIRs;
	bool HLA_C;
	int KIRspecificity;
	int KIRLoci;
	bool onlyOneInitialKIR;
	int KIRGeneType;
	int InitialKIRGeneType ;
	int MHCLoci;
	int sizeMHCPool;
	bool similarMHCsInThePool;
	unsigned int initHostPop;
	int mutationType;
	int mutationTypeHost;
	int mutationTypeVirus;

	bool decoyStealsMHC; //parameter to change the decoy type: either it steals one random MHC or it is a new one!
	bool specificMhcDownregulation;

	bool noExtraViruses; //parameter meaning that there can only be uniwue infections within an individual
	bool extraWTviruses; //here one host can be infected with several WT viruses
	bool extraMHCdownregulatingViruses; //here one host can be infected with several MHC down viruses

protected:
	double birthRate;
	double deathRate;
	double infectionRate;
	double escapeRate;

	double timeStep;
	double timeEnd;
	double simulationTime;
	double outfileRate;
	double backupRate;
	double populationSizeRate;
	int maxHostPop;

	int contactRate;
	double timeIntroducingInfection;
	double timeSecondVirus;
	string secondVirusName;

	Virus secondVirus;
	Virus firstVirus;

	double deltaVirus;
	double timeInfection;
	double downregulationRate;
	double decoyRate;

	double transmissionRateAcute;
	double transmissionRateChronic;

	vector<int> shuffledHosts;
	vector<int> virtualHosts;

    vector<double> deathRates;
    vector<double> birthRates;

	//int acute_infected;
	//int chronic_infected;
	//int immune;

	int wt_susceptible_mhc_susceptible;
	int wt_susceptible_mhc_chronic;
	int wt_susceptible_mhc_immune;

	int wt_chronic_mhc_susceptible;
	int wt_chronic_mhc_chronic;
	int wt_chronic_mhc_immune;

    int wt_immune_mhc_susceptible;
    int wt_immune_mhc_chronic;
    int wt_immune_mhc_immune;

    int best_protection; //p = 0.95
    int high_protection; // p = 0.8
    int normal_protection; // p = 0.7
    int medium_protection; // p = 0.45
    int low_protection; // p = 0.25
    int zero_protection; //p=0

	int simpleInfection;
	int doubleInfection;
	int tripleInfection;

	int maxNumberOfInfectionsPerHost;
	
	bool isFileOpen;
	ofstream populationSize;
	fstream genesFile;
	fstream parameterFile;
	fstream backupFile;
	fstream notBornChildren;
	fstream dyingHosts;
	fstream mapFile;

	bool invasionAnalysis;
	double timeInvasion;
};

#endif /* WORLD_H_ */

/*
 * new member variable: sizeMHCpool -> it's initialized in LoadParameterFile()
 */

