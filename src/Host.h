/*
 * Host.h
 *
 *  Created on: Apr 21, 2011
 *      Author: paola
 */

#ifndef HOST_H_
#define HOST_H_
#include "Genes.h"
#define TWO 2.0
//#define LOCI_KIR 5
//#define LOCI_MHC 1
/*const double YEAR = 31207680.0;
const double MONTH = 2600640.0;
const double WEEK = 604800.0;*/
const double YEAR = 52.0;
const double MONTH = 4.0;
const double WEEK = 1.0;

class Virus{
public:
	Virus();
	virtual ~Virus(){};

	void SetViralParameters(double _downregulation, double _decoy, double _viralLoad, double _lifeTimeVirus, int virus_type, bool only_acute, bool specific_mhc, int _mhcID, bool decoy_type);
	enum type{wildType, downregulating};

	Virus& Copy(Virus& rhsVirus);//works
	bool operator ==(Virus& rhs);

	double GetLifeTimeVirus()const{return lifeTimeVirus;}
	double GetViralLoad()const{return viralLoad;}
	double GetOriginalViralLoad()const{return originalViralLoad;}
	void SetViralLoad(double vl){viralLoad=vl;};
	double GetRateDownregulation()const{return mutationRateDownregulation;}
	double GetRateDecoy()const{return mutationRateDecoy;}
	int GetVirusType()const{return virusType;}

	void BuildViralMolecule(int mhcID, int mutationType);
	void DownregulateMHC(int mhcID);
	bool IsDownregulatingMHC();
	bool IsExpressingMolecules();
	bool IsWildType();
	bool IsOnlyAcute()const{return onlyAcute;};
	bool IsDownregulatingMHCSpecifically(){return specificMHC;};

	void SaveBackupVirus(fstream& backupFile);//works
	string RestoreVirus(stringstream& svline);//works
	void PrintParametersVirus();
	void SaveParametersVirus(fstream& outfile);
	Gene viralMolecule;
	Gene downregulatedMHC;

protected:
	type virusType;
	double mutationRateDownregulation;
	double mutationRateDecoy;
	double viralLoad; //as the increase of the intrinsic death rate
	double originalViralLoad;
	double lifeTimeVirus; // time that a virus can live in one's organism
	bool onlyAcute;
	bool specificMHC;
	bool isDecoyAnMHC;

};

class Infection{
public:
	Infection();
	Infection(Virus& _virus);
	virtual ~Infection(){};

	enum state{incubating, acute, chronic, memory, cleared};
	enum protection{best, high, normal, medium, low, zero, dummy}; //dummy protection means that the host has not been ifnected yet.. so it doesn't really count!
	//corresponds to (0.95, 0.8, 0.7, 0.45, 0.25, 0, n.a)
	void SetInfectionParameters(state _type, protection _protectionwt, double inf_time, double immune_time, double clearance_time); //works
	void ResetInfection(double simulationTime); //ok
	void SetInfectionType(double simulationTime); //ok
	double GetInfectionTime(){return infectionTime;};
	double GetClearanceTime(){return clearanceTime;};
	int GetInfectionType(){return infectionType;};
	double GetViralLoad(){return pathogen.GetViralLoad();};
	bool GetDeadFlag(){return deadFlagForHost;}
	//Virus& GetVirus(){return pathogen;}
	Infection& Copy (Infection& rhsInfection);

	void SetProtectionLevel(const string protectionlevel);

	bool IsPathogenNew(Virus& _newVirus);
	bool IsCured();
	bool IsAcute();
	bool IsChronic();
	bool IsImmune();
	bool IsIncubating();

	bool IsProtectionBest();
	bool IsProtectionHigh();
	bool IsProtectionNormal();
	bool IsProtectionMedium();
	bool IsProtectionLow();
	bool IsProtectionZero();

	void TransmitInfection(Virus& nastyInfection, double simulationTime);

	void SaveParametersInfection(fstream& outfile); //works
	void SaveBackupInfection(fstream & backup);
	string RestoreInfection(stringstream& siline);
	void PrintParameters();
	Virus pathogen;

protected:

	bool deadFlagForHost;
    double infectionTime;
    double clearanceTime;
    int immunityTime;
    double ageInfection;
    double ageClearance;
    state infectionType;
    protection protectionLevel;
};


class Host {
public:
	Host(){}; //for creating memory
	Host(int loci_kir, int loci_mhc, double _mutationRate, bool _tuning,
			int numberOfExtraKirs, Map& kirMap, MHCGenePool& mhcPool, bool hla); //for initialization of the population //works
	
	Host(int loci_kir, int loci_mhc, vector<Gene>& mhcGenesParent, GenePool& mhcPool, bool dist,
			vector<KIRGene>& kirGenesMother, vector<KIRGene>& kirGenesFather, double _mutationRate,
			bool _tuning,int numberOfExtraKirs, Map& kirMap, int mutationType, int gene_type, double simulationTime, double time_invation,bool invasion_analysis);
	virtual ~Host(){};
	void SetDead(){dead = true;}
	bool IsDead()const{return dead;}

	void InitializeHostParameters(double mutationRate, bool _tuning, int kirloci,int loci_mhc);
	void SetHostParameters(bool t, double mut, int inftyp, double inftim, double viraldeathm ,double clrtim);
	void EducateKIRs(); //works with activating receptors also!
	void MutateGenes(int mutationType, KIRGene& kir_hap2, Map& kirMap, GenePool& mhcPool, int gene_type); //works
	void MutateGenesForMutualInvasion(int mutationType, KIRGene& kir_hap2, Map& kirMap, GenePool& mhcPool, double simulationTime, double time_invasion, int gene_type); //works

	void ExpressKIRs(int numberOfExpressedKirs); //ok
	double GetAge()const{return age;}
	void SetAge(double number){age = number;}
	void CountFunctionalInhibitoryKIRs();
	void CountFunctionalActivatingKIRs();
	int CountExpressedKIRs(); //ok
	Host& Copy(Host& rightHandSideHost);//works

	void ResetKIRs();

	//double GetIntrinsicDeathRate(const vector<double>& rates);
	//double GetAgeDependentBirthRate(vector<double>& rates);

	void InfectWith(Virus& nastyVirus, double simulationTime, int maxNumberInfections); //works
	void InfectWithMoreDownregulating(Virus& nastyVirus, double simulationTime, int maxNumberInfections);
	void InfectWithMoreWildType(Virus& nastyVirus, double simulationTime, int maxNumberInfections);


	int IsInfected();
	void ClearInfection(double simulationTime, Infection& _infection);
	void ClearDecoyWithInhibitoryOnly(int inhibiting_signal, double simulationTime, Infection& _infection);
	void ClearDecoyWithActivatingOnly(int activating_signal, double simulationTime, Infection& _infection);
	void ClearDecoyWithActivatingAndInhibitory(int inhibiting_signal, int activating_signal, double simulationTime, Infection& _infection);
	
	Virus& GetAcuteInfection(Virus& virus); //works
	Virus& GetChronicInfection(Virus& virus); //works
	int GetMainInfectionType(); //works
	int CountInfections();//ok
	
	double GetMutationRate(){return mutationRateHost;};
	double GetViralDeathRate(){return viralDeathRate;};


	void UpdateParameters(double timeStep, double simulationTime);//works
	void SaveGenes(fstream& outfile);//works
	void SaveParameters(fstream& outfile);//works
	void SaveAgeDyingHost(fstream& outfile);

	void SaveBackupHost(fstream& file);//works
	void RestoreHost(const string& sline);//works
	
	void PrintParametersHost();

	vector<Gene> mhcGenes;
	vector<KIRGene> kirGenes;
	list<Infection> infections;
	enum state{susceptible, acute, chronic, immune};

	void Set_Host_ID(const unsigned long int _id){host_ID = _id;}
	unsigned long int Get_Host_ID(){return host_ID;}
protected:

	int lociKIR;
	int lociMHC;
	double age;
	bool dead;
	bool tuning;
	state mainInfectionType;
	double mutationRateHost;

	//double intrinsicDeathRate;
	//double ageDependentBirthrate;
	double viralDeathRate;
    int inhibitoryKIRs; //number of FUNTIONAL inhibitory / activating receptors
    int activatingKIRs;
	unsigned long int host_ID;
	int numberOfInfections;
};

#endif /* HOST_H_ */

/*
 * 06.06.2014: crete fix hot fix
 *
 * 23.08.2012: removed int specificity from all function arguments -> all specificity values were taken from genes.GetSpecificity() ....
 *
 * 06.09.2012: changed Host constructors. New argument: kirMap. Changed the mutation-> now it picks a new gene
 *
 * 18.10.2012: changed Host constructors: introduced a new type of mutation-> to PointMutation() and MutateSpecificity()
 */
