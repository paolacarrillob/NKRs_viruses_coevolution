/*
 * World.cpp
 *
 *  Created on: Apr 26, 2011
 *      Author: paola
 */

#include "World.h"
#include "backup.h"

World::World() {
}

void World :: LoadParameterFile(const string& fileName)
{
	double t_outfile;
	double t_backup;
	double t_popFile;
	bool onlyAcuteInfection;
	char buff[256];
	
	kaBackup paramFile(fileName,false);
	paramFile.Load("Maximal population size", &maxHostPop);
	paramFile.Load("Initial population size", &initHostPop);
	paramFile.Load("Time Step", &timeStep);
	paramFile.Load("End of Simulation", &timeEnd);
	paramFile.Load("Time introducing infection", &timeIntroducingInfection);
	paramFile.Load("Time second infection", &timeSecondVirus);
	//paramFile.Load("Time third infection", &timeThirdVirus);
	paramFile.Load("Time outfile", &t_outfile);
	//paramFile.Load("Time recording", &timeRecording);
	paramFile.Load("Time backup", &t_backup);
	paramFile.Load("Population Time", &t_popFile);
	paramFile.Load("Mutation rate host", &mutationRate);
	paramFile.Load("Mutation type host", &mutationTypeHost);
	paramFile.Load("Mutation type virus", &mutationTypeVirus);
	paramFile.Load("Contacts per week", &contactRate);
	paramFile.Load("KIR Loci", &KIRLoci);
	paramFile.Load("only one initial KIR",&onlyOneInitialKIR);
	paramFile.Load("KIR type", &KIRGeneType);
	paramFile.Load("Initial KIR type", &InitialKIRGeneType);
	paramFile.Load("MHC Loci", &MHCLoci);
	paramFile.Load("HLA-C alleles distribution", &HLA_C);
	paramFile.Load("Number HLA-C alleles",&sizeMHCPool);
	paramFile.Load("Similar MHCs",&similarMHCsInThePool);
	paramFile.Load("MHC-KIR specificity", &KIRspecificity);
	paramFile.Load("Tuning", &education);
	paramFile.Load("Extra Number KIRs to be expressed", &expressionExtraKIRs);
	paramFile.Load("Viral load", &deltaVirus);
	paramFile.Load("Life time of the virus", &timeInfection);
	paramFile.Load("MHC dowregulation rate", &downregulationRate);
	paramFile.Load("Decoy rate", &decoyRate);
	paramFile.Load("Transmission rate acute infection", &transmissionRateAcute);
	paramFile.Load("Transmission rate chronic infection", &transmissionRateChronic);
	paramFile.Load("Type of mutation", &mutationType);
	paramFile.Load("Decoy Steals MHC", &decoyStealsMHC);
	paramFile.Load("Type of infection", &onlyAcuteInfection);
	paramFile.Load("Maximal Number of Infections", &maxNumberOfInfectionsPerHost);
	paramFile.LoadString("Second Virus", buff, 256);
	secondVirusName.assign(buff);
	paramFile.Load("MHC downregulation type", &specificMhcDownregulation);
	//paramFile.Load("Decoy type", &random_decoy);

	paramFile.Load("Invasion analysis", &invasionAnalysis);
	paramFile.Load("Time for invasion", &timeInvasion);

	noExtraViruses = true; //parameter meaning that there can only be uniwue infections within an individual
	extraWTviruses = false; //here one host can be infected with several WT viruses
	extraMHCdownregulatingViruses = false;

	timeStep= timeStep*WEEK;
	simulationTime = 0.0;
	timeEnd = timeEnd*YEAR;
	timeInvasion = timeInvasion*YEAR;

	timeIntroducingInfection = timeIntroducingInfection*YEAR;
	timeSecondVirus = timeSecondVirus *YEAR;

	outfileRate = 1.0/(t_outfile*YEAR);
	backupRate = 1.0/(t_backup*YEAR);
	populationSizeRate = 1.0/(t_popFile*YEAR);

	birthRate = 0.5*timeStep/YEAR; // every birth event will happen once every four years
	deathRate = timeStep/YEAR; // every death event will happen once every year
	infectionRate = 52.0*timeStep/YEAR; //every infection event will happen every week (pop concerts removed!)
	escapeRate = 0.0001*timeStep/YEAR;

	nastyVirus.SetViralParameters(downregulationRate, decoyRate,deltaVirus,timeInfection, 0, onlyAcuteInfection, specificMhcDownregulation, 0, decoyStealsMHC);
	downregulatingVirus.SetViralParameters(downregulationRate, decoyRate,deltaVirus,timeInfection, 1, onlyAcuteInfection,specificMhcDownregulation,0, decoyStealsMHC);

	wt_susceptible_mhc_susceptible = 0;
	wt_susceptible_mhc_chronic = 0;
	wt_susceptible_mhc_immune = 0;

	wt_chronic_mhc_susceptible = 0;
	wt_chronic_mhc_chronic = 0;
	wt_chronic_mhc_immune = 0;

    wt_immune_mhc_susceptible = 0;
    wt_immune_mhc_chronic = 0;
    wt_immune_mhc_immune = 0;
	/*downregulating = 0;
	downregulating_immune = 0;
	decoy = 0;
	decoy_immune = 0;
	wildtype = 0;
	wildtype_immune = 0;

	simpleInfection = 0;
	doubleInfection = 0;
	tripleInfection = 0;*/

	WriteInfo();
}

void World :: WriteInfo()
{
	cout << "Beginning simulation with: \n";
	cout << "End Simulation after " << timeEnd/YEAR << " years \n";
	cout << "Maximal population size " << maxHostPop << " individuals \n";
	cout << "Initial population size " << initHostPop << " individuals \n";
	cout << "Introducing virus after " << timeIntroducingInfection/YEAR << " years \n";
	cout << "Mutation rate host " << mutationRate << "\n";
	cout << "Mutation type host " << mutationTypeHost << "\n";
	cout << "KIR loci " << KIRLoci << "\n";
	cout << "only one initial KIR" << onlyOneInitialKIR << "\n";
	cout << "KIR Gene Type" << KIRGeneType << "\n";
	cout << "Initial KIR Gene Type" << InitialKIRGeneType << "\n";
	cout << "MHC loci " << MHCLoci << "\n";
	cout << "Contacts per week: " << contactRate << "\n";
	cout << "HLA-C alleles distribution " << HLA_C<<"\n";
	cout << "Number of HLA alleles " << sizeMHCPool<<"\n";
	cout << "MHC-KIR specificity " << KIRspecificity<<"\n"; //this is only for the initial pool, basically, the maximum specificity is set with this
	cout << "Tuning "<< education << "\n";
	cout << "Expression extra KIRs "<< expressionExtraKIRs << "\n";
	cout << "Viral load: " << deltaVirus << "\n";
	cout << "Life time virus: " << timeInfection << "\n";
	cout << "MHC dowregulation rate " << downregulationRate <<"\n";
	cout << "Decoy rate "<< decoyRate << "\n";
	cout << "Tranmission rate acute infection "<< transmissionRateAcute << "\n";
	cout << "Type of mutation " << mutationType << "\n";
	cout << "Decoy Steals MHC " << decoyStealsMHC << "\n";
	cout << "Maximal Number of Infections " << maxNumberOfInfectionsPerHost << "\n";
	cout << "MHC downregulation type" << specificMhcDownregulation <<"\n";
	
}

// initialize host population
bool World::Initialize()
{
	//initialize MHcPool
	if(similarMHCsInThePool)
		MHCPool.FillMHCGenePoolWithSimilarMHCs(sizeMHCPool);
	else
		MHCPool.FillMHCGenePool(sizeMHCPool);

	MHCPool.WriteOutGenes();

	KIRGene dummy_gene(KIRspecificity);
	if(onlyOneInitialKIR) //if I want to initialize only with the same KIR...
	{
		KIRGene dummy_gene(KIRspecificity); //create only one functional KIR
		dummy_gene.SetGeneType(InitialKIRGeneType);
		KIRGenesMap.FillMap(MHCPool, dummy_gene);
		//create another functional KIR (with the opposite type than the first one, but with the same specificity

		KIRGene next_dummy_gene(KIRspecificity);
		int otherGeneType = 1-InitialKIRGeneType;
		next_dummy_gene.SetGeneType(otherGeneType);
		KIRGenesMap.FillMap(MHCPool, next_dummy_gene);

		for(int i = 2; i<KIRLoci; i++) //and fill the rest with pseudo genes!
		{
			KIRGene dummy_gene(15);
			dummy_gene.SetGeneType(InitialKIRGeneType);
			KIRGenesMap.FillMap(MHCPool, dummy_gene);
		}
	}
	else //otherwise...
	{
		//initialize the map of KIR with the first KIRloci genes to initialize the population
		for (int i =  0; i <KIRLoci; i++)
		{
			KIRGene dummy_gene(KIRspecificity);
			dummy_gene.SetGeneType(InitialKIRGeneType);
			if (!KIRGenesMap.IsGeneInMap(dummy_gene))
				KIRGenesMap.FillMap(MHCPool, dummy_gene);
		}
	}

	//initialize the population with the MHC genes of the pools
	for (unsigned int i = 0; i< initHostPop; i++)
	{
		Host dummyhost(KIRLoci, MHCLoci, mutationRate, education, expressionExtraKIRs, KIRGenesMap, MHCPool, HLA_C);
		dummyhost.Set_Host_ID(i);
		hosts.push_back(dummyhost);
		//dummyhost.PrintParametersHost();
	}

	CreateBirthAndDeathRates();
	return true;
}

/* This function creates a table upon initialization that calculates the birth/death rates for all possible ages*/
void World ::CreateBirthAndDeathRates()
{
	try
	{
		for(unsigned long int j=0; j <150*52; j++)
		{
			double i = j/52.0;
			double birth = -1/(1+ exp(i-20)) + 1/(1+exp(i-45));
			birthRates.push_back(birth);

			double death = exp(0.1*i-10.5)+ exp(-0.4*i-8);
			deathRates.push_back(death);
			//cout << i << " "<< birth <<" "<< death << endl;
		}
		cout <<birthRates.size() << "|" << deathRates.size()<<endl;
	}
	catch (...) {
		cout << "caught something" << endl;
	}
}

/*EVENT functions*/
/*Birth function: creates a child with the haplotype of one parent and another randomly chosen host*/
bool World::Birth(int index,unsigned long int next_id)//, Host& baby_host)
{
	if(hosts.size() <=1)
		return false;
	//check if host's age allows him to become a parent
	//cout << "do i get stuck here??? birth event!"<< endl;
	double ageDependentBirth = 0.0;
	try
	{
		//ageDependentBirth = hosts.at(index).GetAgeDependentBirthRate(birthRates);
		ageDependentBirth = GetAgeDependentBirthRate(hosts.at(index).GetAge());
	}
	catch (...)
	{
		cout << "caught something!!!"<< endl;
		cout << "index:" <<index << endl;
		cout << "hosts:" <<hosts.size() << endl;		
	}
	//cout << "do i get stuck here??? birth event!"<< endl;
	if(RandomNumberDouble()< birthRate*ageDependentBirth*(1-(hosts.size()/(maxHostPop*0.99753))))
	{
		int randomindex = RandomNumber(0,shuffledHosts.size()-1);
		//cout << "do i get stuck here??? birth event!"<< endl;
		//check if the potential parent is himself
		while(randomindex == index)
		{
			randomindex = RandomNumber(0,shuffledHosts.size()-1);
		//	cout << "do i get stuck here??? birth event!"<< endl;
		}

		//choose randomly which parent is going to "donate" his/her mhc molecule
		int parent;
		if((RandomNumberDouble()<0.5))
			parent = index;
		else
			parent = randomindex;
		//cout << "do i get stuck here??? birth event!"<< endl;
		Host testHost(KIRLoci,MHCLoci, hosts.at(parent).mhcGenes ,MHCPool, HLA_C,/*KIRPool, */hosts.at(index).kirGenes,hosts.at(randomindex).kirGenes, mutationRate,education,expressionExtraKIRs, KIRGenesMap, mutationTypeHost, KIRGeneType, simulationTime, timeInvasion, invasionAnalysis);
		//cout << "do i get stuck here??? birth event!"<< endl;
		testHost.Set_Host_ID(next_id);
		hosts.push_back(testHost);
		//baby_host.Copy(testHost);
		//cout << "do i get stuck here??? birth event!"<< endl;
		return true;
	}
	return false;
}

/*Death function: according to an age-dependent rate, a host will be removed from the population*/
bool World::Death(int index)
{
	//cout << "do i get stuck here??? death event!"<< endl;
	double intrinsicDeath = 0;
	vector<Host>::iterator it = hosts.begin() + index;
	//double intrinsicDeathRate = it->GetIntrinsicDeathRate(deathRates);
	double intrinsicDeathRate = GetIntrinsicDeathRate(it->GetAge(),it->GetViralDeathRate());
	intrinsicDeath = deathRate*intrinsicDeathRate;
	double r = RandomNumberDouble();
	if(r<intrinsicDeath)
	{
		//cout << "r , intrinsicdeathrate, intrinsicdeath:" << r << " , " << intrinsicDeathRate << " , " << intrinsicDeath << (r<intrinsicDeath?" ====> will DIE":"") << endl; 		
		it->SetDead();
		return true;
	}
	return false;
}
/*Infect function: upon contact between two hosts, virus can spread*/

void World::Infect(int index)
{
	//cout << "do i get stuck here??? infect event!"<< endl;
	if(hosts.size() <=1)
		return;
	if(RandomNumberDouble()<infectionRate)
	{
		for(int i=0; i<contactRate; i++)
		{
			//pick random partner that is NOT yourself!
			int randomindex = RandomNumber(0,shuffledHosts.size()-1);
			while (randomindex == index)
				randomindex=RandomNumber(0,shuffledHosts.size()-1);
			//check whether partner is infectious

			//if(hosts.at(randomindex).IsSusceptible() || hosts.at(randomindex).IsImmune())
			if(!hosts.at(randomindex).IsInfected())
				continue;
			//transmit the virus according to the type of infection
			else
			{
				int infectionState = hosts.at(randomindex).GetMainInfectionType();
				Virus stupidVirus;
				switch (infectionState)
				{
					case 0: break;
					case 1:
						if(RandomNumberDouble()<transmissionRateAcute)
						{
							//get the acute virus
							Virus dummy;
							//stupidVirus.Copy(hosts.at(randomindex).GetAcuteInfection(dummy));
							hosts.at(randomindex).GetAcuteInfection(dummy);
							stupidVirus.Copy(dummy);
							//stupidVirus.PrintParametersVirus(); cout <<"Infect acute"<<endl;
							hosts.at(index).InfectWith(stupidVirus, simulationTime, maxNumberOfInfectionsPerHost);
						}
						break;
					case 2:
						if(RandomNumberDouble()<transmissionRateChronic)
						{
							// get the chronic virus
							Virus dummy;
							//stupidVirus.Copy(hosts.at(randomindex).GetChronicInfection(dummy));
							hosts.at(randomindex).GetChronicInfection(dummy);
							stupidVirus.Copy(dummy);
							hosts.at(index).InfectWith(stupidVirus, simulationTime, maxNumberOfInfectionsPerHost);
						}
						break;
					case 3: break;
					case 4: break;
					default: cout <<"ERROR!!!!!! shouldn't happen!" <<endl; exit(1);
				}
			}
		}
	}
}

void World::EscapeOnlyDecoy(int index)
{
	int hap = RandomNumberSmallInt(0,MHCLoci*2 - 1);
	int mhcID = hosts.at(index).mhcGenes.at(hap).GetGeneID(); //we pass the MHC of the host to the virus here, because it
															//	will only be set as the "id" of the viral molecules if the virus has that property
															// the function BuildViralMolecule() takes care of that
	list<Infection>::iterator it;
	for(it = hosts.at(index).infections.begin(); it!= hosts.at(index).infections.end(); it++)
	{
			it->pathogen.BuildViralMolecule(mhcID, mutationTypeVirus); //then mutate it!
	}
}

/*SIMULATION functions*/

void World::ShuffleHosts()
{
	//virtualHosts.clear();
	shuffledHosts.clear();

	for(unsigned int i=0; i<hosts.size(); i++)
	{
		shuffledHosts.push_back(i);
	}

   for(unsigned int i = 0 ; i < hosts.size() ; i ++ )
   {
      int j = RandomNumber(i,hosts.size()-1);
      // swap a[i] and a[j]
      int t = shuffledHosts.at(j);
      shuffledHosts.at(j) = shuffledHosts.at(i);
      shuffledHosts.at(i) = t;
   }
}

void World::Simulate()
{
	const string populationFile("PopulationSize.log");
	populationSize.open(populationFile.c_str());
	//populationSize <<"#time\tpopSize\tsimpleInfections\tdoubleInfections\ttripleInfections\twt\twt_i\tmhc_down\tmhc_down_i\tdecoy\tdecoy_i\n";
	populationSize <<"#time\tpopSize\twt_s_mhc_s\twt_s_mhc_c\twt_s_mhc_i\twt_c_mhc_s\twt_c_mhc_c\twt_c_mhc_i\twt_i_mhc_s\twt_i_mhc_c\twt_i_mhc_i\tbest\thigh\tnormal\tmedium\tlow\tzero\n";
	double lastPopulationOutfileTime = 0.0;
	double lastOutfileTime = 0.0;
	double lastBackupTime = 0.0;
	double lastAcuteInfectionTime = 0.0;
	double lastCountingProtectionLevels = 0.0;
	double rateSavingProtectionLevels = 1.0/(100*YEAR);
	isFileOpen = false;
	SaveMap();
	cout << "simulation Time: "<<simulationTime <<"\n"<<endl;
	unsigned long int id_counter = initHostPop;

	DetermineSecondVirus(); //assign the second and third types of the virus!
	while(simulationTime <=timeEnd)
	{
		// printing out the backup files
		if(floor((simulationTime-lastBackupTime)*backupRate)>0)
		{
			cout <<"\tSaving Backup\n"<<endl;
			SaveBackupFile();
			//SaveMap(); i don't really need to save the map always... it will save tons of files i never use anyway!
			lastBackupTime = simulationTime;
		}

		int number_babies = 0;
		int number_dead_people = 0;

		ShuffleHosts();
		vector<int>::iterator shuffledHostsit;

		//introduce infections
		IntroduceVirus(firstVirus, timeIntroducingInfection, mutationTypeVirus);
		AddMoreViruses();

		//reset the number of protection types after one year of gathering information... otherwise i don't have anything!
		if(floor((simulationTime-lastCountingProtectionLevels)*rateSavingProtectionLevels)>0)
		{
			best_protection = 0; //p = 0.95
			high_protection = 0; // p = 0.8
			normal_protection = 0; // p = 0.7
			medium_protection = 0; // p = 0.45
			low_protection = 0; // p = 0.25
			zero_protection = 0; //p=0
			lastCountingProtectionLevels = simulationTime;
		}

		for(shuffledHostsit = shuffledHosts.begin(); shuffledHostsit!=shuffledHosts.end(); shuffledHostsit++)
		{
			//let event happen for every random host
			int index = *shuffledHostsit;
			Host babyHost;
			if(Birth(index,id_counter))
			{
				number_babies++;
				id_counter++;
			}
			Infect(index);

			EscapeOnlyDecoy(index); //allow for the decoy(viral molecules) to mutate within the hosts!

			if(Death(index))
			{
				number_dead_people++;
			}
			
			//clear the infection
			list<Infection>::iterator inf;
			for(inf = hosts.at(index).infections.begin(); inf!=hosts.at(index).infections.end(); inf++) //check for every infection within one host, whether it's time to be cleared
			{
				if(inf->IsAcute())
				{
					if((simulationTime - inf->GetInfectionTime()) == (1.0 + 4.0 *timeInfection)*WEEK)
					{
						hosts.at(index).ClearInfection(simulationTime,(*inf));
						CountProtectionType((*inf));
						inf->SetProtectionLevel("dummy");//reset the protection level
					}
				}
			}
			hosts.at(index).UpdateParameters(timeStep,simulationTime);
		}
		RemoveDeadHosts_HappyNewYear();

		//printing out the population size not every week but with a rate (otherwise the file will turn huge!!!!)
		bool timeToPrintPopulationSize = floor((simulationTime-lastPopulationOutfileTime)*populationSizeRate)>0;
		if(simulationTime == 0.0 || timeToPrintPopulationSize)
		{
			//cout << "what's going on?????" <<endl;
			SavePopulationSize();
			lastPopulationOutfileTime = simulationTime;
		}
		if(hosts.size() == 0)
		{
			cout << "Oooops pao is angry and by the way host size is:" << hosts.size() << endl;
			break;
		}
		TrackInfectedIndividuals();

		// printing out the gene files
		bool timeToPrintOut = floor((simulationTime-lastOutfileTime)*outfileRate)>0;
		if(simulationTime == 0.0 || timeToPrintOut)
		{
			cout <<"\tPrinting parameters\n";
			cout<< "Time: " << simulationTime/YEAR <<endl;
			SaveGenes();
			SaveParameters();
			populationSize.flush();
			if(populationSize.bad())
				cout << "something bad happened \n";
			lastOutfileTime = simulationTime;
		}

		//re introducing acute infection
		if(nastyVirus.IsOnlyAcute())
		{
			if(floor((simulationTime -lastAcuteInfectionTime)*(1.0/(1.0*YEAR)))>0)
			{
				cout << "reintroducing the infection"<<endl;
				IntroduceVirus(nastyVirus,timeIntroducingInfection, mutationTypeVirus);
				lastAcuteInfectionTime = simulationTime;
			}
		}

		simulationTime+=timeStep;
	}
	populationSize.close();
	SaveMap();
	WriteInfo();
}
void World :: CountProtectionType(Infection& _infection)
{
	if(_infection.IsProtectionBest())
		best_protection++; //p = 0.95
	if(_infection.IsProtectionHigh())
		high_protection++; // p = 0.8
    if(_infection.IsProtectionNormal())
    	normal_protection++; // p = 0.7
    if(_infection.IsProtectionMedium())
    	medium_protection++; // p = 0.45
    if(_infection.IsProtectionLow())
    	low_protection++; // p = 0.25
    if(_infection.IsProtectionZero())
    	zero_protection++; //p=0
    //cout << simulationTime<<": "<<best_protection << "|" << high_protection << "|"<<normal_protection << "|"<<medium_protection << "|"<<low_protection << "|"<<zero_protection <<"\n";
}

void World :: IntroduceVirus(Virus& secondVirus, double timeToIntroduceTheVirus, int _mutationTypeVirus)
{
	if(floor(timeToIntroduceTheVirus -simulationTime)>0 && floor(timeToIntroduceTheVirus -simulationTime)<3.0*WEEK)
	{
		cout <<"\t Introducing the infection"<<endl;
		for(int i= 0; i<0.05* hosts.size(); i++)
		{
			int dummy_mhc = RandomNumber(0,65000); //just pick any random number in case we need to mutate ... not the best design.. i know
			/*
			 * the secondVirus has only one random viral molecules. When we choose now 5% of the hosts
			 * to be infected, the second virus will mutate with a small probability... hence, it shouldn't be a big change..
			 * we are basically starting with one viral molecule
			 */
			secondVirus.BuildViralMolecule(dummy_mhc,_mutationTypeVirus);

			int randomindex = RandomNumber(0,hosts.size()-1);
			hosts.at(randomindex).InfectWith(secondVirus, simulationTime, maxNumberOfInfectionsPerHost);
			/*the InfectWith function will take care that there are not super infections of the same viral type within one host!*/
		}
	}
}

void World ::DetermineSecondVirus()
{
	string wildTypeVirusString("wild type");
	string mhcDownVirusString("mhc down");

	if(secondVirusName.compare(mhcDownVirusString) == 0)
	{
		secondVirus.Copy(downregulatingVirus);
		firstVirus.Copy(nastyVirus);

	}

	if(secondVirusName.compare(wildTypeVirusString) == 0)
	{
		secondVirus.Copy(nastyVirus);
		firstVirus.Copy(downregulatingVirus);
	}
	cout << "first virus" <<endl;
	firstVirus.PrintParametersVirus();
	cout <<endl;

	cout << "second virus" <<endl;
	secondVirus.PrintParametersVirus();
	cout <<endl;

}

/*This function inoculates the population with additional viruses, depending on how
 * many infections per Host are allowed!*/
void World::AddMoreViruses()
{
	switch(maxNumberOfInfectionsPerHost)
	{
		case 1:break; //if only one, then return immediately
		case 2: //if two, then add the second virus!

				IntroduceVirus(secondVirus,timeSecondVirus,mutationTypeVirus);

		break;
		//default: cout <<"ERROR!!!!!! an impossible number if infections! shouldn't happen! AddMoreViruses()" <<endl; exit(1);
	}//*/
}

void World::TrackInfectedIndividuals()
{
	wt_susceptible_mhc_susceptible = 0;
	wt_susceptible_mhc_chronic = 0;
	wt_susceptible_mhc_immune = 0;

	wt_chronic_mhc_susceptible = 0;
	wt_chronic_mhc_chronic = 0;
	wt_chronic_mhc_immune = 0;

    wt_immune_mhc_susceptible = 0;
    wt_immune_mhc_chronic = 0;
    wt_immune_mhc_immune = 0;

	vector<Host>::iterator it_host = hosts.begin();
	//check what kind of infections the hosts have
	list<Infection>:: iterator it;

	/*THIS HAS BEEN WRITTEN TO THE SPECIAL CASE WITH MAXIMALLY TWO INFECTIONS PER HOST!!!!!!!!!!!!!!!!!!*/
	while(it_host != hosts.end())
	{
		// check how many infections each host has
		int chronic_and_acute_infections= it_host->CountInfections();
		int total_infections = it_host->infections.size();

		if(total_infections == 0) //there are no infections: susceptible to both infections!
			wt_susceptible_mhc_susceptible ++;

		if(chronic_and_acute_infections == 0) //if there are no acute/chronic ifnections
		{
			if(total_infections == 1)
			{	//if there is only one infection check which type of virus it is
				it = it_host->infections.begin();

				if(it->pathogen.IsWildType())
				{
					if(it->IsImmune()) //then check what kind of infection it is
					{
						wt_immune_mhc_susceptible++;
					}
				}
				if(it->pathogen.IsDownregulatingMHC())
				{
					if(it->IsImmune())
					{
						wt_susceptible_mhc_immune++;
					}
				}
			}
			if (total_infections == 2) //if there are two infections.. then it means that both states are immune
			{
				wt_immune_mhc_immune++;
			}
		}
		if(chronic_and_acute_infections == 1) //check what type the chronic infection is and what is happening to the other
		{
			if(total_infections == 1)
			{
				it = it_host->infections.begin();
				if(it->pathogen.IsWildType())
				{
					wt_chronic_mhc_susceptible++;
				}
				if(it->pathogen.IsDownregulatingMHC())
				{
					wt_susceptible_mhc_chronic ++;
				}
			}
			if(total_infections == 2)
			{
				it = it_host->infections.begin();
				if(it->pathogen.IsWildType())
				{
					if(it->IsImmune())
						wt_immune_mhc_chronic++;
					else
						wt_chronic_mhc_immune++;
				}
				if(it->pathogen.IsDownregulatingMHC())
				{
					if(it->IsImmune())
						wt_chronic_mhc_immune++;
					else
						wt_immune_mhc_chronic++;
				}
			}
		}
		if(chronic_and_acute_infections == 2)
		{
				//both infections are chronic
				wt_chronic_mhc_chronic ++;
		}

		it_host++;
	}
}

void World::RemoveDeadHosts_HappyNewYear()
{	
//	int deathcount=0;
	vector<Host>::iterator it = hosts.begin();
	while(it!=hosts.end())
	{
		if(it->IsDead())
		{
			it = hosts.erase(it);
//			deathcount++;
			continue;
		}
		it++;
	}
}

/*Functions which do the OUTPUT files*/

/*this function saves the genes of each host -> to keep track of MHC and KIR diversity*/
void World:: SaveGenes()
{
	stringstream ss;
	ss << simulationTime/YEAR;
	string s(ss.str());
	s+=".Genes.log";
	
	genesFile.open(s.c_str(), ios::out);

	genesFile << "#Host Index\t";
	for(unsigned int i = 0; i< hosts.at(0).mhcGenes.size(); i++)
	{
		genesFile << "Mhc"<< i+1 <<"\t";
	}

	for(unsigned int i = 0; i< hosts.at(0).kirGenes.size(); i++)
	{
		genesFile << "Kir"<< i+1 <<"\t";
	}

	genesFile << "\tfunctional Kirs\texpressed Kirs\t host_id \n";
	int index = 0;
	vector<Host>::iterator hostIt;
	for(hostIt = hosts.begin(); hostIt != hosts.end(); hostIt ++)
	{
		genesFile << index << "\t";
		hostIt->SaveGenes(genesFile);
		index ++;
	}

	genesFile.close();
}

/*This function keeps track (and saves) of the population size*/
void World::SavePopulationSize()
{
	//populationSize <<simulationTime/YEAR <<"\t"<< hosts.size() <<"\t"<<simpleInfection << "\t"<< doubleInfection<< "\t"<< tripleInfection<< "\t"<<wildtype<<"\t"<<wildtype_immune << "\t"<<downregulating<< "\t"<<downregulating_immune<<"\t"<< decoy <<"\t"<<decoy_immune<<"\n";
	populationSize <<simulationTime/YEAR <<"\t"<< hosts.size() <<"\t" <<wt_susceptible_mhc_susceptible<<"\t"<<wt_susceptible_mhc_chronic<< "\t"<<wt_susceptible_mhc_immune<< "\t"<<wt_chronic_mhc_susceptible<<"\t"<< wt_chronic_mhc_chronic<<"\t"<<wt_chronic_mhc_immune<<"\t"<<wt_immune_mhc_susceptible<<"\t"<<wt_immune_mhc_chronic<<"\t"<<wt_immune_mhc_immune<<"\t";
	populationSize <<  best_protection <<"\t" << high_protection <<"\t" << normal_protection <<"\t" <<medium_protection << "\t" << low_protection << "\t"<< zero_protection <<"\n";
}

/*This function keeps track (ans saves) several parameters of the host and virus*/
void World::SaveParameters()
{
	stringstream ss;
	ss << simulationTime/YEAR;
	string s(ss.str());
	s+=".log";
	parameterFile.open(s.c_str(), ios::out);
	parameterFile << "#HostIndex\t Age\t totalInfections\t infectionTime\t infectionType\t virusType\t viralLoad\t originalViralLoad\t decoyID\t onlyAcute\t mhcID_down\t specificMHCdown\t Host_ID\n";
	vector<Host>::iterator hostIt;
	int index = 0;
	for(hostIt = hosts.begin(); hostIt != hosts.end(); hostIt ++)
	{
		parameterFile << index <<"\t";
		hostIt->SaveParameters(parameterFile);
		index ++;
	}

	parameterFile.close();
}

void World::SaveMap()
{
	stringstream ss;
	ss << simulationTime/YEAR;
	string s(ss.str());
	s+=".Map.log";
	mapFile.open(s.c_str(), ios::out);

	map< pair <int,int>, pair <int, int> > ::iterator it;
	for(it = KIRGenesMap.GetMap().begin(); it != KIRGenesMap.GetMap().end(); it ++)
	{
		mapFile << (*it).first.first << "|" <<(*it).first.second << "|" << (*it).second.first <<"|"<< (*it).second.second <<endl;
	}
	mapFile.close();
}

bool World::SaveAgeDyingHosts(double lastOutfileTime, double lastStopOutfileTime)
{

	bool timeToPrintOut = floor((simulationTime-lastOutfileTime)*outfileRate);
	bool timeToStopOut = floor((simulationTime-lastStopOutfileTime)*outfileRate);

	if(timeToPrintOut>0)
	{
		isFileOpen = true;
		stringstream ss;
		ss << simulationTime/YEAR;
		string s(ss.str());
		s+=".Age.log";
		dyingHosts.open(s.c_str(), ios::out);
		dyingHosts << "#Age \t KIR genes\t virus type\t viral load\t decoy ID\t onlyAcute?\t infection type\t infection time\t clearance time\n";
	}

	if(isFileOpen)
	{
		vector<Host>::iterator hostIt;
		//int index = 0;
		for(hostIt = hosts.begin(); hostIt != hosts.end(); hostIt ++)
		{
			if (hostIt->IsDead())
			{
				//write in a file their age, gene content, and infection type
				hostIt->SaveAgeDyingHost(dyingHosts);
				dyingHosts << simulationTime/YEAR <<"\n";
			}
			//index++;
		}
		if(timeToStopOut>0)
		{
			isFileOpen = false;
			dyingHosts.close();
			return true;
		}
		else
			return false;
	}
	return false;
}


/*Functions which do the BACKUP files*/

void World::SaveBackupFile()
{
	stringstream ss;
	ss << simulationTime/YEAR;
	string s(ss.str());
	s+=".Backup.data";
	backupFile.open(s.c_str(), ios::out);

	try
	{
		backupFile.is_open();
	}
	catch (...)
	{
		cout << "Backup: Could not Open File" <<endl;
	}

	cout <<"Saving backup file\n"<<endl;

	/*if(!backupFile.is_open())
	{
		throw OussException("Backup: Could not Open File");
	}
	else
	{
		cout <<"Saving backup file\n"<<endl;
	}*/

	backupFile << simulationTime<< "\t "<< MHCPool.GetPoolSize() << "\t";
	//save mhc gene pool
	for(unsigned int i = 0; i<MHCPool.GetPoolSize(); i++)
	{
		backupFile << MHCPool.GetGenes().at(i) << "\t";
	}
	backupFile << "\n";

	//save the KIRMap
	KIRGenesMap.SaveBackupMap(backupFile);
	backupFile << "\n";
	vector<Host>::iterator hostIt;
	for(hostIt = hosts.begin(); hostIt != hosts.end(); hostIt ++)
	{
		hostIt->SaveBackupHost(backupFile);
	}
	backupFile.close();
}

void World ::LoadBackupFile(const string& fileName)
{
	CreateBirthAndDeathRates();
	cout << "in the loading function now..."<<endl;
	int l = 0;
	int poolSize = 0;
	int mapSize = 0;
	string sline;
	vector<Host>::iterator hostIt;
	cout <<"opening the back up file...."<<endl;
	backupFile.open(fileName.c_str(), ios::in);
	if(!backupFile.is_open())
	{
		throw OussException("Backup: Could not Open File");
	}
	cout << "back up file open... reading it now..."<<endl;
	while(!backupFile.eof())
	{
		if(l==0)
		{
			cout <<"getting the sim time and the pool size...."<<endl;
			getline(backupFile,sline);
			stringstream ssline(sline);
			ssline >>simulationTime;
			ssline >>poolSize;
			//restore gene pools
			cout <<"done! now getting the MHC genes ...."<<endl;
			for(int i = 0; i<poolSize; i++)
			{
				int mhc_gene;
				ssline>>mhc_gene;
				MHCPool.GetGenes().push_back(mhc_gene);
			}
			// get the second line with the information of the Map
			cout <<"done! now getting the Map with KIR genes ...."<<endl;
			string secondline;
			string mline;
			getline(backupFile, secondline);
			stringstream msline (secondline);
			msline >> mapSize;
			for(int i= 0; i<mapSize; i++)
			{
				mline = KIRGenesMap.RestoreMap(msline);
				msline.str() = mline;
			}
		}

		else
		{
			cout <<"done! and now getting the infos from the Hosts...."<<endl;
			string sline;
			getline(backupFile, sline);
			if(sline.size()!=0)
			{
				//cout <<"sline1:  "<<endl;
				//cout<<sline <<endl;
				Host tempHost;
				tempHost.RestoreHost(sline);
				//cout <<"restored: "<<endl;
				//tempHost.PrintParametersHost();sleep(10);
				hosts.push_back(tempHost);

			}
		}
		l++;
	}
	backupFile.close();
}

double World :: GetIntrinsicDeathRate(const double age, const double viralDeathRate)const
{
	unsigned int a = round(age*52.0);
	if(a >= deathRates.size())
		cout << a << "|" << deathRates.size() << "age" << "ERROR in the deathrate!!!!" <<endl;

	double intrinsicDeathRate = deathRates.at(a) + viralDeathRate;
	// 	cout << "a, age, rate at, intrinsicdeathRates: "<< a << ", " << age << " , " << rates.at(a) << " , "<<intrinsicDeathRate << endl;
	/*
	 * CHANGE SOMETHING HERE!!! i NEED THIS FUNCTION:
	 * (n+1)*delta + k, where n is the number of chronic infections per host and k is the acute viral load
	 *
	 */

	return intrinsicDeathRate;
}

double World :: GetAgeDependentBirthRate(const double age)const
{
	unsigned int a = round(age*52.0);
	if(a >= birthRates.size())
		cout << a << "|" << birthRates.size() << "age" << "ERROR in the birth rate!!!!" <<endl;
	double ageDependentBirthrate = birthRates.at(a);
	return ageDependentBirthrate;
}
