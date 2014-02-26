//============================================================================
// Name        : CoevolutionKIR.cpp
// Author      : P.Carrillo-Bustamante
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <string.h>
#include "World.h"
#include "Exceptions.h"
using namespace std;

int main(int argc, char*argv[])
{
/*	//testing the distribution of Lc,NKR<->MHC
	vector<int> Lcs;

	MHCGenePool mhc;
	//mhc.FillMHCGenePool(30);
	mhc.FillMHCGenePoolWithSimilarMHCs(30);
	double size = mhc.GetPoolSize();
	for(unsigned int j = 0; j <mhc.GetPoolSize(); j++)
	{
		Gene dummy_mhc;
		dummy_mhc.SetGeneID(mhc.GetGenes().at(j));

		for(int i = 0; i<10000; i++)
		{
			KIRGene dummy(2);
			int dummy_Lc = dummy.BindMolecule(dummy_mhc);
			Lcs.push_back(dummy_Lc);
		}


	}

	for(int i = 0; i<=16; i++)
	{
		double mycount = count(Lcs.begin(), Lcs.end(), i);
		mycount = mycount/10000;
		double percentage_mhc = mycount/size;

		cout << i << " " << mycount <<" "<<percentage_mhc<<endl;
	}


	exit(-1);*/

	cerr << "Testing if stderr is redirected to stdoutput" << endl;
	if (argc<2)
	{
		cerr << "Usage: "<<argv[0] << " <Parameter file> <-b Loading Backup file> \n CoevolutionKIR simulates the evolution of the complex KIR system. For correct usage you have to indicate a parameter file and a backup file name. If you are not loading any backup, just give any random name"<< endl;
		exit(-1);
	}
	// delete the existing host_file
	char buffer[512];
	//Ouss: Use boost library here for boost:filesystem to delete files... etc
	cout <<"deleting old files ..."<<endl;
	//sprintf(buffer, "rm *.txt");
	sprintf(buffer, "rm *.log");
	system(buffer);

	string parameterFile(argv[1]);
	string backupFile;
	bool loadingBackup=false;
	for(int i=2; i< argc;i++)
	{
		if(strcmp(argv[i],"-b")==0 && argc > i)
		{
			loadingBackup = true;
			backupFile.assign(argv[++i]);
		}
	}

	World theWorld;
	try
	{
		theWorld.LoadParameterFile(parameterFile);
	}
	catch (OussException& e ) {
		cout << e.GetErrorData() << endl;
		exit(-1);
	}
	catch (...)
	{
		cout << "unknown exception thrown" <<endl;
		exit(-1);
	}
	
	//theWorld.LoadParameterFile(parameterFile);
	cout <<"welcome, hello" << endl;
	if(!loadingBackup)
	{
		// initialize host population
		theWorld.Initialize();
	}
	else
	{
		cout << "\n Loading data from backup file: "  << backupFile << endl;
		theWorld.LoadBackupFile(backupFile);
	}

	theWorld.Simulate();
	cout << "bye bye \n";

	return 0;
}


/*
 *24.10.2012 fixed bug in LoadBackupFile! -> fixed it again (07.10.2013)... learned about ss.clear() and getline :-)
 *
 */
