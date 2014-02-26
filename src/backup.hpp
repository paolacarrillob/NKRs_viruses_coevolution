#include "backup.h"

kaBackup::kaBackup(const std::string& sFilename, bool bWrite):bFileAttached(false)
{
	try
	{
		DetachFile();
		AttachFile(sFilename,bWrite);
	}
	catch(...)
	{
		throw;
	}
}

void kaBackup::AttachFile(const std::string& sFilename, bool bWrite)
{
	if(sFilename.empty())
	{
		//throw kaBaseException("Backup: File name is mandatory");
		throw OussException("Backup: File name is mandatory");
	}
	string temp(sFilename);
	bWriteMode = bWrite;
	if(bWrite)
	{
		fBackup.open(sFilename.c_str(), ios::out);
		//fBackup = fopen(sFilename.c_str(),"w");
		if(!fBackup.is_open())
		{
			//throw kaBaseException("Backup: Could not Open File");
			throw OussException("Backup: Could not Open File");
		}
		sFile = temp;
		bFileAttached = true;
	}
	else
	{
		//fBackup = fopen(sFilename.c_str(),"r");
		fBackup.open(sFilename.c_str(), ios::in);
		if(!fBackup.is_open())
		{
			//throw kaBaseException("Backup: Could not Open File");
			throw OussException("Backup: Could not Open File");
		}
		sFile = temp;
		bFileAttached = true;
	}
	sBlock.clear();
}


kaBackup::~kaBackup()
{
	//FreeBlock();
	DetachFile();
}

void kaBackup::DetachFile()
{
	if(bFileAttached)
	{
		//fclose(fBackup);
		fBackup.close();
		bFileAttached = false;
	}
}

void kaBackup::RoutineFileCheck()
{
	if(!bFileAttached)
	{
		//throw kaBaseException("Backup: File Detached!");
		throw OussException("Backup: File Detached!");
	}
	if(!fBackup.is_open())
	{
		//throw kaBaseException("Backup: Corrupt File Pointer");
		throw OussException("Backup: Corrupt File Pointer");
	}
	//if(ferror(fBackup))
	if(fBackup.bad())
	{
		//throw kaBaseException("Backup: File Error !!! FILE Error :-O");
		throw OussException("Backup: File Error !!! FILE Error :-O");
	}
}

void kaBackup::RoutineWriteCheck()
{
	if(!bWriteMode)
	{
		//throw kaBaseException("Backup: Trying to write in Read Mode");
		throw OussException("Backup: Trying to write in Read Mode");
	}
	RoutineFileCheck();
}
void kaBackup::RoutineReadCheck()
{
	if(bWriteMode)
	{
		//throw kaBaseException("Backup: Trying to read in Write Mode");
		throw OussException("Backup: Trying to read in Write Mode");
	}
	RoutineFileCheck();
}


void kaBackup::SaveString(const std::string& sName,const std::string& sBuff)
{
	WriteSection(sName);
	WriteString(sBuff);
}

template<class V>
void kaBackup::Save(const std::string& sName,V bVal)
{
	WriteSection(sName);
	Write<V>(bVal);
}

void kaBackup::SaveComment(const std::string& comm,char c)
{
	char sBuff[512];
	sprintf(sBuff,"%c%s",c,comm.c_str());
	WriteString(sBuff);
}

void kaBackup::WriteSection(const std::string& sName)
{
	//#To do to check if the Section Exists, if it does, append the new value directly after!
	RoutineWriteCheck();
	//fprintf(fBackup,"[%s]\n",sName);
	fBackup << "[" << sName << "]\n";
}

template <>
void kaBackup::Write(bool bVal)
{
	RoutineWriteCheck();
	//fprintf(fBackup,"%d\n",(bVal?1:0));
	fBackup << int((bVal?1:0)) << "\n";
}

template <>
void kaBackup::Write(double dVal)
{
	RoutineWriteCheck();
	//fprintf(fBackup,"%e\n",dVal);
	fBackup << scientific << dVal << "\n";
}

template <>
void kaBackup::Write(unsigned int iVal)
{
	RoutineWriteCheck();
	//fprintf(fBackup,"%d\n",iVal);
	fBackup << iVal <<"\n";
}

template <>
void kaBackup::Write(float fVal)
{
	RoutineWriteCheck();
	//fprintf(fBackup,"%f\n",fVal);
	fBackup << fVal << "\n";
}

template <>
void kaBackup::Write(int iVal)
{
	RoutineWriteCheck();
	//fprintf(fBackup,"%d\n",iVal);
	fBackup << iVal << "\n";
}


void kaBackup::WriteString(const std::string& sBuff)
{
	RoutineWriteCheck();
	//fprintf(fBackup,"%s\n",sBuff);
	fBackup << sBuff << "\n";
}

template <>
void kaBackup::Write(unsigned long int iVal)
{
	RoutineWriteCheck();
	//fprintf(fBackup,"%d\n",iVal);
	fBackup << iVal << "\n";
}

template <>
void kaBackup::Write(long int iVal)
{
	RoutineWriteCheck();
	//fprintf(fBackup,"%d\n",iVal);
	fBackup << iVal << "\n";
}

template <>
void kaBackup::Write(unsigned char cVal)
{
	RoutineWriteCheck();
	//fprintf(fBackup,"%d\n",cVal);
	fBackup << cVal << "\n";
}

void kaBackup::LoadString(const std::string& sName,char* sVal,int iValLength)
{
	RoutineReadCheck();
	FindSection(sName);
	if(!IgnoreComment())//Do the Actual Reading to sBlock
	{
		return;
	}
	ReadString(sVal,iValLength);
}
void kaBackup::LoadBlock(const std::string& sName)
{
	//FreeBlock();
	sBlock.clear();
	RoutineReadCheck();
	FindSection(sName);
	ReadBlock();
}


template<class V>
void kaBackup::Load(const std::string& sName,V* refVal)
{
	RoutineReadCheck();
	FindSection(sName);
	if(!IgnoreComment())//Do the Actual Reading to sBlock
	{
		return;
	}
	Read(refVal);
}

bool kaBackup::FindNextSectionAndReadBlock()
{
	if(!FindNextSection())
		return false;
	sBlock.clear();
	ReadBlock();
	return true;
}

bool kaBackup::FindNextSection()
{
	sSection.clear();
	char sLine[512]="";
	//while (!feof(fBackup))//Search until the end of file is reached...
	while (!fBackup.eof())//Search until the end of file is reached...
	{
		//fgets(sLine,512,fBackup);//read one line from file
		fBackup.getline(sLine,512);//read one line from file
		if(sLine[0]=='[')
		//when [ is found
		{
			//read the hole line and remove the newline char...
			if(sLine[strlen(sLine)-1]=='\n')
				sLine[strlen(sLine)-1]='\0';//remove the newline char
										//fscanf(fBackup,"%s",sLine);
			if(sLine[strlen(sLine)-1]!=']')
			{
				if(strlen(sLine)<512)
					sLine[strlen(sLine)] = ']';
				if(strlen(sLine)<512)
					sLine[strlen(sLine)] = '\0';
			}
			//assign the sSection to the value...
			sSection.assign(sLine);
			//return true;
			return true;
		}
		//iterate
	}
	return false;
}

void kaBackup::FindSection(const std::string& sName)
{
	sSection.clear();
	char sLine[512]="";
	//char sBuff[512]="";
	string sBuff;
	//sprintf(sBuff,"[%s]",sName);//fill the buffer with the section name enclosed with the middle braces [].
	sBuff.clear();
	sBuff = "[";
	sBuff+=sName;
	sBuff+="]";

	//if(feof(fBackup))//if the file pointer points to the end of the file rewind the file to start searching from the beginning
	if(fBackup.eof())//if the file pointer points to the end of the file rewind the file to start searching from the beginning
	{
		//rewind(fBackup);
		fBackup.clear();
		fBackup.seekg(0, ios_base::beg);
	}
	//The position the search started from
	//long int init_pos;
	int init_pos;
	//init_pos = ftell(fBackup);
	init_pos = fBackup.tellg();
	//while (!feof(fBackup))//Search until the end of file is reached...
	while (!fBackup.eof())//Search until the end of file is reached...
	{
		//fgets(sLine,512,fBackup);//read one line from file
		fBackup.getline(sLine,512);//read one line from file
		if(sLine[strlen(sLine)-1]=='\n')//if the last char is a newline char
			sLine[strlen(sLine)-1]='\0';//replace it with the null char
		//fscanf(fBackup,"%s",sLine);
		if(strcmp(sBuff.c_str(),sLine)==0)//compare it with the section name
		{
			sSection.assign(sName);
			return;//the section is found the pointer points to the line right after the section.
		}
	}

	fBackup.clear();//if the section was not found rewind the file and search until the initial position of the search.
	fBackup.seekg(0, ios_base::beg);
	//rewind(fBackup);//if the section was not found rewind the file and search until the initial position of the search.
	int curr_pos;
	//long int curr_pos;
	//curr_pos = ftell(fBackup);
	curr_pos = fBackup.tellg();
	while (curr_pos != init_pos)
	{
		//fgets(sLine,512,fBackup);
		fBackup.getline(sLine,512);

		if(sLine[strlen(sLine)-1]=='\n')//if the last char is a newline char
			sLine[strlen(sLine)-1]='\0';//replace it with the null char

		if(strcmp(sBuff.c_str(),sLine)==0)
		{
			sSection.assign(sName);
			return;//the section is found the pointer points to the line right after the section.
		}
		//curr_pos = ftell(fBackup);
		curr_pos = fBackup.tellg();
	}
	//throw kaBaseException("Backup: Section Not Found in opened file");//an exception is thrown if the section was not found.
	throw OussException("Backup: Section Not Found in opened file");//an exception is thrown if the section was not found.
}

bool kaBackup::IgnoreComment()
{
	//if(feof(fBackup))
	if(fBackup.eof())
	{
		//throw kaBaseException("Backup: EOF right after the section.");
		throw OussException("Backup: EOF right after the section.");
	}
	char sLine[512]="";
	do
	{
		//fgets(sLine,512,fBackup);
		fBackup.getline(sLine,512);
	}
	//while((sLine[0]=='#' || (sLine[0]=='/'&&sLine[1]=='/') ) && !feof(fBackup));
	while((sLine[0]=='#' || (sLine[0]=='/'&&sLine[1]=='/') ) && !fBackup.eof());
	//if(feof(fBackup))
	if(fBackup.eof())
	{
		fBackup.clear();
		//throw kaBaseException("Backup: EOF reached and the value was not found.");
		throw OussException("Backup: EOF reached and the value was not found.");
	}
	if(sLine[0]=='[')
	{
		//throw kaBaseException("Backup: Next Section reached and the value was not found.");
		throw OussException("Backup: Next Section reached and the value was not found.");
	}
	if(sLine[strlen(sLine)-1]=='\n')//if the last char is a newline char
		sLine[strlen(sLine)-1]='\0';//replace it with the null char
	sBlock.clear();
	sBlock.assign(sLine);
	return true;
}

template<>
void kaBackup::Read(unsigned int* iVal)
{
	RoutineReadCheck();
	*iVal = atoi(sBlock.c_str());
}

template<>
void kaBackup::Read(bool* bVal)
{
	RoutineReadCheck();
	int i=0;
	i = atoi(sBlock.c_str());
	*bVal = (i?true:false);
}//ok

template<>
void kaBackup::Read(double* dVal)
{
	RoutineReadCheck();
	//double d=0;
	//sscanf(sBlock.c_str(),"%f",&d);
	//cerr << sBlock.c_str() << endl;
	*dVal = atof(sBlock.c_str());
	//d = strtod(sBlock.c_str(),NULL);
	//*dVal = double(d);
}//ok

template<>
void kaBackup::Read(float* fVal)
{
	RoutineReadCheck();
	//float d=0;
	//sscanf(sBlock.c_str(),"%f",&d);
	//*fVal = d;
	*fVal = atof(sBlock.c_str());
}//ok

template<>
void kaBackup::Read(unsigned char* cVal)
{
	RoutineReadCheck();
	//unsigned char c=0;
	//sscanf(sBlock.c_str(),"%c",&c);
	*cVal = sBlock[0];
}

template<>
void kaBackup::Read(int* iVal)
{
	RoutineReadCheck();
	//int i=0;
	//sscanf(sBlock.c_str(),"%d",&i);
	//*iVal = i;
	*iVal = atoi(sBlock.c_str());
}//ok

template<>
void kaBackup::Read(long int* iVal)
{
	RoutineReadCheck();
	//unsigned long int i=0;
	//sscanf(sBlock.c_str(),"%d",&i);
	//*iVal = i;
	*iVal = atol(sBlock.c_str());
}//ok

template<>
void kaBackup::Read(unsigned long int* iVal)
{
	RoutineReadCheck();
	//unsigned long int i=0;
	//sscanf(sBlock.c_str(),"%d",&i);
	//*iVal = i;
	*iVal = atol(sBlock.c_str());
}//ok

void kaBackup::ReadString(char* sVal,int iValLength)
{
	RoutineReadCheck();
	if(sBlock.size()>iValLength)
	{
		//throw kaBaseException("Backup: Buff Size too small");
		throw OussException("Backup: Buff Size too small");
	}
	strcpy(sVal,sBlock.c_str());
//	sVal[strlen(sVal)-1]='\0';//remove the newline char
}

void kaBackup::ReadBlock()
{
	RoutineReadCheck();
	char sBuff[512] = "";
	bool bNoReasonToStop = true;
	//while(bNoReasonToStop && !feof(fBackup))
	int pos=0;
	while(bNoReasonToStop && !fBackup.eof())
	{
		//fgets(sBuff,512,fBackup);
		pos = fBackup.tellg();
		fBackup.getline(sBuff,512);
		//add a new line
		if(strlen(sBuff)<512)
			strcat(sBuff,"\n");
		else
		{
			cerr << "hehehe" << endl;
			//throw kaBaseException("Backup: Buff Size too small");
			throw OussException("Backup: Buff Size too small");
		}

		if(sBuff[0]=='#')
		{
			continue;
		}
		if(sBuff[0]=='/' && sBuff[1]=='/' )
		{
			continue;
		}
		if(sBuff[0]=='\n')
		{
			continue;
		}
		if(sBuff[0]=='[')
		{
			bNoReasonToStop = false;
			if(fBackup.eof())
				fBackup.clear();
			fBackup.seekg(pos,ios_base::beg);
			break;
		}
		//add the line to the buffer
		sBlock.append(sBuff);
	}
	return;
}

void kaBackup::PrintBlock()
{
	if(!sBlock.empty())
	{
		cerr << sBlock << endl;
	}
}

template<class V>
bool kaBackup::ExtractIndexedParameter(int* pIndex, V* refVal)
{
	if(sBlock.empty())
	   return false;
	stringstream ssBlock(sBlock);
	ssBlock >> *pIndex >> *refVal;
	//find the new line char and tellg()
	ssBlock.ignore(512,'\n');
	int length = ssBlock.tellg();

	if(length < 0)
	{
		//no end of line found... maybe at the end of the file...
		//if(feof(fBackup))
		if(fBackup.eof())
		{
			//erase the hole block.
			sBlock.erase(0);
			return true;
		}
		//the first part of the Block could be extracted, the second not!
		cerr << " Could not extract Indexed Parameter at "<< sSection << " !!!" << endl;
		return false;
	}

	sBlock.erase(0,length);
	return true;
}
