/*
 * backup.h
 *
 *  Created on: Nov 30, 2010
 *      Author: paola
 */

#ifndef BACKUP_H_
#define BACKUP_H_
//#include <kaExceptions.h>
#include "Exceptions.h"
#include <fstream>
#include <string.h>
#include <stdio.h>
#include <sstream>

class kaBackup
{
public:
	kaBackup(const std::string& sFilename,bool bWrite);//ok
	virtual ~kaBackup();//ok

	template<class V> void Save(const std::string& sName,V Var);//!< Saves Var in Section sName
	void SaveString(const std::string& sName,const std::string& sBuff);//!< Saves sBuff String in Section sName (doesn't search for the string first)
	template<class V> void Load(const std::string& sName,V* refVal);//ok

	#warning to replace the buffer here with somthing else!
	void LoadString(const std::string& sName,char* sBuff,int iBufflength);//ok

	void LoadBlock(const std::string& sName);//!< Loads The hole Section Block into dynamically allocated memory
	template<class V> bool ExtractIndexedParameter(int* pIndex,V* refVal);

	bool IsBlockEmpty(){return sBlock.empty();}
	void PrintBlock();

	void SaveComment(const std::string& buff,char c='#');
	void SaveLineFeed(){WriteString("");}

	bool CheckWriteMode(){return bWriteMode;}
	bool CheckFileAttached(){return bFileAttached;}

	void DetachFile();
	void ReAttachFile(bool _bWriteMode){if(!bFileAttached)AttachFile(sFile.c_str(),_bWriteMode);}

	std::string GetFileName(){return sFile;}
	void AttachFile(const std::string& sFilename, bool bWrite);

	bool FindNextSection();//<! Search for the next section in the file.\n
	//<! If no section was found and eof is reached sSection is cleared and false is returned.\n
	//<!	If a section is found, sSection is set to the found section and the function returns true, the file pointer points to the line right after the section;

	bool FindNextSectionAndReadBlock();

	std::string GetSection(){return sSection;}
	std::string GetBlock(){return sBlock;}
protected:
	void RoutineFileCheck();//ok
	void RoutineWriteCheck();//ok
	void RoutineReadCheck();//ok

	bool IgnoreComment();

	void WriteSection(const std::string& sName);//bug fixed ok
	void FindSection(const std::string& sName);//ok

	//template<class V> void Write(V Val){throw kaBaseException("Backup: No suitable handler for passed datatype!!!");}//ok
	template<class V> void Write(V Val){throw OussException("Backup: No suitable handler for passed datatype!!!");}//ok
	//template<class V> void Read(V* refVal){throw kaBaseException("Backup: No suitable handler for passed datatype!!!");}//ok
	template<class V> void Read(V* refVal){throw OussException("Backup: No suitable handler for passed datatype!!!");}//ok

	void WriteString(const std::string& sBuff);//ok

	#warning replace the char* buffer with string& or something... or replace the hole function
	void ReadString(char* sBuff,int iValLength);//ok
	void ReadBlock();//ok

	//FILE* fBackup;
	fstream fBackup;
	bool bWriteMode;
	bool bFileAttached;

	std::string sBlock;
	std::string sSection;
	std::string sFile;
};

#include "backup.hpp"
#endif /* BACKUP_H_ */
