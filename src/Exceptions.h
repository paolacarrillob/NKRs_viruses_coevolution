/*
 * Exceptions.h
 *
 *  Created on: Nov 29, 2010
 *      Author: paola
 */

#ifndef EXCEPTIONS_H_
#define EXCEPTIONS_H_

/*
 *  Exceptions.h
 *  PluginFramework
 *
 *  Created by Oussama Jarrousse on 29/07/2009.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include <exception>
#include <iostream>
#include <string>

using namespace std;

class OussException: public exception
{
protected:
	string errorData;
public:
	//! Constructor with arguments similiar to printf
	OussException(std::string _errorData)
	{
		errorData = _errorData;
	}

	virtual ~OussException() throw() {};

	//! Get exception text
	string & GetErrorData(){return errorData;}

	friend ostream & operator<<(ostream &s, const OussException & e);
};

//! Output operator for class kaBaseException. The errstr is printed.
inline ostream & operator<<(ostream &s, const OussException & e)
{
	return s << e.errorData;
}

class OussExceptionEmptyString: public OussException
{
public:
	OussExceptionEmptyString(std::string _errorData):OussException(_errorData){};
	virtual ~OussExceptionEmptyString()throw(){};
};

class OussExceptionFileNotFound: public OussException
{
public:
	OussExceptionFileNotFound(std::string _errorData):OussException(_errorData){};
	virtual ~OussExceptionFileNotFound()throw(){};
};

#endif /* EXCEPTIONS_H_ */
