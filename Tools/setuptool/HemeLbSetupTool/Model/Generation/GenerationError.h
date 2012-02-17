#ifndef HEMELBSETUPTOOL_GENERATIONERROR_H
#define HEMELBSETUPTOOL_GENERATIONERROR_H
#include <exception>
#include <string>

struct GenerationError: public std::exception {
	virtual const char* what() const throw () {
		return "GenerationError";
	}
};

struct GenerationErrorMessage: public GenerationError {
	GenerationErrorMessage(const std::string errorMessage) :
		msg(errorMessage) {
	}
	~GenerationErrorMessage() throw () {
	}
	virtual const char* what() const throw () {
		return msg.c_str();
	}

	const std::string msg;
};

#endif
