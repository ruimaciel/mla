#ifndef MLA_EXCEPTION
#define MLA_EXCEPTION

#include <string>
#include <stdexcept>


/**
  * class LAException
  * This class is intended to be used by all <mla routines
  **/
class LAException
	: public std::runtime_error
{
public:
	explicit LAException(const std::string &message): std::runtime_error(std::string("LA: "+message)) {}
};



#endif
